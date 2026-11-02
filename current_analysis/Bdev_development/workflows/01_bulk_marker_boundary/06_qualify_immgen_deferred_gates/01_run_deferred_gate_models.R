#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(jsonlite)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop("Missing argument: ", flag)
  args[[index + 1L]]
}

config_path <- arg_value("--config")
output_root <- arg_value("--output-root")
config <- yaml::read_yaml(config_path)
if (!identical(config$schema_version, "bdev.immgen_deferred_gate_qualification.v1")) stop("Unexpected config schema")
config_dir <- dirname(normalizePath(config_path))
contract <- yaml::read_yaml(file.path(config_dir, config$inputs$statistical_contract))
input_root <- config$inputs$core_qualification_root
input_path <- function(name) file.path(input_root, config$inputs[[name]]$path)
output_path <- function(name) file.path(output_root, config$outputs$tables[[name]])

counts_frame <- read.delim(input_path("synthetic_counts"), check.names = FALSE, stringsAsFactors = FALSE)
gene_symbols <- counts_frame[[1L]]
counts <- as.matrix(counts_frame[, -1L, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- gene_symbols
metadata <- read.delim(input_path("synthetic_metadata"), check.names = FALSE, stringsAsFactors = FALSE)
plan <- read.delim(input_path("contrast_plan"), check.names = FALSE, stringsAsFactors = FALSE)
selected_sets <- read.delim(input_path("synthetic_selected_genesets"), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(colnames(counts), metadata$library_id), ncol(counts) == 126L)

families <- c("superlineage", "broad_lineage", "b_branch", "b_local_stage")
group_column <- c(superlineage = "superlineage", broad_lineage = "table_s1_lineage", b_branch = "b_branch", b_local_stage = "population")

write_tsv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(value, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_tsv_gz <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- gzfile(path, open = "wt", compression = 9)
  write.table(value, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  close(connection)
}

aggregate_by_population <- function(input_counts, input_metadata) {
  populations <- unique(input_metadata$population)
  aggregated <- vapply(populations, function(population) rowSums(input_counts[, input_metadata$population == population, drop = FALSE]), numeric(nrow(input_counts)))
  rownames(aggregated) <- rownames(input_counts)
  population_metadata <- do.call(rbind, lapply(populations, function(population) input_metadata[match(population, input_metadata$population), , drop = FALSE]))
  population_metadata$library_id <- population_metadata$population
  list(counts = aggregated, metadata = population_metadata)
}

family_subset <- function(input_counts, input_metadata, family) {
  selected <- input_metadata$role == "discovery"
  if (family %in% c("b_branch", "b_local_stage")) selected <- selected & input_metadata$table_s1_lineage == "B"
  if (family == "b_local_stage") {
    local_populations <- unique(c(plan$target[plan$family == family], unlist(lapply(plan$weights_json[plan$family == family], function(value) names(jsonlite::fromJSON(value))))))
    selected <- selected & input_metadata$population %in% local_populations
  }
  list(counts = input_counts[, selected, drop = FALSE], metadata = input_metadata[selected, , drop = FALSE])
}

normalized_logcpm <- function(input_counts) {
  y <- DGEList(counts = input_counts)
  y <- calcNormFactors(y, method = "TMM")
  cpm(y, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)
}

group_means <- function(logcpm, groups) {
  levels <- unique(groups)
  result <- vapply(levels, function(level) rowMeans(logcpm[, groups == level, drop = FALSE]), numeric(nrow(logcpm)))
  rownames(result) <- rownames(logcpm)
  result
}

effect_threshold <- function(family, test_type) {
  family_contract <- contract$contrast_families[[family]]
  if (family == "superlineage") return(family_contract$minimum_global_log2_fold_change)
  if (family == "broad_lineage") return(if (test_type == "global") family_contract$minimum_global_log2_fold_change else family_contract$minimum_each_competitor_log2_fold_change)
  if (family == "b_branch") return(if (test_type == "global") family_contract$minimum_global_log2_fold_change else family_contract$minimum_each_competitor_log2_fold_change)
  family_contract$minimum_each_neighbor_log2_fold_change
}

fdr_threshold <- function(family, test_type) {
  family_contract <- contract$contrast_families[[family]]
  if (family == "superlineage") return(family_contract$maximum_global_fdr)
  if (family %in% c("broad_lineage", "b_branch")) return(if (test_type == "global") family_contract$maximum_global_fdr else family_contract$maximum_each_competitor_fdr)
  family_contract$maximum_each_neighbor_fdr
}

contrast_effect <- function(means, weights_json) {
  weights <- jsonlite::fromJSON(weights_json)
  missing <- setdiff(names(weights), colnames(means))
  if (length(missing)) stop("Missing contrast groups: ", paste(missing, collapse = ","))
  as.numeric(means[, names(weights), drop = FALSE] %*% as.numeric(weights))
}

# Deterministic, balanced one-library-per-population schedule.
discovery_metadata <- metadata[metadata$role == "discovery", , drop = FALSE]
populations <- unique(discovery_metadata$population)
iterations <- config$balanced_resampling$iterations
set.seed(config$balanced_resampling$seed)
permutations <- lapply(populations, function(population) sample(sort(discovery_metadata$library_id[discovery_metadata$population == population])))
names(permutations) <- populations
schedule_rows <- vector("list", iterations * length(populations))
row_index <- 1L
for (iteration in seq_len(iterations)) {
  for (population in populations) {
    choices <- permutations[[population]]
    library_id <- choices[((iteration - 1L) %% length(choices)) + 1L]
    schedule_rows[[row_index]] <- data.frame(iteration = iteration, population = population, library_id = library_id, stringsAsFactors = FALSE)
    row_index <- row_index + 1L
  }
}
schedule <- do.call(rbind, schedule_rows)
write_tsv(schedule, output_path("resample_schedule"))

family_targets <- lapply(families, function(family) unique(plan$target[plan$family == family]))
names(family_targets) <- families
stability_counts <- list()
for (family in families) for (target in family_targets[[family]]) stability_counts[[paste(family, target, sep = "|")]] <- integer(nrow(counts))

for (iteration in seq_len(iterations)) {
  selected_ids <- schedule$library_id[schedule$iteration == iteration]
  iteration_metadata <- metadata[match(selected_ids, metadata$library_id), , drop = FALSE]
  iteration_counts <- counts[, selected_ids, drop = FALSE]
  for (family in families) {
    subset <- family_subset(iteration_counts, iteration_metadata, family)
    logcpm <- normalized_logcpm(subset$counts)
    means <- group_means(logcpm, subset$metadata[[group_column[[family]]]])
    family_plan <- plan[plan$family == family, , drop = FALSE]
    for (target in unique(family_plan$target)) {
      target_plan <- family_plan[family_plan$target == target, , drop = FALSE]
      passed <- rep(TRUE, nrow(counts))
      for (test_index in seq_len(nrow(target_plan))) {
        effects <- contrast_effect(means, target_plan$weights_json[[test_index]])
        passed <- passed & effects >= effect_threshold(family, target_plan$test_type[[test_index]])
      }
      key <- paste(family, target, sep = "|")
      stability_counts[[key]] <- stability_counts[[key]] + as.integer(passed)
    }
  }
}

balanced_rows <- do.call(rbind, lapply(names(stability_counts), function(key) {
  parts <- strsplit(key, "|", fixed = TRUE)[[1L]]
  data.frame(family = parts[[1L]], target = parts[[2L]], gene_symbol = gene_symbols,
             passing_iterations = stability_counts[[key]], iterations = iterations,
             stability_fraction = stability_counts[[key]] / iterations, stringsAsFactors = FALSE)
}))
write_tsv_gz(balanced_rows, output_path("balanced_stability"))

fit_target_gate <- function(input_counts, input_metadata, family, target) {
  subset <- family_subset(input_counts, input_metadata, family)
  family_counts <- subset$counts
  family_metadata <- subset$metadata
  if (family != "b_local_stage") {
    aggregated <- aggregate_by_population(family_counts, family_metadata)
    family_counts <- aggregated$counts
    family_metadata <- aggregated$metadata
  }
  groups <- factor(family_metadata[[group_column[[family]]]])
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) stop("Invalid holdout design for ", family, "/", target)
  y_all <- DGEList(counts = family_counts)
  keep <- filterByExpr(y_all, design = design, min.count = contract$count_processing$expression_filter$min_count,
                       min.total.count = contract$count_processing$expression_filter$min_total_count)
  y <- y_all[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y_all$samples$norm.factors <- y$samples$norm.factors
  y <- estimateDisp(y, design = design)
  fit <- glmQLFit(y, design = design, robust = TRUE)
  normalized_cpm <- cpm(y_all, normalized.lib.sizes = TRUE)
  target_plan <- plan[plan$family == family & plan$target == target, , drop = FALSE]
  passed <- keep
  support_values <- rep(1, nrow(y_all))
  for (test_index in seq_len(nrow(target_plan))) {
    weights <- jsonlite::fromJSON(target_plan$weights_json[[test_index]])
    contrast <- setNames(rep(0, ncol(design)), colnames(design))
    contrast[names(weights)] <- as.numeric(weights)
    test <- glmTreat(fit, contrast = contrast, lfc = contract$count_processing$model$treat_log2_fold_change)
    table <- topTags(test, n = Inf, sort.by = "none")$table
    table <- table[match(rownames(y), rownames(table)), , drop = FALSE]
    full_logfc <- rep(NA_real_, nrow(y_all)); full_logfc[keep] <- table$logFC
    full_fdr <- rep(NA_real_, nrow(y_all)); full_fdr[keep] <- table$FDR
    passed <- passed & !is.na(full_logfc) & full_logfc >= effect_threshold(family, target_plan$test_type[[test_index]]) &
      !is.na(full_fdr) & full_fdr <= fdr_threshold(family, target_plan$test_type[[test_index]])
    competitor_groups <- names(weights)[as.numeric(weights) < 0]
    target_samples <- groups == target
    competitor_samples <- groups %in% competitor_groups
    expression_cutoff <- if (family == "b_local_stage" && sum(target_samples) == 1L) contract$contrast_families$b_local_stage$singleton_target_minimum_cpm else 1
    expression_fraction <- rowMeans(normalized_cpm[, target_samples, drop = FALSE] >= expression_cutoff)
    competitor_median <- apply(normalized_cpm[, competitor_samples, drop = FALSE], 1L, median)
    direction_fraction <- rowMeans(normalized_cpm[, target_samples, drop = FALSE] > competitor_median)
    support_values <- pmin(support_values, expression_fraction, direction_fraction)
  }
  if (family == "superlineage") passed <- passed & support_values >= contract$contrast_families$superlineage$minimum_target_population_support
  if (family == "broad_lineage") passed <- passed & support_values >= contract$contrast_families$broad_lineage$minimum_target_population_support
  if (family == "b_branch") passed <- passed & support_values >= contract$contrast_families$b_branch$minimum_target_population_support
  if (family == "b_local_stage") passed <- passed & support_values >= contract$contrast_families$b_local_stage$minimum_target_library_fraction
  setNames(passed, rownames(y_all))
}

selected_gene_map <- split(selected_sets$gene_symbol, paste(selected_sets$family, selected_sets$target, sep = "|"))
population_loo_rows <- list()
population_loo_index <- 1L
for (family in c("superlineage", "broad_lineage", "b_branch")) {
  column <- group_column[[family]]
  subset <- family_subset(counts, metadata, family)
  for (target in family_targets[[family]]) {
    target_populations <- unique(subset$metadata$population[subset$metadata[[column]] == target])
    genes <- selected_gene_map[[paste(family, target, sep = "|")]]
    if (is.null(genes)) genes <- character()
    if (length(target_populations) < contract$stability$leave_one_population_out$minimum_target_populations) {
      for (gene in genes) {
        population_loo_rows[[population_loo_index]] <- data.frame(family = family, target = target, gene_symbol = gene, holdouts = 0L, passing_holdouts = NA_integer_, holdout_fraction = NA_real_, status = "not_assessable_low_replication", stringsAsFactors = FALSE)
        population_loo_index <- population_loo_index + 1L
      }
      next
    }
    pass_counts <- setNames(integer(length(genes)), genes)
    for (held_out_population in target_populations) {
      retained <- metadata$population != held_out_population
      passed <- fit_target_gate(counts[, retained, drop = FALSE], metadata[retained, , drop = FALSE], family, target)
      pass_counts <- pass_counts + as.integer(passed[genes])
    }
    for (gene in genes) {
      population_loo_rows[[population_loo_index]] <- data.frame(family = family, target = target, gene_symbol = gene, holdouts = length(target_populations), passing_holdouts = pass_counts[[gene]], holdout_fraction = pass_counts[[gene]] / length(target_populations), status = "assessed", stringsAsFactors = FALSE)
      population_loo_index <- population_loo_index + 1L
    }
  }
}
population_loo <- do.call(rbind, population_loo_rows)
write_tsv(population_loo, output_path("population_loo"))

library_loo_rows <- list()
library_loo_index <- 1L
local_subset <- family_subset(counts, metadata, "b_local_stage")
for (target in family_targets$b_local_stage) {
  target_libraries <- local_subset$metadata$library_id[local_subset$metadata$population == target]
  genes <- selected_gene_map[[paste("b_local_stage", target, sep = "|")]]
  if (is.null(genes)) genes <- character()
  if (length(target_libraries) < contract$stability$leave_one_library_out$minimum_target_libraries) {
    for (gene in genes) {
      library_loo_rows[[library_loo_index]] <- data.frame(family = "b_local_stage", target = target, gene_symbol = gene, holdouts = 0L, passing_holdouts = NA_integer_, holdout_fraction = NA_real_, status = "not_assessable_low_replication", stringsAsFactors = FALSE)
      library_loo_index <- library_loo_index + 1L
    }
    next
  }
  pass_counts <- setNames(integer(length(genes)), genes)
  for (held_out_library in target_libraries) {
    retained <- metadata$library_id != held_out_library
    passed <- fit_target_gate(counts[, retained, drop = FALSE], metadata[retained, , drop = FALSE], "b_local_stage", target)
    pass_counts <- pass_counts + as.integer(passed[genes])
  }
  for (gene in genes) {
    library_loo_rows[[library_loo_index]] <- data.frame(family = "b_local_stage", target = target, gene_symbol = gene, holdouts = length(target_libraries), passing_holdouts = pass_counts[[gene]], holdout_fraction = pass_counts[[gene]] / length(target_libraries), status = "assessed", stringsAsFactors = FALSE)
    library_loo_index <- library_loo_index + 1L
  }
}
library_loo <- do.call(rbind, library_loo_rows)
write_tsv(library_loo, output_path("library_loo"))

# Training-standardized broad scores for discovery LOPO and validation-only populations.
population_data <- aggregate_by_population(counts, metadata)
population_logcpm <- normalized_logcpm(population_data$counts)
population_meta <- population_data$metadata
broad_sets <- split(selected_sets$gene_symbol[selected_sets$family == "broad_lineage"], selected_sets$target[selected_sets$family == "broad_lineage"])
score_sets <- function(sample_index, training_indices) {
  vapply(names(broad_sets), function(set_name) {
    genes <- broad_sets[[set_name]]
    train <- population_logcpm[genes, training_indices, drop = FALSE]
    centers <- rowMeans(train)
    scales <- apply(train, 1L, sd)
    scales[!is.finite(scales) | scales == 0] <- 1
    mean((population_logcpm[genes, sample_index] - centers) / scales)
  }, numeric(1L))
}

lopo_scores <- list(); lopo_predictions <- list(); score_index <- 1L; prediction_index <- 1L
discovery_indices <- which(population_meta$role == "discovery")
for (sample_index in discovery_indices) {
  training <- setdiff(discovery_indices, sample_index)
  scores <- score_sets(sample_index, training)
  for (set_name in names(scores)) {
    lopo_scores[[score_index]] <- data.frame(population = population_meta$population[[sample_index]], expected_lineage = population_meta$table_s1_lineage[[sample_index]], gene_set = set_name, score = scores[[set_name]], stringsAsFactors = FALSE)
    score_index <- score_index + 1L
  }
  predicted <- names(scores)[which.max(scores)]
  lopo_predictions[[prediction_index]] <- data.frame(population = population_meta$population[[sample_index]], expected_lineage = population_meta$table_s1_lineage[[sample_index]], predicted_lineage = predicted, correct = predicted == population_meta$table_s1_lineage[[sample_index]], stringsAsFactors = FALSE)
  prediction_index <- prediction_index + 1L
}
write_tsv(do.call(rbind, lopo_scores), output_path("broad_lopo_scores"))
write_tsv(do.call(rbind, lopo_predictions), output_path("broad_lopo_predictions"))

heldout_scores <- list(); heldout_predictions <- list(); score_index <- 1L; prediction_index <- 1L
heldout_indices <- which(population_meta$role == "held_out_validation")
for (sample_index in heldout_indices) {
  scores <- score_sets(sample_index, discovery_indices)
  for (set_name in names(scores)) {
    heldout_scores[[score_index]] <- data.frame(population = population_meta$population[[sample_index]], expected_lineage = population_meta$table_s1_lineage[[sample_index]], gene_set = set_name, score = scores[[set_name]], stringsAsFactors = FALSE)
    score_index <- score_index + 1L
  }
  predicted <- names(scores)[which.max(scores)]
  heldout_predictions[[prediction_index]] <- data.frame(population = population_meta$population[[sample_index]], expected_lineage = population_meta$table_s1_lineage[[sample_index]], predicted_lineage = predicted, correct = predicted == population_meta$table_s1_lineage[[sample_index]], stringsAsFactors = FALSE)
  prediction_index <- prediction_index + 1L
}
write_tsv(do.call(rbind, heldout_scores), output_path("heldout_scores"))
write_tsv(do.call(rbind, heldout_predictions), output_path("heldout_predictions"))

auc_rank <- function(values, positive) {
  n_pos <- sum(positive); n_neg <- sum(!positive)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  ranks <- rank(values, ties.method = "average")
  (sum(ranks[positive]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

b_metadata <- metadata[metadata$role == "discovery" & metadata$table_s1_lineage == "B", , drop = FALSE]
b_logcpm <- normalized_logcpm(counts[, b_metadata$library_id, drop = FALSE])
local_rows <- list(); local_index <- 1L
local_plan <- plan[plan$family == "b_local_stage", , drop = FALSE]
for (target in unique(local_plan$target)) {
  genes <- selected_gene_map[[paste("b_local_stage", target, sep = "|")]]
  if (is.null(genes) || !length(genes)) next
  gene_values <- b_logcpm[genes, , drop = FALSE]
  centers <- rowMeans(gene_values); scales <- apply(gene_values, 1L, sd); scales[!is.finite(scales) | scales == 0] <- 1
  sample_scores <- colMeans((gene_values - centers) / scales)
  neighbors <- unique(unlist(lapply(local_plan$weights_json[local_plan$target == target], function(value) names(jsonlite::fromJSON(value))[as.numeric(jsonlite::fromJSON(value)) < 0])))
  for (neighbor in neighbors) {
    included <- b_metadata$population %in% c(target, neighbor)
    positive <- b_metadata$population[included] == target
    n_target <- sum(positive); n_neighbor <- sum(!positive)
    local_rows[[local_index]] <- data.frame(target = target, neighbor = neighbor, target_libraries = n_target, neighbor_libraries = n_neighbor,
      auroc = auc_rank(sample_scores[included], positive), status = if (n_target >= 2L && n_neighbor >= 2L) "assessed" else "insufficient_replication", stringsAsFactors = FALSE)
    local_index <- local_index + 1L
  }
}
local_auroc <- do.call(rbind, local_rows)
write_tsv(local_auroc, output_path("local_neighbor_auroc"))

summary <- list(
  schema_version = "bdev.immgen_deferred_gate_model_summary.v1",
  status = "model_outputs_complete_pending_validation",
  synthetic_only = TRUE,
  immgen_expression_values_read = FALSE,
  iterations = iterations,
  resample_schedule_rows = nrow(schedule),
  stability_rows = nrow(balanced_rows),
  population_loo_rows = nrow(population_loo),
  library_loo_rows = nrow(library_loo),
  broad_lopo_populations = nrow(do.call(rbind, lopo_predictions)),
  heldout_populations = nrow(do.call(rbind, heldout_predictions)),
  local_neighbor_tests = nrow(local_auroc),
  importable_gmt_created = FALSE,
  scworkbench_invoked = FALSE
)
report_path <- file.path(output_root, config$outputs$reports$model_summary)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
writeLines(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), report_path)
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
