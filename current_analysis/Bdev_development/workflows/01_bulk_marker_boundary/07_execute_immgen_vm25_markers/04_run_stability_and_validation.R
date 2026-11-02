#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(jsonlite)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) { index <- match(flag, args); if (is.na(index) || index == length(args)) stop("Missing argument: ", flag); args[[index + 1L]] }
config_path <- arg_value("--config")
output_root <- arg_value("--output-root")
config <- yaml::read_yaml(config_path)
if (!identical(config$schema_version, "bdev.immgen_vm25_marker_execution.v1")) stop("Unexpected execution schema")
config_dir <- dirname(normalizePath(config_path))
contract <- yaml::read_yaml(file.path(config_dir, config$inputs$statistical_contract$path))
table_path <- function(name) file.path(output_root, config$outputs$tables[[name]])

counts_frame <- read.delim(table_path("selected_counts"), check.names = FALSE, stringsAsFactors = FALSE)
gene_symbols <- counts_frame[[1L]]
counts <- as.matrix(counts_frame[, -1L, drop = FALSE]); storage.mode(counts) <- "integer"; rownames(counts) <- gene_symbols
metadata <- read.delim(table_path("metadata"), check.names = FALSE, stringsAsFactors = FALSE)
plan <- read.delim(table_path("contrast_plan"), check.names = FALSE, stringsAsFactors = FALSE)
candidates <- read.delim(table_path("preliminary_candidates"), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(colnames(counts), metadata$library_id))

families <- c("superlineage", "broad_lineage", "b_branch", "b_local_stage")
group_column <- c(superlineage = "superlineage", broad_lineage = "table_s1_lineage", b_branch = "b_branch", b_local_stage = "population")
write_tsv <- function(value, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); write.table(value, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "") }
write_tsv_gz <- function(value, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); connection <- gzfile(path, "wt", compression = 9); write.table(value, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = ""); close(connection) }

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
    family_plan <- plan[plan$family == family, , drop = FALSE]
    local_populations <- unique(c(family_plan$target, unlist(lapply(family_plan$weights_json, function(value) names(jsonlite::fromJSON(value))))))
    selected <- selected & input_metadata$population %in% local_populations
  }
  list(counts = input_counts[, selected, drop = FALSE], metadata = input_metadata[selected, , drop = FALSE])
}

effect_threshold <- function(family, test_type) {
  threshold <- contract$contrast_families[[family]]
  if (family == "superlineage") return(threshold$minimum_global_log2_fold_change)
  if (family == "broad_lineage") return(if (test_type == "global") threshold$minimum_global_log2_fold_change else threshold$minimum_each_competitor_log2_fold_change)
  if (family == "b_branch") return(if (test_type == "global") threshold$minimum_global_log2_fold_change else threshold$minimum_each_competitor_log2_fold_change)
  threshold$minimum_each_neighbor_log2_fold_change
}
fdr_threshold <- function(family, test_type) {
  threshold <- contract$contrast_families[[family]]
  if (family == "superlineage") return(threshold$maximum_global_fdr)
  if (family %in% c("broad_lineage", "b_branch")) return(if (test_type == "global") threshold$maximum_global_fdr else threshold$maximum_each_competitor_fdr)
  threshold$maximum_each_neighbor_fdr
}

candidate_map <- split(candidates$gene_symbol, paste(candidates$family, candidates$target, sep = "|"))
candidate_indices <- lapply(candidate_map, function(genes) match(genes, gene_symbols))

# Seeded cycling produces one library per population per iteration with maximum imbalance one.
discovery_metadata <- metadata[metadata$role == "discovery", , drop = FALSE]
populations <- unique(discovery_metadata$population)
iterations <- config$execution$balanced_resampling_iterations
set.seed(config$execution$seed)
permutations <- lapply(populations, function(population) sample(sort(discovery_metadata$library_id[discovery_metadata$population == population])))
names(permutations) <- populations
schedule <- do.call(rbind, lapply(seq_len(iterations), function(iteration) do.call(rbind, lapply(populations, function(population) {
  choices <- permutations[[population]]
  data.frame(iteration = iteration, population = population, library_id = choices[((iteration - 1L) %% length(choices)) + 1L], stringsAsFactors = FALSE)
}))))
write_tsv(schedule, table_path("resample_schedule"))

stability_counts <- lapply(candidate_map, function(genes) integer(length(genes)))
for (iteration in seq_len(iterations)) {
  selected_ids <- schedule$library_id[schedule$iteration == iteration]
  iteration_metadata <- metadata[match(selected_ids, metadata$library_id), , drop = FALSE]
  iteration_counts <- counts[, selected_ids, drop = FALSE]
  for (family in families) {
    subset <- family_subset(iteration_counts, iteration_metadata, family)
    y <- calcNormFactors(DGEList(counts = subset$counts), method = "TMM")
    family_keys <- names(candidate_map)[startsWith(names(candidate_map), paste0(family, "|"))]
    family_genes <- unique(unlist(candidate_map[family_keys], use.names = FALSE))
    family_index <- match(family_genes, rownames(subset$counts))
    logcpm <- cpm(y[family_index, , keep.lib.sizes = TRUE], log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)
    rownames(logcpm) <- family_genes
    groups <- subset$metadata[[group_column[[family]]]]
    means <- vapply(unique(groups), function(group) rowMeans(logcpm[, groups == group, drop = FALSE]), numeric(length(family_genes)))
    rownames(means) <- family_genes
    family_plan <- plan[plan$family == family, , drop = FALSE]
    for (key in family_keys) {
      target <- sub("^[^|]+\\|", "", key)
      target_plan <- family_plan[family_plan$target == target, , drop = FALSE]
      genes <- candidate_map[[key]]
      passed <- rep(TRUE, length(genes))
      for (test_index in seq_len(nrow(target_plan))) {
        weights <- jsonlite::fromJSON(target_plan$weights_json[[test_index]])
        effects <- as.numeric(means[genes, names(weights), drop = FALSE] %*% as.numeric(weights))
        passed <- passed & effects >= effect_threshold(family, target_plan$test_type[[test_index]])
      }
      stability_counts[[key]] <- stability_counts[[key]] + as.integer(passed)
    }
  }
  if (iteration %% 25L == 0L) cat("completed balanced iteration", iteration, "of", iterations, "\n")
}
balanced <- do.call(rbind, lapply(names(candidate_map), function(key) {
  parts <- strsplit(key, "|", fixed = TRUE)[[1L]]
  data.frame(family = parts[[1L]], target = parts[[2L]], gene_symbol = candidate_map[[key]], passing_iterations = stability_counts[[key]], iterations = iterations,
             stability_fraction = stability_counts[[key]] / iterations, stringsAsFactors = FALSE)
}))
write_tsv_gz(balanced, table_path("balanced_stability"))

fit_target_gate <- function(input_counts, input_metadata, family, target, genes) {
  subset <- family_subset(input_counts, input_metadata, family)
  family_counts <- subset$counts; family_metadata <- subset$metadata
  source_counts <- family_counts; source_metadata <- family_metadata
  source_y <- calcNormFactors(DGEList(counts = source_counts), method = "TMM")
  source_gene_index <- match(genes, rownames(source_counts))
  source_cpm <- cpm(source_y[source_gene_index, , keep.lib.sizes = TRUE], normalized.lib.sizes = TRUE)
  rownames(source_cpm) <- genes
  if (family != "b_local_stage") { aggregated <- aggregate_by_population(family_counts, family_metadata); family_counts <- aggregated$counts; family_metadata <- aggregated$metadata }
  groups <- factor(family_metadata[[group_column[[family]]]])
  design <- model.matrix(~ 0 + groups); colnames(design) <- levels(groups)
  if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) stop("Invalid holdout design: ", family, "/", target)
  y_all <- DGEList(counts = family_counts)
  keep <- filterByExpr(y_all, design = design, min.count = contract$count_processing$expression_filter$min_count, min.total.count = contract$count_processing$expression_filter$min_total_count)
  y <- y_all[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y, method = "TMM"); y_all$samples$norm.factors <- y$samples$norm.factors
  y <- estimateDisp(y, design = design); fit <- glmQLFit(y, design = design, robust = TRUE)
  normalized_cpm <- cpm(y_all, normalized.lib.sizes = TRUE)
  gene_index <- match(genes, rownames(y_all)); passed <- keep[gene_index]; support <- rep(1, length(genes))
  target_plan <- plan[plan$family == family & plan$target == target, , drop = FALSE]
  for (test_index in seq_len(nrow(target_plan))) {
    weights <- jsonlite::fromJSON(target_plan$weights_json[[test_index]])
    contrast <- setNames(rep(0, ncol(design)), colnames(design)); contrast[names(weights)] <- as.numeric(weights)
    test <- glmTreat(fit, contrast = contrast, lfc = contract$count_processing$model$treat_log2_fold_change)
    table <- topTags(test, n = Inf, sort.by = "none")$table
    tested_index <- match(genes, rownames(table)); logfc <- table$logFC[tested_index]; fdr <- table$FDR[tested_index]
    passed <- passed & !is.na(logfc) & logfc >= effect_threshold(family, target_plan$test_type[[test_index]]) & !is.na(fdr) & fdr <= fdr_threshold(family, target_plan$test_type[[test_index]])
    competitor_groups <- names(weights)[as.numeric(weights) < 0]; target_samples <- groups == target; competitor_samples <- groups %in% competitor_groups
    cutoff <- if (family == "b_local_stage" && sum(target_samples) == 1L) contract$contrast_families$b_local_stage$singleton_target_minimum_cpm else 1
    expression <- rowMeans(normalized_cpm[gene_index, target_samples, drop = FALSE] >= cutoff)
    if (family != "b_local_stage") {
      target_populations <- family_metadata$population[target_samples]
      library_support <- vapply(target_populations, function(population) apply(source_cpm[, source_metadata$population == population, drop = FALSE], 1L, max) >= 1, logical(length(genes)))
      if (is.vector(library_support)) library_support <- matrix(library_support, ncol = 1L)
      expression <- rowMeans((normalized_cpm[gene_index, target_samples, drop = FALSE] >= 1) & library_support)
    }
    competitor_median <- apply(normalized_cpm[gene_index, competitor_samples, drop = FALSE], 1L, median)
    direction <- rowMeans(normalized_cpm[gene_index, target_samples, drop = FALSE] > competitor_median)
    support <- pmin(support, expression, direction)
  }
  threshold <- contract$contrast_families[[family]]
  if (family == "superlineage") passed <- passed & support >= threshold$minimum_target_population_support
  if (family %in% c("broad_lineage", "b_branch")) passed <- passed & support >= threshold$minimum_target_population_support
  if (family == "b_local_stage") passed <- passed & support >= threshold$minimum_target_library_fraction
  setNames(passed, genes)
}

population_rows <- list(); population_index <- 1L
for (family in c("superlineage", "broad_lineage", "b_branch")) {
  subset <- family_subset(counts, metadata, family); column <- group_column[[family]]
  keys <- names(candidate_map)[startsWith(names(candidate_map), paste0(family, "|"))]
  for (key in keys) {
    target <- sub("^[^|]+\\|", "", key); genes <- candidate_map[[key]]
    target_populations <- unique(subset$metadata$population[subset$metadata[[column]] == target])
    if (length(target_populations) < contract$stability$leave_one_population_out$minimum_target_populations) {
      population_rows[[population_index]] <- data.frame(family = family, target = target, gene_symbol = genes, holdouts = 0L, passing_holdouts = NA_integer_, holdout_fraction = NA_real_, status = "not_assessable_low_replication", stringsAsFactors = FALSE)
      population_index <- population_index + 1L; next
    }
    pass_counts <- setNames(integer(length(genes)), genes)
    for (held_out_population in target_populations) {
      retained <- metadata$population != held_out_population
      pass_counts <- pass_counts + as.integer(fit_target_gate(counts[, retained, drop = FALSE], metadata[retained, , drop = FALSE], family, target, genes))
    }
    population_rows[[population_index]] <- data.frame(family = family, target = target, gene_symbol = genes, holdouts = length(target_populations), passing_holdouts = pass_counts,
      holdout_fraction = pass_counts / length(target_populations), status = "assessed", stringsAsFactors = FALSE)
    population_index <- population_index + 1L
    cat("completed population LOO", family, target, "\n")
  }
}
population_loo <- do.call(rbind, population_rows); write_tsv_gz(population_loo, table_path("population_loo"))

library_rows <- list(); library_index <- 1L; local_subset <- family_subset(counts, metadata, "b_local_stage")
for (key in names(candidate_map)[startsWith(names(candidate_map), "b_local_stage|")]) {
  target <- sub("^[^|]+\\|", "", key); genes <- candidate_map[[key]]
  target_libraries <- local_subset$metadata$library_id[local_subset$metadata$population == target]
  if (length(target_libraries) < contract$stability$leave_one_library_out$minimum_target_libraries) {
    library_rows[[library_index]] <- data.frame(family = "b_local_stage", target = target, gene_symbol = genes, holdouts = 0L, passing_holdouts = NA_integer_, holdout_fraction = NA_real_, status = "not_assessable_low_replication", stringsAsFactors = FALSE)
    library_index <- library_index + 1L; next
  }
  pass_counts <- setNames(integer(length(genes)), genes)
  for (held_out_library in target_libraries) {
    retained <- metadata$library_id != held_out_library
    pass_counts <- pass_counts + as.integer(fit_target_gate(counts[, retained, drop = FALSE], metadata[retained, , drop = FALSE], "b_local_stage", target, genes))
  }
  library_rows[[library_index]] <- data.frame(family = "b_local_stage", target = target, gene_symbol = genes, holdouts = length(target_libraries), passing_holdouts = pass_counts,
    holdout_fraction = pass_counts / length(target_libraries), status = "assessed", stringsAsFactors = FALSE)
  library_index <- library_index + 1L
}
library_loo <- if (length(library_rows)) do.call(rbind, library_rows) else data.frame(family=character(),target=character(),gene_symbol=character(),holdouts=integer(),passing_holdouts=integer(),holdout_fraction=numeric(),status=character())
write_tsv(library_loo, table_path("library_loo"))

# Freeze validation sets after stability/LOO gates, retaining low-replication candidates but never accepting them.
balanced_key <- setNames(balanced$stability_fraction, paste(balanced$family, balanced$target, balanced$gene_symbol, sep = "|"))
pop_key <- setNames(population_loo$holdout_fraction, paste(population_loo$family, population_loo$target, population_loo$gene_symbol, sep = "|"))
pop_status <- setNames(population_loo$status, paste(population_loo$family, population_loo$target, population_loo$gene_symbol, sep = "|"))
lib_key <- setNames(library_loo$holdout_fraction, paste(library_loo$family, library_loo$target, library_loo$gene_symbol, sep = "|"))
lib_status <- setNames(library_loo$status, paste(library_loo$family, library_loo$target, library_loo$gene_symbol, sep = "|"))
validation_rows <- list(); validation_index <- 1L
for (key in names(candidate_map)) {
  parts <- strsplit(key, "|", fixed = TRUE)[[1L]]; family <- parts[[1L]]; target <- parts[[2L]]
  target_candidates <- candidates[candidates$family == family & candidates$target == target, , drop = FALSE]
  row_keys <- paste(family, target, target_candidates$gene_symbol, sep = "|")
  stable <- balanced_key[row_keys] >= contract$contrast_families[[family]]$minimum_stability_fraction
  if (family == "b_local_stage") {
    holdout_ok <- lib_status[row_keys] == "not_assessable_low_replication" | (!is.na(lib_key[row_keys]) & lib_key[row_keys] >= 0.80)
  } else {
    holdout_ok <- pop_status[row_keys] == "not_assessable_low_replication" | (!is.na(pop_key[row_keys]) & pop_key[row_keys] >= 0.80)
  }
  eligible <- target_candidates[stable & holdout_ok, , drop = FALSE]
  threshold <- contract$contrast_families[[family]]
  if (nrow(eligible) >= threshold$minimum_genes) eligible <- eligible[seq_len(min(nrow(eligible), threshold$maximum_genes)), , drop = FALSE] else eligible <- eligible[0, , drop = FALSE]
  if (nrow(eligible)) {
    eligible$validation_rank <- seq_len(nrow(eligible)); eligible$stability_fraction <- balanced_key[paste(family, target, eligible$gene_symbol, sep = "|")]
    validation_rows[[validation_index]] <- eligible; validation_index <- validation_index + 1L
  }
}
validation_sets <- if (length(validation_rows)) do.call(rbind, validation_rows) else candidates[0, , drop = FALSE]
write_tsv(validation_sets, table_path("validation_gene_sets"))

# Training-standardized broad population scores.
population_data <- aggregate_by_population(counts, metadata); population_y <- calcNormFactors(DGEList(counts = population_data$counts), method = "TMM")
population_logcpm <- cpm(population_y, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE); population_meta <- population_data$metadata
broad_sets <- split(validation_sets$gene_symbol[validation_sets$family == "broad_lineage"], validation_sets$target[validation_sets$family == "broad_lineage"])
score_sets <- function(sample_index, training_indices) vapply(names(broad_sets), function(set_name) {
  genes <- broad_sets[[set_name]]; train <- population_logcpm[genes, training_indices, drop = FALSE]; centers <- rowMeans(train); scales <- apply(train, 1L, sd); scales[!is.finite(scales) | scales == 0] <- 1
  mean((population_logcpm[genes, sample_index] - centers) / scales)
}, numeric(1L))
score_population_group <- function(indices, training_function) {
  score_rows <- list(); prediction_rows <- list(); si <- 1L
  for (sample_index in indices) {
    scores <- score_sets(sample_index, training_function(sample_index)); expected <- population_meta$table_s1_lineage[[sample_index]]
    score_rows[[si]] <- data.frame(population = population_meta$population[[sample_index]], expected_lineage = expected, gene_set = names(scores), score = as.numeric(scores), stringsAsFactors = FALSE)
    prediction_rows[[si]] <- data.frame(population = population_meta$population[[sample_index]], expected_lineage = expected, predicted_lineage = if (length(scores)) names(scores)[which.max(scores)] else "", correct = length(scores) && names(scores)[which.max(scores)] == expected, stringsAsFactors = FALSE); si <- si + 1L
  }
  list(scores = do.call(rbind, score_rows), predictions = do.call(rbind, prediction_rows))
}
discovery_indices <- which(population_meta$role == "discovery")
lopo <- score_population_group(discovery_indices, function(index) setdiff(discovery_indices, index)); write_tsv(lopo$scores, table_path("broad_lopo_scores")); write_tsv(lopo$predictions, table_path("broad_lopo_predictions"))
heldout_indices <- which(population_meta$role == "held_out_validation")
heldout <- score_population_group(heldout_indices, function(index) discovery_indices); write_tsv(heldout$scores, table_path("heldout_scores")); write_tsv(heldout$predictions, table_path("heldout_predictions"))

auc_rank <- function(values, positive) { n_pos <- sum(positive); n_neg <- sum(!positive); if (!n_pos || !n_neg) return(NA_real_); ranks <- rank(values, ties.method = "average"); (sum(ranks[positive]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg) }
b_metadata <- metadata[metadata$role == "discovery" & metadata$table_s1_lineage == "B", , drop = FALSE]
b_y <- calcNormFactors(DGEList(counts = counts[, b_metadata$library_id, drop = FALSE]), method = "TMM")
b_logcpm <- cpm(b_y, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)
local_rows <- list(); local_index <- 1L; local_plan <- plan[plan$family == "b_local_stage", , drop = FALSE]
for (target in unique(local_plan$target)) {
  genes <- validation_sets$gene_symbol[validation_sets$family == "b_local_stage" & validation_sets$target == target]
  if (!length(genes)) next
  values <- b_logcpm[genes, , drop = FALSE]; centers <- rowMeans(values); scales <- apply(values, 1L, sd); scales[!is.finite(scales) | scales == 0] <- 1; scores <- colMeans((values - centers) / scales)
  target_plan <- local_plan[local_plan$target == target, , drop = FALSE]
  neighbors <- unique(unlist(lapply(target_plan$weights_json, function(value) { weights <- jsonlite::fromJSON(value); names(weights)[as.numeric(weights) < 0] })))
  for (neighbor in neighbors) {
    included <- b_metadata$population %in% c(target, neighbor); positive <- b_metadata$population[included] == target
    local_rows[[local_index]] <- data.frame(target = target, neighbor = neighbor, target_libraries = sum(positive), neighbor_libraries = sum(!positive), auroc = auc_rank(scores[included], positive),
      status = if (sum(positive) >= 2L && sum(!positive) >= 2L) "assessed" else "insufficient_replication", stringsAsFactors = FALSE); local_index <- local_index + 1L
  }
}
local_auroc <- if (length(local_rows)) do.call(rbind, local_rows) else data.frame(target=character(),neighbor=character(),target_libraries=integer(),neighbor_libraries=integer(),auroc=numeric(),status=character())
write_tsv(local_auroc, table_path("local_neighbor_auroc"))
cat(jsonlite::toJSON(list(status = "stability_and_validation_complete", preliminary_candidates = nrow(candidates), validation_gene_rows = nrow(validation_sets), iterations = iterations), auto_unbox = TRUE, pretty = TRUE), "\n")
