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
if (!identical(config$schema_version, "bdev.immgen_vm25_marker_execution.v1")) stop("Unexpected execution schema")
config_dir <- dirname(normalizePath(config_path))
contract <- yaml::read_yaml(file.path(config_dir, config$inputs$statistical_contract$path))
table_path <- function(name) file.path(output_root, config$outputs$tables[[name]])

counts_frame <- read.delim(table_path("selected_counts"), check.names = FALSE, stringsAsFactors = FALSE)
gene_symbols <- counts_frame[[1L]]
counts <- as.matrix(counts_frame[, -1L, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- gene_symbols
metadata <- read.delim(table_path("metadata"), check.names = FALSE, stringsAsFactors = FALSE)
plan <- read.delim(table_path("contrast_plan"), check.names = FALSE, stringsAsFactors = FALSE)
annotation <- read.delim(config$inputs$gene_annotation$path, check.names = FALSE, stringsAsFactors = FALSE)
annotation <- annotation[match(gene_symbols, annotation$gene_symbol), , drop = FALSE]
stopifnot(identical(colnames(counts), metadata$library_id), identical(gene_symbols, annotation$gene_symbol), nrow(counts) == 55011L)

aggregate_by_population <- function(input_counts, input_metadata) {
  populations <- unique(input_metadata$population)
  aggregated <- vapply(populations, function(population) rowSums(input_counts[, input_metadata$population == population, drop = FALSE]), numeric(nrow(input_counts)))
  rownames(aggregated) <- rownames(input_counts)
  population_metadata <- do.call(rbind, lapply(populations, function(population) input_metadata[match(population, input_metadata$population), , drop = FALSE]))
  population_metadata$library_id <- population_metadata$population
  list(counts = aggregated, metadata = population_metadata)
}

auc_rank <- function(values, positive) {
  n_pos <- sum(positive); n_neg <- sum(!positive)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  ranks <- rank(values, ties.method = "average")
  (sum(ranks[positive]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

weights_for <- function(text, levels) {
  parsed <- jsonlite::fromJSON(text)
  result <- setNames(rep(0, length(levels)), levels)
  if (length(setdiff(names(parsed), levels))) stop("Contrast contains absent groups")
  result[names(parsed)] <- as.numeric(parsed)
  result
}

statistics_path <- table_path("complete_statistics")
dir.create(dirname(statistics_path), recursive = TRUE, showWarnings = FALSE)
statistics_connection <- gzfile(statistics_path, open = "wt", compression = 9)
header_written <- FALSE
design_rows <- list()
total_rows <- 0L

families <- c("superlineage", "broad_lineage", "b_branch", "b_local_stage")
for (family in families) {
  family_plan <- plan[plan$family == family, , drop = FALSE]
  scope <- unique(family_plan$sample_scope)
  group_column <- unique(family_plan$group_column)
  stopifnot(length(scope) == 1L, length(group_column) == 1L)
  selected <- metadata$role == "discovery"
  if (scope %in% c("b_discovery", "b_local_discovery")) selected <- selected & metadata$table_s1_lineage == "B"
  if (scope == "b_local_discovery") {
    local_populations <- unique(c(family_plan$target, unlist(lapply(family_plan$weights_json, function(value) names(jsonlite::fromJSON(value))))))
    selected <- selected & metadata$population %in% local_populations
  }
  source_counts <- counts[, selected, drop = FALSE]
  source_metadata <- metadata[selected, , drop = FALSE]
  library_y <- calcNormFactors(DGEList(counts = source_counts), method = "TMM")
  library_cpm <- cpm(library_y, normalized.lib.sizes = TRUE)
  family_counts <- source_counts
  family_metadata <- source_metadata
  if (family != "b_local_stage") {
    aggregated <- aggregate_by_population(source_counts, source_metadata)
    family_counts <- aggregated$counts
    family_metadata <- aggregated$metadata
  }
  groups <- factor(family_metadata[[group_column]])
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient design for ", family)
  y_all <- DGEList(counts = family_counts)
  keep <- filterByExpr(y_all, design = design, min.count = contract$count_processing$expression_filter$min_count,
                       min.total.count = contract$count_processing$expression_filter$min_total_count)
  y <- y_all[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y_all$samples$norm.factors <- y$samples$norm.factors
  y <- estimateDisp(y, design = design)
  fit <- glmQLFit(y, design = design, robust = TRUE)
  normalized_cpm <- cpm(y_all, normalized.lib.sizes = TRUE)
  design_rows[[family]] <- data.frame(family = family, unit = if (family == "b_local_stage") "rna_library" else "population_profile",
    samples = ncol(family_counts), groups = nlevels(groups), design_columns = ncol(design), design_rank = qr(design)$rank,
    residual_df = nrow(design) - qr(design)$rank, filtered_in = sum(keep), filtered_out = sum(!keep), stringsAsFactors = FALSE)

  for (test_index in seq_len(nrow(family_plan))) {
    plan_row <- family_plan[test_index, , drop = FALSE]
    contrast <- weights_for(plan_row$weights_json, levels(groups))
    test <- glmTreat(fit, contrast = contrast, lfc = contract$count_processing$model$treat_log2_fold_change)
    table <- topTags(test, n = Inf, sort.by = "none")$table
    table <- table[match(rownames(y), rownames(table)), , drop = FALSE]
    logFC <- logCPM <- PValue <- FDR <- rep(NA_real_, nrow(y_all))
    logFC[keep] <- table$logFC; logCPM[keep] <- table$logCPM; PValue[keep] <- table$PValue; FDR[keep] <- table$FDR
    parsed_weights <- jsonlite::fromJSON(plan_row$weights_json)
    competitor_groups <- names(parsed_weights)[as.numeric(parsed_weights) < 0]
    target_samples <- groups == plan_row$target
    competitor_samples <- groups %in% competitor_groups
    expression_cutoff <- if (family == "b_local_stage" && sum(target_samples) == 1L) contract$contrast_families$b_local_stage$singleton_target_minimum_cpm else 1
    target_expression_fraction <- rowMeans(normalized_cpm[, target_samples, drop = FALSE] >= expression_cutoff)
    competitor_median <- apply(normalized_cpm[, competitor_samples, drop = FALSE], 1L, median)
    target_direction_fraction <- rowMeans(normalized_cpm[, target_samples, drop = FALSE] > competitor_median)
    if (family != "b_local_stage") {
      target_populations <- family_metadata$population[target_samples]
      library_support <- vapply(target_populations, function(population) apply(library_cpm[, source_metadata$population == population, drop = FALSE], 1L, max) >= 1, logical(nrow(y_all)))
      if (is.vector(library_support)) library_support <- matrix(library_support, ncol = 1L)
      target_expression_fraction <- rowMeans((normalized_cpm[, target_samples, drop = FALSE] >= 1) & library_support)
    }
    auc <- vapply(seq_len(nrow(normalized_cpm)), function(gene_index) auc_rank(normalized_cpm[gene_index, target_samples | competitor_samples], target_samples[target_samples | competitor_samples]), numeric(1L))
    complete <- data.frame(
      family = family, target = plan_row$target, test_id = plan_row$test_id, test_type = plan_row$test_type,
      competitor_groups = paste(competitor_groups, collapse = ";"), gene_symbol = gene_symbols, filter_pass = keep,
      logFC = logFC, logCPM = logCPM, PValue = PValue, FDR = FDR,
      target_expression_fraction = target_expression_fraction, target_direction_fraction = target_direction_fraction,
      target_support_fraction = pmin(target_expression_fraction, target_direction_fraction),
      target_versus_competitor_auroc = auc,
      mapping_status = annotation$mapping_status, ensembl_gene_id = annotation$ensembl_gene_id,
      chromosome = annotation$chromosome, gene_biotype = annotation$gene_biotype,
      state_flag = annotation$state_flag, mandatory_identity_exclusion = annotation$mandatory_identity_exclusion,
      core_identity_eligible = annotation$core_identity_eligible, exclusion_reasons = annotation$exclusion_reasons,
      stringsAsFactors = FALSE
    )
    write.table(complete, statistics_connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "", col.names = !header_written)
    header_written <- TRUE
    total_rows <- total_rows + nrow(complete)
    rm(complete, table, test)
    gc(FALSE)
  }
}
close(statistics_connection)
design_audit <- do.call(rbind, design_rows)
write.table(design_audit, table_path("design_audit"), sep = "\t", quote = FALSE, row.names = FALSE)
summary <- list(schema_version = "bdev.immgen_vm25_model_summary.v1", status = "primary_models_complete",
  genes = nrow(counts), libraries = ncol(counts), contrast_tests = nrow(plan), complete_statistic_rows = total_rows,
  designs = design_audit, immgen_expression_values_read = TRUE, scworkbench_invoked = FALSE)
report_path <- file.path(output_root, config$outputs$reports$model_summary)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
writeLines(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), report_path)
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
