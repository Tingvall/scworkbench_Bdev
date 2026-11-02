#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(jsonlite)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing argument: ", flag)
  args[[i + 1L]]
}
config_path <- normalizePath(arg_value("--config"))
config <- yaml::read_yaml(config_path)
if (!identical(config$schema_version, "bdev.immgen_replication_adaptive_parse_transfer.v1")) stop("Unexpected config schema")
config_dir <- dirname(config_path)
resolve_path <- function(value) if (startsWith(value, "/")) value else normalizePath(file.path(config_dir, value))
input_path <- function(name) resolve_path(config$inputs[[name]]$path)
output_root <- config$outputs$root
table_path <- function(name) file.path(output_root, config$outputs$tables[[name]])
report_path <- function(name) file.path(output_root, config$outputs$reports[[name]])

contract <- yaml::read_yaml(input_path("statistical_contract"))
counts_frame <- read.delim(input_path("selected_counts"), check.names = FALSE, stringsAsFactors = FALSE)
gene_symbols <- counts_frame[[1L]]
counts <- as.matrix(counts_frame[, -1L, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- gene_symbols
metadata <- read.delim(input_path("metadata"), check.names = FALSE, stringsAsFactors = FALSE)
plan <- read.delim(input_path("contrast_plan"), check.names = FALSE, stringsAsFactors = FALSE)
annotation <- read.delim(input_path("gene_annotation"), check.names = FALSE, stringsAsFactors = FALSE)
annotation <- annotation[match(gene_symbols, annotation$gene_symbol), , drop = FALSE]
annotation$core_identity_eligible <- tolower(as.character(annotation$core_identity_eligible)) == "true"
stopifnot(identical(colnames(counts), metadata$library_id), identical(gene_symbols, annotation$gene_symbol))

write_tsv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(value, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
write_tsv_gz <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, "wt", compression = 9)
  write.table(value, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  close(con)
}
aggregate_by_population <- function(input_counts, input_metadata) {
  populations <- unique(input_metadata$population)
  aggregated <- vapply(populations, function(population) {
    rowSums(input_counts[, input_metadata$population == population, drop = FALSE])
  }, numeric(nrow(input_counts)))
  rownames(aggregated) <- rownames(input_counts)
  population_metadata <- do.call(rbind, lapply(populations, function(population) {
    input_metadata[match(population, input_metadata$population), , drop = FALSE]
  }))
  population_metadata$library_id <- population_metadata$population
  list(counts = aggregated, metadata = population_metadata)
}

families <- c("superlineage", "broad_lineage", "b_branch", "b_local_stage")
group_column <- c(superlineage = "superlineage", broad_lineage = "table_s1_lineage", b_branch = "b_branch", b_local_stage = "population")
base_selection <- function(family) {
  selected <- metadata$role == "discovery"
  if (family %in% c("b_branch", "b_local_stage")) selected <- selected & metadata$table_s1_lineage == "B"
  if (family == "b_local_stage") selected <- selected & metadata$population %in% unique(plan$target[plan$family == family])
  selected
}
target_units <- function(family, md, target) {
  if (family == "b_local_stage") md$library_id[md$population == target] else unique(md$population[md[[group_column[[family]]]] == target])
}

plan_rows <- list()
plan_i <- 1L
for (family in families) {
  md <- metadata[base_selection(family), , drop = FALSE]
  for (target in unique(plan$target[plan$family == family])) {
    units <- target_units(family, md, target)
    n_units <- length(units)
    if (n_units >= 3L) {
      plan_rows[[plan_i]] <- data.frame(family = family, target = target, independent_units = n_units, evidence_tier = "high_threefold", fold = NA_integer_, heldout_unit = "", model_required = FALSE, stringsAsFactors = FALSE)
      plan_i <- plan_i + 1L
    } else if (n_units == 2L) {
      for (fold in 1:2) {
        plan_rows[[plan_i]] <- data.frame(family = family, target = target, independent_units = n_units, evidence_tier = "moderate_twofold", fold = fold, heldout_unit = units[[fold]], model_required = TRUE, stringsAsFactors = FALSE)
        plan_i <- plan_i + 1L
      }
    } else if (n_units == 1L) {
      plan_rows[[plan_i]] <- data.frame(family = family, target = target, independent_units = n_units, evidence_tier = "provisional_singleton", fold = 0L, heldout_unit = "", model_required = TRUE, stringsAsFactors = FALSE)
      plan_i <- plan_i + 1L
    } else stop("Target has no discovery unit: ", family, " / ", target)
  }
}
adaptive_plan <- do.call(rbind, plan_rows)
write_tsv(adaptive_plan, table_path("adaptive_plan"))

model_plan <- adaptive_plan[adaptive_plan$model_required, , drop = FALSE]
ranking_rows <- list()
audit_rows <- list()
ranking_i <- 1L
audit_i <- 1L
for (row_i in seq_len(nrow(model_plan))) {
  item <- model_plan[row_i, , drop = FALSE]
  family <- item$family
  target <- item$target
  selected <- base_selection(family)
  if (item$evidence_tier == "moderate_twofold") {
    if (family == "b_local_stage") selected <- selected & metadata$library_id != item$heldout_unit else selected <- selected & metadata$population != item$heldout_unit
  }
  source_counts <- counts[, selected, drop = FALSE]
  source_metadata <- metadata[selected, , drop = FALSE]
  if (family != "b_local_stage") {
    aggregated <- aggregate_by_population(source_counts, source_metadata)
    source_counts <- aggregated$counts
    source_metadata <- aggregated$metadata
  }
  groups <- factor(source_metadata[[group_column[[family]]]])
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  if (!target %in% levels(groups)) stop("Target absent: ", family, " / ", target, " fold ", item$fold)
  if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) stop("Invalid design: ", family, " / ", target, " fold ", item$fold)
  y_all <- DGEList(counts = source_counts)
  keep <- filterByExpr(y_all, design = design, min.count = contract$count_processing$expression_filter$min_count, min.total.count = contract$count_processing$expression_filter$min_total_count)
  y <- y_all[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design = design)
  fit <- glmQLFit(y, design = design, robust = TRUE)
  contrast <- setNames(rep(-1 / (nlevels(groups) - 1), nlevels(groups)), levels(groups))
  contrast[[target]] <- 1
  test <- glmQLFTest(fit, contrast = contrast)
  tab <- topTags(test, n = Inf, sort.by = "none")$table
  eligible <- tab$logFC > config$method$minimum_oriented_log2_fold_change
  eligible <- eligible & annotation$core_identity_eligible[match(rownames(tab), gene_symbols)]
  ranked <- data.frame(gene_symbol = rownames(tab)[eligible], logFC = tab$logFC[eligible], F = tab$F[eligible], FDR = tab$FDR[eligible], stringsAsFactors = FALSE)
  if (nrow(ranked)) {
    ranked$strict_significant <- ranked$FDR < config$method$strict_fdr
    ranked <- ranked[order(-as.integer(ranked$strict_significant), -ranked$F, -ranked$logFC, ranked$FDR, ranked$gene_symbol), , drop = FALSE]
    ranked <- head(ranked, config$method$fold_top_genes)
    ranked$fold_rank <- seq_len(nrow(ranked))
    ranked <- cbind(data.frame(family = family, target = target, evidence_tier = item$evidence_tier, fold = item$fold, heldout_unit = item$heldout_unit, stringsAsFactors = FALSE), ranked)
    ranking_rows[[ranking_i]] <- ranked
    ranking_i <- ranking_i + 1L
  }
  audit_rows[[audit_i]] <- data.frame(family = family, target = target, evidence_tier = item$evidence_tier, fold = item$fold, heldout_unit = item$heldout_unit, training_profiles = ncol(source_counts), training_target_units = sum(groups == target), groups = nlevels(groups), residual_df = nrow(design) - qr(design)$rank, filter_pass_genes = sum(keep), positive_effect_eligible_genes = sum(eligible), top_genes = min(sum(eligible), config$method$fold_top_genes), stringsAsFactors = FALSE)
  audit_i <- audit_i + 1L
  cat("completed", row_i, "of", nrow(model_plan), family, target, item$evidence_tier, "fold", item$fold, "\n")
}
rankings <- do.call(rbind, ranking_rows)
audit <- do.call(rbind, audit_rows)
write_tsv_gz(rankings, table_path("low_replication_rankings"))
write_tsv(audit, table_path("model_audit"))

candidate_rows <- list()
candidate_i <- 1L
low_targets <- unique(model_plan[, c("family", "target", "independent_units", "evidence_tier")])
for (i in seq_len(nrow(low_targets))) {
  item <- low_targets[i, , drop = FALSE]
  local <- rankings[rankings$family == item$family & rankings$target == item$target, , drop = FALSE]
  folds <- sort(unique(local$fold))
  fold_lists <- lapply(folds, function(fold) local$gene_symbol[local$fold == fold])
  candidate <- if (length(fold_lists) == 1L) fold_lists[[1L]] else Reduce(intersect, fold_lists)
  rank_matrix <- vapply(folds, function(fold) match(candidate, fold_lists[[which(folds == fold)]]), integer(length(candidate)))
  if (is.vector(rank_matrix)) rank_matrix <- matrix(rank_matrix, ncol = 1L)
  frame <- data.frame(family = item$family, target = item$target, independent_units = item$independent_units, evidence_tier = item$evidence_tier, gene_symbol = candidate, mean_fold_rank = rowMeans(rank_matrix), max_fold_rank = apply(rank_matrix, 1L, max), stringsAsFactors = FALSE)
  frame <- frame[order(frame$mean_fold_rank, frame$max_fold_rank, frame$gene_symbol), , drop = FALSE]
  frame$candidate_rank <- seq_len(nrow(frame))
  frame$construction_pass <- nrow(frame) >= config$method$minimum_candidate_genes
  candidate_rows[[candidate_i]] <- frame
  candidate_i <- candidate_i + 1L
}
candidates <- do.call(rbind, candidate_rows)
write_tsv_gz(candidates, table_path("low_replication_candidates"))

candidate_summary <- unique(candidates[, c("family", "target", "evidence_tier", "construction_pass")])
summary <- list(schema_version = "bdev.immgen_replication_adaptive_model_summary.v1", status = "complete", targets = length(unique(paste(adaptive_plan$family, adaptive_plan$target))), high_threefold_targets = length(unique(paste(adaptive_plan$family[adaptive_plan$evidence_tier == "high_threefold"], adaptive_plan$target[adaptive_plan$evidence_tier == "high_threefold"]))), moderate_twofold_targets = sum(candidate_summary$evidence_tier == "moderate_twofold"), provisional_singleton_targets = sum(candidate_summary$evidence_tier == "provisional_singleton"), low_replication_models = nrow(model_plan), low_replication_candidates_passing_minimum = sum(candidate_summary$construction_pass), single_cell_expression_read = FALSE, scworkbench_invoked = FALSE)
dir.create(dirname(report_path("model_summary")), recursive = TRUE, showWarnings = FALSE)
writeLines(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), report_path("model_summary"))
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
