#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(jsonlite)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) { i <- match(flag, args); if (is.na(i) || i == length(args)) stop("Missing argument: ", flag); args[[i + 1L]] }
config_path <- normalizePath(arg_value("--config"))
config <- yaml::read_yaml(config_path)
if (!identical(config$schema_version, "bdev.immgen_positive_refinement_nested_qualification.v1")) stop("Unexpected config schema")
config_dir <- dirname(config_path)
resolve_path <- function(value) if (startsWith(value, "/")) value else normalizePath(file.path(config_dir, value))
input_path <- function(name) resolve_path(config$inputs[[name]]$path)
output_root <- config$outputs$root
table_path <- function(name) file.path(output_root, config$outputs$tables[[name]])
report_path <- function(name) file.path(output_root, config$outputs$reports[[name]])

contract <- yaml::read_yaml(input_path("statistical_contract"))
counts_frame <- read.delim(input_path("selected_counts"), check.names = FALSE, stringsAsFactors = FALSE)
gene_symbols <- counts_frame[[1L]]
counts <- as.matrix(counts_frame[, -1L, drop = FALSE]); storage.mode(counts) <- "integer"; rownames(counts) <- gene_symbols
metadata <- read.delim(input_path("metadata"), check.names = FALSE, stringsAsFactors = FALSE)
plan <- read.delim(input_path("contrast_plan"), check.names = FALSE, stringsAsFactors = FALSE)
annotation <- read.delim(input_path("gene_annotation"), check.names = FALSE, stringsAsFactors = FALSE)
annotation <- annotation[match(gene_symbols, annotation$gene_symbol), , drop = FALSE]
annotation$core_identity_eligible <- tolower(as.character(annotation$core_identity_eligible)) == "true"
stopifnot(identical(colnames(counts), metadata$library_id), identical(gene_symbols, annotation$gene_symbol))

write_tsv <- function(value, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); write.table(value, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "") }
write_tsv_gz <- function(value, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); con <- gzfile(path, "wt", compression = 9); write.table(value, con, sep = "\t", quote = FALSE, row.names = FALSE, na = ""); close(con) }
aggregate_by_population <- function(input_counts, input_metadata) {
  populations <- unique(input_metadata$population)
  aggregated <- vapply(populations, function(population) rowSums(input_counts[, input_metadata$population == population, drop = FALSE]), numeric(nrow(input_counts)))
  rownames(aggregated) <- rownames(input_counts)
  population_metadata <- do.call(rbind, lapply(populations, function(population) input_metadata[match(population, input_metadata$population), , drop = FALSE]))
  population_metadata$library_id <- population_metadata$population
  list(counts = aggregated, metadata = population_metadata)
}
weights_for <- function(text, levels) {
  parsed <- jsonlite::fromJSON(text); result <- setNames(rep(0, length(levels)), levels)
  if (length(setdiff(names(parsed), levels))) stop("Contrast contains absent groups")
  result[names(parsed)] <- as.numeric(parsed); result
}
effect_threshold <- function(family, test_type) {
  x <- contract$contrast_families[[family]]
  if (family == "superlineage") return(x$minimum_global_log2_fold_change)
  if (family %in% c("broad_lineage", "b_branch")) return(if (test_type == "global") x$minimum_global_log2_fold_change else x$minimum_each_competitor_log2_fold_change)
  x$minimum_each_neighbor_log2_fold_change
}
fdr_threshold <- function(family, test_type) {
  x <- contract$contrast_families[[family]]
  if (family == "superlineage") return(x$maximum_global_fdr)
  if (family %in% c("broad_lineage", "b_branch")) return(if (test_type == "global") x$maximum_global_fdr else x$maximum_each_competitor_fdr)
  x$maximum_each_neighbor_fdr
}
slug <- function(x) toupper(gsub("(^_+|_+$)", "", gsub("[^A-Za-z0-9]+", "_", x)))

families <- c("superlineage", "broad_lineage", "b_branch", "b_local_stage")
group_column <- c(superlineage = "superlineage", broad_lineage = "table_s1_lineage", b_branch = "b_branch", b_local_stage = "population")
base_selection <- function(family) {
  selected <- metadata$role == "discovery"
  if (family %in% c("b_branch", "b_local_stage")) selected <- selected & metadata$table_s1_lineage == "B"
  if (family == "b_local_stage") {
    fp <- plan[plan$family == family, , drop = FALSE]
    local_populations <- unique(c(fp$target, unlist(lapply(fp$weights_json, function(v) names(jsonlite::fromJSON(v))))))
    selected <- selected & metadata$population %in% local_populations
  }
  selected
}

fold_plan_rows <- list(); fold_i <- 1L
for (family in families) {
  selected <- base_selection(family); md <- metadata[selected, , drop = FALSE]
  targets <- unique(plan$target[plan$family == family])
  for (target in targets) {
    units <- if (family == "b_local_stage") md$library_id[md$population == target] else unique(md$population[md[[group_column[[family]]]] == target])
    eligible <- length(units) >= config$qualification$minimum_independent_target_units
    if (eligible) {
      for (unit_index in seq_along(units)) {
        unit <- units[[unit_index]]
        fold_plan_rows[[fold_i]] <- data.frame(family = family, target = target, heldout_unit = unit, fold_id = paste(slug(family), slug(target), sprintf("FOLD%02d", unit_index), slug(unit), sep = "__"), target_units = length(units), qualification_eligible = TRUE, exclusion_reason = "", stringsAsFactors = FALSE); fold_i <- fold_i + 1L
      }
    } else {
      fold_plan_rows[[fold_i]] <- data.frame(family = family, target = target, heldout_unit = "", fold_id = "", target_units = length(units), qualification_eligible = FALSE, exclusion_reason = "fewer_than_three_independent_target_units", stringsAsFactors = FALSE); fold_i <- fold_i + 1L
    }
  }
}
fold_plan <- do.call(rbind, fold_plan_rows); write_tsv(fold_plan, table_path("fold_plan"))

pool_rows <- list(); pool_i <- 1L; audit_rows <- list(); audit_i <- 1L
eligible_folds <- fold_plan[fold_plan$qualification_eligible, , drop = FALSE]
for (fold_index in seq_len(nrow(eligible_folds))) {
  fold <- eligible_folds[fold_index, , drop = FALSE]; family <- fold$family; target <- fold$target; heldout <- fold$heldout_unit
  selected <- base_selection(family)
  if (family == "b_local_stage") selected <- selected & metadata$library_id != heldout else selected <- selected & metadata$population != heldout
  source_counts <- counts[, selected, drop = FALSE]; source_metadata <- metadata[selected, , drop = FALSE]
  library_y <- calcNormFactors(DGEList(counts = source_counts), method = "TMM")
  library_cpm <- cpm(library_y, normalized.lib.sizes = TRUE)
  family_counts <- source_counts; family_metadata <- source_metadata
  if (family != "b_local_stage") { aggregated <- aggregate_by_population(source_counts, source_metadata); family_counts <- aggregated$counts; family_metadata <- aggregated$metadata }
  groups <- factor(family_metadata[[group_column[[family]]]])
  design <- model.matrix(~ 0 + groups); colnames(design) <- levels(groups)
  if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) stop("Invalid fold design: ", fold$fold_id)
  y_all <- DGEList(counts = family_counts)
  keep <- filterByExpr(y_all, design = design, min.count = contract$count_processing$expression_filter$min_count, min.total.count = contract$count_processing$expression_filter$min_total_count)
  y <- y_all[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y, method = "TMM"); y_all$samples$norm.factors <- y$samples$norm.factors
  y <- estimateDisp(y, design = design); fit <- glmQLFit(y, design = design, robust = TRUE)
  normalized_cpm <- cpm(y_all, normalized.lib.sizes = TRUE)
  target_plan <- plan[plan$family == family & plan$target == target, , drop = FALSE]
  effect_pass <- rep(TRUE, nrow(y_all)); fdr_passes <- integer(nrow(y_all)); tests <- integer(nrow(y_all))
  worst_fdr_ratio <- numeric(nrow(y_all)); worst_effect_margin <- rep(Inf, nrow(y_all)); support <- rep(1, nrow(y_all)); expression <- rep(1, nrow(y_all))
  for (test_index in seq_len(nrow(target_plan))) {
    pr <- target_plan[test_index, , drop = FALSE]; contrast <- weights_for(pr$weights_json, levels(groups))
    test <- glmTreat(fit, contrast = contrast, lfc = contract$count_processing$model$treat_log2_fold_change)
    tab <- topTags(test, n = Inf, sort.by = "none")$table
    logfc <- rep(-Inf, nrow(y_all)); fdr <- rep(1, nrow(y_all)); idx <- match(rownames(y), rownames(tab)); logfc[keep] <- tab$logFC[idx]; fdr[keep] <- tab$FDR[idx]
    lfc_cut <- effect_threshold(family, pr$test_type); fdr_cut <- fdr_threshold(family, pr$test_type)
    effect_pass <- effect_pass & logfc >= lfc_cut; fdr_passes <- fdr_passes + as.integer(fdr <= fdr_cut); tests <- tests + 1L
    worst_fdr_ratio <- pmax(worst_fdr_ratio, fdr / fdr_cut); worst_effect_margin <- pmin(worst_effect_margin, logfc - lfc_cut)
    parsed <- jsonlite::fromJSON(pr$weights_json); competitors <- names(parsed)[as.numeric(parsed) < 0]
    target_samples <- groups == target; competitor_samples <- groups %in% competitors
    cutoff <- if (family == "b_local_stage" && sum(target_samples) == 1L) contract$contrast_families$b_local_stage$singleton_target_minimum_cpm else 1
    test_expression <- rowMeans(normalized_cpm[, target_samples, drop = FALSE] >= cutoff)
    competitor_median <- apply(normalized_cpm[, competitor_samples, drop = FALSE], 1L, median)
    direction <- rowMeans(normalized_cpm[, target_samples, drop = FALSE] > competitor_median)
    if (family != "b_local_stage") {
      target_populations <- family_metadata$population[target_samples]
      library_support <- vapply(target_populations, function(population) apply(library_cpm[, source_metadata$population == population, drop = FALSE], 1L, max) >= 1, logical(nrow(y_all)))
      if (is.vector(library_support)) library_support <- matrix(library_support, ncol = 1L)
      test_expression <- rowMeans((normalized_cpm[, target_samples, drop = FALSE] >= 1) & library_support)
    }
    expression <- pmin(expression, test_expression); support <- pmin(support, test_expression, direction)
  }
  threshold <- contract$contrast_families[[family]]
  support_ok <- if (family == "b_local_stage") expression >= threshold$minimum_target_library_fraction else support >= threshold$minimum_target_population_support
  eligible <- keep & annotation$core_identity_eligible & effect_pass & support_ok
  frame <- data.frame(gene_symbol = gene_symbols[eligible], required_tests = tests[eligible], fdr_passing_tests = fdr_passes[eligible], fdr_pass_fraction = fdr_passes[eligible] / tests[eligible], worst_fdr_ratio = worst_fdr_ratio[eligible], target_support_fraction = support[eligible], target_expression_fraction = expression[eligible], worst_effect_margin = worst_effect_margin[eligible], stringsAsFactors = FALSE)
  if (nrow(frame)) frame <- frame[order(-frame$fdr_pass_fraction, frame$worst_fdr_ratio, -frame$target_support_fraction, -frame$worst_effect_margin, frame$gene_symbol), , drop = FALSE]
  frame <- head(frame, config$qualification$candidate_pool_maximum_genes); frame$pool_rank <- seq_len(nrow(frame))
  if (nrow(frame)) {
    frame <- cbind(data.frame(family = family, target = target, heldout_unit = heldout, fold_id = fold$fold_id, stringsAsFactors = FALSE), frame)
    pool_rows[[pool_i]] <- frame; pool_i <- pool_i + 1L
  }
  audit_rows[[audit_i]] <- data.frame(family = family, target = target, heldout_unit = heldout, fold_id = fold$fold_id, training_profiles = ncol(family_counts), training_target_units = sum(groups == target), groups = nlevels(groups), residual_df = nrow(design) - qr(design)$rank, filter_pass_genes = sum(keep), positive_eligible_genes = sum(eligible), pool_genes = nrow(frame), stringsAsFactors = FALSE); audit_i <- audit_i + 1L
  cat("completed", fold_index, "of", nrow(eligible_folds), fold$fold_id, "pool", nrow(frame), "\n")
}
pool_table <- do.call(rbind, pool_rows); audit_table <- do.call(rbind, audit_rows)
write_tsv_gz(pool_table, table_path("fold_candidate_pools")); write_tsv(audit_table, table_path("fold_model_audit"))
summary <- list(schema_version = "bdev.immgen_positive_nested_model_summary.v1", status = "complete", eligible_targets = length(unique(paste(eligible_folds$family, eligible_folds$target))), excluded_targets = sum(!fold_plan$qualification_eligible), folds = nrow(eligible_folds), folds_below_family_minimum = sum(vapply(seq_len(nrow(audit_table)), function(i) audit_table$pool_genes[[i]] < contract$contrast_families[[audit_table$family[[i]]]]$minimum_genes, logical(1L))), positive_only = TRUE, scworkbench_invoked = FALSE)
dir.create(dirname(report_path("model_summary")), recursive = TRUE, showWarnings = FALSE); writeLines(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), report_path("model_summary")); cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
