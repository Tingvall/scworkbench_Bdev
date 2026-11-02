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
if (!identical(config$schema_version, "bdev.immgen_all_population_parse.v1")) stop("Unexpected config schema")
resolve_path <- function(value) if (startsWith(value, "/")) value else normalizePath(file.path(dirname(config_path), value))
input_path <- function(name) resolve_path(config$inputs[[name]]$path)
intermediate_root <- config$outputs$intermediate_root
result_root <- config$outputs$results_root
intermediate_table <- function(name) file.path(intermediate_root, config$outputs$intermediate_tables[[name]])
report_path <- function(name) file.path(result_root, config$outputs$reports[[name]])

contract <- yaml::read_yaml(input_path("statistical_contract"))
counts_frame <- read.delim(input_path("selected_counts"), check.names = FALSE, stringsAsFactors = FALSE)
gene_symbols <- counts_frame[[1L]]
counts <- as.matrix(counts_frame[, -1L, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- gene_symbols
metadata <- read.delim(input_path("metadata"), check.names = FALSE, stringsAsFactors = FALSE)
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

population_order <- unique(metadata$population)
population_md <- metadata[match(population_order, metadata$population), , drop = FALSE]
if (length(population_order) != config$target_scope$expected_populations) stop("Unexpected population count")

evaluation_selection <- function(target) {
  target_row <- population_md[population_md$population == target, , drop = FALSE]
  lineage <- target_row$table_s1_lineage[[1L]]
  target_role <- target_row$role[[1L]]
  selected <- metadata$table_s1_lineage == lineage & metadata$role == config$target_scope$discovery_role
  if (target_role == config$target_scope$held_out_role) selected <- selected | metadata$population == target
  selected
}

plan_rows <- list()
plan_i <- 1L
for (target in population_order) {
  selected <- evaluation_selection(target)
  target_rows <- metadata[selected & metadata$population == target, , drop = FALSE]
  units <- target_rows$library_id
  n_units <- length(units)
  if (n_units >= 3L) {
    tier <- "replicated_threefold"
    n_folds <- config$method$crossfit_folds_for_three_or_more_units
    assignment <- ((seq_along(units) - 1L) %% n_folds) + 1L
    for (fold in seq_len(n_folds)) {
      heldout <- units[assignment == fold]
      plan_rows[[plan_i]] <- data.frame(
        target = target, source_role = target_rows$role[[1L]], lineage = target_rows$table_s1_lineage[[1L]],
        independent_units = n_units, evidence_tier = tier, fold = fold,
        heldout_units = paste(heldout, collapse = ";"), stringsAsFactors = FALSE
      )
      plan_i <- plan_i + 1L
    }
  } else if (n_units == 2L) {
    for (fold in 1:2) {
      plan_rows[[plan_i]] <- data.frame(
        target = target, source_role = target_rows$role[[1L]], lineage = target_rows$table_s1_lineage[[1L]],
        independent_units = n_units, evidence_tier = "moderate_twofold", fold = fold,
        heldout_units = units[[fold]], stringsAsFactors = FALSE
      )
      plan_i <- plan_i + 1L
    }
  } else if (n_units == 1L) {
    plan_rows[[plan_i]] <- data.frame(
      target = target, source_role = target_rows$role[[1L]], lineage = target_rows$table_s1_lineage[[1L]],
      independent_units = n_units, evidence_tier = "provisional_singleton", fold = 0L,
      heldout_units = "", stringsAsFactors = FALSE
    )
    plan_i <- plan_i + 1L
  } else stop("Target has no units: ", target)
}
fold_plan <- do.call(rbind, plan_rows)
write_tsv(fold_plan, intermediate_table("fold_plan"))

fit_one <- function(row_i) {
  item <- fold_plan[row_i, , drop = FALSE]
  target <- item$target[[1L]]
  selected <- evaluation_selection(target)
  heldout <- if (nzchar(item$heldout_units[[1L]])) strsplit(item$heldout_units[[1L]], ";", fixed = TRUE)[[1L]] else character()
  if (length(heldout)) selected <- selected & !metadata$library_id %in% heldout
  source_counts <- counts[, selected, drop = FALSE]
  source_metadata <- metadata[selected, , drop = FALSE]
  groups <- factor(source_metadata$population)
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  if (!target %in% levels(groups)) stop("Target absent after holdout: ", target, " fold ", item$fold)
  if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) stop("Invalid design: ", target, " fold ", item$fold)

  y_all <- DGEList(counts = source_counts)
  keep <- filterByExpr(
    y_all, design = design,
    min.count = contract$count_processing$expression_filter$min_count,
    min.total.count = contract$count_processing$expression_filter$min_total_count
  )
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
  ranked <- data.frame(
    gene_symbol = rownames(tab)[eligible], logFC = tab$logFC[eligible],
    F = tab$F[eligible], FDR = tab$FDR[eligible], stringsAsFactors = FALSE
  )
  if (nrow(ranked)) {
    ranked$strict_significant <- ranked$FDR < config$method$strict_fdr
    ranked <- ranked[order(-as.integer(ranked$strict_significant), -ranked$F, -ranked$logFC, ranked$FDR, ranked$gene_symbol), , drop = FALSE]
    ranked <- head(ranked, config$method$fold_top_genes)
    ranked$fold_rank <- seq_len(nrow(ranked))
    ranked <- cbind(
      data.frame(
        target = target, source_role = item$source_role, lineage = item$lineage,
        independent_units = item$independent_units, evidence_tier = item$evidence_tier,
        fold = item$fold, heldout_units = item$heldout_units, stringsAsFactors = FALSE
      ),
      ranked
    )
  } else {
    ranked <- NULL
  }
  audit <- data.frame(
    target = target, source_role = item$source_role, lineage = item$lineage,
    independent_units = item$independent_units, evidence_tier = item$evidence_tier,
    fold = item$fold, heldout_units = item$heldout_units,
    training_libraries = ncol(source_counts), training_target_libraries = sum(groups == target),
    populations = nlevels(groups), residual_df = nrow(design) - qr(design)$rank,
    filter_pass_genes = sum(keep), positive_effect_eligible_genes = sum(eligible),
    top_genes = min(sum(eligible), config$method$fold_top_genes), stringsAsFactors = FALSE
  )
  cat("completed", row_i, "of", nrow(fold_plan), target, item$evidence_tier, "fold", item$fold, "\n")
  list(ranking = ranked, audit = audit)
}

threads <- as.integer(config$runtime$threads)
fit_results <- parallel::mclapply(seq_len(nrow(fold_plan)), fit_one, mc.cores = threads, mc.preschedule = FALSE)
ranking_rows <- lapply(fit_results, `[[`, "ranking")
ranking_rows <- ranking_rows[!vapply(ranking_rows, is.null, logical(1L))]
if (!length(ranking_rows)) stop("No fold rankings were generated")
rankings <- do.call(rbind, ranking_rows)
audit <- do.call(rbind, lapply(fit_results, `[[`, "audit"))
write_tsv_gz(rankings, intermediate_table("fold_rankings"))
write_tsv(audit, intermediate_table("model_audit"))

summary <- list(
  schema_version = "bdev.immgen_all_population_parse_model_summary.v1",
  status = "complete",
  targets = length(population_order),
  discovery_targets = sum(population_md$role == config$target_scope$discovery_role),
  held_out_derived_targets = sum(population_md$role == config$target_scope$held_out_role),
  model_fits = nrow(fold_plan),
  evidence_tier_counts = as.list(table(unique(fold_plan[, c("target", "evidence_tier")])$evidence_tier)),
  single_cell_expression_read = FALSE,
  scworkbench_invoked = FALSE
)
dir.create(dirname(report_path("model_summary")), recursive = TRUE, showWarnings = FALSE)
writeLines(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), report_path("model_summary"))
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
