#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 54_reviewed_cytokine_permutation_validation.R",
  "<bundle_dir> <inventory.tsv> <membership.tsv> <output_dir> <n_permutations>"))

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
inventory_path <- normalizePath(args[[2L]], mustWork = TRUE)
membership_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
n_permutations <- as.integer(args[[5L]])
stopifnot(is.finite(n_permutations), n_permutations > 0L)
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
})

read_tsv <- function(path) read.delim(path, check.names = FALSE,
  na.strings = "NA", stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt"); on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ..., "\n")
  flush.console()
}

inventory <- read_tsv(inventory_path)
members <- read_tsv(membership_path)
required_inventory <- c("signature_id", "cytokine", "direction", "final_size",
  "final_overall_pass")
required_members <- c("signature_id", "cytokine", "direction", "final_rank",
  "feature_id")
stopifnot(all(required_inventory %in% names(inventory)),
  all(required_members %in% names(members)), nrow(inventory) == 38L,
  !anyDuplicated(inventory$signature_id),
  setequal(unique(members$signature_id), inventory$signature_id),
  all(tolower(as.character(inventory$final_overall_pass)) %in%
    c("true", "t", "1")))
member_counts <- table(members$signature_id)
stopifnot(all(member_counts[inventory$signature_id] == inventory$final_size),
  all(inventory$final_size >= 20L))

namespace <- "reviewed_cytokine_38_trajectory_permutation_v1"
import_table <- data.frame(feature_set = members$signature_id,
  set_name = paste(members$cytokine, toupper(members$direction)),
  feature_id = members$feature_id,
  description = paste("Reviewed competitor-aware cytokine signature",
    members$cytokine, members$direction), stringsAsFactors = FALSE)
import_csv <- file.path(stage_dir, "reviewed_cytokine_38_memberships.csv")
write.csv(import_table, import_csv, row.names = FALSE)
write_tsv(inventory, file.path(stage_dir, "reviewed_cytokine_38_inventory.tsv"))

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "pt"))
stopifnot(identical(as.character(metadata$cell_id), cells),
  length(cells) == 26723L, sum(is.finite(metadata$pt)) == 12642L,
  sum(is.na(metadata$pt)) == 14081L)
use <- is.finite(metadata$pt)
pt_cells <- cells[use]
pt <- as.numeric(metadata$pt[use])
pt_family <- wb_prim_cell_vector_family(rep("bdev_pt", length(pt)),
  pt_cells, pt, meta = list(source_bundle = "bdev_reference",
    source_field = "pt", source_artifact_id = "cell_vector_20260914_162517_7q4y8xbb",
    authoritative_values = "exact Bdev metadata values", direction = "as_is",
    fitted_model = FALSE))
stopifnot(identical(pt_family$value, pt))

session <- wb_session_create(bundle, workspace_root = file.path(stage_dir, "session"),
  session_id = "reviewed_cytokine_38_permutation_validation")
log_event("Importing 38 reviewed cytokine signatures")
imported <- wb_session_run_op(session, bundle, "feature_set_import",
  inputs = list(matrix = "RNA_log1p_cp10k", cell_set = "all", min_cells = 10L,
    source_path = import_csv, source_format = "table", namespace = namespace,
    collection = "reviewed_cytokine_signatures", species = "Mus musculus",
    organism = "Mus musculus", gene_id_type = "gene_symbol",
    set_id_col = "feature_set", set_name_col = "set_name",
    gene_col = "feature_id", description_col = "description",
    split_by_collection = FALSE, overwrite = FALSE, strict = TRUE,
    created_by = "scworkbench_Bdev trajectory permutation validation",
    source_id = "all_cytokine_refined_ucell_v1_20260925",
    source_version = "2026-09-25", origin = "reviewed Parse signatures",
    is_readonly = TRUE),
  notes = "[reviewed_cytokine_38_import] Exact memberships for the 38 reviewed passing cytokine signatures; staged only.")
session <- imported$session
import_manifest <- imported$run$outputs$import_manifest$table
stopifnot(nrow(import_manifest) == 38L,
  all(import_manifest$status == "imported"),
  setequal(import_manifest$feature_set, inventory$signature_id))
write_tsv(imported$run$outputs$import_summary$table,
  file.path(stage_dir, "import_summary.tsv"))
write_tsv(import_manifest, file.path(stage_dir, "import_manifest.tsv"))

log_event("Running overlap QC")
qc <- wb_session_run_op(session, bundle, "feature_set_overlap_qc",
  inputs = list(matrix = "RNA_log1p_cp10k", cell_set = "all",
    namespace = namespace, feature_sets = inventory$signature_id,
    min_set_size = 5L, max_set_size = 5000L,
    min_detected_fraction = 0, deduplicate_sets = FALSE),
  notes = "[reviewed_cytokine_38_overlap_qc] Matrix overlap QC before the trajectory screen.")
session <- qc$session
write_tsv(qc$run$outputs$overlap_qc$table,
  file.path(stage_dir, "overlap_qc.tsv"))
write_tsv(qc$run$outputs$namespace_summary$table,
  file.path(stage_dir, "namespace_summary.tsv"))

log_event("Running native scWorkbench validation with", n_permutations,
  "permutations per signature")
native <- wb_session_run_op(session, bundle, "trajectory_geneset_dynamics",
  inputs = list(matrix = "RNA_log1p_cp10k", pseudotime = pt_family,
    pseudotime_key = "bdev_pt", membership = pt_cells, cell_set = "all",
    namespace = namespace, feature_sets = inventory$signature_id,
    score_method = "mean", pseudotime_direction = "as_is", n_bins = 30L,
    summary_fun = "mean", min_cells_per_bin = 30L,
    min_set_size = 5L, max_set_size = 5000L, top_n = 38L,
    top_driver_genes = 0L, min_amplitude = 0, min_bins_observed = 5L,
    min_features_detected = 1L, min_detected_fraction = 0,
    smooth_method = "none", bootstrap_n = 0L, bootstrap_seed = 1729L,
    n_permutations = n_permutations, permutation_seed = 1729L,
    p_adjust_method = "BH"),
  notes = paste0("[reviewed_cytokine_38_native_permutation] Exact Bdev pt; ",
    "RNA_log1p_cp10k mean member expression; 30 equal-width bins; ",
    n_permutations, " native permutations; BH within the 38-set collection."))
session <- native$session
native_curves <- native$run$outputs$signature_table$table
native_results <- native$run$outputs$dynamics_table$table
native_fit <- native$run$outputs$fit_summary$table
stopifnot(nrow(native_results) == 38L,
  setequal(native_results$feature_set, inventory$signature_id),
  all(native_results$n_permutations == n_permutations))
write_tsv_gz(native_curves, file.path(stage_dir, "native_trajectory_curves.tsv.gz"))
write_tsv(native_results, file.path(stage_dir, "native_dynamics.tsv"))
write_tsv(native_fit, file.path(stage_dir, "native_fit_summary.tsv"))

# Optimized equivalent of scWorkbench's mean-score, equal-width-bin permutation
# loop. The RNG call order, bin summaries, and dynamic-score formula are kept
# identical and are checked term by term against the native result above.
view <- wb_session_view(session)
set_order <- wb_feature_sets(view, matrix_id = "RNA_log1p_cp10k",
  namespace = namespace)$feature_set
set_order <- intersect(unique(as.character(set_order)), inventory$signature_id)
stopifnot(length(set_order) == 38L)
bins <- as.integer(cut(pt, breaks = 30L, include.lowest = TRUE, labels = FALSE))
bin_counts <- tabulate(bins, nbins = 30L)
retained_bins <- which(bin_counts >= 30L)
pt_by_bin <- vapply(retained_bins, function(b) mean(pt[bins == b]), numeric(1))
qn <- max(1L, ceiling(length(retained_bins) * 0.2))

dynamic_scores <- function(curves) {
  if (is.null(dim(curves))) curves <- matrix(curves, ncol = 1L)
  amp <- apply(curves, 2L, function(y) max(y) - min(y))
  sdv <- apply(curves, 2L, stats::sd)
  delta <- colMeans(curves[(nrow(curves) - qn + 1L):nrow(curves), , drop = FALSE]) -
    colMeans(curves[seq_len(qn), , drop = FALSE])
  rho <- apply(curves, 2L, function(y)
    suppressWarnings(stats::cor(pt_by_bin, y, method = "spearman")))
  amp + sdv + abs(delta) + ifelse(is.finite(rho), abs(rho) * sdv, 0)
}

log_event("Running vectorized equivalent permutations for validation")
set.seed(1729L)
optimized_rows <- vector("list", length(set_order))
optimized_curve_rows <- vector("list", length(set_order))
for (ii in seq_along(set_order)) {
  signature_id <- set_order[[ii]]
  genes <- members$feature_id[members$signature_id == signature_id]
  expression <- as.matrix(wb_matrix_read_features(bundle,
    "RNA_log1p_cp10k", genes, pt_cells, missing = "NA"))
  values <- colMeans(expression, na.rm = TRUE)
  observed_curve <- vapply(retained_bins,
    function(b) mean(values[bins == b]), numeric(1))
  observed_dynamic <- dynamic_scores(observed_curve)[[1L]]
  native_row <- native_results[native_results$feature_set == signature_id, , drop = FALSE]
  stopifnot(nrow(native_row) == 1L,
    isTRUE(all.equal(observed_dynamic, native_row$dynamic_score,
      tolerance = 1e-12)))

  null <- numeric(n_permutations)
  batch_size <- 100L
  starts <- seq.int(1L, n_permutations, by = batch_size)
  for (start in starts) {
    take <- min(batch_size, n_permutations - start + 1L)
    permuted <- replicate(take, sample(values, length(values), replace = FALSE))
    if (is.null(dim(permuted))) permuted <- matrix(permuted, ncol = 1L)
    perm_curves <- vapply(retained_bins,
      function(b) colMeans(permuted[bins == b, , drop = FALSE]),
      numeric(ncol(permuted)))
    perm_curves <- t(perm_curves)
    null[start:(start + take - 1L)] <- dynamic_scores(perm_curves)
  }
  p_dynamic <- (1 + sum(null >= observed_dynamic)) / (1 + n_permutations)
  null_mean <- mean(null)
  null_sd <- stats::sd(null)
  optimized_rows[[ii]] <- data.frame(feature_set = signature_id,
    dynamic_score = observed_dynamic, p_dynamic = p_dynamic,
    null_mean_dynamic_score = null_mean,
    null_sd_dynamic_score = null_sd,
    z_dynamic_score = (observed_dynamic - null_mean) / null_sd,
    n_permutations = n_permutations, stringsAsFactors = FALSE)
  optimized_curve_rows[[ii]] <- data.frame(feature_set = signature_id,
    bin = retained_bins, pseudotime = pt_by_bin, value = observed_curve,
    stringsAsFactors = FALSE)
  if (ii %% 5L == 0L || ii == length(set_order))
    log_event("Optimized validation completed", ii, "of", length(set_order))
}
optimized <- do.call(rbind, optimized_rows)
optimized$q_dynamic <- p.adjust(optimized$p_dynamic, method = "BH")
optimized$rank <- rank(-optimized$dynamic_score, ties.method = "first")
optimized <- optimized[order(optimized$rank), , drop = FALSE]
optimized_curves <- do.call(rbind, optimized_curve_rows)

compare_cols <- c("dynamic_score", "p_dynamic", "q_dynamic",
  "null_mean_dynamic_score", "null_sd_dynamic_score", "z_dynamic_score")
validation <- merge(native_results[, c("feature_set", compare_cols)],
  optimized[, c("feature_set", compare_cols)], by = "feature_set",
  suffixes = c("_native", "_optimized"), sort = FALSE)
for (column in compare_cols) {
  validation[[paste0(column, "_abs_diff")]] <- abs(
    validation[[paste0(column, "_native")]] -
      validation[[paste0(column, "_optimized")]])
}
curve_validation <- merge(native_curves[, c("feature_set", "bin", "pseudotime", "value")],
  optimized_curves, by = c("feature_set", "bin"), suffixes = c("_native", "_optimized"))
curve_validation$pseudotime_abs_diff <- abs(curve_validation$pseudotime_native -
  curve_validation$pseudotime_optimized)
curve_validation$value_abs_diff <- abs(curve_validation$value_native -
  curve_validation$value_optimized)
validation_pass <- all(vapply(compare_cols, function(column)
  max(validation[[paste0(column, "_abs_diff")]], na.rm = TRUE) <= 1e-12,
  logical(1))) && max(curve_validation$pseudotime_abs_diff) <= 1e-12 &&
  max(curve_validation$value_abs_diff) <= 1e-12
write_tsv(optimized, file.path(stage_dir, "optimized_dynamics.tsv"))
write_tsv(validation, file.path(stage_dir, "native_vs_optimized_validation.tsv"))
write_tsv_gz(curve_validation, file.path(stage_dir,
  "native_vs_optimized_curve_validation.tsv.gz"))
if (!validation_pass) stop("Optimized permutation validation did not match native scWorkbench")

significant <- native_results[is.finite(native_results$q_dynamic) &
  native_results$q_dynamic < 0.05 & native_results$passes_filters, , drop = FALSE]
significant <- significant[order(significant$q_dynamic,
  -significant$dynamic_score, significant$feature_set), , drop = FALSE]
significant <- merge(significant,
  inventory[, c("signature_id", "cytokine", "direction", "final_size")],
  by.x = "feature_set", by.y = "signature_id", all.x = TRUE, sort = FALSE)
write_tsv(significant, file.path(stage_dir, "significant_cytokine_gene_sets.tsv"))

completion <- list(status = "complete", collection = "reviewed cytokine signatures",
  requested_sets = 38L, tested_sets = nrow(native_results),
  significant_q_lt_0_05 = nrow(significant), n_permutations = n_permutations,
  p_adjust_method = "BH within the 38-set cytokine collection",
  native_validation = list(passed = validation_pass,
    max_abs_differences = setNames(lapply(compare_cols, function(column)
      max(validation[[paste0(column, "_abs_diff")]], na.rm = TRUE)), compare_cols),
    max_curve_value_abs_difference = max(curve_validation$value_abs_diff),
    max_curve_pseudotime_abs_difference = max(curve_validation$pseudotime_abs_diff)),
  matrix = "RNA_log1p_cp10k", score_method = "mean", n_bins = 30L,
  summary_fun = "mean", min_cells_per_bin = 30L, smooth_method = "none",
  pseudotime_source = "exact Bdev metadata pt", defined_pt_cells = length(pt),
  missing_pt_cells = sum(!use), promoted = FALSE,
  bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = wb_manifest(bundle)$current_state$fingerprint,
  bundle_unchanged = identical(fingerprint_before,
    wb_manifest(bundle)$current_state$fingerprint))
stopifnot(completion$bundle_unchanged)
write_json(completion, file.path(stage_dir, "completion.json"))

preview <- head(significant, 38L)
rows <- if (nrow(preview)) paste0("| ", preview$cytokine, " | ",
  preview$direction, " | ", preview$curve_class, " | ",
  signif(preview$q_dynamic, 4), " | ", signif(preview$dynamic_score, 4),
  " | ", signif(preview$peak_pseudotime, 5), " |") else
  "No signatures passed q_dynamic < 0.05."
report <- c("# Reviewed cytokine trajectory permutation analysis", "",
  paste0("All 38 reviewed cytokine signatures were tested with ",
    n_permutations, " permutations each."),
  paste0("Significant at within-collection BH q < 0.05: ", nrow(significant), "."),
  "The optimized permutation implementation matched native scWorkbench term by term within 1e-12 for curves, dynamic scores, null summaries, p-values, q-values, and z-scores.",
  "", "| Cytokine | Direction | Curve class | q_dynamic | Dynamic score | Peak pt |",
  "|---|---|---|---:|---:|---:|", rows, "",
  "Significance is for overall trajectory dynamics; peak pt localizes the maximum and is not a bin-specific p-value.")
writeLines(report, file.path(stage_dir, "REPORT.md"))

wb_close_bundle(view)
wb_close_bundle(bundle)
file.rename(stage_dir, output_dir)
log_event("Completed validated 38-signature cytokine analysis at", output_dir)
