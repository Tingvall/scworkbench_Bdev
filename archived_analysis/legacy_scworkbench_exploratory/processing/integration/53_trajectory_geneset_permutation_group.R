#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L || length(args) > 5L) stop(paste(
  "Usage: 53_trajectory_geneset_permutation_group.R",
  "<bundle_dir> <output_dir> <cytokines|bp|mf|cc> <n_permutations>",
  "[max_sets_for_pilot]"))

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
group <- tolower(args[[3L]])
n_permutations <- as.integer(args[[4L]])
max_sets <- if (length(args) == 5L) as.integer(args[[5L]]) else NA_integer_
stopifnot(group %in% c("cytokines", "bp", "mf", "cc"),
  is.finite(n_permutations), n_permutations > 0L,
  is.na(max_sets) || (is.finite(max_sets) && max_sets > 0L))

if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
  if (group != "cytokines") library(msigdbr)
})

write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ..., "\n")
  flush.console()
}

group_spec <- list(
  cytokines = list(namespace = "parse_cytokine_signatures",
    label = "Parse cytokine signatures", subcollection = NA_character_),
  bp = list(namespace = "msigdb_m5_go_bp_2026_1_mm",
    label = "GO Biological Process", subcollection = "GO:BP"),
  mf = list(namespace = "msigdb_m5_go_mf_2026_1_mm",
    label = "GO Molecular Function", subcollection = "GO:MF"),
  cc = list(namespace = "msigdb_m5_go_cc_2026_1_mm",
    label = "GO Cellular Component", subcollection = "GO:CC")
)[[group]]

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "pt"))
stopifnot(identical(as.character(metadata$cell_id), cells),
  length(cells) == 26723L, sum(is.finite(metadata$pt)) == 12642L,
  sum(is.na(metadata$pt)) == 14081L)
pt_use <- is.finite(metadata$pt)
pt_family <- wb_prim_cell_vector_family(rep("bdev_pt", sum(pt_use)),
  cells[pt_use], as.numeric(metadata$pt[pt_use]),
  meta = list(source_bundle = "bdev_reference", source_field = "pt",
    source_artifact_id = "cell_vector_20260914_162517_7q4y8xbb",
    authoritative_values = "exact Bdev metadata values",
    missing_policy = "preserve in metadata; exclude from trajectory operation",
    direction = "as_is", fitted_model = FALSE))
stopifnot(identical(pt_family$value, as.numeric(metadata$pt[pt_use])))

session_root <- file.path(stage_dir, "session")
session <- wb_session_create(bundle, workspace_root = session_root,
  session_id = paste0("trajectory_permutation_", group))

feature_subset <- NULL
source_version <- NA_character_
source_rows <- NA_integer_
if (group != "cytokines") {
  log_event("Retrieving mouse-native MSigDB", group_spec$subcollection)
  go <- msigdbr::msigdbr(db_species = "MM", species = "Mus musculus",
    collection = "M5", subcollection = group_spec$subcollection)
  required <- c("gene_symbol", "gs_name", "gs_description", "db_version",
    "gs_collection", "gs_subcollection")
  stopifnot(all(required %in% names(go)), nrow(go) > 0L,
    length(unique(go$db_version)) == 1L,
    identical(unique(go$gs_collection), "M5"),
    identical(unique(go$gs_subcollection), group_spec$subcollection))
  source_version <- unique(go$db_version)
  go <- unique(go[, required])
  go <- go[order(go$gs_name, go$gene_symbol), , drop = FALSE]
  if (is.finite(max_sets)) {
    keep_sets <- head(unique(go$gs_name), max_sets)
    go <- go[go$gs_name %in% keep_sets, , drop = FALSE]
  }
  source_rows <- nrow(go)
  feature_subset <- sort(unique(go$gs_name))
  import_table <- data.frame(feature_set = go$gs_name,
    set_name = go$gs_name, feature_id = go$gene_symbol,
    description = go$gs_description, stringsAsFactors = FALSE)
  import_path <- file.path(stage_dir, paste0(group, "_memberships.csv"))
  write.csv(import_table, import_path, row.names = FALSE)
  write_json(list(database = "Mouse MSigDB", db_version = source_version,
    collection = "M5", subcollection = group_spec$subcollection,
    species = "Mus musculus", msigdbr_version = as.character(packageVersion("msigdbr")),
    membership_rows = nrow(import_table), source_sets = length(feature_subset)),
    file.path(stage_dir, paste0(group, "_source.json")))

  log_event("Importing", length(feature_subset), group_spec$label, "sets")
  imported <- wb_session_run_op(session, bundle, "feature_set_import",
    inputs = list(matrix = "RNA_log1p_cp10k", cell_set = "all",
      min_cells = 10L, source_path = import_path, source_format = "table",
      namespace = group_spec$namespace, collection = "M5",
      subcollection = group_spec$subcollection, species = "Mus musculus",
      organism = "Mus musculus", gene_id_type = "gene_symbol",
      set_id_col = "feature_set", set_name_col = "set_name",
      gene_col = "feature_id", description_col = "description",
      split_by_collection = FALSE, overwrite = FALSE, strict = FALSE,
      created_by = "scworkbench_Bdev trajectory permutation screen",
      source_id = "MSigDB mouse-native M5",
      source_version = source_version, origin = "MSigDB via msigdbr",
      is_readonly = TRUE),
    notes = paste0("[", group, "_import] Mouse-native MSigDB ",
      source_version, " ", group_spec$subcollection,
      "; staged only, not promoted."))
  session <- imported$session
  write_tsv(imported$run$outputs$import_summary$table,
    file.path(stage_dir, "import_summary.tsv"))
  write_tsv_gz(imported$run$outputs$import_manifest$table,
    file.path(stage_dir, "import_manifest.tsv.gz"))
  write_tsv_gz(imported$run$outputs$rejected_sets$table,
    file.path(stage_dir, "rejected_sets.tsv.gz"))

  log_event("Running feature-set overlap QC")
  checked <- wb_session_run_op(session, bundle, "feature_set_overlap_qc",
    inputs = list(matrix = "RNA_log1p_cp10k", cell_set = "all",
      namespace = group_spec$namespace, feature_sets = feature_subset,
      min_set_size = 5L, max_set_size = 5000L,
      min_detected_fraction = 0, deduplicate_sets = FALSE),
    notes = paste0("[", group, "_overlap_qc] Pre-screen overlap QC for ",
      group_spec$label, "; no expression or pseudotime changes."))
  session <- checked$session
  write_tsv_gz(checked$run$outputs$overlap_qc$table,
    file.path(stage_dir, "overlap_qc.tsv.gz"))
  write_tsv(checked$run$outputs$namespace_summary$table,
    file.path(stage_dir, "namespace_summary.tsv"))
}

if (group == "cytokines") {
  available <- wb_feature_sets(wb_session_view(session),
    matrix_id = "RNA_log1p_cp10k", namespace = group_spec$namespace)
  stopifnot(nrow(available) == 23L)
  feature_subset <- sort(as.character(available$feature_set))
  if (is.finite(max_sets)) feature_subset <- head(feature_subset, max_sets)
  source_rows <- length(feature_subset)
  source_version <- "accepted bdev_integrated namespace"
}

log_event("Running trajectory_geneset_dynamics for", group_spec$label,
  "with", length(feature_subset), "sets and", n_permutations, "permutations")
dynamics <- wb_session_run_op(session, bundle, "trajectory_geneset_dynamics",
  inputs = list(matrix = "RNA_log1p_cp10k", pseudotime = pt_family,
    pseudotime_key = "bdev_pt", membership = pt_family$cell_id,
    cell_set = "all", namespace = group_spec$namespace,
    feature_sets = feature_subset, score_method = "mean",
    pseudotime_direction = "as_is", n_bins = 30L,
    summary_fun = "mean", min_cells_per_bin = 30L,
    min_set_size = 5L, max_set_size = 5000L, top_n = 100L,
    top_driver_genes = 0L, min_amplitude = 0,
    min_bins_observed = 5L, min_features_detected = 1L,
    min_detected_fraction = 0, smooth_method = "none",
    bootstrap_n = 0L, bootstrap_seed = 1729L,
    n_permutations = n_permutations, permutation_seed = 1729L,
    p_adjust_method = "BH"),
  notes = paste0("[", group, "_trajectory_permutation] Exact Bdev metadata pt; ",
    "RNA_log1p_cp10k mean member expression; 30 equal-width bins; mean bin summaries; ",
    "no smoothing or bootstrap; ", n_permutations,
    " score-to-pseudotime permutations; BH correction within ", group, "."))
session <- dynamics$session

curves <- dynamics$run$outputs$signature_table$table
results <- dynamics$run$outputs$dynamics_table$table
top_sets <- dynamics$run$outputs$top_sets$table
fit <- dynamics$run$outputs$fit_summary$table
results$analysis_group <- group
results$collection_label <- group_spec$label
significant <- results[is.finite(results$q_dynamic) & results$q_dynamic < 0.05 &
  results$passes_filters, , drop = FALSE]
significant <- significant[order(significant$q_dynamic,
  -significant$dynamic_score, significant$feature_set), , drop = FALSE]

write_tsv_gz(curves, file.path(stage_dir, "trajectory_curves.tsv.gz"))
write_tsv_gz(results, file.path(stage_dir, "all_gene_sets.tsv.gz"))
write_tsv(significant, file.path(stage_dir, "significant_gene_sets.tsv"))
write_tsv(top_sets, file.path(stage_dir, "top_100_gene_sets.tsv"))
write_tsv(fit, file.path(stage_dir, "fit_summary.tsv"))

class_counts <- if (nrow(significant)) {
  z <- as.data.frame(table(significant$curve_class), stringsAsFactors = FALSE)
  names(z) <- c("curve_class", "n_significant")
  z[z$n_significant > 0L, , drop = FALSE]
} else data.frame(curve_class = character(0), n_significant = integer(0))
write_tsv(class_counts, file.path(stage_dir, "significant_curve_class_counts.tsv"))

summary <- list(status = "complete", group = group,
  collection_label = group_spec$label, namespace = group_spec$namespace,
  source_version = source_version, source_membership_rows = source_rows,
  requested_sets = length(feature_subset), tested_sets = nrow(results),
  significant_q_lt_0_05 = nrow(significant), n_permutations = n_permutations,
  matrix = "RNA_log1p_cp10k", score_method = "mean", n_bins = 30L,
  summary_fun = "mean", min_cells_per_bin = 30L,
  smooth_method = "none", bootstrap_n = 0L,
  pseudotime_source = "exact Bdev metadata pt",
  defined_pseudotime_cells = sum(pt_use), missing_pseudotime_cells = sum(!pt_use),
  p_adjust_method = "BH within collection", q_threshold = 0.05,
  promoted = FALSE, session_workspace = session$workspace,
  bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = wb_manifest(bundle)$current_state$fingerprint,
  bundle_unchanged = identical(fingerprint_before,
    wb_manifest(bundle)$current_state$fingerprint))
stopifnot(isTRUE(summary$bundle_unchanged))
write_json(summary, file.path(stage_dir, "completion.json"))

top_preview <- head(significant, 25L)
preview_lines <- if (nrow(top_preview)) paste0("| ", top_preview$feature_set,
  " | ", top_preview$curve_class, " | ", signif(top_preview$q_dynamic, 4),
  " | ", signif(top_preview$dynamic_score, 4), " | ",
  signif(top_preview$peak_pseudotime, 5), " |") else
  "No gene sets passed `q_dynamic < 0.05` and `passes_filters`."
report <- c(paste0("# ", group_spec$label,
  " trajectory permutation screen"), "",
  paste0("- Tested sets: ", nrow(results)),
  paste0("- Significant sets (`q_dynamic < 0.05`): ", nrow(significant)),
  paste0("- Permutations per set: ", n_permutations),
  "- Multiple-testing correction: Benjamini-Hochberg within this collection",
  "- Matrix/score: RNA_log1p_cp10k, mean member-gene expression",
  "- Trajectory: exact existing Bdev `pt`; 30 equal-width bins; no smoothing",
  "- Interpretation: significance tests overall trajectory dynamics; peak pseudotime localizes the maximum but is not a bin-specific p-value.",
  "", "## Top significant sets", "",
  if (nrow(top_preview))
    c("| Gene set | Curve class | q_dynamic | Dynamic score | Peak pt |",
      "|---|---:|---:|---:|---:|", preview_lines) else preview_lines,
  "", "Full results are in `significant_gene_sets.tsv` and `all_gene_sets.tsv.gz`.")
writeLines(report, file.path(stage_dir, "REPORT.md"))

wb_close_bundle(wb_session_view(session))
wb_close_bundle(bundle)
file.rename(stage_dir, output_dir)
log_event("Completed", group_spec$label, "screen at", output_dir)
