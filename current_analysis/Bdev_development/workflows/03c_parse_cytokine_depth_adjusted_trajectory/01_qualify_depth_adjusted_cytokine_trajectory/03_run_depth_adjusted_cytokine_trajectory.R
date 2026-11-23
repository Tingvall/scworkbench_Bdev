#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config"))
  stop("Usage: 03_run_depth_adjusted_cytokine_trajectory.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbRNA)
  library(scwbGenePrograms)
  library(scwbTrajectory)
  library(Matrix)
  library(yaml)
  library(jsonlite)
  library(digest)
})

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_depth_adjusted_cytokine_trajectory.log"),
  split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
write_tsv <- function(x, name) write.table(x, file.path(dirs[["tables"]], name),
  sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
write_json <- function(x, name) jsonlite::write_json(x,
  file.path(dirs[["manifests"]], name), pretty = TRUE, auto_unbox = TRUE,
  digits = NA, null = "null", na = "null")
run_ref <- function(run, output) wb_ref_run(run$run_id, output)
artifact_ref <- function(promotion, output) {
  rows <- promotion$promoted$status == "ok" & promotion$promoted$output == output
  stopifnot(sum(rows) == 1L)
  wb_ref_artifact(promotion$promoted$artifact_id[rows], output)
}
read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  out <- lapply(fields, function(x) unique(x[seq.int(3L, length(x))]))
  names(out) <- vapply(fields, `[[`, character(1L), 1L)
  out
}
package_identity <- function(package, expected) {
  observed <- as.character(packageVersion(package))
  stopifnot(identical(observed, as.character(expected)))
  observed
}

versions <- list(
  scworkbench = package_identity("scworkbench", cfg$runtime$scworkbench$version),
  scwbRNA = package_identity("scwbRNA", cfg$runtime$scwbRNA$version),
  scwbGenePrograms = package_identity("scwbGenePrograms", cfg$runtime$scwbGenePrograms$version),
  scwbTrajectory = package_identity("scwbTrajectory", cfg$runtime$scwbTrajectory$version),
  fgsea = package_identity("fgsea", cfg$runtime$fgsea),
  BiocParallel = package_identity("BiocParallel", cfg$runtime$BiocParallel),
  slingshot = package_identity("slingshot", cfg$runtime$slingshot))
stopifnot(identical(digest::digest(cfg$runtime$sif, algo = "sha256", file = TRUE),
  cfg$runtime$sif_sha256))

wb_operation_pack_load("scwbRNA")
wb_operation_pack_load("scwbGenePrograms")
wb_operation_pack_load("scwbTrajectory")
pack_records <- wb_operation_packs()
manifest_by_package <- stats::setNames(as.character(pack_records$manifest_sha256),
  as.character(pack_records$package))
stopifnot(identical(unname(manifest_by_package[["scwbRNA"]]),
    cfg$runtime$scwbRNA$manifest_sha256),
  identical(unname(manifest_by_package[["scwbGenePrograms"]]),
    cfg$runtime$scwbGenePrograms$manifest_sha256),
  identical(unname(manifest_by_package[["scwbTrajectory"]]),
    cfg$runtime$scwbTrajectory$manifest_sha256))

operations <- c(cfg$normalization$operation, "metadata_cell_vector",
  "feature_set_import", "scwbgeneprograms::feature_set_overlap_qc",
  cfg$score_generation$operation, cfg$methods$dynamics$operation,
  cfg$methods$intervals$operation, cfg$methods$gsea$operation)
for (operation in operations)
  write_json(wb_op_interface(operation), paste0("interface_",
    gsub("[^A-Za-z0-9]+", "_", operation), ".json"))

source_manifest <- file.path(cfg$inputs$source_bundle$path, "manifest.json")
stopifnot(identical(digest::digest(source_manifest, algo = "sha256", file = TRUE),
  cfg$inputs$source_bundle$manifest_sha256))
gmt_path <- file.path(dirs[["tables"]], "01_parse38_cytokine_gene_sets.gmt")
stopifnot(identical(digest::digest(gmt_path, algo = "sha256", file = TRUE),
  cfg$inputs$gene_sets$sha256))
gene_sets <- read_gmt(gmt_path)
stopifnot(length(gene_sets) == cfg$inputs$gene_sets$expected_sets,
  sum(lengths(gene_sets)) == cfg$inputs$gene_sets$expected_memberships,
  !anyDuplicated(names(gene_sets)))

log_event("Opening isolated 03c execution clone")
bundle <- wb_open_bundle(cfg$storage$execution_clone, validate = FALSE,
  read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
all_cells <- as.character(wb_cells(bundle, "all"))
stopifnot(length(all_cells) == cfg$inputs$source_bundle$expected_cells)
metadata <- wb_meta(bundle, cells = all_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$sample_field,
  cfg$inputs$depth_field,
  cfg$inputs$preserved_pseudotime_field))
stopifnot(identical(as.character(metadata$cell_id), all_cells),
  all(is.finite(as.numeric(metadata[[cfg$inputs$depth_field]]))),
  all(as.numeric(metadata[[cfg$inputs$depth_field]]) > 0))
preserved_pt_before <- as.numeric(metadata[[cfg$inputs$preserved_pseudotime_field]])
preserved_pt_digest <- digest::digest(serialize(list(cell_id = metadata$cell_id,
  value = preserved_pt_before), NULL, version = 3L), algo = "sha256",
  serialize = FALSE)

log_event("Loading the existing authenticated 03b fitted L4")
l4_input <- cfg$inputs$authenticated_l4
stopifnot(
  identical(digest::digest(l4_input$pseudotime_path, algo = "sha256", file = TRUE),
    l4_input$pseudotime_sha256),
  identical(digest::digest(l4_input$model_path, algo = "sha256", file = TRUE),
    l4_input$model_sha256),
  identical(digest::digest(l4_input$fixture_path, algo = "sha256", file = TRUE),
    l4_input$fixture_sha256))
pseudotime_input <- readRDS(l4_input$pseudotime_path)
stopifnot(inherits(pseudotime_input, "wb_primitive"),
  identical(pseudotime_input$kind, "cell_vector_family"))
saveRDS(readRDS(l4_input$model_path),
  file.path(dirs[["objects"]], "01_authenticated_l4_model.rds"), version = 3L)

family <- pseudotime_input
selected <- family$key == cfg$fit$lineage & is.finite(as.numeric(family$value))
l4 <- data.frame(cell_id = as.character(family$cell_id[selected]),
  pseudotime = as.numeric(family$value[selected]), stringsAsFactors = FALSE)
stopifnot(nrow(l4) > 0L, !anyDuplicated(l4$cell_id))
l4$normalized_pseudotime <- (l4$pseudotime - min(l4$pseudotime)) /
  diff(range(l4$pseudotime))
l4$bin <- as.integer(cut(l4$normalized_pseudotime,
  breaks = seq(0, 1, length.out = cfg$fixture$n_bins + 1L),
  include.lowest = TRUE, labels = FALSE))
fixture <- read.delim(l4_input$fixture_path, check.names = FALSE,
  quote = "\"", na.strings = "NA")
fixture_cells <- as.character(fixture$cell_id)
stopifnot(length(fixture_cells) <= cfg$fixture$maximum_cells,
  length(fixture_cells) >= cfg$fixture$n_bins * 20L,
  !anyDuplicated(fixture_cells), all(fixture_cells %in% all_cells),
  identical(fixture_cells, all_cells[all_cells %in% fixture_cells]))
l4_match <- match(fixture_cells, l4$cell_id)
stopifnot(!anyNA(l4_match),
  isTRUE(all.equal(as.numeric(fixture$pseudotime), l4$pseudotime[l4_match],
    tolerance = 1e-12)),
  isTRUE(all.equal(as.numeric(fixture$normalized_pseudotime),
    l4$normalized_pseudotime[l4_match], tolerance = 1e-12)),
  identical(as.integer(fixture$bin), l4$bin[l4_match]))
expected_samples <- as.character(metadata[[cfg$inputs$sample_field]][
  match(fixture_cells, metadata$cell_id)])
stopifnot(identical(as.character(fixture$sample_id), expected_samples))
fixture$depth <- as.numeric(metadata[[cfg$inputs$depth_field]][
  match(fixture$cell_id, metadata$cell_id)])
stopifnot(!anyNA(fixture$sample_id), !anyNA(fixture$depth),
  length(unique(fixture$bin)) == cfg$fixture$n_bins)
write_tsv(fixture, "02_l4_fixture_cells.tsv")
sample_audit <- as.data.frame(table(bin = fixture$bin,
  sample_id = fixture$sample_id), stringsAsFactors = FALSE)
sample_audit <- sample_audit[sample_audit$Freq > 0L, , drop = FALSE]
names(sample_audit)[3L] <- "cells"
write_tsv(sample_audit, "03_l4_bin_sample_composition.tsv")

log_event("Materializing authenticated nCount_RNA depth vector")
depth_session <- wb_session_create(bundle, workspace_root = cfg$storage$sessions,
  session_id = "materialize_rna_depth")
depth_run <- wb_session_run_op(depth_session, bundle, "metadata_cell_vector",
  inputs = list(source_column = cfg$inputs$depth_field,
    vector_name = "03c_nCount_RNA_depth", value_type = "continuous",
    missing_policy = "error", cell_set = "all", min_cells = 20L),
  notes = "Authenticated RNA depth covariate for 03c.")
depth_promotion <- wb_session_promote_current(bundle, depth_run$session,
  outputs = "vector", name_prefix = "03c_depth",
  description = "Authenticated nCount_RNA depth vector for 03c.",
  created_by = "scworkbench_Bdev development qualification", fail_fast = TRUE)
depth_ref <- artifact_ref(depth_promotion, "vector")

log_event("Auditing RNA-count library size against nCount_RNA")
count_depth <- stats::setNames(numeric(length(fixture_cells)), fixture_cells)
iterator <- wb_slice_chunks(bundle, cfg$normalization$source_matrix,
  cells = fixture_cells, axis = "features", chunk_size = 1000L)
while (iterator$has_next()) {
  chunk <- iterator$next_chunk()
  count_depth <- count_depth + as.numeric(Matrix::colSums(chunk$value))
}
metadata_depth <- stats::setNames(fixture$depth, fixture$cell_id)[fixture_cells]
depth_delta <- count_depth - metadata_depth
depth_audit <- data.frame(cell_id = fixture_cells,
  matrix_library_size = as.numeric(count_depth),
  metadata_depth = as.numeric(metadata_depth), delta = as.numeric(depth_delta),
  stringsAsFactors = FALSE)
write_tsv(depth_audit, "04_depth_library_size_audit.tsv")
if (isTRUE(cfg$normalization$require_depth_match))
  stopifnot(max(abs(depth_delta)) <= cfg$normalization$depth_match_tolerance)

log_event("Creating authenticated RNA_log1p_cp10k matrix")
normalization_session <- wb_session_create(bundle,
  workspace_root = cfg$storage$sessions, session_id = "normalize_rna_log1p_cp10k")
normalization <- wb_session_run_op(normalization_session, bundle,
  cfg$normalization$operation, inputs = list(
    matrix = cfg$normalization$source_matrix, cell_ids = fixture_cells,
    min_cells = 20L, output_matrix_id = cfg$normalization$output_matrix,
    scale_factor = as.numeric(cfg$normalization$scale_factor),
    chunk_size = as.integer(cfg$normalization$chunk_size), overwrite = FALSE),
  notes = "Authenticated log1p(counts/library_size*10000) matrix for 03c.")
write_tsv(normalization$run$outputs$summary$table,
  "05_RNA_log1p_cp10k_normalization_summary.tsv")
normalization_promotion <- wb_session_promote_current(bundle,
  normalization$session, outputs = "summary", name_prefix = "03c_log1p_cp10k",
  description = "Fixture-scoped RNA log1p CP10K normalization for 03c.",
  created_by = "scworkbench_Bdev development qualification", fail_fast = TRUE)
stopifnot(identical(as.character(wb_cells_with_matrix(bundle,
  cfg$normalization$output_matrix)), fixture_cells))

import_one <- function(matrix_id, table_name) {
  existing <- wb_feature_sets(bundle, matrix_id = matrix_id,
    namespace = cfg$gene_set_import$namespace)
  if (nrow(existing) == cfg$inputs$gene_sets$expected_sets)
    return(list(reused = TRUE, imported = nrow(existing)))
  stopifnot(nrow(existing) == 0L)
  session <- wb_session_create(bundle, workspace_root = cfg$storage$sessions,
    session_id = paste0("import_parse38_", matrix_id))
  result <- wb_session_run_op(session, bundle, "feature_set_import", inputs = list(
    matrix = matrix_id, cell_ids = fixture_cells, min_cells = 20L,
    source_path = gmt_path, source_format = "gmt",
    namespace = cfg$gene_set_import$namespace,
    collection = cfg$gene_set_import$collection,
    species = cfg$gene_set_import$species,
    gene_id_type = cfg$gene_set_import$gene_id_type,
    overwrite = TRUE, strict = FALSE,
    created_by = "scworkbench_Bdev development qualification",
    source_id = cfg$run_identity, source_version = "20261103",
    origin = "reconstructed_archived_parse38_method", is_readonly = TRUE),
    notes = paste("03c Parse38 import for", matrix_id))
  imported <- result$run$outputs$import_manifest$table
  stopifnot(nrow(imported) == cfg$inputs$gene_sets$expected_sets,
    all(imported$status == "imported"))
  wb_session_promote_current(bundle, result$session,
    outputs = c("import_summary", "import_manifest", "rejected_sets"),
    name_prefix = paste0("03c_parse38_", matrix_id),
    description = paste("03c Parse38 registry for", matrix_id),
    created_by = "scworkbench_Bdev development qualification", fail_fast = TRUE)
  write_tsv(imported, table_name)
  list(reused = FALSE, imported = nrow(imported))
}

imports <- list(
  log1p_cp10k = import_one(cfg$score_generation$matrix,
    "06_RNA_log1p_cp10k_feature_set_import.tsv"),
  SCT_data = import_one(cfg$methods$gsea$matrix,
    "07_SCT_data_feature_set_import.tsv"))

analysis_session <- wb_session_create(bundle, workspace_root = cfg$storage$sessions,
  session_id = "depth_adjusted_cytokine_trajectory_qualification")
log_event("Auditing Parse38 overlap")
log_overlap <- wb_session_run_op(analysis_session, bundle,
  "scwbgeneprograms::feature_set_overlap_qc", inputs = list(
    matrix = cfg$score_generation$matrix, cell_ids = fixture_cells,
    min_cells = 20L, namespace = cfg$gene_set_import$namespace,
    min_set_size = as.integer(cfg$score_generation$min_set_size),
    min_detected_fraction = 0, deduplicate_sets = TRUE),
  notes = "Parse38 log1p-CP10K overlap audit.")
analysis_session <- log_overlap$session
sct_overlap <- wb_session_run_op(analysis_session, bundle,
  "scwbgeneprograms::feature_set_overlap_qc", inputs = list(
    matrix = cfg$methods$gsea$matrix, cell_ids = fixture_cells,
    min_cells = 20L, namespace = cfg$gene_set_import$namespace,
    min_set_size = as.integer(cfg$methods$gsea$min_set_size),
    min_detected_fraction = 0, deduplicate_sets = TRUE),
  notes = "Parse38 SCT-data overlap audit.")
analysis_session <- sct_overlap$session
write_tsv(log_overlap$run$outputs$overlap_qc$table,
  "08_RNA_log1p_cp10k_overlap_qc.tsv")
write_tsv(sct_overlap$run$outputs$overlap_qc$table,
  "09_SCT_data_overlap_qc.tsv")
stopifnot(sum(log_overlap$run$outputs$overlap_qc$table$screening_status == "include") == 38L,
  sum(sct_overlap$run$outputs$overlap_qc$table$screening_status == "include") == 38L)

pseudotime_direction <- if (identical(cfg$fit$orientation, "native"))
  "as_is" else cfg$fit$orientation
stopifnot(pseudotime_direction %in% c("as_is", "invert"))
common <- list(cell_ids = fixture_cells, min_cells = 20L,
  pseudotime_key = cfg$fit$lineage,
  pseudotime_direction = pseudotime_direction)

log_event("Computing df=8 gene-level depth-adjusted scores")
score <- wb_session_run_op(analysis_session, bundle,
  cfg$score_generation$operation, inputs = c(list(
    matrix = cfg$score_generation$matrix, pseudotime = pseudotime_input,
    depth = depth_ref, namespace = cfg$gene_set_import$namespace,
    feature_sets = names(gene_sets),
    score_method = cfg$score_generation$score_method,
    adjust_depth = isTRUE(cfg$score_generation$adjust_depth),
    depth_field = cfg$inputs$depth_field,
    adjustment_spline_df = as.integer(cfg$score_generation$adjustment_spline_df),
    min_set_size = as.integer(cfg$score_generation$min_set_size),
    max_set_size = as.integer(cfg$score_generation$max_set_size),
    chunk_size = as.integer(cfg$score_generation$chunk_size)), common),
  notes = "Archived-style per-gene depth adjustment on RNA_log1p_cp10k.")
analysis_session <- score$session
score_primitive <- score$run$outputs$adjusted_scores
score_table <- data.frame(signal_key = score_primitive$key,
  cell_id = score_primitive$cell_id, value = score_primitive$value,
  stringsAsFactors = FALSE)
write_tsv(score_table, "10_depth_adjusted_scores.tsv")
write_tsv(score$run$outputs$adjustment_diagnostics$table,
  "11_depth_adjustment_diagnostics.tsv")
write_tsv(score$run$outputs$fit_summary$table,
  "12_score_fit_summary.tsv")
stopifnot(isTRUE(score_primitive$meta$depth_adjusted),
  identical(as.integer(score_primitive$meta$adjustment_spline_df), 8L),
  length(unique(score_table$signal_key)) == 38L)

signal_ref <- run_ref(score$run, "adjusted_scores")
log_event("Running descriptive dynamics on adjusted scores")
dynamics <- wb_session_run_op(analysis_session, bundle,
  cfg$methods$dynamics$operation, inputs = c(list(
    signal = signal_ref, pseudotime = pseudotime_input,
    n_bins = as.integer(cfg$methods$dynamics$n_bins),
    summary_fun = cfg$methods$dynamics$summary_fun,
    smooth_method = cfg$methods$dynamics$smooth_method,
    bootstrap_n = as.integer(cfg$methods$dynamics$bootstrap_n),
    min_cells_per_bin = as.integer(cfg$methods$dynamics$min_cells_per_bin),
    min_bins_observed = as.integer(cfg$methods$dynamics$min_bins_observed)), common),
  notes = "03c adjusted score dynamics on the common L4 scale.")
analysis_session <- dynamics$session

log_event("Running permutation-qualified intervals on the same adjusted scores")
intervals <- wb_session_run_op(analysis_session, bundle,
  cfg$methods$intervals$operation, inputs = c(list(
    signal = signal_ref, pseudotime = pseudotime_input,
    window_width_fraction = as.numeric(cfg$methods$intervals$window_width_fraction),
    window_step_fraction = as.numeric(cfg$methods$intervals$window_step_fraction),
    summary_fun = cfg$methods$intervals$summary_fun,
    permutation_mode = cfg$methods$intervals$permutation_mode,
    n_permutations = as.integer(cfg$methods$intervals$n_permutations),
    seed = as.integer(cfg$methods$intervals$seed),
    min_cells_per_window = as.integer(cfg$methods$intervals$min_cells_per_window),
    min_samples_per_window = as.integer(cfg$methods$intervals$min_samples_per_window),
    global_q_cutoff = as.numeric(cfg$methods$intervals$global_q_cutoff),
    local_q_cutoff = as.numeric(cfg$methods$intervals$local_q_cutoff),
    merge_gap_fraction = as.numeric(cfg$methods$intervals$merge_gap_fraction),
    adjust_depth = TRUE, depth_field = cfg$inputs$depth_field,
    adjustment_spline_df = as.integer(cfg$score_generation$adjustment_spline_df)), common),
  notes = "03c intervals reuse already adjusted scores without double adjustment.")
analysis_session <- intervals$session

log_event("Running SCT-based GSEA with depth adjustment disabled")
gsea <- wb_session_run_op(analysis_session, bundle,
  cfg$methods$gsea$operation, inputs = c(list(
    matrix = cfg$methods$gsea$matrix, pseudotime = pseudotime_input,
    namespace = cfg$gene_set_import$namespace,
    adjust_depth = isTRUE(cfg$methods$gsea$adjust_depth),
    analysis_scope = cfg$methods$gsea$analysis_scope,
    global_rank_method = cfg$methods$gsea$global_rank_method,
    window_width_fraction = as.numeric(cfg$methods$gsea$window_width_fraction),
    window_step_fraction = as.numeric(cfg$methods$gsea$window_step_fraction),
    permutation_mode = cfg$methods$gsea$permutation_mode,
    n_permutations = as.integer(cfg$methods$gsea$n_permutations),
    seed = as.integer(cfg$methods$gsea$seed),
    min_cells_per_window = as.integer(cfg$methods$gsea$min_cells_per_window),
    min_samples_per_window = as.integer(cfg$methods$gsea$min_samples_per_window),
    min_set_size = as.integer(cfg$methods$gsea$min_set_size),
    max_set_size = as.integer(cfg$methods$gsea$max_set_size),
    score_type = cfg$methods$gsea$score_type,
    gsea_param = as.numeric(cfg$methods$gsea$gsea_param),
    sample_size = as.integer(cfg$methods$gsea$sample_size),
    n_perm_simple = as.integer(cfg$methods$gsea$n_perm_simple),
    eps = as.numeric(cfg$methods$gsea$eps),
    p_adjust_method = cfg$methods$gsea$p_adjust_method,
    global_q_cutoff = as.numeric(cfg$methods$gsea$global_q_cutoff),
    local_q_cutoff = as.numeric(cfg$methods$gsea$local_q_cutoff),
    merge_gap_fraction = as.numeric(cfg$methods$gsea$merge_gap_fraction),
    consensus_leading_edge_fraction = as.numeric(
      cfg$methods$gsea$consensus_leading_edge_fraction),
    chunk_size = as.integer(cfg$methods$gsea$chunk_size)), common),
  notes = "03c SCT-data GSEA; no additional depth regression.")

for (entry in list(
  list(run = dynamics$run, prefix = "dynamics", start = 13L),
  list(run = intervals$run, prefix = "intervals", start = 17L))) {
  for (ii in seq_along(entry$run$outputs))
    write_tsv(entry$run$outputs[[ii]]$table, sprintf("%02d_%s_%s.tsv",
      entry$start + ii - 1L, entry$prefix, names(entry$run$outputs)[[ii]]))
}
gsea_order <- c("global_enrichment", "window_enrichment",
  "enrichment_intervals", "leading_edge", "ranked_features", "fit_summary")
for (ii in seq_along(gsea_order))
  write_tsv(gsea$run$outputs[[gsea_order[[ii]]]]$table,
    sprintf("%02d_gsea_%s.tsv", 23L + ii - 1L, gsea_order[[ii]]))

dynamics_table <- dynamics$run$outputs$dynamics_table$table
interval_global <- intervals$run$outputs$global_test$table
gsea_global <- gsea$run$outputs$global_enrichment$table
comparison <- merge(dynamics_table,
  interval_global[, c("signal_key", "p_global", "q_global",
    "global_significant", "observed_signal_range"), drop = FALSE],
  by = "signal_key", all = TRUE, sort = FALSE)
comparison$feature_set <- sub(paste0("^", cfg$fit$lineage, "::"), "",
  comparison$signal_key)
gsea_summary <- gsea_global[, c("feature_set", "normalized_enrichment_score",
  "p_value", "q_within_ranking", "q_global_family", "supported"), drop = FALSE]
names(gsea_summary)[-1L] <- paste0("gsea_", names(gsea_summary)[-1L])
comparison <- merge(comparison, gsea_summary, by = "feature_set",
  all = TRUE, sort = FALSE)
write_tsv(comparison, "29_method_comparison.tsv")

interval_adjustment <- intervals$run$outputs$adjustment_summary$table
interval_fit <- intervals$run$outputs$fit_summary$table
gsea_fit <- gsea$run$outputs$fit_summary$table
stopifnot(all(interval_adjustment$depth_adjustment_status == "already_adjusted"),
  identical(as.character(interval_fit$permutation_mode_used), "unblocked"),
  all(as.character(gsea_fit$permutation_mode) == "unblocked"),
  all(!gsea_fit$adjust_depth),
  identical(as.numeric(interval_fit$window_width_fraction), 0.05),
  identical(as.numeric(interval_fit$window_step_fraction), 0.025))

metadata_after <- wb_meta(bundle, cells = all_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$preserved_pseudotime_field))
preserved_pt_after_digest <- digest::digest(serialize(list(
  cell_id = metadata_after$cell_id,
  value = as.numeric(metadata_after[[cfg$inputs$preserved_pseudotime_field]])),
  NULL, version = 3L), algo = "sha256", serialize = FALSE)
stopifnot(identical(preserved_pt_digest, preserved_pt_after_digest))

diagnostics <- score$run$outputs$adjustment_diagnostics$table
report <- list(
  schema_version = "bdev.parse38_depth_adjusted_trajectory_summary.v1",
  status = "complete_development_qualification",
  run_identity = cfg$run_identity,
  fixture = list(cells = length(fixture_cells), bins = length(unique(fixture$bin)),
    samples = length(unique(fixture$sample_id)), l4_available_cells = nrow(l4)),
  normalization = list(matrix = cfg$normalization$output_matrix,
    formula = "log1p(RNA_counts / per-cell library size * 10000)",
    maximum_depth_delta = max(abs(depth_delta))),
  score = list(feature_sets = length(unique(score_table$signal_key)),
    adjust_depth = TRUE, depth_field = cfg$inputs$depth_field,
    adjustment_spline_df = cfg$score_generation$adjustment_spline_df,
    genes_adjusted = nrow(diagnostics),
    genes_standardized = sum(diagnostics$standardized)),
  methods = list(
    fitted_trajectory = list(
      source_run_identity = l4_input$source_run_identity,
      lineage = cfg$fit$lineage, reused_without_refit = TRUE,
      pseudotime_sha256 = l4_input$pseudotime_sha256,
      model_sha256 = l4_input$model_sha256,
      fixture_sha256 = l4_input$fixture_sha256),
    dynamics = list(run_id = dynamics$run$run_id,
      bins = cfg$methods$dynamics$n_bins),
    intervals = list(run_id = intervals$run$run_id,
      score_adjustment_status = unique(interval_adjustment$depth_adjustment_status),
      permutations = cfg$methods$intervals$n_permutations,
      supported_intervals = nrow(intervals$run$outputs$signal_intervals$table)),
    gsea = list(run_id = gsea$run$run_id, matrix = cfg$methods$gsea$matrix,
      adjust_depth = FALSE, permutations = cfg$methods$gsea$n_permutations,
      supported_intervals = nrow(gsea$run$outputs$enrichment_intervals$table))),
  preservation = list(archive_mutated = FALSE,
    preserved_pt_digest_before = preserved_pt_digest,
    preserved_pt_digest_after = preserved_pt_after_digest,
    preserved_pt_unchanged = TRUE),
  runtime = c(versions, list(pack_manifest_sha256 = as.list(manifest_by_package))))
jsonlite::write_json(report,
  file.path(dirs[["reports"]], "01_depth_adjusted_qualification_summary.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
writeLines(c(
  "# 03c depth-adjusted cytokine trajectory qualification",
  "",
  paste0("Status: complete development qualification on ", length(fixture_cells), " cells."),
  "",
  "The score dynamics and interval analyses use the same reusable cell-level scores built from RNA_log1p_cp10k after per-gene adjustment for nCount_RNA with ns(L4, df=8).",
  "Signal intervals recognized those scores as already adjusted and did not regress depth a second time.",
  "GSEA uses SCT_data with adjust_depth=FALSE.",
  "All analyses use native fitted L4; figures reverse only the display axis to LSK (0) -> CD19 (1).",
  "This remains a bounded development qualification and does not unlock Bdev implementation."
), file.path(dirs[["reports"]], "02_depth_adjusted_qualification_summary.md"))

log_event("03c analysis completed")
