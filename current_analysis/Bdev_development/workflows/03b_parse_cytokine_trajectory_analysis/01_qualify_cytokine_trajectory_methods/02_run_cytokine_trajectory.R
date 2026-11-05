#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config"))
  stop("Usage: 02_run_cytokine_trajectory.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbTrajectory)
  library(yaml)
  library(jsonlite)
  library(digest)
})

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports", "04_logs", "06_manifests")),
  c("tables", "figures", "reports", "logs", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_cytokine_trajectory.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
run_ref <- function(run, output) wb_ref_run(run$run_id, output)
read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  out <- lapply(fields, function(x) unique(x[seq.int(3L, length(x))]))
  names(out) <- vapply(fields, `[[`, character(1), 1L)
  out
}
save_outputs <- function(run, prefix, start_index) {
  output_names <- names(run$outputs)
  for (ii in seq_along(output_names)) {
    name <- output_names[[ii]]
    primitive <- run$outputs[[name]]
    stopifnot(inherits(primitive, "wb_primitive"), identical(primitive$kind, "obs_table"))
    write_tsv(primitive$table, file.path(dirs[["tables"]], sprintf("%02d_%s_%s.tsv",
      start_index + ii - 1L, prefix, name)))
  }
  start_index + length(output_names)
}
package_identity <- function(package, expected) {
  observed <- as.character(packageVersion(package))
  stopifnot(identical(observed, as.character(expected)))
  observed
}

versions <- list(
  scworkbench = package_identity("scworkbench", cfg$runtime$scworkbench$version),
  scwbTrajectory = package_identity("scwbTrajectory", cfg$runtime$scwbTrajectory$version),
  scwbGenePrograms = package_identity("scwbGenePrograms", cfg$runtime$scwbGenePrograms$version),
  fgsea = package_identity("fgsea", cfg$runtime$fgsea),
  BiocParallel = package_identity("BiocParallel", cfg$runtime$BiocParallel),
  slingshot = package_identity("slingshot", cfg$runtime$slingshot))
stopifnot(identical(digest::digest(cfg$runtime$sif, algo = "sha256", file = TRUE),
  cfg$runtime$sif_sha256))

wb_operation_pack_load("scwbGenePrograms")
wb_operation_pack_load("scwbTrajectory")
pack_records <- wb_operation_packs()
observed_manifests <- stats::setNames(as.character(pack_records$manifest_sha256),
  as.character(pack_records$package))
stopifnot(cfg$runtime$scwbGenePrograms$manifest_sha256 %in% observed_manifests,
  cfg$runtime$scwbTrajectory$manifest_sha256 %in% observed_manifests)

operations <- c("feature_set_import", "scwbgeneprograms::feature_set_overlap_qc",
  "scwbgeneprograms::ucell", "scwbtrajectory::slingshot_fit",
  cfg$methods$dynamics$operation, cfg$methods$intervals$operation,
  cfg$methods$gsea$operation)
for (operation in operations)
  write_json(wb_op_interface(operation), file.path(dirs[["manifests"]],
    paste0("interface_", gsub("[^A-Za-z0-9]+", "_", operation), ".json")))

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

log_event("Opening isolated 03b execution clone")
bundle <- wb_open_bundle(cfg$storage$execution_clone, validate = FALSE,
  read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
all_cells <- as.character(wb_cells(bundle, "all"))
stopifnot(length(all_cells) == cfg$inputs$source_bundle$expected_cells)
metadata <- wb_meta(bundle, cells = all_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$cluster_field,
  cfg$inputs$sample_field, cfg$inputs$preserved_pseudotime_field))
stopifnot(identical(as.character(metadata$cell_id), all_cells))
preserved_pt_before <- as.numeric(metadata[[cfg$inputs$preserved_pseudotime_field]])
preserved_pt_digest <- digest::digest(serialize(list(cell_id = metadata$cell_id,
  value = preserved_pt_before), NULL, version = 3L), algo = "sha256", serialize = FALSE)

import_one <- function(matrix_id) {
  existing <- wb_feature_sets(bundle, matrix_id = matrix_id,
    namespace = cfg$gene_set_import$namespace)
  if (nrow(existing) == cfg$inputs$gene_sets$expected_sets)
    return(list(reused = TRUE, run_id = "", imported = nrow(existing)))
  stopifnot(nrow(existing) == 0L)
  session <- wb_session_create(bundle, workspace_root = cfg$storage$sessions,
    session_id = paste0("import_parse38_", matrix_id))
  result <- wb_session_run_op(session, bundle, "feature_set_import", inputs = list(
    matrix = matrix_id, cell_ids = all_cells, min_cells = 10L,
    source_path = gmt_path, source_format = "gmt",
    namespace = cfg$gene_set_import$namespace,
    collection = cfg$gene_set_import$collection,
    species = cfg$gene_set_import$species,
    gene_id_type = cfg$gene_set_import$gene_id_type,
    overwrite = TRUE, strict = FALSE,
    created_by = "scworkbench_Bdev development qualification",
    source_id = cfg$run_identity, source_version = "20261103",
    origin = "reconstructed_archived_parse38_method", is_readonly = TRUE),
    notes = paste("Development-only Parse38 import for", matrix_id))
  imported <- result$run$outputs$import_manifest$table
  stopifnot(nrow(imported) == cfg$inputs$gene_sets$expected_sets,
    all(imported$status == "imported"))
  wb_session_promote_current(bundle, result$session,
    outputs = c("import_summary", "import_manifest", "rejected_sets"),
    name_prefix = paste0("parse38_", matrix_id),
    description = paste("Development-only Parse38 registry for", matrix_id),
    created_by = "scworkbench_Bdev development qualification", fail_fast = TRUE)
  write_tsv(imported, file.path(dirs[["tables"]],
    paste0(if (identical(matrix_id, "RNA_counts")) "02" else "03",
      "_", matrix_id, "_feature_set_import.tsv")))
  list(reused = FALSE, run_id = result$run$run_id, imported = nrow(imported))
}

imports <- list(RNA_counts = import_one(cfg$signal_generation$matrix),
  SCT_data = import_one(cfg$methods$gsea$matrix))

session <- wb_session_create(bundle, workspace_root = cfg$storage$sessions,
  session_id = "parse38_cytokine_trajectory_qualification")
clusters <- wb_prim_cell_vector(stats::setNames(
  as.character(metadata[[cfg$inputs$cluster_field]]), metadata$cell_id),
  value_type = "categorical")
log_event("Fitting authenticated Slingshot lineage family")
fit <- wb_session_run_op(session, bundle, "scwbtrajectory::slingshot_fit",
  inputs = list(fit_space = cfg$inputs$fit_space, clusters = clusters,
    start.clus = cfg$fit$start_cluster,
    exclude_clusters = unlist(cfg$fit$exclude_clusters, use.names = FALSE),
    stretch = as.numeric(cfg$fit$stretch), omega = isTRUE(cfg$fit$omega),
    omega_scale = as.numeric(cfg$fit$omega_scale)),
  notes = "Development-only refit of the approved Bdev Slingshot recipe for 03b qualification.")
session <- fit$session
family <- fit$run$outputs$pseudotime
selected <- family$key == cfg$fit$lineage & is.finite(as.numeric(family$value))
l4 <- data.frame(cell_id = as.character(family$cell_id[selected]),
  pseudotime = as.numeric(family$value[selected]), stringsAsFactors = FALSE)
stopifnot(nrow(l4) > 0L, !anyDuplicated(l4$cell_id))
l4$normalized_pseudotime <- (l4$pseudotime - min(l4$pseudotime)) /
  diff(range(l4$pseudotime))
l4$bin <- as.integer(cut(l4$normalized_pseudotime,
  breaks = seq(0, 1, length.out = cfg$fixture$n_bins + 1L),
  include.lowest = TRUE, labels = FALSE))
bundle_order <- stats::setNames(seq_along(all_cells), all_cells)
l4$bundle_order <- unname(bundle_order[l4$cell_id])
choose_cells <- function(rows, n) {
  rows <- rows[order(rows$bundle_order), , drop = FALSE]
  if (nrow(rows) <= n) return(rows$cell_id)
  rows$cell_id[unique(round(seq(1, nrow(rows), length.out = n)))]
}
chosen <- unlist(lapply(split(l4, l4$bin), choose_cells,
  n = as.integer(cfg$fixture$cells_per_bin)), use.names = FALSE)
fixture_cells <- all_cells[all_cells %in% chosen]
stopifnot(length(fixture_cells) <= cfg$fixture$maximum_cells,
  length(fixture_cells) >= cfg$fixture$n_bins * 20L)
fixture <- l4[match(fixture_cells, l4$cell_id), , drop = FALSE]
fixture$sample_id <- as.character(metadata[[cfg$inputs$sample_field]][
  match(fixture$cell_id, metadata$cell_id)])
stopifnot(!anyNA(fixture$sample_id), length(unique(fixture$bin)) == cfg$fixture$n_bins)
write_tsv(fixture, file.path(dirs[["tables"]], "04_l4_fixture_cells.tsv"))
sample_audit <- as.data.frame(table(bin = fixture$bin,
  sample_id = fixture$sample_id), stringsAsFactors = FALSE)
sample_audit <- sample_audit[sample_audit$Freq > 0L, , drop = FALSE]
names(sample_audit)[3L] <- "cells"
write_tsv(sample_audit, file.path(dirs[["tables"]],
  "05_l4_bin_sample_composition.tsv"))

log_event("Auditing Parse38 overlap on the bounded L4 fixture")
rna_overlap <- wb_session_run_op(session, bundle,
  "scwbgeneprograms::feature_set_overlap_qc", inputs = list(
    matrix = cfg$signal_generation$matrix, cell_ids = fixture_cells,
    min_cells = 20L, namespace = cfg$gene_set_import$namespace,
    min_set_size = 15L, min_detected_fraction = 0,
    deduplicate_sets = TRUE), notes = "Parse38 RNA-count overlap audit.")
session <- rna_overlap$session
sct_overlap <- wb_session_run_op(session, bundle,
  "scwbgeneprograms::feature_set_overlap_qc", inputs = list(
    matrix = cfg$methods$gsea$matrix, cell_ids = fixture_cells,
    min_cells = 20L, namespace = cfg$gene_set_import$namespace,
    min_set_size = as.integer(cfg$methods$gsea$min_set_size),
    min_detected_fraction = 0, deduplicate_sets = TRUE),
  notes = "Parse38 SCT-data overlap audit.")
session <- sct_overlap$session
write_tsv(rna_overlap$run$outputs$overlap_qc$table,
  file.path(dirs[["tables"]], "06_RNA_counts_overlap_qc.tsv"))
write_tsv(sct_overlap$run$outputs$overlap_qc$table,
  file.path(dirs[["tables"]], "07_SCT_data_overlap_qc.tsv"))
stopifnot(sum(rna_overlap$run$outputs$overlap_qc$table$screening_status == "include") == 38L,
  sum(sct_overlap$run$outputs$overlap_qc$table$screening_status == "include") == 38L)

log_event("Generating 38 UCell cytokine signals")
score <- wb_session_run_op(session, bundle, "scwbgeneprograms::ucell",
  inputs = list(matrix = cfg$signal_generation$matrix,
    cell_ids = fixture_cells, min_cells = 20L, features = gene_sets,
    name = "PARSE38_CYTOKINE", max_rank = as.integer(cfg$signal_generation$max_rank),
    ncores = as.integer(cfg$signal_generation$ncores),
    parallel_backend = cfg$signal_generation$parallel_backend,
    chunk_size = as.integer(cfg$signal_generation$chunk_size),
    missing_features_policy = cfg$signal_generation$missing_features_policy),
  notes = "Raw UCell scores for the frozen 38-set Parse cytokine collection.")
session <- score$session
write_tsv(score$run$outputs$signature_summary$table,
  file.path(dirs[["tables"]], "08_ucell_signature_summary.tsv"))

pseudotime_direction <- if (identical(cfg$fit$orientation, "native")) {
  "as_is"
} else {
  cfg$fit$orientation
}
stopifnot(pseudotime_direction %in% c("as_is", "invert"))

common <- list(cell_ids = fixture_cells, min_cells = 20L,
  pseudotime_key = cfg$fit$lineage,
  pseudotime_direction = pseudotime_direction)
signal_ref <- run_ref(score$run, "signature_scores")
pseudotime_ref <- run_ref(fit$run, "pseudotime")

log_event("Running descriptive signal dynamics")
dynamics_inputs <- c(list(signal = signal_ref, pseudotime = pseudotime_ref),
  common, list(n_bins = as.integer(cfg$methods$dynamics$n_bins),
    summary_fun = cfg$methods$dynamics$summary_fun,
    smooth_method = cfg$methods$dynamics$smooth_method,
    bootstrap_n = as.integer(cfg$methods$dynamics$bootstrap_n),
    min_cells_per_bin = as.integer(cfg$methods$dynamics$min_cells_per_bin),
    min_bins_observed = as.integer(cfg$methods$dynamics$min_bins_observed)))
dynamics <- wb_session_run_op(session, bundle, cfg$methods$dynamics$operation,
  inputs = dynamics_inputs,
  notes = "Parse38 UCell dynamics on the common five-percent L4 scale.")
session <- dynamics$session

log_event("Running unblocked permutation-qualified signal intervals")
interval_inputs <- c(list(signal = signal_ref, pseudotime = pseudotime_ref),
  common, list(
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
    merge_gap_fraction = as.numeric(cfg$methods$intervals$merge_gap_fraction)))
intervals <- wb_session_run_op(session, bundle, cfg$methods$intervals$operation,
  inputs = interval_inputs,
  notes = "Unblocked Parse38 signal intervals on five-percent L4 windows.")
session <- intervals$session

log_event("Running unblocked global and local preranked GSEA")
gsea_inputs <- c(list(matrix = cfg$methods$gsea$matrix,
    pseudotime = pseudotime_ref), common, list(
      namespace = cfg$gene_set_import$namespace,
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
      chunk_size = as.integer(cfg$methods$gsea$chunk_size)))
gsea <- wb_session_run_op(session, bundle, cfg$methods$gsea$operation,
  inputs = gsea_inputs,
  notes = "Unblocked Parse38 preranked GSEA on five-percent L4 windows.")

next_index <- save_outputs(dynamics$run, "dynamics", 9L)
next_index <- save_outputs(intervals$run, "intervals", next_index)
next_index <- save_outputs(gsea$run, "gsea", next_index)

dynamics_table <- dynamics$run$outputs$dynamics_table$table
interval_global <- intervals$run$outputs$global_test$table
gsea_global <- gsea$run$outputs$global_enrichment$table
comparison <- merge(dynamics_table,
  interval_global[, c("signal_key", "p_global", "q_global",
    "global_significant", "observed_signal_range"), drop = FALSE],
  by = "signal_key", all = TRUE, sort = FALSE)
gsea_summary <- gsea_global[, c("feature_set", "normalized_enrichment_score",
  "p_value", "q_within_ranking", "q_global_family", "supported"), drop = FALSE]
names(gsea_summary)[names(gsea_summary) == "feature_set"] <- "signal_key"
names(gsea_summary)[-1L] <- paste0("gsea_", names(gsea_summary)[-1L])
comparison <- merge(comparison, gsea_summary, by = "signal_key",
  all = TRUE, sort = FALSE)
write_tsv(comparison, file.path(dirs[["tables"]],
  sprintf("%02d_method_comparison.tsv", next_index)))

fit_summary <- intervals$run$outputs$fit_summary$table
gsea_fit <- gsea$run$outputs$fit_summary$table
stopifnot(identical(as.character(fit_summary$permutation_mode_used), "unblocked"),
  all(as.character(gsea_fit$permutation_mode) == "unblocked"),
  identical(as.numeric(fit_summary$window_width_fraction), 0.05),
  identical(as.numeric(fit_summary$window_step_fraction), 0.025),
  all(as.numeric(gsea_fit$window_start_fraction[is.finite(gsea_fit$window_start_fraction)]) >= 0))

metadata_after <- wb_meta(bundle, cells = all_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$preserved_pseudotime_field))
preserved_pt_after <- as.numeric(metadata_after[[cfg$inputs$preserved_pseudotime_field]])
preserved_pt_after_digest <- digest::digest(serialize(list(cell_id = metadata_after$cell_id,
  value = preserved_pt_after), NULL, version = 3L), algo = "sha256", serialize = FALSE)
stopifnot(identical(preserved_pt_digest, preserved_pt_after_digest))

report <- list(
  schema_version = "bdev.parse38_trajectory_qualification_summary.v1",
  status = "complete_development_qualification",
  run_identity = cfg$run_identity,
  fixture = list(cells = length(fixture_cells), bins = length(unique(fixture$bin)),
    samples = length(unique(fixture$sample_id)), l4_available_cells = nrow(l4)),
  gene_sets = list(expected = 38L,
    rna_included = sum(rna_overlap$run$outputs$overlap_qc$table$screening_status == "include"),
    sct_included = sum(sct_overlap$run$outputs$overlap_qc$table$screening_status == "include")),
  methods = list(
    dynamics = list(run_id = dynamics$run$run_id,
      bins = cfg$methods$dynamics$n_bins),
    intervals = list(run_id = intervals$run$run_id,
      permutation_mode = unique(fit_summary$permutation_mode_used),
      permutations = cfg$methods$intervals$n_permutations,
      supported_intervals = nrow(intervals$run$outputs$signal_intervals$table)),
    gsea = list(run_id = gsea$run$run_id,
      permutation_mode = unique(gsea_fit$permutation_mode),
      permutations = cfg$methods$gsea$n_permutations,
      global_results = nrow(gsea_global),
      local_results = nrow(gsea$run$outputs$window_enrichment$table),
      enrichment_intervals = nrow(gsea$run$outputs$enrichment_intervals$table))),
  common_scale = list(width_fraction = 0.05, overlapping_step = 0.025),
  versions = versions,
  source_bundle_mutated = FALSE,
  execution_clone_only = TRUE,
  outputs_promoted = FALSE,
  preserved_pt_unchanged = TRUE,
  implementation_unlocked = FALSE)
write_json(report, file.path(dirs[["reports"]],
  "01_qualification_summary.json"))
write_json(list(
  schema_version = "bdev.parse38_trajectory_session_receipt.v1",
  run_identity = cfg$run_identity,
  import_runs = imports,
  slingshot_run_id = fit$run$run_id,
  ucell_run_id = score$run$run_id,
  dynamics_run_id = dynamics$run$run_id,
  intervals_run_id = intervals$run$run_id,
  gsea_run_id = gsea$run$run_id,
  staged_outputs = names(gsea$run$outputs),
  promoted = FALSE), file.path(dirs[["manifests"]],
    "02_session_receipt.json"))
log_event("03b bounded qualification complete")
