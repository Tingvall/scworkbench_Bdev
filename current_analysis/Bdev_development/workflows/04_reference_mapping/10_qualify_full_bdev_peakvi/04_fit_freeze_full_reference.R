#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 04_fit_freeze_full_reference.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbLatent)
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
sink(file.path(dirs[["logs"]], "03_full_reference_freeze.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
promoted_ref <- function(promoted, output, kind) {
  selected <- promoted$promoted$status == "ok" & promoted$promoted$output == output
  stopifnot(sum(selected) == 1L)
  list(ref = wb_ref_artifact(promoted$promoted$artifact_id[selected][[1L]], kind),
    artifact_id = promoted$promoted$artifact_id[selected][[1L]])
}

qualification_path <- file.path(dirs[["reports"]],
  "02_candidate_qualification_summary.json")
stopifnot(file.exists(qualification_path))
qualification <- jsonlite::read_json(qualification_path, simplifyVector = TRUE)
stopifnot(identical(qualification$status, "passed"),
  identical(qualification$selected_candidate, "hybrid20k"),
  isFALSE(qualification$external_query_use), isFALSE(qualification$real_immgen_use))

manifest <- system.file("scworkbench", "operation-pack.json", package = "scwbLatent")
stopifnot(nzchar(manifest), identical(sha256(manifest),
  cfg$runtime$scwbLatent$manifest_sha256),
  identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256))
wb_operation_pack_activate(manifest, lib.loc = .libPaths())

source_manifest <- file.path(cfg$inputs$source_bundle$path, "manifest.json")
source_manifest_before <- sha256(source_manifest)
stopifnot(identical(source_manifest_before,
  cfg$inputs$source_bundle$manifest_sha256))
peak_manifest_path <- file.path(dirs[["tables"]], "02_hybrid20k_peak_manifest.tsv")
peak_manifest <- read.delim(peak_manifest_path, check.names = FALSE,
  quote = "\"", na.strings = "NA")
stopifnot(nrow(peak_manifest) == cfg$feature_selection$target_features,
  sum(peak_manifest$hybrid_selection_route == "cluster_10000") ==
    cfg$feature_selection$cluster_features,
  sum(peak_manifest$hybrid_selection_route == "overall_10000") ==
    cfg$feature_selection$overall_features,
  !anyDuplicated(peak_manifest$feature_id),
  all(diff(peak_manifest$source_order) > 0L))

source <- wb_open_bundle(cfg$inputs$source_bundle$path, validate = FALSE,
  read_only = TRUE)
on.exit(try(wb_close_bundle(source), silent = TRUE), add = TRUE)
cells <- as.character(wb_cells(source, "all"))
stopifnot(length(cells) == cfg$inputs$source_bundle$expected_cells)
metadata <- wb_meta(source, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$cluster_field,
  cfg$inputs$sample_field, cfg$inputs$depth_field))
pt_family <- readRDS(cfg$inputs$authenticated_l4$pseudotime_path)
pt_selected <- as.character(pt_family$key) == cfg$inputs$authenticated_l4$key
pt <- stats::setNames(as.numeric(pt_family$value[pt_selected]),
  as.character(pt_family$cell_id[pt_selected]))
weight_family <- readRDS(cfg$inputs$authenticated_l4$weights_path)
weight_selected <- as.character(weight_family$key) == cfg$inputs$authenticated_l4$key
l4_weight <- stats::setNames(as.numeric(weight_family$value[weight_selected]),
  as.character(weight_family$cell_id[weight_selected]))
stopifnot(!anyDuplicated(names(l4_weight)), all(names(l4_weight) %in% cells))
metadata$l4_weight <- 0
metadata$l4_weight[match(names(l4_weight), cells)] <- unname(l4_weight)
metadata$native_l4 <- unname(pt[cells])
metadata$l4_weight_axis_status <- ifelse(cells %in% names(l4_weight),
  "authenticated_weight", "excluded_from_fitted_lineage_family")

log_event("Reading selected 20,000-peak matrix for all", length(cells), "Bdev cells")
counts <- wb_slice(source, cfg$inputs$count_matrix,
  features = as.character(peak_manifest$feature_id), cells = cells)
stopifnot(inherits(counts, "sparseMatrix"),
  identical(rownames(counts), as.character(peak_manifest$feature_id)),
  identical(colnames(counts), cells))
wb_close_bundle(source); source <- NULL

if (dir.exists(cfg$storage$frozen_bundle))
  stop("Frozen bundle already exists; preserve it and use a new run identity")
bundle <- wb_init_bundle(cfg$storage$frozen_bundle, overwrite = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
wb_define_cells(bundle, cells)
wb_attach_matrix(bundle, "Bdev_ATAC_hybrid20k_counts", counts,
  layer = "counts", feature_type = "peak", cell_ids = cells,
  description = "Frozen Bdev-only 50/50 cluster and overall-variable 20k PeakVI axis",
  source = "active_03c_Bdev_ATAC_counts_Bdev_only_frozen_feature_recipe",
  created_by = "scworkbench_Bdev development")
wb_attach_metadata(bundle, metadata, cell_ids = cells,
  source = "active_03c_Bdev_metadata_plus_authenticated_native_L4",
  created_by = "scworkbench_Bdev development")
rm(counts); invisible(gc())

log_event("Fitting frozen PeakVI model on all Bdev cells")
fit_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_frozen_full_bdev_fit"))
fit <- wb_session_run_op(fit_session, bundle, "scwblatent::peakvi_fit",
  inputs = list(matrix = "Bdev_ATAC_hybrid20k_counts",
    n_latent = as.integer(cfg$peakvi$n_latent),
    n_layers = as.integer(cfg$peakvi$n_layers),
    n_hidden = as.integer(cfg$peakvi$n_hidden),
    dropout = as.numeric(cfg$peakvi$dropout),
    max_epochs = as.integer(cfg$peakvi$max_epochs),
    lr = as.numeric(cfg$peakvi$learning_rate),
    early_stopping = isTRUE(cfg$peakvi$early_stopping),
    seed = as.integer(cfg$peakvi$seed),
    chunk_size = as.integer(cfg$peakvi$chunk_size), batch_key = NULL,
    use_gpu = FALSE, cell_set = "all", min_cells = 1L),
  notes = "Frozen full-Bdev hybrid20k PeakVI reference; external queries absent")
reference_space_id <- as.character(fit$run$outputs$space$space_id %||%
  fit$run$outputs$space$allocated_id)[[1L]]
fit_summary <- fit$run$outputs$fit_summary$table
stopifnot(nrow(fit_summary) == 1L,
  as.integer(fit_summary$n_cells) == cfg$inputs$source_bundle$expected_cells,
  as.integer(fit_summary$n_features) == cfg$feature_selection$target_features,
  as.integer(fit_summary$n_latent) == cfg$peakvi$n_latent)
fit$session <- wb_session_resume(bundle, fit$session$durable$workspace)
promoted <- wb_session_promote_current(bundle, fit$session,
  outputs = c("model", "fit_summary"),
  name_prefix = "frozen_full_bdev_hybrid20k_peakvi",
  description = "Frozen Bdev-only whole-reference PeakVI model",
  created_by = "scworkbench_Bdev development", fail_fast = TRUE)
model_record <- promoted_ref(promoted, "model", "model")
summary_record <- promoted_ref(promoted, "fit_summary", "obs_table")
promoted$session <- wb_session_discard(promoted$session)
wb_close_bundle(bundle); bundle <- NULL

bundle_manifest <- file.path(cfg$storage$frozen_bundle, "manifest.json")
stopifnot(file.exists(bundle_manifest),
  identical(sha256(source_manifest), source_manifest_before))
feature_sha <- sha256(peak_manifest_path)
freeze <- list(
  schema_version = "bdev.full_peakvi_frozen_reference.v1",
  status = "frozen_development_reference",
  run_identity = cfg$run_identity,
  model = list(engine = "peakvi", candidate = "hybrid20k",
    artifact_id = model_record$artifact_id,
    fit_summary_artifact_id = summary_record$artifact_id,
    reference_space_id = reference_space_id,
    cells = as.integer(fit_summary$n_cells), features = as.integer(fit_summary$n_features),
    latent_dimensions = as.integer(fit_summary$n_latent),
    latent_sha256 = as.character(fit_summary$latent_sha256)),
  features = list(path = peak_manifest_path, sha256 = feature_sha,
    cluster_selected = cfg$feature_selection$cluster_features,
    overall_selected = cfg$feature_selection$overall_features,
    ordering = "source_ATAC_counts_feature_order"),
  trajectory_mapping = qualification$frozen_l4_rule,
  runtime = list(sif = cfg$runtime$sif, sif_sha256 = cfg$runtime$sif_sha256,
    scworkbench_commit = cfg$runtime$scworkbench$source_commit,
    scwbLatent_commit = cfg$runtime$scwbLatent$source_commit,
    scwbLatent_manifest_sha256 = cfg$runtime$scwbLatent$manifest_sha256),
  frozen_bundle = list(path = cfg$storage$frozen_bundle,
    manifest_sha256 = sha256(bundle_manifest)),
  source_bundle_manifest_sha256 = source_manifest_before,
  external_query_use_during_design = FALSE,
  real_immgen_use = FALSE,
  implementation_unlocked = FALSE,
  sif_built = FALSE)
write_json(freeze, file.path(dirs[["reports"]], "03_frozen_reference_manifest.json"))
write_json(list(schema_version = "bdev.full_peakvi_freeze_execution.v1",
  status = "complete", run_identity = cfg$run_identity,
  config_path = cfg_path, config_sha256 = sha256(cfg_path),
  qualification_sha256 = sha256(qualification_path),
  frozen_reference_manifest_sha256 = sha256(file.path(dirs[["reports"]],
    "03_frozen_reference_manifest.json")),
  source_manifest_unchanged = TRUE, external_query_use = FALSE,
  real_immgen_use = FALSE),
  file.path(dirs[["manifests"]], "03_freeze_execution_receipt.json"))
log_event("Frozen full-Bdev hybrid20k PeakVI reference complete")
