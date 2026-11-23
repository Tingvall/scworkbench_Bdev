#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 03_map_and_plot_fixed_umap.R --config <config.yaml>", call. = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbLatent)
  library(scwbRepresentation)
  library(scwbReferenceMap)
  library(BPCells)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(slingshot)
  library(yaml)
  library(jsonlite)
  library(digest)
})

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
atomic_save_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), ".tmp-"),
    tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(value, temporary, version = 3L)
  stopifnot(file.rename(temporary, path))
  invisible(path)
}
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}
space_matrix <- function(x) {
  if (is.data.frame(x) && "cell_id" %in% names(x)) {
    ids <- as.character(x$cell_id)
    x <- as.matrix(x[, setdiff(names(x), "cell_id"), drop = FALSE])
    rownames(x) <- ids
  } else x <- as.matrix(x)
  x
}
promoted_id <- function(x, output) {
  take <- x$promoted$status == "ok" & x$promoted$output == output
  stopifnot(sum(take) == 1L)
  as.character(x$promoted$artifact_id[take][[1L]])
}
promoted_ref <- function(x, output, kind) {
  wb_ref_artifact(promoted_id(x, output), kind)
}
merge_lists <- function(left, right) {
  result <- left
  for (name in names(right)) {
    if (is.list(right[[name]]) && is.list(result[[name]]) &&
        !is.null(names(right[[name]])) && !is.null(names(result[[name]])))
      result[[name]] <- merge_lists(result[[name]], right[[name]])
    else result[[name]] <- right[[name]]
  }
  result
}
read_config <- function(path) {
  raw <- yaml::read_yaml(path)
  if (is.null(raw$extends)) return(raw)
  parent <- normalizePath(file.path(dirname(path), raw$extends), mustWork = TRUE)
  raw$extends <- NULL
  merge_lists(read_config(parent), raw)
}

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- read_config(cfg_path)
root <- cfg$storage$results_root
invisible(lapply(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(root, "04_logs", "02_map_and_plot_fixed_umap.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(
  sha256(cfg$runtime$base_sif) == cfg$runtime$base_sif_sha256,
  sha256(file.path(cfg$inputs$frozen_model_bundle$path, "manifest.json")) ==
    cfg$inputs$frozen_model_bundle$manifest_sha256,
  sha256(file.path(cfg$inputs$fixed_bdev_bundle$path, "manifest.json")) ==
    cfg$inputs$fixed_bdev_bundle$manifest_sha256,
  sha256(cfg$inputs$authenticated_l4$model_path) ==
    cfg$inputs$authenticated_l4$model_sha256,
  sha256(cfg$inputs$authenticated_l4$pseudotime_path) ==
    cfg$inputs$authenticated_l4$pseudotime_sha256,
  as.character(packageVersion("scworkbench")) == "0.1.2",
  as.character(packageVersion("scwbLatent")) == "0.1.3",
  as.character(packageVersion("scwbRepresentation")) == "0.1.1",
  as.character(packageVersion("scwbReferenceMap")) == "0.1.0")

for (package in c("scwbLatent", "scwbRepresentation", "scwbReferenceMap")) {
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  namespace <- switch(package, scwbLatent = "scwblatent",
    scwbRepresentation = "scwbrepresentation",
    scwbReferenceMap = "scwbreferencemap")
  try(wb_operation_pack_unload(namespace), silent = TRUE)
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
}
required <- c("scwblatent::peakvi_project",
  "scwblatent::model_bound_latent_union",
  "scwbrepresentation::knn_neighbors",
  "scwbreferencemap::reference_map_calibrate",
  "scwbreferencemap::reference_map_transfer")
stopifnot(all(required %in% wb_list_ops()))

stopifnot(sha256(cfg$inputs$prepared_cells$path) ==
  cfg$inputs$prepared_cells$sha256,
  sha256(cfg$inputs$parent_quantification_audit$path) ==
    cfg$inputs$parent_quantification_audit$sha256)
cells <- as.data.table(read.delim(gzfile(cfg$inputs$prepared_cells$path),
  check.names = FALSE, quote = "\"",
  na.strings = "NA"))
query_matrix <- BPCells::open_matrix_dir(cfg$storage$quantified_matrix)
stopifnot(nrow(cells) == cfg$scope$expected_query_cells,
  ncol(query_matrix) == cfg$scope$expected_query_cells,
  identical(as.character(colnames(query_matrix)), as.character(cells$cell_id)))
slots <- sprintf("external_ko_query_%05d", seq_len(nrow(cells)))
slot_map <- copy(cells)
slot_map[, slot_id := slots]
setcolorder(slot_map, c("slot_id", setdiff(names(slot_map), "slot_id")))
fwrite(slot_map, file.path(root, "01_tables", "04_full_query_slot_map.tsv.gz"),
  sep = "\t", quote = TRUE, na = "NA")

log_event("Reading frozen model and reference space")
source <- wb_open_bundle(cfg$inputs$frozen_model_bundle$path,
  validate = FALSE, read_only = TRUE)
source_cells <- as.character(wb_cells(source, "all"))
reference_cells <- source_cells[!grepl("^external_query_[0-9]+$", source_cells)]
stopifnot(length(reference_cells) ==
  cfg$inputs$frozen_model_bundle$expected_reference_cells)
reference_space <- space_matrix(wb_space(source,
  cfg$inputs$frozen_model_bundle$reference_space, cells = reference_cells))
reference_meta <- wb_meta(source, cells = reference_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "orig.ident", "native_l4"))
source_model <- wb_art_load(source,
  cfg$inputs$frozen_model_bundle$model_artifact_id, expected_kind = "model_ref")
source_model_payload <- file.path(source_model$model_tree_root,
  source_model$model_primary_path)
stopifnot(sha256(source_model_payload) ==
  cfg$inputs$frozen_model_bundle$model_payload_sha256,
  identical(rownames(reference_space), reference_cells),
  identical(as.character(reference_meta$cell_id), reference_cells))

wb_close_bundle(source)
source <- NULL
all_cells <- c(source_cells, slots)
query_meta <- data.frame(
  cell_id = slots,
  orig.ident = as.character(slot_map$library_id),
  native_l4 = rep(NA_real_, nrow(slot_map)),
  cell_role = rep("flt3lko_il7ko_query", nrow(slot_map)),
  stringsAsFactors = FALSE)
if (dir.exists(cfg$storage$mapping_bundle)) {
  log_event("Resuming verified frozen-bundle clone")
  bundle <- wb_open_bundle(cfg$storage$mapping_bundle, validate = FALSE,
    read_only = FALSE)
} else {
  log_event("Reflink-cloning frozen bundle with its durable model artifact")
  dir.create(dirname(cfg$storage$mapping_bundle), recursive = TRUE,
    showWarnings = FALSE)
  status <- system2("cp", c("--recursive", "--reflink=auto", "--",
    cfg$inputs$frozen_model_bundle$path, cfg$storage$mapping_bundle))
  stopifnot(status == 0L)
  bundle <- wb_open_bundle(cfg$storage$mapping_bundle, validate = FALSE,
    read_only = FALSE)
}
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
current_cells <- as.character(wb_cells(bundle, "all"))
if (identical(current_cells, source_cells)) {
  log_event("Ingesting all query cells into the frozen-bundle clone")
  wb_ingest_apply(bundle, set_name = "flt3lko_il7ko_full_queries",
    cell_ids = slots, metadata = query_meta,
    source = "checksummed_full_KO_query_metadata",
    created_by = "Bdev development capability 05 full descriptive mapping",
    apply_mode = "atomic", emit_report_artifact = TRUE)
  current_cells <- as.character(wb_cells(bundle, "all"))
}
stopifnot(identical(current_cells, all_cells))
query_sparse <- as(query_matrix, "dgCMatrix")
colnames(query_sparse) <- slots
query_matrix_id <- "hvp_depth20k_flt3lko_il7ko_full_counts"
if (!(query_matrix_id %in% wb_matrices(bundle)$matrix_id))
  wb_attach_matrix(bundle, query_matrix_id, query_sparse, layer = "counts",
    feature_type = "peak", cell_ids = slots,
    description = "All 23,252 QC-selected FLT3LKO/IL7KO cells on exact frozen HVP-depth20k axis",
    source = "reciprocal_mm10_mm39_liftover_and_read_only_fragment_quantification",
    created_by = "Bdev development capability 05 full descriptive mapping")
rm(query_sparse); invisible(gc())
model <- wb_art_load(bundle,
  cfg$inputs$frozen_model_bundle$model_artifact_id, expected_kind = "model_ref")
model_ref <- wb_ref_artifact(cfg$inputs$frozen_model_bundle$model_artifact_id,
  "model")
model_payload <- file.path(model$model_tree_root, model$model_primary_path)
stopifnot(sha256(model_payload) ==
  cfg$inputs$frozen_model_bundle$model_payload_sha256)

query_space_id <- "frozen_hvp_depth20k_flt3lko_il7ko_full"
artifacts <- wb_art_list(bundle)
projection_take <- artifacts$source_op == "scwblatent::peakvi_project" &
  artifacts$type == "obs_table" & artifacts$lifecycle_state == "current" &
  grepl("^flt3lko_il7ko_full_projection - ", artifacts$name)
if (query_space_id %in% wb_spaces(bundle)$space_id &&
    sum(projection_take) == 1L) {
  log_event("Recovering already-published 23,252-cell projection")
  projection_id <- as.character(artifacts$artifact_id[projection_take][[1L]])
  projection_summary <- wb_art_load(bundle, projection_id,
    expected_kind = "obs_table")$table
  stopifnot(projection_summary$transformation[[1L]] == "transfer_without_refit",
    projection_summary$n_query_cells[[1L]] == length(slots))
  projection_ref <- wb_ref_artifact(projection_id, "obs_table")
} else {
  stopifnot(!(query_space_id %in% wb_spaces(bundle)$space_id),
    sum(projection_take) == 0L)
  log_event("Projecting", length(slots), "query cells without model refit")
  project_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_peakvi_project_retry2"))
  project <- tryCatch(wb_session_run_op(project_session, bundle,
      "scwblatent::peakvi_project", inputs = list(
        model = model_ref, query_matrix = query_matrix_id,
        space_id_out = query_space_id,
        chunk_size = as.integer(cfg$mapping$peakvi_chunk_size), batch_key = NULL,
        use_gpu = isTRUE(cfg$mapping$use_gpu),
        pickle_execution_trust = cfg$mapping$pickle_execution_trust,
        cell_set = "set:flt3lko_il7ko_full_queries", min_cells = 1L),
      notes = "Full descriptive FLT3LKO/IL7KO projection through frozen HVP-depth20k model"),
    error = function(error) {
      diagnostic <- list(status = "failed_before_promotion",
        message = conditionMessage(error), code = error$code,
        phase = error$phase, expected = error$expected,
        observed = error$observed)
      write_json(diagnostic, file.path(root, "03_reports",
        "02_peakvi_projection_failure_diagnostic.json"))
      stop(error)
    })
  stopifnot(project$run$outputs$projection_summary$table$transformation[[1L]] ==
    "transfer_without_refit")
  project$session <- wb_session_resume(bundle, project$session$durable$workspace)
  project_promoted <- wb_session_promote_current(bundle, project$session,
    outputs = "projection_summary", name_prefix = "flt3lko_il7ko_full_projection",
    description = "All-query projection without model refit",
    created_by = "Bdev development capability 05 full descriptive mapping",
    fail_fast = TRUE)
  projection_ref <- promoted_ref(project_promoted, "projection_summary", "obs_table")
  project_promoted$session <- wb_session_discard(project_promoted$session)
}

joint_space_id <- "hvp_depth20k_bdev_flt3lko_il7ko_full_union"
union_cell_set <- "set:flt3lko_il7ko_full_union"
if (!(union_cell_set %in% wb_cell_sets(bundle))) {
  union_set_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_union_cell_set"))
  union_set <- wb_session_run_op(union_set_session, bundle,
    "cell_set_derive", inputs = list(
      cell_ids = c(reference_cells, slots), target_cell_set = union_cell_set,
      role = "analysis_universe", status = "accepted",
      source = "exact_frozen_reference_and_full_query_union",
      created_by = "Bdev development capability 05 full descriptive mapping",
      description = "Exact cells in the authenticated full mapping union",
      min_cells = 1L),
    notes = "Provenance-bound exact full reference/query union cell set")
  union_set$session <- wb_session_resume(bundle,
    union_set$session$durable$workspace)
  union_set_promoted <- wb_session_promote_current(bundle, union_set$session,
    outputs = c("manifest", "rule_table"),
    name_prefix = "flt3lko_il7ko_full_union_cell_set",
    description = "Exact full reference/query union cell set",
    created_by = "Bdev development capability 05 full descriptive mapping",
    fail_fast = TRUE)
  union_set_promoted$session <- wb_session_discard(
    union_set_promoted$session)
}
stopifnot(identical(as.character(wb_cells(bundle, union_cell_set)),
  c(reference_cells, slots)))
if (joint_space_id %in% wb_spaces(bundle)$space_id) {
  log_event("Recovering already-published frozen-model reference/query union")
} else {
  log_event("Authenticating frozen-model reference/query union")
  union_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_latent_union_retry2"))
  union <- wb_session_run_op(union_session, bundle,
    "scwblatent::model_bound_latent_union", inputs = list(
      reference_space = cfg$inputs$frozen_model_bundle$reference_space,
      query_space = query_space_id, model = model_ref,
      projection_summary = projection_ref, space_id_out = joint_space_id,
      chunk_size = 1024L, cell_set = union_cell_set, min_cells = 1L),
    notes = "Authenticated frozen-model Bdev/full-KO latent union")
  union$session <- wb_session_resume(bundle, union$session$durable$workspace)
  union_promoted <- wb_session_promote_current(bundle, union$session,
    outputs = "union_summary", name_prefix = "flt3lko_il7ko_full_union",
    description = "Authenticated frozen-model full-query latent union",
    created_by = "Bdev development capability 05 full descriptive mapping",
    fail_fast = TRUE)
  union_promoted$session <- wb_session_discard(union_promoted$session)
}

reference_l4 <- wb_prim_cell_vector_family(rep("L4", length(reference_cells)),
  reference_cells, as.numeric(reference_meta$native_l4), value_type = "continuous")
authorizer <- Matrix::sparseMatrix(i = rep.int(1L, length(reference_cells)),
  j = seq_along(reference_cells), x = 1,
  dims = c(1L, length(reference_cells)),
  dimnames = list("mapping_authorizer_peak", reference_cells))
if (!("hvp_depth20k_full_mapping_authorizer" %in%
    wb_matrices(bundle)$matrix_id))
  wb_attach_matrix(bundle, "hvp_depth20k_full_mapping_authorizer", authorizer,
    layer = "counts", feature_type = "peak", cell_ids = reference_cells,
    description = "One-feature sparse authorizer for reference mapping",
    source = "frozen_reference_cell_axis",
    created_by = "Bdev development capability 05 full descriptive mapping")

calibration_prefix <- if (is.null(cfg$mapping$calibration_name_prefix))
  "flt3lko_il7ko_full_calibration" else
    as.character(cfg$mapping$calibration_name_prefix)
current_calibration <- artifacts$source_op ==
  "scwbreferencemap::reference_map_calibrate" &
  artifacts$lifecycle_state == "current" &
  startsWith(artifacts$name, paste0(calibration_prefix, " - "))
calibration_model_take <- current_calibration & artifacts$type == "model_ref"
threshold_take <- current_calibration & artifacts$type == "obs_table" &
  grepl("\\(thresholds\\)$", artifacts$name)
metrics_take <- current_calibration & artifacts$type == "obs_table" &
  grepl("\\(calibration_metrics\\)$", artifacts$name)
if (sum(calibration_model_take) == 1L && sum(threshold_take) == 1L &&
    sum(metrics_take) == 1L) {
  log_event("Recovering already-published Bdev reference calibration")
  calibration_ref <- wb_ref_artifact(
    artifacts$artifact_id[calibration_model_take][[1L]], "model")
  thresholds <- wb_art_load(bundle,
    artifacts$artifact_id[threshold_take][[1L]],
    expected_kind = "obs_table")$table
  cal_metrics <- wb_art_load(bundle,
    artifacts$artifact_id[metrics_take][[1L]],
    expected_kind = "obs_table")$table
} else {
  log_event("Calibrating Bdev reference support")
  cal_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_calibrate_checkpointed_",
      format(Sys.time(), "%Y%m%dT%H%M%S")))
  cal <- wb_session_run_op(cal_session, bundle,
    "scwbreferencemap::reference_map_calibrate", inputs = list(
      reference_matrix = "hvp_depth20k_full_mapping_authorizer",
      reference_space = cfg$inputs$frozen_model_bundle$reference_space,
      reference_continuous = reference_l4,
      holdout_group_field = "orig.ident",
      k_values = as.integer(unlist(cfg$mapping$k_values)),
      primary_k = as.integer(cfg$mapping$primary_k),
      support_quantiles = as.list(c(0.95, 0.99)),
      primary_continuous_key = "L4",
      space_chunk_size = as.integer(cfg$mapping$space_chunk_size),
      query_chunk_size = as.integer(cfg$mapping$query_chunk_size)),
    notes = "Bdev-only support and native-L4 calibration for full KO mapping")
  thresholds <- cal$run$outputs$thresholds$table
  cal_metrics <- cal$run$outputs$calibration_metrics$table
  cal$session <- wb_session_resume(bundle, cal$session$durable$workspace)
  cal_promoted <- wb_session_promote_current(bundle, cal$session,
    outputs = names(cal$run$outputs),
    name_prefix = calibration_prefix,
    description = "Bdev-only calibration for full KO mapping",
    created_by = "Bdev development capability 05 full descriptive mapping",
    fail_fast = TRUE)
  calibration_ref <- promoted_ref(cal_promoted, "calibration", "model")
  cal_promoted$session <- wb_session_discard(cal_promoted$session)
}
supported_threshold <- thresholds$threshold[thresholds$role == "supported_upper"][[1L]]
out_threshold <- thresholds$threshold[thresholds$role == "borderline_upper"][[1L]]
if (!is.null(cfg$inputs$baseline_multi_k_calibration)) {
  baseline_cfg <- cfg$inputs$baseline_multi_k_calibration
  stopifnot(sha256(file.path(baseline_cfg$bundle_path, "manifest.json")) ==
    baseline_cfg$bundle_manifest_sha256)
  baseline_bundle <- wb_open_bundle(baseline_cfg$bundle_path,
    validate = FALSE, read_only = TRUE)
  baseline_thresholds <- wb_art_load(baseline_bundle,
    baseline_cfg$threshold_artifact_id,
    expected_kind = "obs_table")$table
  wb_close_bundle(baseline_bundle)
  threshold_comparison <- merge(
    baseline_thresholds[, c("role", "threshold")],
    thresholds[, c("role", "threshold")], by = "role",
    suffixes = c("_multi_k", "_primary_k_only"), sort = FALSE)
  stopifnot(nrow(threshold_comparison) == nrow(thresholds),
    isTRUE(all.equal(threshold_comparison$threshold_multi_k,
      threshold_comparison$threshold_primary_k_only,
      tolerance = 1e-12, check.attributes = FALSE)))
  fwrite(threshold_comparison, file.path(root, "01_tables",
    "10_primary_k_calibration_threshold_comparison.tsv"), sep = "\t",
    quote = TRUE, na = "NA")
}

batch_size <- as.integer(cfg$mapping$transfer_batch_size)
stopifnot(length(batch_size) == 1L, is.finite(batch_size), batch_size >= 1L)
batches <- split(slots, ceiling(seq_along(slots) / batch_size))
neighbor_batches <- vector("list", length(batches))
prediction_batches <- vector("list", length(batches))
continuous_batches <- vector("list", length(batches))
query_coordinates <- space_matrix(wb_space(bundle, query_space_id,
  cells = slots))
stopifnot(identical(rownames(query_coordinates), slots))
available_spaces <- wb_spaces(bundle)$space_id
checkpoint_root <- file.path(root, "05_objects",
  "01_bounded_transfer_checkpoints")
dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)
batch_contract <- digest::digest(list(
  schema_version = "bdev.bounded_reference_transfer_checkpoint.v1",
  run_identity = cfg$run_identity,
  joint_space_id = joint_space_id,
  query_space_id = query_space_id,
  reference_space_id = cfg$inputs$frozen_model_bundle$reference_space,
  model_payload_sha256 = cfg$inputs$frozen_model_bundle$model_payload_sha256,
  k_values = as.integer(unlist(cfg$mapping$k_values)),
  primary_k = as.integer(cfg$mapping$primary_k),
  maximum_k = as.integer(cfg$mapping$maximum_k),
  supported_threshold = as.numeric(supported_threshold),
  out_of_reference_threshold = as.numeric(out_threshold),
  batch_size = batch_size), algo = "sha256")
validate_checkpoint <- function(value, batch_index, batch_cells) {
  stopifnot(is.list(value),
    identical(value$schema_version,
      "bdev.bounded_reference_transfer_checkpoint.v1"),
    identical(value$batch_contract_sha256, batch_contract),
    identical(as.integer(value$batch_index), as.integer(batch_index)),
    identical(as.character(value$query_cells), as.character(batch_cells)),
    is.data.frame(value$neighbors), is.data.frame(value$prediction),
    is.data.frame(value$continuous),
    nrow(value$neighbors) == length(batch_cells) *
      as.integer(cfg$mapping$maximum_k),
    identical(unique(as.character(value$neighbors$query_cell_id)),
      as.character(batch_cells)),
    identical(as.character(value$prediction$query_cell_id),
      as.character(batch_cells)),
    all(value$continuous$query_cell_id %in% batch_cells))
  invisible(TRUE)
}
log_event("Finding exact neighbors and transferring support/native L4 in",
  length(batches), "bounded query batches")
for (batch_index in seq_along(batches)) {
  batch_cells <- batches[[batch_index]]
  checkpoint_path <- file.path(checkpoint_root,
    sprintf("batch_%03d.rds", batch_index))
  if (file.exists(checkpoint_path)) {
    checkpoint <- readRDS(checkpoint_path)
    validate_checkpoint(checkpoint, batch_index, batch_cells)
    neighbor_batches[[batch_index]] <- as.data.table(checkpoint$neighbors)
    prediction_batches[[batch_index]] <- as.data.table(checkpoint$prediction)
    continuous_batches[[batch_index]] <- as.data.table(checkpoint$continuous)
    log_event("Recovered validated batch", batch_index, "of", length(batches),
      "from", checkpoint_path)
    next
  }
  log_event("Mapping batch", batch_index, "of", length(batches),
    "with", length(batch_cells), "cells")
  batch_space_id <- sprintf("%s_batch_%03d", query_space_id, batch_index)
  batch_coordinates <- query_coordinates[batch_cells, , drop = FALSE]
  if (!(batch_space_id %in% available_spaces)) {
    batch_coordinates_frame <- data.frame(cell_id = batch_cells,
      batch_coordinates, check.names = FALSE, row.names = NULL)
    wb_attach_space(bundle, batch_space_id, batch_coordinates_frame,
      cell_ids = batch_cells,
      source = paste0("provenance_recorded_subset_of:", query_space_id),
      created_by = "Bdev development capability 05 bounded transfer",
      emit_report_artifact = FALSE)
    available_spaces <- c(available_spaces, batch_space_id)
  } else {
    existing_batch_coordinates <- space_matrix(wb_space(bundle,
      batch_space_id, cells = batch_cells))
    stopifnot(identical(rownames(existing_batch_coordinates), batch_cells),
      identical(unname(existing_batch_coordinates),
        unname(batch_coordinates)))
  }
  batch_session <- wb_session_create(bundle, cfg$storage$sessions,
    sprintf("%s_bounded_transfer_%03d_%s_%s", cfg$run_identity,
      batch_index, format(Sys.time(), "%Y%m%dT%H%M%S"), Sys.getpid()))
  batch_knn <- wb_session_run_op(batch_session, bundle,
    "scwbrepresentation::knn_neighbors", inputs = list(
      space = joint_space_id, query_cells = batch_cells,
      candidate_cells = reference_cells,
      k = as.integer(cfg$mapping$maximum_k), include_query = FALSE,
      space_chunk_size = as.integer(cfg$mapping$space_chunk_size)),
    notes = sprintf("Bounded full-query neighbor batch %d/%d",
      batch_index, length(batches)))
  batch_transfer <- wb_session_run_op(batch_knn$session, bundle,
    "scwbreferencemap::reference_map_transfer", inputs = list(
      in_space = batch_space_id,
      reference_space = cfg$inputs$frozen_model_bundle$reference_space,
      neighbors = batch_knn$run$outputs$neighbors,
      calibration = calibration_ref,
      support_output_key = "flt3lko_il7ko_full_support",
      annotation_output_key = "unused_no_categorical_transfer",
      continuous_output_keys = "mapped_native_L4",
      supported_threshold = as.numeric(supported_threshold),
      out_of_reference_threshold = as.numeric(out_threshold),
      cell_ids = batch_cells, min_cells = 1L,
      space_chunk_size = as.integer(cfg$mapping$space_chunk_size)),
    notes = sprintf("Bounded support/native-L4 transfer batch %d/%d",
      batch_index, length(batches)))
  neighbor_batches[[batch_index]] <- as.data.table(
    batch_knn$run$outputs$neighbors$table)
  prediction_batches[[batch_index]] <- as.data.table(
    batch_transfer$run$outputs$per_query_summary$table)
  continuous_batches[[batch_index]] <- as.data.table(
    batch_transfer$run$outputs$continuous_estimates$table)
  checkpoint <- list(
    schema_version = "bdev.bounded_reference_transfer_checkpoint.v1",
    batch_contract_sha256 = batch_contract,
    batch_index = as.integer(batch_index),
    query_cells = as.character(batch_cells),
    neighbors = as.data.frame(neighbor_batches[[batch_index]]),
    prediction = as.data.frame(prediction_batches[[batch_index]]),
    continuous = as.data.frame(continuous_batches[[batch_index]]),
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  validate_checkpoint(checkpoint, batch_index, batch_cells)
  atomic_save_rds(checkpoint, checkpoint_path)
  log_event("Checkpointed batch", batch_index, "of", length(batches),
    "sha256", sha256(checkpoint_path))
  batch_transfer$session <- wb_session_discard(batch_transfer$session)
}
neighbors <- rbindlist(neighbor_batches, use.names = TRUE)
prediction <- rbindlist(prediction_batches, use.names = TRUE)
continuous <- rbindlist(continuous_batches, use.names = TRUE)
stopifnot(identical(unique(neighbors$query_cell_id), slots),
  identical(prediction$query_cell_id, slots))
primary <- continuous[k == cfg$mapping$primary_k & reference_key == "L4",
  .(query_cell_id, weighted_q25, raw_native_l4 = raw_estimate,
    weighted_q75, iqr, cross_k_sensitivity)]
outcomes <- merge(slot_map, prediction, by.x = "slot_id",
  by.y = "query_cell_id", sort = FALSE)
outcomes <- merge(outcomes, primary, by.x = "slot_id",
  by.y = "query_cell_id", sort = FALSE)
outcomes <- outcomes[match(slots, slot_id)]
outcomes[, accepted_native_l4 := fifelse(support_call == "out_of_reference",
  NA_real_, raw_native_l4)]
stopifnot(nrow(outcomes) == cfg$scope$expected_query_cells,
  all(is.na(outcomes$accepted_native_l4[
    outcomes$support_call == "out_of_reference"])))
fwrite(outcomes, file.path(root, "01_tables",
  "05_full_support_aware_native_l4.tsv.gz"), sep = "\t", quote = TRUE, na = "NA")
fwrite(thresholds, file.path(root, "01_tables", "06_support_thresholds.tsv"),
  sep = "\t", quote = TRUE, na = "NA")
fwrite(cal_metrics, file.path(root, "01_tables", "07_calibration_metrics.tsv"),
  sep = "\t", quote = TRUE, na = "NA")

log_event("Computing placements on unchanged ATAC UMAP")
fixed <- wb_open_bundle(cfg$inputs$fixed_bdev_bundle$path,
  validate = FALSE, read_only = TRUE)
atac_umap <- space_matrix(wb_space(fixed, "umap.atac", cells = reference_cells))
wb_close_bundle(fixed)
umap <- data.table(cell_id = rownames(atac_umap),
  UMAP_1 = atac_umap[, 1L], UMAP_2 = atac_umap[, 2L])
primary_neighbors <- neighbors[neighbor_rank <= cfg$mapping$primary_k]
kth <- primary_neighbors[neighbor_rank == cfg$mapping$primary_k,
  .(query_cell_id, kth_distance = distance)]
primary_neighbors <- merge(primary_neighbors, kth, by = "query_cell_id",
  sort = FALSE)
primary_neighbors <- merge(primary_neighbors, umap, by = "cell_id", sort = FALSE)
primary_neighbors[, weight := fifelse(kth_distance == 0,
  as.numeric(distance == 0), exp(-(distance / kth_distance)^2))]
placement <- primary_neighbors[, .(
  UMAP_1 = sum(weight * UMAP_1) / sum(weight),
  UMAP_2 = sum(weight * UMAP_2) / sum(weight),
  neighbor_spread = sqrt(sum(weight * ((UMAP_1 - sum(weight * UMAP_1) /
    sum(weight))^2 + (UMAP_2 - sum(weight * UMAP_2) / sum(weight))^2)) /
    sum(weight)),
  nearest_distance = min(distance), kth_distance = max(kth_distance)),
  by = query_cell_id]
placement <- merge(outcomes, placement, by.x = "slot_id",
  by.y = "query_cell_id", sort = FALSE)
placement <- placement[match(slots, slot_id)]
labels <- c(FLT3LKO_BM_WT = "Flt3l WT", FLT3LKO_BM_KO = "Flt3l KO",
  IL7KO_BM_WT = "Il7 WT", IL7KO_BM_KO = "Il7 KO")
placement[, library_label := factor(labels[library_id],
  levels = c("Flt3l WT", "Flt3l KO", "Il7 WT", "Il7 KO"))]
fwrite(placement, file.path(root, "01_tables", "08_full_fixed_umap_placements.tsv.gz"),
  sep = "\t", quote = TRUE, na = "NA")

family <- readRDS(cfg$inputs$authenticated_l4$pseudotime_path)
take <- as.character(family$key) == "L4" & is.finite(as.numeric(family$value)) &
  family$cell_id %in% reference_cells
trajectory <- data.table(cell_id = as.character(family$cell_id[take]),
  native_l4 = as.numeric(family$value[take]))
trajectory <- merge(trajectory, umap, by = "cell_id", sort = FALSE)
slingshot_model <- readRDS(cfg$inputs$authenticated_l4$model_path)
curve_value <- slingCurves(slingshot_model)[["Lineage4"]]
curve_xy <- curve_value$s[curve_value$ord, , drop = FALSE]
curve <- data.frame(UMAP_1 = curve_xy[, 1L], UMAP_2 = curve_xy[, 2L])
l4_limits <- range(c(trajectory$native_l4, placement$accepted_native_l4),
  finite = TRUE)

p <- ggplot() +
  geom_point(data = umap, aes(UMAP_1, UMAP_2), color = "#D4D4D4",
    size = 0.12, alpha = 0.28) +
  geom_path(data = curve, aes(UMAP_1, UMAP_2), color = "#111111",
    linewidth = 0.55, arrow = grid::arrow(type = "closed",
      length = grid::unit(0.09, "inches"))) +
  geom_point(data = placement[is.finite(accepted_native_l4)],
    aes(UMAP_1, UMAP_2, color = accepted_native_l4),
    size = 0.38, alpha = 0.52) +
  geom_point(data = placement[!is.finite(accepted_native_l4)],
    aes(UMAP_1, UMAP_2), color = "#CC0077", shape = 4,
    size = 0.45, alpha = 0.48, stroke = 0.35) +
  scale_color_viridis_c(option = "inferno", limits = l4_limits,
    name = "Transferred native L4") +
  facet_wrap(~library_label, ncol = 2) +
  coord_equal() +
  labs(title = "All QC-selected FLT3LKO and IL7KO cells on the fixed Bdev ATAC UMAP",
    subtitle = "Frozen HVP-depth20k projection; fixed 30-neighbor UMAP placement; magenta × lacks accepted native L4",
    x = "ATAC UMAP 1", y = "ATAC UMAP 2",
    caption = "The 26,723-cell Bdev UMAP, PeakVI model, and native Slingshot L4 curve were not refit.") +
  theme_classic(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"), legend.position = "right")
ggsave(file.path(root, "02_figures", "01_all_flt3lko_il7ko_fixed_atac_umap.png"),
  p, width = 13.5, height = 10.2, dpi = 260, bg = "white")
ggsave(file.path(root, "02_figures", "01_all_flt3lko_il7ko_fixed_atac_umap.pdf"),
  p, width = 13.5, height = 10.2, device = cairo_pdf, bg = "white")

group_summary <- placement[, .(
  cells = .N,
  supported = sum(support_call == "supported"),
  borderline = sum(support_call == "borderline"),
  out_of_reference = sum(support_call == "out_of_reference"),
  finite_native_l4 = sum(is.finite(accepted_native_l4)),
  median_native_l4 = median(accepted_native_l4, na.rm = TRUE),
  median_neighbor_spread = median(neighbor_spread)),
  by = .(knockout_target, genotype, library_id, library_label)]
fwrite(group_summary, file.path(root, "01_tables", "09_full_mapping_group_summary.tsv"),
  sep = "\t", quote = TRUE, na = "NA")

model_payload_after <- sha256(model_payload)
source_manifest_after <- sha256(file.path(cfg$inputs$frozen_model_bundle$path,
  "manifest.json"))
fixed_manifest_after <- sha256(file.path(cfg$inputs$fixed_bdev_bundle$path,
  "manifest.json"))
stopifnot(model_payload_after == cfg$inputs$frozen_model_bundle$model_payload_sha256,
  source_manifest_after == cfg$inputs$frozen_model_bundle$manifest_sha256,
  fixed_manifest_after == cfg$inputs$fixed_bdev_bundle$manifest_sha256)

support_counts <- as.list(table(factor(outcomes$support_call,
  levels = c("supported", "borderline", "out_of_reference"))))
summary <- list(
  schema_version = "bdev.flt3lko_il7ko_full_mapping_fixed_umap_summary.v1",
  status = "complete", run_identity = cfg$run_identity,
  reference_cells = length(reference_cells), query_cells = nrow(outcomes),
  library_counts = as.list(table(outcomes$library_id)),
  support_counts = support_counts,
  finite_accepted_native_l4 = sum(is.finite(outcomes$accepted_native_l4)),
  fixed_umap_placement = "30_neighbor_distance_kernel_barycenter",
  model_refit = FALSE, umap_refit = FALSE, slingshot_refit = FALSE,
  genotype_contrast = FALSE, paper_analysis = FALSE,
  implementation_unlocked = FALSE, sif_built = FALSE,
  source_manifests_unchanged = TRUE)
write_json(summary, file.path(root, "03_reports", "03_full_mapping_fixed_umap_summary.json"))
write_json(c(summary, list(
  config_path = cfg_path, config_sha256 = sha256(cfg_path),
  base_sif_sha256 = sha256(cfg$runtime$base_sif),
  scworkbench_version = as.character(packageVersion("scworkbench")),
  scwbLatent_version = as.character(packageVersion("scwbLatent")),
  scwbRepresentation_version = as.character(packageVersion("scwbRepresentation")),
  scwbReferenceMap_version = as.character(packageVersion("scwbReferenceMap")),
  frozen_model_payload_sha256 = model_payload_after,
  frozen_model_bundle_manifest_sha256 = source_manifest_after,
  fixed_bdev_bundle_manifest_sha256 = fixed_manifest_after)),
  file.path(root, "06_manifests", "03_full_mapping_execution_receipt.json"))
log_event("Completed all-cell fixed-UMAP mapping")
