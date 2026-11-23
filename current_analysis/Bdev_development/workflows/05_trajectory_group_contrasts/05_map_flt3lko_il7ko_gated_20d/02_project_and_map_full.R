#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_index <- match("--config", args)
if (is.na(config_index) || config_index == length(args)) {
  stop("Usage: 02_project_and_map_full.R --config CONFIG", call. = FALSE)
}
config_path <- normalizePath(args[[config_index + 1L]], mustWork = TRUE)

options(stringsAsFactors = FALSE, warn = 1, width = 180)
suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbLatent)
  library(scwbRepresentation)
  library(scwbReferenceMap)
  library(BPCells)
  library(Matrix)
  library(data.table)
  library(digest)
  library(jsonlite)
  library(yaml)
})

cfg <- yaml::read_yaml(config_path)
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
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
promoted_id <- function(x, output) {
  take <- x$promoted$status == "ok" & x$promoted$output == output
  stopifnot(sum(take) == 1L)
  as.character(x$promoted$artifact_id[take][[1L]])
}
as_space_matrix <- function(x) {
  if (is.data.frame(x)) {
    ids <- as.character(x$cell_id)
    x <- as.matrix(x[, setdiff(names(x), "cell_id"), drop = FALSE])
    rownames(x) <- ids
  } else x <- as.matrix(x)
  x
}

root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")), c("tables", "figures",
  "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_project_and_map_full.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(sha256(cfg$runtime$sif) == cfg$runtime$sif_sha256,
  sha256(file.path(cfg$inputs$qualified_bundle$path, "manifest.json")) ==
    cfg$inputs$qualified_bundle$manifest_sha256,
  sha256(cfg$inputs$qualification_report$path) ==
    cfg$inputs$qualification_report$sha256,
  sha256(cfg$inputs$qualification_inventory$path) ==
    cfg$inputs$qualification_inventory$sha256,
  sha256(cfg$inputs$prepared_cells$path) == cfg$inputs$prepared_cells$sha256,
  sha256(cfg$inputs$quantification_audit$path) ==
    cfg$inputs$quantification_audit$sha256)

for (package in c("scwbLatent", "scwbRepresentation", "scwbReferenceMap")) {
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  stopifnot(nzchar(manifest), sha256(manifest) == cfg$runtime[[package]]$manifest_sha256)
  namespace <- switch(package, scwbLatent = "scwblatent",
    scwbRepresentation = "scwbrepresentation",
    scwbReferenceMap = "scwbreferencemap")
  try(wb_operation_pack_unload(namespace), silent = TRUE)
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
}
required_ops <- c("scwblatent::peakvi_project",
  "scwblatent::model_bound_latent_union",
  "scwbrepresentation::knn_neighbors",
  "scwbreferencemap::reference_map_transfer")
stopifnot(all(required_ops %in% wb_list_ops()))

cells <- as.data.table(read.delim(gzfile(cfg$inputs$prepared_cells$path),
  check.names = FALSE, quote = "\"", na.strings = "NA"))
query_matrix <- BPCells::open_matrix_dir(cfg$inputs$quantified_matrix$path)
stopifnot(nrow(cells) == cfg$scope$expected_query_cells,
  nrow(query_matrix) == cfg$inputs$quantified_matrix$expected_features,
  ncol(query_matrix) == cfg$inputs$quantified_matrix$expected_cells,
  identical(as.character(colnames(query_matrix)), as.character(cells$cell_id)),
  length(unique(cells$library_id)) == 4L)
slots <- sprintf("external_gated20d_query_%05d", seq_len(nrow(cells)))
slot_map <- copy(cells)
slot_map[, slot_id := slots]
setcolorder(slot_map, c("slot_id", setdiff(names(slot_map), "slot_id")))
fwrite(slot_map, file.path(dirs[["tables"]], "01_full_query_slot_map.tsv.gz"),
  sep = "\t", quote = TRUE, na = "NA")

bundle <- wb_open_bundle(cfg$storage$mapping_bundle,
  validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
current_cells <- as.character(wb_cells(bundle, "all"))
full_query_pattern <- "^external_gated20d_query_[0-9]+$"
source_cells <- current_cells[!grepl(full_query_pattern, current_cells)]
reference_cells <- source_cells[!grepl("^external_query_[0-9]+$", source_cells)]
stopifnot(length(reference_cells) == cfg$scope$expected_reference_cells)
all_expected_cells <- c(source_cells, slots)
if (identical(current_cells, source_cells)) {
  log_event("Ingesting", length(slots), "full-query cell identities")
  query_meta <- data.frame(cell_id = slots,
    orig.ident = as.character(slot_map$library_id),
    native_l4 = rep(NA_real_, length(slots)),
    l4_weight = rep(NA_real_, length(slots)),
    cell_role = rep("flt3lko_il7ko_full_gated20d_query", length(slots)),
    stringsAsFactors = FALSE)
  wb_ingest_apply(bundle,
    set_name = sub("^set:", "", cfg$mapping$query_cell_set),
    cell_ids = slots, metadata = query_meta,
    source = "checksummed_existing_full_KO_quantification",
    created_by = "Bdev development full gated 20D mapping",
    apply_mode = "atomic", emit_report_artifact = TRUE)
  current_cells <- as.character(wb_cells(bundle, "all"))
}
stopifnot(identical(current_cells, all_expected_cells),
  identical(as.character(wb_cells(bundle, cfg$mapping$query_cell_set)), slots))

if (!(cfg$mapping$query_matrix_id %in% wb_matrices(bundle)$matrix_id)) {
  log_event("Attaching exact-axis full-query count matrix")
  query_sparse <- as(query_matrix, "dgCMatrix")
  colnames(query_sparse) <- slots
  wb_attach_matrix(bundle, cfg$mapping$query_matrix_id, query_sparse,
    layer = "counts", feature_type = "peak", cell_ids = slots,
    description = "All QC-selected FLT3LKO/IL7KO cells on qualified HVP20k axis",
    source = "existing_reciprocal_liftover_read_only_fragment_quantification",
    created_by = "Bdev development full gated 20D mapping")
  rm(query_sparse); invisible(gc())
}

model <- wb_art_load(bundle, cfg$inputs$qualified_bundle$model_artifact_id,
  expected_kind = "model_ref")
model_payload <- file.path(model$model_tree_root, model$model_primary_path)
stopifnot(sha256(model_payload) ==
  cfg$inputs$qualified_bundle$model_payload_sha256,
  as.integer(model$model_linkage$latent_dimension_count) == 20L)
model_ref <- wb_ref_artifact(cfg$inputs$qualified_bundle$model_artifact_id,
  "model")

artifacts <- wb_art_list(bundle)
projection_take <- artifacts$source_op == "scwblatent::peakvi_project" &
  artifacts$type == "obs_table" & artifacts$lifecycle_state == "current" &
  grepl("^full_gated20d_projection - ", artifacts$name)
if (cfg$mapping$query_space_id %in% wb_spaces(bundle)$space_id &&
    sum(projection_take) == 1L) {
  projection_artifact_id <- as.character(
    artifacts$artifact_id[projection_take][[1L]])
  projection_ref <- wb_ref_artifact(projection_artifact_id, "obs_table")
  log_event("Recovered existing full 20D query projection")
} else {
  stopifnot(!(cfg$mapping$query_space_id %in% wb_spaces(bundle)$space_id),
    sum(projection_take) == 0L)
  log_event("Projecting", length(slots), "query cells without model refit")
  project_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_project"))
  project <- wb_session_run_op(project_session, bundle,
    "scwblatent::peakvi_project", inputs = list(model = model_ref,
      query_matrix = cfg$mapping$query_matrix_id,
      space_id_out = cfg$mapping$query_space_id,
      chunk_size = as.integer(cfg$mapping$peakvi_chunk_size), batch_key = NULL,
      use_gpu = isTRUE(cfg$mapping$use_gpu),
      pickle_execution_trust = cfg$mapping$pickle_execution_trust,
      cell_set = cfg$mapping$query_cell_set, min_cells = 1L),
    notes = "Full KO projection through qualified 20D HVP20k model")
  projection_summary <- project$run$outputs$projection_summary$table
  stopifnot(nrow(projection_summary) == 1L,
    projection_summary$transformation[[1L]] == "transfer_without_refit",
    as.integer(projection_summary$n_latent[[1L]]) == 20L,
    as.integer(projection_summary$n_query_cells[[1L]]) == length(slots))
  project$session <- wb_session_resume(bundle,
    project$session$durable$workspace)
  project_promoted <- wb_session_promote_current(bundle, project$session,
    outputs = "projection_summary", name_prefix = "full_gated20d_projection",
    description = "All-query projection through qualified 20D model",
    created_by = "Bdev development full gated 20D mapping", fail_fast = TRUE)
  projection_artifact_id <- promoted_id(project_promoted, "projection_summary")
  projection_ref <- wb_ref_artifact(projection_artifact_id, "obs_table")
  project_promoted$session <- wb_session_discard(project_promoted$session)
}
stopifnot(sha256(model_payload) ==
  cfg$inputs$qualified_bundle$model_payload_sha256)

if (!(cfg$mapping$union_cell_set %in% wb_cell_sets(bundle))) {
  union_set_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_union_cell_set"))
  union_set <- wb_session_run_op(union_set_session, bundle,
    "cell_set_derive", inputs = list(cell_ids = c(reference_cells, slots),
      target_cell_set = cfg$mapping$union_cell_set,
      role = "analysis_universe", status = "accepted",
      source = "qualified_20d_reference_and_full_query_union",
      created_by = "Bdev development full gated 20D mapping",
      description = "Exact qualified reference and full-query union",
      min_cells = 1L), notes = "Create authenticated full gated20d union")
  union_set$session <- wb_session_resume(bundle,
    union_set$session$durable$workspace)
  union_promoted <- wb_session_promote_current(bundle, union_set$session,
    outputs = c("manifest", "rule_table"),
    name_prefix = "full_gated20d_union_cell_set",
    description = "Exact qualified reference and full-query union",
    created_by = "Bdev development full gated 20D mapping", fail_fast = TRUE)
  union_promoted$session <- wb_session_discard(union_promoted$session)
}
stopifnot(identical(as.character(wb_cells(bundle, cfg$mapping$union_cell_set)),
  c(reference_cells, slots)))

if (!(cfg$mapping$joint_space_id %in% wb_spaces(bundle)$space_id)) {
  log_event("Authenticating qualified-model reference/query union")
  union_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_latent_union"))
  union <- wb_session_run_op(union_session, bundle,
    "scwblatent::model_bound_latent_union", inputs = list(
      reference_space = cfg$inputs$qualified_bundle$reference_space,
      query_space = cfg$mapping$query_space_id, model = model_ref,
      projection_summary = projection_ref,
      space_id_out = cfg$mapping$joint_space_id, chunk_size = 1024L,
      cell_set = cfg$mapping$union_cell_set, min_cells = 1L),
    notes = "Authenticated qualified 20D Bdev/full-KO latent union")
  union$session <- wb_session_resume(bundle, union$session$durable$workspace)
  union_promoted <- wb_session_promote_current(bundle, union$session,
    outputs = "union_summary", name_prefix = "full_gated20d_union",
    description = "Authenticated qualified 20D Bdev/full-KO union",
    created_by = "Bdev development full gated 20D mapping", fail_fast = TRUE)
  union_artifact_id <- promoted_id(union_promoted, "union_summary")
  union_promoted$session <- wb_session_discard(union_promoted$session)
} else {
  artifacts <- wb_art_list(bundle)
  take <- artifacts$source_op == "scwblatent::model_bound_latent_union" &
    artifacts$type == "obs_table" & artifacts$lifecycle_state == "current" &
    grepl("^full_gated20d_union - ", artifacts$name)
  stopifnot(sum(take) == 1L)
  union_artifact_id <- as.character(artifacts$artifact_id[take][[1L]])
}

artifacts <- wb_art_list(bundle)
calibration <- wb_art_load(bundle,
  cfg$inputs$qualified_bundle$calibration_artifact_id,
  expected_kind = "model_ref")
stopifnot(as.character(calibration$model_space_id) ==
  cfg$inputs$qualified_bundle$reference_space)
calibration_ref <- wb_ref_artifact(
  cfg$inputs$qualified_bundle$calibration_artifact_id, "model")
calibration_prefix <- "gated_transfer_balanced_fold_calibration"
threshold_take <- artifacts$source_op ==
  "scwbreferencemap::reference_map_calibrate" &
  artifacts$type == "obs_table" & artifacts$lifecycle_state == "current" &
  startsWith(artifacts$name, paste0(calibration_prefix, " - ")) &
  grepl("\\(thresholds\\)$", artifacts$name)
stopifnot(sum(threshold_take) == 1L)
thresholds <- wb_art_load(bundle,
  artifacts$artifact_id[threshold_take][[1L]], expected_kind = "obs_table")$table
supported_threshold <- as.numeric(thresholds$threshold[
  thresholds$role == "supported_upper"][[1L]])
out_threshold <- as.numeric(thresholds$threshold[
  thresholds$role == "borderline_upper"][[1L]])

query_coordinates <- as_space_matrix(wb_space(bundle,
  cfg$mapping$query_space_id, cells = slots))
stopifnot(identical(rownames(query_coordinates), slots),
  ncol(query_coordinates) == 20L)
batch_size <- as.integer(cfg$mapping$transfer_batch_size)
batches <- split(slots, ceiling(seq_along(slots) / batch_size))
checkpoint_root <- file.path(dirs[["objects"]],
  "01_bounded_transfer_checkpoints")
dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)
batch_contract <- digest::digest(list(
  schema_version = "bdev.full_gated20d_batch.v1",
  run_identity = cfg$run_identity,
  model_payload_sha256 = cfg$inputs$qualified_bundle$model_payload_sha256,
  calibration_artifact_id = cfg$inputs$qualified_bundle$calibration_artifact_id,
  primary_k = as.integer(cfg$mapping$primary_k),
  supported_threshold = supported_threshold,
  out_of_reference_threshold = out_threshold,
  minimum_annotation_confidence =
    as.numeric(cfg$mapping$minimum_annotation_confidence),
  batch_size = batch_size), algo = "sha256")

validate_checkpoint <- function(value, batch_index, batch_cells) {
  stopifnot(is.list(value),
    identical(value$schema_version, "bdev.full_gated20d_batch.v1"),
    identical(value$batch_contract_sha256, batch_contract),
    identical(as.integer(value$batch_index), as.integer(batch_index)),
    identical(as.character(value$query_cells), as.character(batch_cells)),
    is.data.frame(value$neighbors), is.data.frame(value$prediction),
    is.data.frame(value$categorical), is.data.frame(value$continuous),
    is.data.frame(value$acceptance),
    nrow(value$neighbors) == length(batch_cells) *
      as.integer(cfg$mapping$primary_k),
    identical(as.character(value$prediction$query_cell_id), batch_cells),
    identical(as.character(value$acceptance$query_cell_id), batch_cells))
  invisible(TRUE)
}

neighbor_batches <- prediction_batches <- categorical_batches <-
  continuous_batches <- acceptance_batches <- vector("list", length(batches))
available_spaces <- wb_spaces(bundle)$space_id
log_event("Mapping", length(slots), "queries in", length(batches),
  "bounded batches")

# Register every bounded query space before creating the durable execution
# session. This keeps the exact same scientific batch contract while avoiding
# one multi-gigabyte bundle snapshot per batch.
for (batch_index in seq_along(batches)) {
  batch_cells <- batches[[batch_index]]
  batch_space_id <- sprintf("%s_batch_%03d", cfg$mapping$query_space_id,
    batch_index)
  if (!(batch_space_id %in% available_spaces)) {
    batch_coordinates <- query_coordinates[batch_cells, , drop = FALSE]
    frame <- data.frame(cell_id = batch_cells, batch_coordinates,
      check.names = FALSE, row.names = NULL)
    wb_attach_space(bundle, batch_space_id, frame, cell_ids = batch_cells,
      source = paste0("provenance_recorded_subset_of:",
        cfg$mapping$query_space_id),
      created_by = "Bdev development full gated 20D mapping",
      emit_report_artifact = FALSE)
    available_spaces <- c(available_spaces, batch_space_id)
  }
}
checkpoint_paths <- file.path(checkpoint_root,
  sprintf("batch_%03d.rds", seq_along(batches)))
missing_batches <- which(!file.exists(checkpoint_paths))
shared_session <- if (length(missing_batches)) wb_session_create(bundle,
  cfg$storage$sessions, sprintf("%s_transfer_resume_%s", cfg$run_identity,
    format(Sys.time(), "%Y%m%dT%H%M%S"))) else NULL
for (batch_index in seq_along(batches)) {
  batch_cells <- batches[[batch_index]]
  checkpoint_path <- file.path(checkpoint_root,
    sprintf("batch_%03d.rds", batch_index))
  if (file.exists(checkpoint_path)) {
    checkpoint <- readRDS(checkpoint_path)
    validate_checkpoint(checkpoint, batch_index, batch_cells)
    neighbor_batches[[batch_index]] <- as.data.table(checkpoint$neighbors)
    prediction_batches[[batch_index]] <- as.data.table(checkpoint$prediction)
    categorical_batches[[batch_index]] <- as.data.table(checkpoint$categorical)
    continuous_batches[[batch_index]] <- as.data.table(checkpoint$continuous)
    acceptance_batches[[batch_index]] <- as.data.table(checkpoint$acceptance)
    log_event("Recovered batch", batch_index, "of", length(batches))
    next
  }
  log_event("Mapping batch", batch_index, "of", length(batches),
    "with", length(batch_cells), "cells")
  batch_space_id <- sprintf("%s_batch_%03d", cfg$mapping$query_space_id,
    batch_index)
  batch_knn <- wb_session_run_op(shared_session, bundle,
    "scwbrepresentation::knn_neighbors", inputs = list(
      space = cfg$mapping$joint_space_id, query_cells = batch_cells,
      candidate_cells = reference_cells,
      k = as.integer(cfg$mapping$primary_k), include_query = FALSE,
      space_chunk_size = as.integer(cfg$mapping$space_chunk_size)),
    notes = sprintf("Full gated20d neighbor batch %d/%d",
      batch_index, length(batches)))
  batch_transfer <- wb_session_run_op(batch_knn$session, bundle,
    "scwbreferencemap::reference_map_transfer", inputs = list(
      in_space = batch_space_id,
      reference_space = cfg$inputs$qualified_bundle$reference_space,
      neighbors = batch_knn$run$outputs$neighbors,
      calibration = calibration_ref,
      support_output_key = "full_gated20d_support",
      annotation_output_key = "full_gated20d_l4_domain",
      continuous_output_keys = cfg$mapping$continuous_output_key,
      continuous_acceptance = list(mode = "support_and_annotation",
        accepted_annotation_values = cfg$mapping$accepted_annotation_value,
        minimum_annotation_confidence =
          as.numeric(cfg$mapping$minimum_annotation_confidence)),
      supported_threshold = supported_threshold,
      out_of_reference_threshold = out_threshold,
      cell_ids = batch_cells, min_cells = 1L,
      space_chunk_size = as.integer(cfg$mapping$space_chunk_size)),
    notes = sprintf("Full gated20d transfer batch %d/%d",
      batch_index, length(batches)))
  checkpoint <- list(
    schema_version = "bdev.full_gated20d_batch.v1",
    batch_contract_sha256 = batch_contract,
    batch_index = as.integer(batch_index),
    query_cells = as.character(batch_cells),
    neighbors = as.data.frame(batch_knn$run$outputs$neighbors$table),
    prediction = as.data.frame(
      batch_transfer$run$outputs$per_query_summary$table),
    categorical = as.data.frame(
      batch_transfer$run$outputs$categorical_estimates$table),
    continuous = as.data.frame(
      batch_transfer$run$outputs$continuous_estimates$table),
    acceptance = as.data.frame(
      batch_transfer$run$outputs$continuous_acceptance$table),
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  validate_checkpoint(checkpoint, batch_index, batch_cells)
  atomic_save_rds(checkpoint, checkpoint_path)
  neighbor_batches[[batch_index]] <- as.data.table(checkpoint$neighbors)
  prediction_batches[[batch_index]] <- as.data.table(checkpoint$prediction)
  categorical_batches[[batch_index]] <- as.data.table(checkpoint$categorical)
  continuous_batches[[batch_index]] <- as.data.table(checkpoint$continuous)
  acceptance_batches[[batch_index]] <- as.data.table(checkpoint$acceptance)
  shared_session <- batch_transfer$session
  log_event("Checkpointed batch", batch_index, "sha256",
    sha256(checkpoint_path))
}
if (!is.null(shared_session)) shared_session <- wb_session_discard(shared_session)

neighbors <- rbindlist(neighbor_batches)
prediction <- rbindlist(prediction_batches)
categorical <- rbindlist(categorical_batches)
continuous <- rbindlist(continuous_batches)
acceptance <- rbindlist(acceptance_batches)
neighbors <- neighbors[order(match(query_cell_id, slots), neighbor_rank)]
prediction <- prediction[match(slots, query_cell_id)]
acceptance <- acceptance[match(slots, query_cell_id)]
primary_categorical <- categorical[k == as.integer(cfg$mapping$primary_k)]
primary_categorical <- primary_categorical[match(slots, query_cell_id)]
primary_continuous <- continuous[k == as.integer(cfg$mapping$primary_k) &
  reference_key == "L4"]
primary_continuous <- primary_continuous[match(slots, query_cell_id)]
stopifnot(nrow(prediction) == length(slots),
  identical(as.character(prediction$query_cell_id), slots),
  nrow(acceptance) == length(slots),
  identical(as.character(acceptance$query_cell_id), slots),
  identical(as.character(primary_categorical$query_cell_id), slots),
  identical(as.character(primary_continuous$query_cell_id), slots),
  all(acceptance$accepted == is.finite(acceptance$accepted_value)),
  all(is.na(acceptance$accepted_value[!acceptance$accepted])))
expected_accept <- prediction$support_call != "out_of_reference" &
  primary_categorical$raw_prediction == cfg$mapping$accepted_annotation_value &
  is.finite(primary_categorical$confidence) &
  primary_categorical$confidence >=
    as.numeric(cfg$mapping$minimum_annotation_confidence) &
  is.finite(primary_continuous$raw_estimate)
stopifnot(identical(as.logical(acceptance$accepted), expected_accept))

result <- data.table(slot_id = slots,
  support_call = as.character(prediction$support_call),
  distance_statistic = as.numeric(prediction$distance_statistic),
  predicted_domain = as.character(primary_categorical$raw_prediction),
  domain_confidence = as.numeric(primary_categorical$confidence),
  raw_native_l4 = as.numeric(primary_continuous$raw_estimate),
  accepted_native_l4 = as.numeric(acceptance$accepted_value),
  accepted = as.logical(acceptance$accepted),
  acceptance_reason = as.character(acceptance$acceptance_reason))
result <- merge(slot_map, result, by = "slot_id", sort = FALSE)
result <- result[match(slots, slot_id)]
stopifnot(nrow(result) == cfg$scope$expected_query_cells)
group_summary <- result[, .(
  cells = .N,
  supported = sum(support_call == "supported"),
  borderline = sum(support_call == "borderline"),
  out_of_reference = sum(support_call == "out_of_reference"),
  accepted_native_l4 = sum(accepted),
  withheld_native_l4 = sum(!accepted)),
  by = .(library_id, knockout_target, genotype)]

fwrite(result, file.path(dirs[["tables"]],
  "02_full_gated_20d_mapping.tsv.gz"), sep = "\t", quote = TRUE, na = "NA")
fwrite(group_summary, file.path(dirs[["tables"]],
  "03_full_gated_20d_group_counts.tsv"), sep = "\t", quote = TRUE, na = "NA")
fwrite(thresholds, file.path(dirs[["tables"]],
  "04_qualified_support_thresholds.tsv"), sep = "\t", quote = TRUE,
  na = "NA")

support_counts <- as.list(table(factor(result$support_call,
  levels = c("supported", "borderline", "out_of_reference"))))
reason_counts <- as.list(table(factor(result$acceptance_reason,
  levels = c("accepted", "out_of_reference", "annotation_not_accepted",
    "annotation_confidence_too_low", "continuous_value_unavailable"))))
report <- list(
  schema_version = "bdev.flt3lko_il7ko_full_gated_20d_report.v1",
  status = "passed", gate_passed = TRUE, run_identity = cfg$run_identity,
  reference_cells = length(reference_cells), query_cells = nrow(result),
  selected_features = 20000L, latent_dimensions = 20L,
  model_artifact_id = cfg$inputs$qualified_bundle$model_artifact_id,
  model_payload_sha256 = sha256(model_payload),
  projection_artifact_id = projection_artifact_id,
  union_artifact_id = union_artifact_id,
  calibration_artifact_id = cfg$inputs$qualified_bundle$calibration_artifact_id,
  primary_k = as.integer(cfg$mapping$primary_k),
  minimum_annotation_confidence =
    as.numeric(cfg$mapping$minimum_annotation_confidence),
  support_counts = support_counts,
  acceptance_reason_counts = reason_counts,
  accepted_query_cells = sum(result$accepted),
  raw_estimate_rows = nrow(result),
  checkpoint_files = length(list.files(checkpoint_root,
    pattern = "^batch_[0-9]+\\.rds$")),
  qualified_source_manifest_sha256 = sha256(file.path(
    cfg$inputs$qualified_bundle$path, "manifest.json")),
  mapping_bundle_manifest_sha256 = sha256(file.path(
    cfg$storage$mapping_bundle, "manifest.json")),
  model_refit = FALSE, genotype_contrast_performed = FALSE,
  paper_analysis = FALSE, sif_built = FALSE)
write_json(report, file.path(dirs[["reports"]],
  "01_full_gated_20d_mapping.json"))
write_json(list(
  schema_version = "bdev.flt3lko_il7ko_full_gated_20d_execution.v1",
  config_path = config_path, config_sha256 = sha256(config_path),
  scworkbench = as.character(packageVersion("scworkbench")),
  scwbLatent = as.character(packageVersion("scwbLatent")),
  scwbRepresentation = as.character(packageVersion("scwbRepresentation")),
  scwbReferenceMap = as.character(packageVersion("scwbReferenceMap")),
  qualified_source_manifest_unchanged = sha256(file.path(
    cfg$inputs$qualified_bundle$path, "manifest.json")) ==
    cfg$inputs$qualified_bundle$manifest_sha256,
  source_archives_modified = FALSE), file.path(dirs[["manifests"]],
  "02_execution_receipt.json"))
log_event("Full gated 20D mapping passed")
