#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config"))
  stop("Usage: 03_run_reference_mapping.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbRepresentation)
  library(scwbReferenceMap)
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
sink(file.path(dirs[["logs"]], "01_reference_mapping_qualification.log"),
  split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
write_tsv <- function(x, name) write.table(x, file.path(dirs[["tables"]], name),
  sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
run_ref <- function(run, output) wb_ref_run(run$run_id, output)
bind_rows <- function(rows) {
  if (!length(rows)) return(data.frame())
  fields <- unique(unlist(lapply(rows, names), use.names = FALSE))
  normalized <- lapply(rows, function(x) {
    for (field in setdiff(fields, names(x))) x[[field]] <- NA
    x[, fields, drop = FALSE]
  })
  out <- do.call(rbind, normalized); rownames(out) <- NULL; out
}
package_identity <- function(package, expected) {
  observed <- as.character(packageVersion(package))
  stopifnot(identical(observed, as.character(expected)))
  observed
}
space_frame <- function(x) data.frame(cell_id = rownames(x),
  as.data.frame(x, check.names = FALSE), check.names = FALSE,
  stringsAsFactors = FALSE)
ensure_space <- function(bundle, space_id, values) {
  existing_ids <- as.character(wb_spaces(bundle)$space_id)
  if (space_id %in% existing_ids) {
    observed <- wb_space(bundle, space_id, cells = rownames(values))
    stopifnot(identical(dim(observed), dim(values)),
      identical(rownames(observed), rownames(values)),
      isTRUE(all.equal(as.numeric(observed), as.numeric(values),
        tolerance = 1e-12, check.attributes = FALSE)))
    return(invisible(FALSE))
  }
  wb_attach_space(bundle, space_id, space_frame(values),
    cell_ids = rownames(values), source = "04_reference_mapping_qualification",
    created_by = "scworkbench_Bdev development qualification")
  invisible(TRUE)
}

stopifnot(identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256))
versions <- list(
  scworkbench = package_identity("scworkbench", cfg$runtime$scworkbench$version),
  scwbRepresentation = package_identity("scwbRepresentation",
    cfg$runtime$scwbRepresentation$version),
  scwbReferenceMap = package_identity("scwbReferenceMap",
    cfg$runtime$scwbReferenceMap$version))
for (package in c("scwbRepresentation", "scwbReferenceMap")) {
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  expected <- cfg$runtime[[package]]$manifest_sha256
  stopifnot(nzchar(manifest), identical(sha256(manifest), expected))
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
}
required_ops <- c("scwbrepresentation::knn_neighbors",
  "scwbreferencemap::reference_map_calibrate",
  "scwbreferencemap::reference_map_transfer")
stopifnot(all(required_ops %in% wb_list_ops()))
for (operation in required_ops)
  write_json(wb_op_interface(operation), file.path(dirs[["manifests"]],
    paste0("interface_", gsub("[^A-Za-z0-9]+", "_", operation), ".json")))

source_manifest <- file.path(cfg$inputs$source_bundle$path, "manifest.json")
source_manifest_before <- sha256(source_manifest)
stopifnot(identical(source_manifest_before,
  cfg$inputs$source_bundle$manifest_sha256),
  identical(sha256(cfg$inputs$authenticated_l4$path),
    cfg$inputs$authenticated_l4$sha256),
  identical(sha256(cfg$inputs$fixture$path), cfg$inputs$fixture$sha256))

fixture <- read.delim(cfg$inputs$fixture$path, check.names = FALSE,
  quote = "\"", na.strings = "NA")
fixture_cells <- as.character(fixture$cell_id)
stopifnot(length(fixture_cells) == cfg$inputs$fixture$expected_cells,
  !anyDuplicated(fixture_cells))

log_event("Reading fixture inputs from the active 03c source bundle")
source_bundle <- wb_open_bundle(cfg$inputs$source_bundle$path, validate = FALSE,
  read_only = TRUE)
all_cells <- as.character(wb_cells(source_bundle, "all"))
stopifnot(length(all_cells) == cfg$inputs$source_bundle$expected_cells,
  all(fixture_cells %in% all_cells),
  identical(fixture_cells, all_cells[all_cells %in% fixture_cells]))
metadata <- wb_meta(source_bundle, cells = fixture_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$sample_field,
  cfg$inputs$annotation_field))
stopifnot(identical(as.character(metadata$cell_id), fixture_cells))
full_lsi <- wb_space(source_bundle, cfg$inputs$space, cells = fixture_cells)
dims <- as.integer(unlist(cfg$inputs$dimensions, use.names = FALSE))
stopifnot(max(dims) <= ncol(full_lsi), min(dims) >= 1L)
qualification_lsi <- full_lsi[, dims, drop = FALSE]
axis_features <- head(wb_feature_ids(source_bundle,
  cfg$inputs$source_axis_matrix), as.integer(cfg$inputs$fixture_axis_features))
axis_matrix <- wb_slice(source_bundle, cfg$inputs$source_axis_matrix,
  features = axis_features, cells = fixture_cells)
wb_close_bundle(source_bundle)

fixture$sample_id <- as.character(metadata[[cfg$inputs$sample_field]])
fixture$annotation <- as.character(metadata[[cfg$inputs$annotation_field]])
samples <- sort(unique(fixture$sample_id), method = "radix")
stopifnot(length(samples) == cfg$inputs$fixture$expected_samples,
  !anyNA(fixture$annotation), all(nzchar(fixture$annotation)))

family <- readRDS(cfg$inputs$authenticated_l4$path)
selected <- as.character(family$key) == cfg$inputs$authenticated_l4$key &
  as.character(family$cell_id) %in% fixture_cells & is.finite(as.numeric(family$value))
l4 <- stats::setNames(as.numeric(family$value[selected]),
  as.character(family$cell_id[selected]))
stopifnot(length(l4) == length(fixture_cells), !anyDuplicated(names(l4)),
  all(fixture_cells %in% names(l4)))
fixture$native_l4 <- unname(l4[fixture_cells])
write_tsv(fixture, "01_fixture_cells_and_truth.tsv")

log_event("Creating fresh fixture-only qualification bundle with LSI dimensions",
  paste(cfg$inputs$dimensions, collapse = ","))
bundle_exists <- dir.exists(cfg$storage$mapping_fixture_bundle)
bundle <- if (bundle_exists) wb_open_bundle(cfg$storage$mapping_fixture_bundle,
  validate = FALSE, read_only = FALSE) else
  wb_init_bundle(cfg$storage$mapping_fixture_bundle, overwrite = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
if (!bundle_exists) {
  wb_define_cells(bundle, fixture_cells)
  wb_attach_matrix(bundle, cfg$inputs$fixture_axis_matrix, axis_matrix,
    layer = "data", feature_type = "motif_activity",
    description = paste("Ten real chromVAR features retained only to authenticate",
      "the fixture cell and orig.ident axes for reference-map calibration."),
    cell_ids = fixture_cells, source = "03c_active_bundle_fixture_extraction",
    created_by = "scworkbench_Bdev development qualification")
  wb_attach_metadata(bundle, metadata, cell_ids = fixture_cells,
    source = "03c_active_bundle_fixture_extraction",
    created_by = "scworkbench_Bdev development qualification")
} else {
  stopifnot(identical(as.character(wb_cells(bundle, "all")), fixture_cells),
    cfg$inputs$fixture_axis_matrix %in% wb_matrices(bundle)$matrix_id)
}
fixture_space_id <- "step4_fixture_lsi_dims2_30"
ensure_space(bundle, fixture_space_id, qualification_lsi)

fold_table <- data.frame(fold = seq_along(samples), held_out_sample = samples,
  n_query = vapply(samples, function(x) sum(fixture$sample_id == x), integer(1)),
  n_reference = vapply(samples, function(x) sum(fixture$sample_id != x), integer(1)),
  stringsAsFactors = FALSE)
write_tsv(fold_table, "02_leave_one_sample_out_folds.tsv")

all_predictions <- list(); fold_metrics <- list(); calibration_metrics <- list()
threshold_rows <- list(); coverage_rows <- list(); grouped_rows <- list()
for (fold_index in seq_along(samples)) {
  held_out <- samples[[fold_index]]
  fold_id <- paste0("fold", fold_index, "_sample", gsub("[^A-Za-z0-9]", "_", held_out))
  query_cells <- fixture$cell_id[fixture$sample_id == held_out]
  reference_cells <- fixture$cell_id[fixture$sample_id != held_out]
  reference_space_id <- paste0("step4_reference_", fold_id)
  query_space_id <- paste0("step4_query_", fold_id)
  log_event("Starting", fold_id, "reference", length(reference_cells),
    "query", length(query_cells))
  ensure_space(bundle, reference_space_id,
    qualification_lsi[reference_cells, , drop = FALSE])
  ensure_space(bundle, query_space_id,
    qualification_lsi[query_cells, , drop = FALSE])

  reference_annotation <- wb_prim_cell_vector(stats::setNames(
    fixture$annotation[match(reference_cells, fixture$cell_id)], reference_cells),
    value_type = "categorical")
  reference_l4 <- wb_prim_cell_vector_family(
    rep(cfg$inputs$authenticated_l4$key, length(reference_cells)),
    reference_cells, unname(l4[reference_cells]), value_type = "continuous")

  calibration_session <- wb_session_create(bundle,
    workspace_root = cfg$storage$sessions,
    session_id = paste0(fold_id, "_calibrate_fresh_bundle_retry3"))
  calibration <- wb_session_run_op(calibration_session, bundle,
    "scwbreferencemap::reference_map_calibrate", inputs = list(
      reference_matrix = cfg$inputs$fixture_axis_matrix,
      reference_space = reference_space_id,
      reference_annotation = reference_annotation,
      reference_continuous = reference_l4,
      holdout_group_field = cfg$inputs$sample_field,
      k_values = as.integer(unlist(cfg$mapping$k_values, use.names = FALSE)),
      primary_k = as.integer(cfg$mapping$primary_k),
      support_quantiles = as.list(as.numeric(unlist(
        cfg$mapping$support_quantiles, use.names = FALSE))),
      primary_continuous_key = cfg$inputs$authenticated_l4$key,
      space_chunk_size = as.integer(cfg$mapping$space_chunk_size),
      query_chunk_size = as.integer(cfg$mapping$query_chunk_size)),
    notes = paste("Step 4 held-out-sample reference calibration", fold_id))
  promoted <- wb_session_promote_current(bundle, calibration$session,
    outputs = names(calibration$run$outputs),
    name_prefix = paste0("step4_retry3_", fold_id),
    description = paste("Step 4 reference-map calibration", fold_id),
    created_by = "scworkbench_Bdev development qualification", fail_fast = TRUE)
  model_row <- promoted$promoted$status == "ok" &
    promoted$promoted$output == "calibration"
  stopifnot(sum(model_row) == 1L)
  calibration_ref <- wb_ref_artifact(
    promoted$promoted$artifact_id[model_row][[1L]], "model")

  knn_session <- wb_session_create(bundle, workspace_root = cfg$storage$sessions,
    session_id = paste0(fold_id, "_knn_fresh_bundle_retry3"))
  knn <- wb_session_run_op(knn_session, bundle,
    "scwbrepresentation::knn_neighbors", inputs = list(
      space = fixture_space_id, query_cells = query_cells,
      candidate_cells = reference_cells,
      k = as.integer(cfg$mapping$neighbor_k), include_query = FALSE,
      backend = cfg$mapping$backend,
      space_chunk_size = as.integer(cfg$mapping$space_chunk_size),
      query_chunk_size = as.integer(cfg$mapping$query_chunk_size)),
    notes = paste("Step 4 exact query-to-reference neighbors", fold_id))

  transfer_session <- wb_session_create(bundle,
    workspace_root = cfg$storage$sessions,
    session_id = paste0(fold_id, "_transfer_fresh_bundle_retry3"))
  calibrated_thresholds <- calibration$run$outputs$thresholds$table
  supported_threshold <- calibrated_thresholds$threshold[
    calibrated_thresholds$role == "supported_upper"][[1L]]
  out_of_reference_threshold <- calibrated_thresholds$threshold[
    calibrated_thresholds$role == "borderline_upper"][[1L]]
  transfer <- wb_session_run_op(transfer_session, bundle,
    "scwbreferencemap::reference_map_transfer", inputs = list(
      in_space = query_space_id, reference_space = reference_space_id,
      neighbors = knn$run$outputs$neighbors, calibration = calibration_ref,
      support_output_key = "step4_map_support",
      annotation_output_key = "step4_mapped_cluster",
      continuous_output_keys = "mapped_native_L4",
      supported_threshold = as.numeric(supported_threshold),
      out_of_reference_threshold = as.numeric(out_of_reference_threshold),
      cell_ids = query_cells, min_cells = 1L,
      space_chunk_size = as.integer(cfg$mapping$space_chunk_size)),
    notes = paste("Step 4 held-out-sample transfer", fold_id))

  prediction <- transfer$run$outputs$per_query_summary$table
  prediction$fold <- fold_index
  prediction$held_out_sample <- held_out
  truth_index <- match(prediction$query_cell_id, fixture$cell_id)
  stopifnot(!anyNA(truth_index), identical(as.character(prediction$query_cell_id),
    query_cells))
  prediction$true_annotation <- fixture$annotation[truth_index]
  prediction$true_native_l4 <- fixture$native_l4[truth_index]
  prediction$absolute_l4_error <- abs(prediction$primary_continuous_value -
    prediction$true_native_l4)
  prediction$annotation_correct <- !is.na(prediction$accepted_annotation) &
    prediction$accepted_annotation == prediction$true_annotation
  all_predictions[[fold_index]] <- prediction

  accepted <- prediction$support_call != "out_of_reference"
  fold_metrics[[fold_index]] <- data.frame(fold = fold_index,
    held_out_sample = held_out, n_query = nrow(prediction),
    n_supported = sum(prediction$support_call == "supported"),
    n_borderline = sum(prediction$support_call == "borderline"),
    n_out_of_reference = sum(prediction$support_call == "out_of_reference"),
    accepted_fraction = mean(accepted),
    annotation_accuracy_accepted = if (any(accepted))
      mean(prediction$annotation_correct[accepted]) else NA_real_,
    l4_mae_accepted = if (any(accepted))
      mean(prediction$absolute_l4_error[accepted], na.rm = TRUE) else NA_real_,
    l4_rmse_accepted = if (any(accepted)) sqrt(mean(
      (prediction$primary_continuous_value[accepted] -
        prediction$true_native_l4[accepted])^2, na.rm = TRUE)) else NA_real_,
    l4_spearman_accepted = if (sum(accepted) >= 3L) suppressWarnings(cor(
      prediction$primary_continuous_value[accepted],
      prediction$true_native_l4[accepted], method = "spearman",
      use = "complete.obs")) else NA_real_, stringsAsFactors = FALSE)

  annotate <- function(table, output) {
    table$fold <- fold_index; table$held_out_sample <- held_out;
    table$output <- output; table
  }
  calibration_metrics[[fold_index]] <- annotate(
    calibration$run$outputs$calibration_metrics$table, "calibration_metrics")
  threshold_rows[[fold_index]] <- annotate(
    calibration$run$outputs$thresholds$table, "thresholds")
  coverage_rows[[fold_index]] <- annotate(
    calibration$run$outputs$coverage_limitations$table, "coverage_limitations")
  grouped_rows[[fold_index]] <- annotate(
    transfer$run$outputs$grouped_summary$table, "grouped_summary")
  log_event("Completed", fold_id)
}

predictions <- bind_rows(all_predictions)
metrics <- bind_rows(fold_metrics)
write_tsv(predictions, "03_held_out_transfer_predictions.tsv")
write_tsv(metrics, "04_fold_metrics.tsv")
write_tsv(bind_rows(calibration_metrics), "05_calibration_metrics.tsv")
write_tsv(bind_rows(threshold_rows), "06_support_thresholds.tsv")
write_tsv(bind_rows(coverage_rows), "07_calibration_coverage.tsv")
write_tsv(bind_rows(grouped_rows), "08_support_grouped_summary.tsv")

source_manifest_after <- sha256(source_manifest)
out_rows <- predictions$support_call == "out_of_reference"
acceptance <- list(
  all_samples_held_out_once = identical(sort(unique(predictions$held_out_sample),
    method = "radix"), samples),
  complete_query_accounting = nrow(predictions) == length(fixture_cells) &&
    !anyDuplicated(predictions$query_cell_id),
  out_of_reference_values_withheld = all(is.na(
    predictions$primary_continuous_value[out_rows])) && all(is.na(
    predictions$accepted_annotation[out_rows])),
  source_manifest_unchanged = identical(source_manifest_before,
    source_manifest_after),
  native_l4_used = identical(cfg$inputs$authenticated_l4$orientation, "native"),
  new_sif_built = FALSE,
  implementation_unlocked = FALSE)
stopifnot(all(unlist(acceptance[c("all_samples_held_out_once",
  "complete_query_accounting", "out_of_reference_values_withheld",
  "source_manifest_unchanged", "native_l4_used")])))

accepted <- predictions$support_call != "out_of_reference"
summary <- list(
  schema_version = "bdev.reference_mapping_qualification_summary.v1",
  status = "complete",
  run_identity = cfg$run_identity,
  scope = cfg$scope,
  runtime_versions = versions,
  fixture_cells = length(fixture_cells),
  samples = as.list(samples),
  lsi_dimensions = as.list(dims),
  folds = nrow(metrics),
  aggregate = list(
    accepted_fraction = mean(accepted),
    annotation_accuracy_accepted = mean(predictions$annotation_correct[accepted]),
    l4_mae_accepted = mean(predictions$absolute_l4_error[accepted], na.rm = TRUE),
    l4_rmse_accepted = sqrt(mean((predictions$primary_continuous_value[accepted] -
      predictions$true_native_l4[accepted])^2, na.rm = TRUE)),
    l4_spearman_accepted = suppressWarnings(cor(
      predictions$primary_continuous_value[accepted],
      predictions$true_native_l4[accepted], method = "spearman",
      use = "complete.obs"))),
  source_bundle_manifest_sha256 = source_manifest_after,
  qualification_bundle_manifest_sha256 = sha256(file.path(
    cfg$storage$mapping_fixture_bundle, "manifest.json")),
  acceptance = acceptance,
  scientific_boundary = list(
    validates = "held-out Bdev LSI reference-map transfer semantics",
    does_not_validate = "bulk ATAC PeakVI projection",
    next_gate = paste("depth-matched Bdev pseudobulk/pseudo-cell projection",
      "through a Bdev-only frozen PeakVI model"),
    infrastructure_limitation = paste("Calibration promotion on the migrated",
      "legacy-derived clone is blocked because artifact_model_ref lacks",
      "linkage_json; qualification therefore uses a fresh fixture-only bundle",
      "with real source data and makes no core modification.")))
write_json(summary, file.path(dirs[["reports"]],
  "01_reference_mapping_qualification_summary.json"))
write_json(list(schema_version = "bdev.reference_mapping_execution_receipt.v1",
  status = "complete", run_identity = cfg$run_identity,
  config_path = cfg_path, config_sha256 = sha256(cfg_path),
  source_manifest_sha256_before = source_manifest_before,
  source_manifest_sha256_after = source_manifest_after,
  archive_accessed = FALSE, sif_built = FALSE),
  file.path(dirs[["manifests"]], "03_execution_receipt.json"))
log_event("Step 4 mapping-semantics qualification complete")
