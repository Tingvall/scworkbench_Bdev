#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 02_select_variance_freeze.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbRNA)
  library(scwbLatent)
  library(Matrix)
  library(yaml)
  library(jsonlite)
  library(digest)
})

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
promoted_id <- function(promoted, output) {
  take <- promoted$promoted$status == "ok" & promoted$promoted$output == output
  stopifnot(sum(take) == 1L)
  as.character(promoted$promoted$artifact_id[take][[1L]])
}

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
stopifnot(is.null(cfg$base_config),
  identical(cfg$scope$external_query_use, "prohibited"),
  identical(cfg$scope$immgen_identity_use, "prohibited"),
  identical(cfg$feature_selection$operation, "scwbrna::hvg_select"),
  identical(cfg$feature_selection$method, "variance"),
  is.null(cfg$feature_selection$batch_key))

root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_select_variance_and_freeze.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256),
  identical(as.character(packageVersion("scworkbench")),
    cfg$runtime$scworkbench$version),
  identical(as.character(packageVersion("scwbRNA")), cfg$runtime$scwbRNA$version),
  identical(as.character(packageVersion("scwbLatent")),
    cfg$runtime$scwbLatent$version))
rna_manifest <- system.file("scworkbench", "operation-pack.json", package = "scwbRNA")
latent_manifest <- system.file("scworkbench", "operation-pack.json",
  package = "scwbLatent")
stopifnot(nzchar(rna_manifest), nzchar(latent_manifest),
  identical(sha256(rna_manifest), cfg$runtime$scwbRNA$manifest_sha256),
  identical(sha256(latent_manifest), cfg$runtime$scwbLatent$manifest_sha256))
wb_operation_pack_activate(rna_manifest, lib.loc = .libPaths())
wb_operation_pack_activate(latent_manifest, lib.loc = .libPaths())
stopifnot(all(c("scwbrna::hvg_select", "scwblatent::peakvi_fit") %in%
  wb_list_ops()))

source_manifest <- file.path(cfg$inputs$source_bundle$path, "manifest.json")
source_manifest_before <- sha256(source_manifest)
runtime_receipt <- file.path(dirs[["manifests"]],
  "01_runtime_and_clone_receipt.json")
stopifnot(file.exists(runtime_receipt),
  identical(source_manifest_before,
    cfg$inputs$source_bundle$manifest_sha256),
  identical(sha256(cfg$inputs$authenticated_l4$model_path),
    cfg$inputs$authenticated_l4$model_sha256),
  identical(sha256(cfg$inputs$authenticated_l4$pseudotime_path),
    cfg$inputs$authenticated_l4$pseudotime_sha256),
  identical(sha256(cfg$inputs$authenticated_l4$weights_path),
    cfg$inputs$authenticated_l4$weights_sha256),
  identical(cfg$inputs$authenticated_l4$orientation, "native"))

selection_bundle <- wb_open_bundle(cfg$storage$selection_clone,
  validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(selection_bundle), silent = TRUE), add = TRUE)
selection_artifacts <- wb_art_list(selection_bundle)
selection_take <- selection_artifacts$source_op == "scwbrna::hvg_select" &
  selection_artifacts$input_matrix_id %in% cfg$inputs$count_matrix &
  selection_artifacts$type == "feature_table" &
  selection_artifacts$lifecycle_state == "current"
summary_take <- selection_artifacts$source_op == "scwbrna::hvg_select" &
  selection_artifacts$input_matrix_id %in% cfg$inputs$count_matrix &
  selection_artifacts$type == "obs_table" &
  selection_artifacts$lifecycle_state == "current"
if (sum(selection_take) == 1L && sum(summary_take) == 1L) {
  selection_id <- as.character(selection_artifacts$artifact_id[selection_take][[1L]])
  selection_summary_id <- as.character(
    selection_artifacts$artifact_id[summary_take][[1L]])
  selection <- wb_art_load(selection_bundle, selection_id,
    expected_kind = "feature_table")$table
  selection_summary <- wb_art_load(selection_bundle, selection_summary_id,
    expected_kind = "obs_table")$table
  log_event("Recovered existing authenticated variance selection")
} else {
  stopifnot(sum(selection_take) == 0L, sum(summary_take) == 0L)
  log_event("Running authenticated ATAC variance selection on raw counts")
  session <- wb_session_create(selection_bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_variance20k_select"))
  selected <- wb_session_run_op(session, selection_bundle,
    "scwbrna::hvg_select",
    inputs = list(matrix = cfg$inputs$count_matrix,
      n_features = as.integer(cfg$feature_selection$n_features),
      method = cfg$feature_selection$method,
      min_mean = NULL, max_mean = NULL, min_dispersion = NULL,
      batch_key = NULL, exclude_features = NULL,
      feature_type = cfg$feature_selection$feature_type,
      seed = as.integer(cfg$feature_selection$seed),
      chunk_size = as.integer(cfg$feature_selection$chunk_size),
      cell_set = "all", min_cells = 1L),
    notes = paste("ATAC-compatible non-depth-adjusted raw-count variance",
      "selection; Bdev only and externally blind"))
  selection <- selected$run$outputs$selected_features$table
  selection_summary <- selected$run$outputs$fit_summary$table
  selected$session <- wb_session_resume(selection_bundle,
    selected$session$durable$workspace)
  promoted <- wb_session_promote_current(selection_bundle, selected$session,
    outputs = c("selected_features", "fit_summary"),
    name_prefix = "atac_variance20k",
    description = "ATAC raw-count non-depth-adjusted variance 20k selection",
    created_by = "scworkbench_Bdev development step 13 blind selection",
    fail_fast = TRUE)
  selection_id <- promoted_id(promoted, "selected_features")
  selection_summary_id <- promoted_id(promoted, "fit_summary")
  promoted$session <- wb_session_discard(promoted$session)
}

stopifnot(nrow(selection_summary) == 1L,
  selection_summary$matrix_id[[1L]] == cfg$inputs$count_matrix,
  selection_summary$method[[1L]] == "variance",
  selection_summary$feature_type[[1L]] == cfg$feature_selection$feature_type,
  as.integer(selection_summary$n_cells[[1L]]) ==
    cfg$inputs$source_bundle$expected_cells,
  as.integer(selection_summary$n_input_features[[1L]]) ==
    cfg$inputs$expected_features,
  as.integer(selection_summary$n_selected_features[[1L]]) ==
    cfg$feature_selection$n_features,
  nrow(selection) == cfg$inputs$expected_features,
  sum(selection$selected) == cfg$feature_selection$n_features,
  !anyDuplicated(selection$feature_id),
  identical(as.integer(selection$rank), seq_len(nrow(selection))))

score_path <- file.path(dirs[["tables"]], "01_variance_peak_scores.tsv.gz")
score_connection <- gzfile(score_path, "wt")
write.table(selection, score_connection, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
close(score_connection)
write_tsv(selection_summary,
  file.path(dirs[["tables"]], "02_variance_selection_summary.tsv"))
wb_close_bundle(selection_bundle); selection_bundle <- NULL

source <- wb_open_bundle(cfg$inputs$source_bundle$path, validate = FALSE,
  read_only = TRUE)
on.exit(try(wb_close_bundle(source), silent = TRUE), add = TRUE)
cells <- as.character(wb_cells(source, "all"))
source_features <- as.character(wb_feature_ids(source, cfg$inputs$count_matrix))
stopifnot(length(cells) == cfg$inputs$source_bundle$expected_cells,
  length(source_features) == cfg$inputs$expected_features,
  !anyDuplicated(cells), !anyDuplicated(source_features),
  setequal(source_features, selection$feature_id))
selected_rank <- as.character(selection$feature_id[selection$selected])
axis <- source_features[source_features %in% selected_rank]
selection_rank <- stats::setNames(selection$rank, selection$feature_id)
selection_score <- stats::setNames(selection$score, selection$feature_id)
axis_table <- data.frame(feature_order = seq_along(axis), feature_id = axis,
  variance_selection_rank = unname(selection_rank[axis]),
  variance_score = unname(selection_score[axis]), stringsAsFactors = FALSE)
stopifnot(length(axis) == cfg$feature_selection$n_features,
  !anyDuplicated(axis), setequal(axis, selected_rank))
axis_path <- file.path(dirs[["tables"]],
  "03_variance20k_source_order_axis.tsv")
write_tsv(axis_table, axis_path)

metadata <- wb_meta(source, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$cluster_field,
  cfg$inputs$sample_field, cfg$inputs$depth_field))
stopifnot(identical(as.character(metadata$cell_id), cells))
pt_family <- readRDS(cfg$inputs$authenticated_l4$pseudotime_path)
pt_take <- as.character(pt_family$key) == cfg$inputs$authenticated_l4$key
pt <- stats::setNames(as.numeric(pt_family$value[pt_take]),
  as.character(pt_family$cell_id[pt_take]))
weight_family <- readRDS(cfg$inputs$authenticated_l4$weights_path)
weight_take <- as.character(weight_family$key) == cfg$inputs$authenticated_l4$key
l4_weight <- stats::setNames(as.numeric(weight_family$value[weight_take]),
  as.character(weight_family$cell_id[weight_take]))
stopifnot(!anyDuplicated(names(l4_weight)), all(names(l4_weight) %in% cells))
metadata$l4_weight <- 0
metadata$l4_weight[match(names(l4_weight), cells)] <- unname(l4_weight)
metadata$native_l4 <- unname(pt[cells])
metadata$cell_role <- "bdev_reference"

slots <- sprintf("external_query_%03d",
  seq_len(cfg$model$external_query_slots))
placeholder <- data.frame(cell_id = slots, stringsAsFactors = FALSE)
placeholder[[cfg$inputs$cluster_field]] <- "external_query_slot"
placeholder[[cfg$inputs$sample_field]] <- "external_query_slot"
placeholder[[cfg$inputs$depth_field]] <- NA_real_
placeholder$l4_weight <- NA_real_
placeholder$native_l4 <- NA_real_
placeholder$cell_role <- "anonymous_external_query_slot"
all_metadata <- rbind(metadata, placeholder)
all_cells <- c(cells, slots)
stopifnot(identical(as.character(all_metadata$cell_id), all_cells))

matrix_id <- "variance20k_bdev_counts"
if (dir.exists(cfg$storage$frozen_bundle)) {
  bundle <- wb_open_bundle(cfg$storage$frozen_bundle, validate = FALSE,
    read_only = FALSE)
  stopifnot(identical(as.character(wb_cells(bundle, "all")), all_cells),
    matrix_id %in% wb_matrices(bundle)$matrix_id,
    identical(as.character(wb_feature_ids(bundle, matrix_id)), axis))
  log_event("Resuming externally blind partial variance20k bundle")
} else {
  bundle <- wb_init_bundle(cfg$storage$frozen_bundle, overwrite = FALSE)
  wb_define_cells(bundle, all_cells)
  wb_attach_metadata(bundle, all_metadata, cell_ids = all_cells,
    source = paste("active_Bdev_metadata_authenticated_native_L4",
      "and_anonymous_slots"),
    created_by = "scworkbench_Bdev development step 13 blind freeze")
  log_event("Attaching source-ordered variance20k full-Bdev counts")
  counts <- wb_slice(source, cfg$inputs$count_matrix,
    features = axis, cells = cells)
  stopifnot(inherits(counts, "sparseMatrix"),
    identical(rownames(counts), axis), identical(colnames(counts), cells))
  wb_attach_matrix(bundle, matrix_id, counts, layer = "counts",
    feature_type = "peak", cell_ids = cells,
    description = "Full-Bdev ATAC variance20k frozen candidate",
    source = "active_03c_Bdev_ATAC_counts_scwbrna_variance20k",
    created_by = "scworkbench_Bdev development step 13 blind freeze")
  rm(counts); invisible(gc())
}
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
wb_close_bundle(source); source <- NULL

artifacts <- wb_art_list(bundle)
model_take <- artifacts$source_op == "scwblatent::peakvi_fit" &
  artifacts$input_matrix_id %in% matrix_id &
  artifacts$type == "model_ref" & artifacts$lifecycle_state == "current"
fit_take <- artifacts$source_op == "scwblatent::peakvi_fit" &
  artifacts$input_matrix_id %in% matrix_id &
  artifacts$type == "obs_table" & artifacts$lifecycle_state == "current"
if (sum(model_take) == 1L && sum(fit_take) == 1L) {
  model_id <- as.character(artifacts$artifact_id[model_take][[1L]])
  fit_summary_id <- as.character(artifacts$artifact_id[fit_take][[1L]])
  model_value <- wb_art_load(bundle, model_id, expected_kind = "model_ref")
  fit_summary <- wb_art_load(bundle, fit_summary_id,
    expected_kind = "obs_table")$table
  reference_space_id <- as.character(model_value$model_space_id)[[1L]]
  log_event("Recovered already-frozen variance20k model")
} else {
  stopifnot(sum(model_take) == 0L, sum(fit_take) == 0L)
  log_event("Fitting externally blind full-Bdev variance20k PeakVI model")
  session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_variance20k_full_bdev_fit"))
  fit <- wb_session_run_op(session, bundle, "scwblatent::peakvi_fit",
    inputs = list(matrix = matrix_id,
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
    notes = paste("Bdev-only external-blind freeze:", cfg$model$candidate))
  reference_space_id <- as.character(fit$run$outputs$space$space_id %||%
    fit$run$outputs$space$allocated_id)[[1L]]
  fit_summary <- fit$run$outputs$fit_summary$table
  fit$session <- wb_session_resume(bundle, fit$session$durable$workspace)
  promoted <- wb_session_promote_current(bundle, fit$session,
    outputs = c("model", "fit_summary"),
    name_prefix = "frozen_full_bdev_variance20k",
    description = "Externally blind frozen full-Bdev variance20k PeakVI model",
    created_by = "scworkbench_Bdev development step 13 blind freeze",
    fail_fast = TRUE)
  model_id <- promoted_id(promoted, "model")
  fit_summary_id <- promoted_id(promoted, "fit_summary")
  promoted$session <- wb_session_discard(promoted$session)
  model_value <- wb_art_load(bundle, model_id, expected_kind = "model_ref")
}

stopifnot(nrow(fit_summary) == 1L,
  as.integer(fit_summary$n_cells) == cfg$inputs$source_bundle$expected_cells,
  as.integer(fit_summary$n_features) == cfg$feature_selection$n_features,
  as.integer(fit_summary$n_latent) == cfg$peakvi$n_latent,
  identical(as.character(model_value$model_space_id), reference_space_id))
model_primary_path <- file.path(model_value$model_tree_root,
  model_value$model_primary_path)
stopifnot(file.exists(model_primary_path))
registry <- data.frame(candidate = cfg$model$candidate,
  features = as.integer(fit_summary$n_features),
  reference_cells = as.integer(fit_summary$n_cells),
  latent_dimensions = as.integer(fit_summary$n_latent),
  latent_sha256 = as.character(fit_summary$latent_sha256),
  reference_space_id = reference_space_id,
  model_artifact_id = model_id,
  fit_summary_artifact_id = fit_summary_id,
  selection_artifact_id = selection_id,
  selection_summary_artifact_id = selection_summary_id,
  model_payload_path = as.character(model_primary_path),
  model_payload_sha256 = sha256(model_primary_path),
  feature_order_sha256 = as.character(
    model_value$model_linkage$feature_order_identity$sha256),
  model_tree_manifest_sha256 = as.character(
    model_value$model_tree_manifest$manifest_identity$sha256),
  stringsAsFactors = FALSE)
registry_path <- file.path(dirs[["tables"]],
  "04_frozen_variance_model_registry.tsv")
write_tsv(registry, registry_path)
wb_close_bundle(bundle); bundle <- NULL

stopifnot(identical(sha256(source_manifest), source_manifest_before))
freeze_receipt <- list(
  schema_version = "bdev.variance20k_peakvi_blind_freeze_receipt.v1",
  status = "variance20k_frozen_external_data_unread",
  run_identity = cfg$run_identity,
  frozen_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  candidate = cfg$model$candidate,
  feature_selection = list(operation = cfg$feature_selection$operation,
    method = cfg$feature_selection$method,
    depth_adjusted = FALSE,
    input_matrix = cfg$inputs$count_matrix,
    selected_features = cfg$feature_selection$n_features,
    selection_artifact_id = selection_id,
    selection_summary_artifact_id = selection_summary_id),
  bdev_cells = cfg$inputs$source_bundle$expected_cells,
  anonymous_external_query_slots = length(slots),
  external_query_use = FALSE,
  immgen_identity_use = FALSE,
  model_registry = list(path = registry_path, sha256 = sha256(registry_path)),
  score_table = list(path = score_path, sha256 = sha256(score_path)),
  feature_axis = list(path = axis_path, sha256 = sha256(axis_path)),
  source_bundle_manifest_sha256 = source_manifest_before,
  config = list(path = cfg_path, sha256 = sha256(cfg_path)),
  runtime = list(sif_sha256 = cfg$runtime$sif_sha256,
    scworkbench_commit = cfg$runtime$scworkbench$source_commit,
    scwbRNA_commit = cfg$runtime$scwbRNA$source_commit,
    scwbRNA_manifest_sha256 = cfg$runtime$scwbRNA$manifest_sha256,
    scwbLatent_commit = cfg$runtime$scwbLatent$source_commit,
    scwbLatent_manifest_sha256 = cfg$runtime$scwbLatent$manifest_sha256),
  current_default_replaced = FALSE,
  implementation_unlocked = FALSE,
  sif_built = FALSE)
freeze_receipt_path <- file.path(dirs[["manifests"]],
  "02_pre_external_freeze_receipt.json")
write_json(freeze_receipt, freeze_receipt_path)
log_event("Variance20k model frozen before external data access")
