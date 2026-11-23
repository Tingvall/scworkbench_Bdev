#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 02_freeze_full_bdev_candidates.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
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
stopifnot(is.null(cfg$base_config), identical(cfg$scope$external_query_use, "prohibited"),
  identical(cfg$scope$immgen_identity_use, "prohibited"))
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_freeze_full_bdev_candidates.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256),
  identical(as.character(packageVersion("scworkbench")), cfg$runtime$scworkbench$version),
  identical(as.character(packageVersion("scwbLatent")), cfg$runtime$scwbLatent$version))
pack_manifest <- system.file("scworkbench", "operation-pack.json", package = "scwbLatent")
stopifnot(nzchar(pack_manifest),
  identical(sha256(pack_manifest), cfg$runtime$scwbLatent$manifest_sha256))
wb_operation_pack_activate(pack_manifest, lib.loc = .libPaths())
stopifnot("scwblatent::peakvi_fit" %in% wb_list_ops())

source_manifest <- file.path(cfg$inputs$source_bundle$path, "manifest.json")
source_manifest_before <- sha256(source_manifest)
stopifnot(identical(source_manifest_before, cfg$inputs$source_bundle$manifest_sha256),
  identical(sha256(cfg$inputs$feature_scores$path), cfg$inputs$feature_scores$sha256),
  identical(sha256(cfg$inputs$authenticated_l4$pseudotime_path),
    cfg$inputs$authenticated_l4$pseudotime_sha256),
  identical(sha256(cfg$inputs$authenticated_l4$weights_path),
    cfg$inputs$authenticated_l4$weights_sha256),
  identical(cfg$inputs$authenticated_l4$orientation, "native"))

source <- wb_open_bundle(cfg$inputs$source_bundle$path, validate = FALSE,
  read_only = TRUE)
on.exit(try(wb_close_bundle(source), silent = TRUE), add = TRUE)
cells <- as.character(wb_cells(source, "all"))
features <- as.character(wb_feature_ids(source, cfg$inputs$count_matrix))
stopifnot(length(cells) == cfg$inputs$source_bundle$expected_cells,
  length(features) == cfg$inputs$expected_features,
  !anyDuplicated(cells), !anyDuplicated(features))
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

score <- read.delim(gzfile(cfg$inputs$feature_scores$path), check.names = FALSE,
  quote = "\"", na.strings = "NA")
stopifnot(nrow(score) == length(features),
  identical(as.character(score$feature_id), features),
  identical(as.integer(score$source_order), seq_along(features)))
candidate_features <- list(
  hybrid20k = features[score$in_hybrid20k],
  cluster20k = features[score$in_cluster20k],
  overall20k = features[score$in_overall20k],
  full191433 = features)
candidates <- unlist(cfg$models$candidates, use.names = FALSE)
stopifnot(identical(names(candidate_features), candidates),
  all(lengths(candidate_features[1:3]) == cfg$models$selected_feature_count),
  length(candidate_features$full191433) == cfg$inputs$expected_features,
  all(vapply(candidate_features, function(x) !anyDuplicated(x), logical(1))))

feature_rows <- do.call(rbind, lapply(candidates, function(candidate)
  data.frame(candidate = candidate, feature_order = seq_along(candidate_features[[candidate]]),
    feature_id = candidate_features[[candidate]], stringsAsFactors = FALSE)))
feature_manifest_path <- file.path(dirs[["tables"]],
  "01_frozen_candidate_feature_axes.tsv.gz")
connection <- gzfile(feature_manifest_path, "wt")
write.table(feature_rows, connection, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
close(connection)

slots <- sprintf("external_query_%03d", seq_len(cfg$models$external_query_slots))
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

if (dir.exists(cfg$storage$frozen_bundle)) {
  bundle <- wb_open_bundle(cfg$storage$frozen_bundle, validate = FALSE,
    read_only = FALSE)
  stopifnot(identical(as.character(wb_cells(bundle, "all")), all_cells),
    all(paste0(candidates, "_bdev_counts") %in% wb_matrices(bundle)$matrix_id))
  log_event("Resuming externally blind partial freeze bundle")
} else {
  bundle <- wb_init_bundle(cfg$storage$frozen_bundle, overwrite = FALSE)
  wb_define_cells(bundle, all_cells)
  wb_attach_metadata(bundle, all_metadata, cell_ids = all_cells,
    source = "active_Bdev_metadata_authenticated_native_L4_and_anonymous_slots",
    created_by = "scworkbench_Bdev development step 11 Bdev-only freeze")
  for (candidate in candidates) {
    log_event("Attaching full-Bdev count matrix for", candidate)
    axis <- candidate_features[[candidate]]
    counts <- wb_slice(source, cfg$inputs$count_matrix, features = axis, cells = cells)
    stopifnot(inherits(counts, "sparseMatrix"), identical(rownames(counts), axis),
      identical(colnames(counts), cells))
    wb_attach_matrix(bundle, paste0(candidate, "_bdev_counts"), counts,
      layer = "counts", feature_type = "peak", cell_ids = cells,
      description = paste("Full-Bdev frozen candidate", candidate),
      source = "active_03c_Bdev_ATAC_counts_Bdev_only_frozen_axis",
      created_by = "scworkbench_Bdev development step 11 Bdev-only freeze")
    rm(counts); invisible(gc())
  }
}
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
wb_close_bundle(source); source <- NULL

model_rows <- list()
for (candidate in candidates) {
  artifacts <- wb_art_list(bundle)
  model_take <- artifacts$source_op == "scwblatent::peakvi_fit" &
    artifacts$input_matrix_id %in% paste0(candidate, "_bdev_counts") &
    artifacts$type == "model_ref" & artifacts$lifecycle_state == "current"
  summary_take <- artifacts$source_op == "scwblatent::peakvi_fit" &
    artifacts$input_matrix_id %in% paste0(candidate, "_bdev_counts") &
    artifacts$type == "obs_table" & artifacts$lifecycle_state == "current"
  if (sum(model_take) == 1L && sum(summary_take) == 1L) {
    model_id <- as.character(artifacts$artifact_id[model_take][[1L]])
    fit_summary_id <- as.character(artifacts$artifact_id[summary_take][[1L]])
    model_value <- wb_art_load(bundle, model_id, expected_kind = "model_ref")
    fit_summary <- wb_art_load(bundle, fit_summary_id,
      expected_kind = "obs_table")$table
    reference_space_id <- as.character(model_value$model_space_id)[[1L]]
    log_event("Recovered already-frozen full-Bdev candidate", candidate)
  } else {
    stopifnot(sum(model_take) == 0L, sum(summary_take) == 0L)
    log_event("Fitting frozen full-Bdev candidate", candidate)
    session <- wb_session_create(bundle, cfg$storage$sessions,
      paste0(cfg$run_identity, "_", candidate, "_full_bdev_fit"))
    fit <- wb_session_run_op(session, bundle, "scwblatent::peakvi_fit",
      inputs = list(matrix = paste0(candidate, "_bdev_counts"),
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
      notes = paste("Bdev-only freeze before external bulk evaluation:", candidate))
    reference_space_id <- as.character(fit$run$outputs$space$space_id %||%
      fit$run$outputs$space$allocated_id)[[1L]]
    fit_summary <- fit$run$outputs$fit_summary$table
    fit$session <- wb_session_resume(bundle, fit$session$durable$workspace)
    promoted <- wb_session_promote_current(bundle, fit$session,
      outputs = c("model", "fit_summary"),
      name_prefix = paste0("frozen_full_bdev_", candidate),
      description = paste("Externally blind frozen full-Bdev", candidate, "PeakVI model"),
      created_by = "scworkbench_Bdev development step 11 Bdev-only freeze",
      fail_fast = TRUE)
    model_id <- promoted_id(promoted, "model")
    fit_summary_id <- promoted_id(promoted, "fit_summary")
    promoted$session <- wb_session_discard(promoted$session)
    model_value <- wb_art_load(bundle, model_id, expected_kind = "model_ref")
  }
  stopifnot(nrow(fit_summary) == 1L,
    as.integer(fit_summary$n_cells) == cfg$inputs$source_bundle$expected_cells,
    as.integer(fit_summary$n_features) == length(candidate_features[[candidate]]),
    as.integer(fit_summary$n_latent) == cfg$peakvi$n_latent)
  model_primary_path <- file.path(model_value$model_tree_root,
    model_value$model_primary_path)
  stopifnot(file.exists(model_primary_path),
    identical(as.character(model_value$model_space_id), reference_space_id))
  model_rows[[candidate]] <- data.frame(candidate = candidate,
    features = as.integer(fit_summary$n_features),
    reference_cells = as.integer(fit_summary$n_cells),
    latent_dimensions = as.integer(fit_summary$n_latent),
    latent_sha256 = as.character(fit_summary$latent_sha256),
    reference_space_id = reference_space_id,
    model_artifact_id = model_id,
    fit_summary_artifact_id = fit_summary_id,
    model_payload_path = as.character(model_primary_path),
    model_payload_sha256 = sha256(model_primary_path),
    feature_order_sha256 = as.character(
      model_value$model_linkage$feature_order_identity$sha256),
    model_tree_manifest_sha256 = as.character(
      model_value$model_tree_manifest$manifest_identity$sha256),
    stringsAsFactors = FALSE)
  log_event("Frozen", candidate, "model payload",
    model_rows[[candidate]]$model_payload_sha256)
  rm(list = intersect(c("fit", "promoted", "model_value"), ls()))
  invisible(gc())
}
registry <- do.call(rbind, model_rows); rownames(registry) <- NULL
registry_path <- file.path(dirs[["tables"]], "02_frozen_model_registry.tsv")
write_tsv(registry, registry_path)
wb_close_bundle(bundle); bundle <- NULL

stopifnot(identical(sha256(source_manifest), source_manifest_before))
freeze_receipt <- list(
  schema_version = "bdev.full_peakvi_four_model_freeze_receipt.v1",
  status = "four_models_frozen_external_data_unread",
  run_identity = cfg$run_identity,
  frozen_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  candidates = candidates,
  bdev_cells_per_model = cfg$inputs$source_bundle$expected_cells,
  anonymous_external_query_slots = length(slots),
  external_query_use = FALSE,
  immgen_identity_use = FALSE,
  model_registry = list(path = registry_path, sha256 = sha256(registry_path)),
  feature_manifest = list(path = feature_manifest_path,
    sha256 = sha256(feature_manifest_path)),
  source_bundle_manifest_sha256 = source_manifest_before,
  config = list(path = cfg_path, sha256 = sha256(cfg_path)),
  runtime = list(sif_sha256 = cfg$runtime$sif_sha256,
    scworkbench_commit = cfg$runtime$scworkbench$source_commit,
    scwbLatent_commit = cfg$runtime$scwbLatent$source_commit,
    scwbLatent_manifest_sha256 = cfg$runtime$scwbLatent$manifest_sha256),
  automatic_model_promotion = FALSE,
  implementation_unlocked = FALSE,
  sif_built = FALSE)
freeze_receipt_path <- file.path(dirs[["manifests"]], "02_pre_external_freeze_receipt.json")
write_json(freeze_receipt, freeze_receipt_path)
log_event("All four full-Bdev models frozen before external data access")
