#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 03_project_immgen_all_candidates.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbLatent)
  library(Matrix)
  library(RANN)
  library(yaml)
  library(jsonlite)
  library(digest)
})

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
read_config <- function(path) {
  current <- yaml::read_yaml(path)
  base <- current$base_config
  current$base_config <- NULL
  if (is.null(base)) return(current)
  utils::modifyList(read_config(normalizePath(file.path(dirname(path), base))),
    current, keep.null = TRUE)
}
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
as_space_matrix <- function(value) {
  if (is.data.frame(value)) {
    ids <- as.character(value$cell_id)
    value <- as.matrix(value[, setdiff(names(value), "cell_id"), drop = FALSE])
    rownames(value) <- ids
  } else value <- as.matrix(value)
  value
}
normalized_entropy <- function(groups, weights) {
  totals <- tapply(weights, groups, sum)
  p <- totals / sum(totals)
  if (length(p) <= 1L) return(0)
  -sum(p * log(p)) / log(length(p))
}
safe_cor <- function(x, y) {
  if (length(x) < 3L || !is.finite(stats::sd(x)) || stats::sd(x) == 0 ||
      !is.finite(stats::sd(y)) || stats::sd(y) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x, y, method = "spearman"))
}
promoted_id <- function(promoted, output) {
  take <- promoted$promoted$status == "ok" & promoted$promoted$output == output
  stopifnot(sum(take) == 1L)
  as.character(promoted$promoted$artifact_id[take][[1L]])
}

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- read_config(cfg_path)
stopifnot(identical(cfg$scope$external_query_use, "evaluation_only_after_freeze"),
  identical(cfg$scope$automatic_model_promotion, "prohibited"))
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "02_immgen_four_model_projection.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

pack_manifest <- system.file("scworkbench", "operation-pack.json", package = "scwbLatent")
stopifnot(nzchar(pack_manifest),
  identical(sha256(pack_manifest), cfg$runtime$scwbLatent$manifest_sha256),
  identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256))
wb_operation_pack_activate(pack_manifest, lib.loc = .libPaths())
stopifnot("scwblatent::peakvi_project" %in% wb_list_ops())

freeze_path <- file.path(dirs[["manifests"]], "02_pre_external_freeze_receipt.json")
registry_path <- file.path(dirs[["tables"]], "02_frozen_model_registry.tsv")
feature_manifest_path <- file.path(dirs[["tables"]],
  "01_frozen_candidate_feature_axes.tsv.gz")
stopifnot(file.exists(freeze_path), file.exists(registry_path),
  file.exists(feature_manifest_path))
freeze <- jsonlite::read_json(freeze_path, simplifyVector = TRUE)
stopifnot(identical(freeze$status, "four_models_frozen_external_data_unread"),
  isFALSE(freeze$external_query_use), isFALSE(freeze$immgen_identity_use),
  identical(sha256(registry_path), freeze$model_registry$sha256),
  identical(sha256(feature_manifest_path), freeze$feature_manifest$sha256))
registry <- read.delim(registry_path, check.names = FALSE, quote = "\"")
candidate_features <- read.delim(gzfile(feature_manifest_path),
  check.names = FALSE, quote = "\"")
candidates <- unlist(cfg$models$candidates, use.names = FALSE)
stopifnot(identical(as.character(registry$candidate), candidates),
  setequal(unique(candidate_features$candidate), candidates))

source_manifest <- file.path(cfg$inputs$source_bundle$path, "manifest.json")
source_manifest_before <- sha256(source_manifest)
query_manifest <- file.path(cfg$external_evaluation$query_bundle$path, "manifest.json")
stopifnot(identical(source_manifest_before, cfg$inputs$source_bundle$manifest_sha256),
  identical(sha256(query_manifest),
    cfg$external_evaluation$query_bundle$manifest_sha256),
  identical(sha256(cfg$external_evaluation$library_manifest$path),
    cfg$external_evaluation$library_manifest$sha256))

libraries <- read.delim(cfg$external_evaluation$library_manifest$path,
  check.names = FALSE, quote = "")
stopifnot(nrow(libraries) == cfg$external_evaluation$query_bundle$expected_libraries,
  !anyDuplicated(libraries$gsm))
slots <- sprintf("external_query_%03d", seq_len(nrow(libraries)))
slot_map <- data.frame(slot_id = slots, gsm = as.character(libraries$gsm),
  population = as.character(libraries$population), lineage = as.character(libraries$lineage),
  stringsAsFactors = FALSE)
write_tsv(slot_map, file.path(dirs[["tables"]], "03_external_query_slot_map.tsv"))

bundle <- wb_open_bundle(cfg$storage$frozen_bundle, validate = FALSE,
  read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
all_cells <- as.character(wb_cells(bundle, "all"))
reference_cells <- all_cells[!grepl("^external_query_[0-9]{3}$", all_cells)]
stopifnot(length(reference_cells) == cfg$inputs$source_bundle$expected_cells,
  identical(tail(all_cells, length(slots)), slots))

model_payload_before <- setNames(character(length(candidates)), candidates)
for (candidate in candidates) {
  row <- registry[registry$candidate == candidate, , drop = FALSE]
  stopifnot(nrow(row) == 1L, file.exists(row$model_payload_path))
  model_payload_before[[candidate]] <- sha256(row$model_payload_path)
  stopifnot(identical(model_payload_before[[candidate]], row$model_payload_sha256))
  model_value <- wb_art_load(bundle, row$model_artifact_id,
    expected_kind = "model_ref")
  stopifnot(identical(as.character(model_value$model_tree_manifest$manifest_identity$sha256),
    as.character(row$model_tree_manifest_sha256)))
}
log_event("Validated all four frozen model payloads before opening ImmGen counts")

query_source <- wb_open_bundle(cfg$external_evaluation$query_bundle$path,
  validate = FALSE, read_only = TRUE)
on.exit(try(wb_close_bundle(query_source), silent = TRUE), add = TRUE)
query_axis <- as.character(wb_feature_ids(query_source,
  cfg$external_evaluation$query_bundle$matrix))
stopifnot(length(query_axis) == cfg$inputs$expected_features)
full_query <- wb_slice(query_source, cfg$external_evaluation$query_bundle$matrix,
  features = query_axis, cells = as.character(libraries$gsm))
stopifnot(inherits(full_query, "sparseMatrix"),
  identical(rownames(full_query), query_axis),
  identical(colnames(full_query), as.character(libraries$gsm)))
colnames(full_query) <- slots
common_total_counts <- Matrix::colSums(full_query)

query_counts_by_candidate <- list()
for (candidate in candidates) {
  axis <- as.character(candidate_features$feature_id[
    candidate_features$candidate == candidate])
  stopifnot(length(axis) == registry$features[registry$candidate == candidate],
    all(axis %in% query_axis))
  counts <- full_query[axis, , drop = FALSE]
  query_counts_by_candidate[[candidate]] <- counts
  wb_attach_matrix(bundle, paste0(candidate, "_immgen_counts"), counts,
    layer = "counts", feature_type = "peak", cell_ids = slots,
    description = paste("Identical frozen ImmGen bulk queries on", candidate, "axis"),
    source = "development_step_03_exact_Bdev_peak_counts_external_evaluation_only",
    created_by = "scworkbench_Bdev development step 11 post-freeze benchmark")
}
wb_close_bundle(query_source); query_source <- NULL
rm(full_query); invisible(gc())

reference_metadata <- wb_meta(bundle, cells = reference_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$cluster_field,
  cfg$inputs$sample_field, cfg$inputs$depth_field, "l4_weight", "native_l4"))
mapping_rows <- list(); neighbor_rows <- list(); metric_rows <- list()
for (candidate in candidates) {
  log_event("Projecting identical ImmGen bulk queries through frozen", candidate)
  row <- registry[registry$candidate == candidate, , drop = FALSE]
  model_ref <- wb_ref_artifact(as.character(row$model_artifact_id), "model")
  query_space_id <- paste0("frozen_", candidate, "_immgen_queries")
  session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_", candidate, "_immgen_project"))
  projected <- wb_session_run_op(session, bundle, "scwblatent::peakvi_project",
    inputs = list(model = model_ref,
      query_matrix = paste0(candidate, "_immgen_counts"),
      space_id_out = query_space_id,
      chunk_size = as.integer(cfg$peakvi$chunk_size), batch_key = NULL,
      use_gpu = FALSE,
      pickle_execution_trust = cfg$peakvi$pickle_execution_trust,
      cell_set = "all", min_cells = 1L),
    notes = paste("Post-freeze identical ImmGen bulk projection:", candidate))
  projection_summary <- projected$run$outputs$projection_summary$table
  stopifnot(nrow(projection_summary) == 1L,
    projection_summary$transformation[[1L]] == "transfer_without_refit",
    as.integer(projection_summary$n_latent) == cfg$peakvi$n_latent)
  projected$session <- wb_session_resume(bundle,
    projected$session$durable$workspace)
  promoted <- wb_session_promote_current(bundle, projected$session,
    outputs = "projection_summary",
    name_prefix = paste0("frozen_", candidate, "_immgen_projection"),
    description = paste("Post-freeze ImmGen projection for", candidate),
    created_by = "scworkbench_Bdev development step 11 post-freeze benchmark",
    fail_fast = TRUE)
  projection_artifact_id <- promoted_id(promoted, "projection_summary")
  promoted$session <- wb_session_discard(promoted$session)

  reference_latent <- as_space_matrix(wb_space(bundle,
    as.character(row$reference_space_id), cells = reference_cells))
  query_latent <- as_space_matrix(wb_space(bundle, query_space_id, cells = slots))
  stopifnot(identical(rownames(reference_latent), reference_cells),
    identical(rownames(query_latent), slots),
    ncol(reference_latent) == cfg$peakvi$n_latent,
    ncol(query_latent) == cfg$peakvi$n_latent)
  k <- min(as.integer(cfg$external_evaluation$neighbor_k), nrow(reference_latent))
  nn <- RANN::nn2(reference_latent, query_latent, k = k)
  reference_nn <- RANN::nn2(reference_latent, reference_latent, k = 2L)$nn.dists[, 2L]
  reference_distance_median <- median(reference_nn[is.finite(reference_nn)])
  ref_cluster <- as.character(reference_metadata[[cfg$inputs$cluster_field]])
  ref_sample <- as.character(reference_metadata[[cfg$inputs$sample_field]])
  ref_l4_weight <- as.numeric(reference_metadata$l4_weight)
  ref_native_l4 <- as.numeric(reference_metadata$native_l4)

  per_query <- vector("list", length(slots))
  per_neighbor <- vector("list", length(slots))
  counts <- query_counts_by_candidate[[candidate]]
  candidate_counts <- Matrix::colSums(counts)
  saturation <- Matrix::colSums(counts != 0) / nrow(counts)
  for (i in seq_along(slots)) {
    index <- nn$nn.idx[i, ]
    distance <- nn$nn.dists[i, ]
    weights <- 1 / pmax(distance, 1e-8)
    cluster_weight <- tapply(weights, ref_cluster[index], sum)
    predicted_cluster <- names(cluster_weight)[order(-cluster_weight,
      names(cluster_weight), method = "radix")][[1L]]
    membership <- weighted.mean(ref_l4_weight[index], weights)
    l4_local <- which(ref_l4_weight[index] >=
      cfg$external_evaluation$l4_weight_positive_threshold &
      is.finite(ref_native_l4[index]))
    provisional_pt <- if (length(l4_local) >=
      cfg$external_evaluation$minimum_l4_neighbors)
      weighted.mean(ref_native_l4[index[l4_local]], weights[l4_local]) else NA_real_
    per_query[[i]] <- data.frame(candidate = candidate, slot_id = slots[[i]],
      gsm = slot_map$gsm[[i]], population = slot_map$population[[i]],
      lineage = slot_map$lineage[[i]], common_total_peak_counts = common_total_counts[[i]],
      candidate_peak_counts = candidate_counts[[i]], peak_saturation = saturation[[i]],
      predicted_bdev_cluster = predicted_cluster,
      nearest_reference_distance = min(distance),
      nearest_distance_over_reference_median = min(distance) / reference_distance_median,
      neighbor_sample_entropy = normalized_entropy(ref_sample[index], weights),
      neighbor_l4_weight = membership,
      provisional_native_l4 = provisional_pt,
      trajectory_output_status = cfg$external_evaluation$trajectory_output_status,
      projection_artifact_id = projection_artifact_id,
      stringsAsFactors = FALSE)
    per_neighbor[[i]] <- data.frame(candidate = candidate, slot_id = slots[[i]],
      gsm = slot_map$gsm[[i]], rank = seq_len(k),
      reference_cell_id = reference_cells[index], distance = distance,
      reference_cluster = ref_cluster[index], reference_sample = ref_sample[index],
      reference_l4_weight = ref_l4_weight[index], stringsAsFactors = FALSE)
  }
  per_query <- do.call(rbind, per_query)
  per_neighbor <- do.call(rbind, per_neighbor)
  mapping_rows[[candidate]] <- per_query
  neighbor_rows[[candidate]] <- per_neighbor

  query_dist <- as.matrix(stats::dist(query_latent))
  pairs <- which(upper.tri(query_dist), arr.ind = TRUE)
  same_population <- slot_map$population[pairs[, 1L]] == slot_map$population[pairs[, 2L]]
  within <- query_dist[pairs[same_population, , drop = FALSE]]
  between <- query_dist[pairs[!same_population, , drop = FALSE]]
  replicate_pairs <- pairs[same_population, , drop = FALSE]
  replicate_agreement <- if (nrow(replicate_pairs)) mean(
    per_query$predicted_bdev_cluster[replicate_pairs[, 1L]] ==
      per_query$predicted_bdev_cluster[replicate_pairs[, 2L]]) else NA_real_
  depth_cor <- apply(query_latent, 2L, function(x)
    safe_cor(x, log1p(common_total_counts)))
  saturation_cor <- apply(query_latent, 2L, function(x)
    safe_cor(x, saturation))
  metric_rows[[candidate]] <- data.frame(candidate = candidate,
    features = nrow(counts), query_libraries = nrow(query_latent),
    within_population_to_between_population_latent_distance_ratio =
      median(within) / median(between),
    replicate_predicted_cluster_agreement = replicate_agreement,
    unique_reference_neighbors = length(unique(per_neighbor$reference_cell_id)),
    maximum_reference_hub_fraction = max(table(per_neighbor$reference_cell_id)) /
      nrow(per_neighbor),
    median_nearest_distance_over_reference_median =
      median(per_query$nearest_distance_over_reference_median),
    median_neighbor_sample_entropy = median(per_query$neighbor_sample_entropy),
    maximum_absolute_depth_latent_spearman = max(abs(depth_cor), na.rm = TRUE),
    maximum_absolute_saturation_latent_spearman =
      max(abs(saturation_cor), na.rm = TRUE),
    median_neighbor_l4_weight = median(per_query$neighbor_l4_weight),
    finite_provisional_l4_libraries = sum(is.finite(per_query$provisional_native_l4)),
    model_payload_sha256_before = model_payload_before[[candidate]],
    stringsAsFactors = FALSE)
  rm(reference_latent, query_latent, nn, reference_nn, counts); invisible(gc())
  log_event("Completed external projection for", candidate)
}

mappings <- do.call(rbind, mapping_rows); rownames(mappings) <- NULL
neighbors <- do.call(rbind, neighbor_rows); rownames(neighbors) <- NULL
metrics <- do.call(rbind, metric_rows); rownames(metrics) <- NULL
write_tsv(mappings, file.path(dirs[["tables"]], "04_immgen_candidate_mappings.tsv"))
neighbor_connection <- gzfile(file.path(dirs[["tables"]],
  "05_immgen_candidate_neighbors.tsv.gz"), "wt")
write.table(neighbors, neighbor_connection, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
close(neighbor_connection)
write_tsv(metrics, file.path(dirs[["tables"]], "06_external_candidate_metrics.tsv"))

population_summary <- do.call(rbind, lapply(split(mappings,
  interaction(mappings$candidate, mappings$population, drop = TRUE)), function(x)
  data.frame(candidate = x$candidate[[1L]], population = x$population[[1L]],
    libraries = nrow(x), predicted_clusters = paste(sort(unique(x$predicted_bdev_cluster)),
      collapse = ";"), median_neighbor_l4_weight = median(x$neighbor_l4_weight),
    median_provisional_native_l4 = if (any(is.finite(x$provisional_native_l4)))
      median(x$provisional_native_l4, na.rm = TRUE) else NA_real_,
    median_support_distance_ratio = median(x$nearest_distance_over_reference_median),
    stringsAsFactors = FALSE)))
write_tsv(population_summary, file.path(dirs[["tables"]],
  "07_immgen_population_candidate_summary.tsv"))

model_payload_after <- setNames(character(length(candidates)), candidates)
for (candidate in candidates) {
  row <- registry[registry$candidate == candidate, , drop = FALSE]
  model_payload_after[[candidate]] <- sha256(row$model_payload_path)
  stopifnot(identical(model_payload_after[[candidate]],
    model_payload_before[[candidate]]))
}
metrics$model_payload_sha256_after <- unname(model_payload_after[metrics$candidate])
write_tsv(metrics, file.path(dirs[["tables"]], "06_external_candidate_metrics.tsv"))

source_manifest_after <- sha256(source_manifest)
query_manifest_after <- sha256(query_manifest)
stopifnot(identical(source_manifest_after, source_manifest_before),
  identical(query_manifest_after, cfg$external_evaluation$query_bundle$manifest_sha256))
finalized <- wb_finalize_recipient_bundle(bundle)
stopifnot(isTRUE(finalized$ready))
wb_close_bundle(bundle); bundle <- NULL

summary <- list(
  schema_version = "bdev.frozen_peakvi_immgen_four_model_benchmark_summary.v1",
  status = "complete_external_comparison_no_automatic_promotion",
  run_identity = cfg$run_identity,
  candidates = candidates,
  full_bdev_cells_per_model = cfg$inputs$source_bundle$expected_cells,
  immgen_libraries = nrow(libraries),
  same_query_matrix_for_all_models = TRUE,
  projection_without_refit = TRUE,
  models_frozen_before_external_data_access = TRUE,
  model_payloads_unchanged = TRUE,
  whole_reference_mapping_first = TRUE,
  trajectory_outputs = "diagnostic_unaccepted_second_stage",
  automatic_model_promotion = FALSE,
  interpretation = paste("Use external metrics to compare the four already-frozen",
    "models. No feature selection, fitting, refitting, or tuning used ImmGen."),
  source_bundle_manifest_unchanged = TRUE,
  query_bundle_manifest_unchanged = TRUE,
  implementation_unlocked = FALSE,
  sif_built = FALSE)
write_json(summary, file.path(dirs[["reports"]],
  "01_external_four_model_benchmark_summary.json"))
write_json(list(
  schema_version = "bdev.frozen_peakvi_immgen_four_model_execution.v1",
  status = "complete", run_identity = cfg$run_identity,
  external_config_path = cfg_path,
  external_config_sha256 = sha256(cfg_path),
  freeze_receipt_path = freeze_path,
  freeze_receipt_sha256 = sha256(freeze_path),
  model_payload_sha256_before = as.list(model_payload_before),
  model_payload_sha256_after = as.list(model_payload_after),
  frozen_bundle_manifest_sha256 = sha256(file.path(cfg$storage$frozen_bundle,
    "manifest.json")),
  source_manifest_unchanged = TRUE,
  query_manifest_unchanged = TRUE,
  automatic_model_promotion = FALSE),
  file.path(dirs[["manifests"]], "03_external_execution_receipt.json"))
log_event("Four-model ImmGen benchmark complete; no candidate auto-promoted")
