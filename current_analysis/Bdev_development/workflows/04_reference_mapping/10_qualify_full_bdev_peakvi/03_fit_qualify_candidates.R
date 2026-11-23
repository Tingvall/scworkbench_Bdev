#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 03_fit_qualify_candidates.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbLatent)
  library(Matrix)
  library(RANN)
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
sink(file.path(dirs[["logs"]], "02_candidate_qualification.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
promoted_ref <- function(promoted, output, kind) {
  selected <- promoted$promoted$status == "ok" & promoted$promoted$output == output
  stopifnot(sum(selected) == 1L)
  wb_ref_artifact(promoted$promoted$artifact_id[selected][[1L]], kind)
}
as_space_matrix <- function(value) {
  if (is.data.frame(value)) {
    ids <- as.character(value$cell_id)
    value <- as.matrix(value[, setdiff(names(value), "cell_id"), drop = FALSE])
    rownames(value) <- ids
  } else value <- as.matrix(value)
  value
}
auc_rank <- function(truth, score) {
  keep <- !is.na(truth) & is.finite(score)
  truth <- truth[keep]; score <- score[keep]
  n1 <- sum(truth); n0 <- sum(!truth)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
calibrate_threshold <- function(truth, score, minimum_specificity) {
  candidates <- sort(unique(c(0, score, 1)), decreasing = TRUE)
  audit <- do.call(rbind, lapply(candidates, function(threshold) {
    call <- score >= threshold
    data.frame(threshold = threshold,
      sensitivity = if (any(truth)) mean(call[truth]) else NA_real_,
      specificity = if (any(!truth)) mean(!call[!truth]) else NA_real_)
  }))
  eligible <- audit[is.finite(audit$specificity) &
    audit$specificity >= minimum_specificity, , drop = FALSE]
  stopifnot(nrow(eligible) > 0L)
  eligible <- eligible[order(-eligible$sensitivity, eligible$threshold), , drop = FALSE]
  list(selected = eligible[1L, , drop = FALSE], audit = audit)
}

stopifnot(identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256),
  identical(as.character(packageVersion("scworkbench")), cfg$runtime$scworkbench$version),
  identical(as.character(packageVersion("scwbLatent")), cfg$runtime$scwbLatent$version))
manifest <- system.file("scworkbench", "operation-pack.json", package = "scwbLatent")
stopifnot(nzchar(manifest), identical(sha256(manifest),
  cfg$runtime$scwbLatent$manifest_sha256))
wb_operation_pack_activate(manifest, lib.loc = .libPaths())
stopifnot(all(c("scwblatent::peakvi_fit", "scwblatent::peakvi_project") %in%
  wb_list_ops()))

source_manifest <- file.path(cfg$inputs$source_bundle$path, "manifest.json")
source_manifest_before <- sha256(source_manifest)
stopifnot(identical(source_manifest_before,
  cfg$inputs$source_bundle$manifest_sha256))

bundle <- wb_open_bundle(cfg$storage$candidate_bundle, validate = FALSE,
  read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
all_cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = all_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$cluster_field,
  cfg$inputs$sample_field, cfg$inputs$depth_field, "partition", "l4_weight",
  "native_l4"))
reference_cells <- as.character(metadata$cell_id[metadata$partition == "reference"])
query_cells <- as.character(metadata$cell_id[metadata$partition != "reference"])
query_type <- ifelse(metadata$partition[match(query_cells, metadata$cell_id)] ==
  "held_out", "held_out_single_cell", "bdev_pseudobulk")

all_predictions <- list(); all_neighbors <- list(); metric_rows <- list()
hybrid_threshold <- NULL
for (candidate in unlist(cfg$qualification$candidates, use.names = FALSE)) {
  log_event("Fitting Bdev-only candidate", candidate)
  fit_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_", candidate, "_fit"))
  fit <- wb_session_run_op(fit_session, bundle, "scwblatent::peakvi_fit",
    inputs = list(matrix = paste0(candidate, "_reference_counts"),
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
    notes = paste("Bdev-only full-reference feature candidate", candidate))
  reference_space_id <- as.character(fit$run$outputs$space$space_id %||%
    fit$run$outputs$space$allocated_id)[[1L]]
  fit_summary <- fit$run$outputs$fit_summary$table
  fit$session <- wb_session_resume(bundle, fit$session$durable$workspace)
  promoted <- wb_session_promote_current(bundle, fit$session,
    outputs = c("model", "fit_summary"),
    name_prefix = paste0("full_bdev_", candidate, "_qualification"),
    description = paste(candidate, "Bdev-only PeakVI qualification model"),
    created_by = "scworkbench_Bdev development", fail_fast = TRUE)
  model_ref <- promoted_ref(promoted, "model", "model")
  promoted$session <- wb_session_discard(promoted$session)

  query_space_id <- paste0("full_bdev_", candidate, "_heldout_projection")
  project_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_", candidate, "_project"))
  project <- wb_session_run_op(project_session, bundle,
    "scwblatent::peakvi_project", inputs = list(model = model_ref,
      query_matrix = paste0(candidate, "_query_counts"),
      space_id_out = query_space_id,
      chunk_size = as.integer(cfg$peakvi$chunk_size), batch_key = NULL,
      use_gpu = FALSE,
      pickle_execution_trust = cfg$peakvi$pickle_execution_trust,
      cell_set = "all", min_cells = 1L),
    notes = paste("Disjoint Bdev qualification projection", candidate))
  projection_summary <- project$run$outputs$projection_summary$table
  project$session <- wb_session_resume(bundle, project$session$durable$workspace)
  projected <- wb_session_promote_current(bundle, project$session,
    outputs = "projection_summary",
    name_prefix = paste0("full_bdev_", candidate, "_projection"),
    description = paste(candidate, "held-out Bdev projection"),
    created_by = "scworkbench_Bdev development", fail_fast = TRUE)
  projected$session <- wb_session_discard(projected$session)

  wb_close_bundle(bundle)
  bundle <- wb_open_bundle(cfg$storage$candidate_bundle, validate = FALSE,
    read_only = FALSE)
  reference_latent <- as_space_matrix(wb_space(bundle, reference_space_id,
    cells = reference_cells))
  query_latent <- as_space_matrix(wb_space(bundle, query_space_id,
    cells = query_cells))
  stopifnot(identical(rownames(reference_latent), reference_cells),
    identical(rownames(query_latent), query_cells),
    ncol(reference_latent) == cfg$peakvi$n_latent,
    ncol(query_latent) == cfg$peakvi$n_latent)

  k <- min(as.integer(cfg$qualification$neighbor_k), nrow(reference_latent))
  nn <- RANN::nn2(reference_latent, query_latent, k = k)
  ref_cluster <- as.character(metadata[[cfg$inputs$cluster_field]][
    match(reference_cells, metadata$cell_id)])
  ref_weight <- as.numeric(metadata$l4_weight[match(reference_cells, metadata$cell_id)])
  ref_pt <- as.numeric(metadata$native_l4[match(reference_cells, metadata$cell_id)])
  truth_cluster <- as.character(metadata[[cfg$inputs$cluster_field]][
    match(query_cells, metadata$cell_id)])
  truth_weight <- as.numeric(metadata$l4_weight[match(query_cells, metadata$cell_id)])
  truth_pt <- as.numeric(metadata$native_l4[match(query_cells, metadata$cell_id)])
  query_depth <- as.numeric(metadata[[cfg$inputs$depth_field]][
    match(query_cells, metadata$cell_id)])

  prediction <- vector("list", length(query_cells))
  neighbor_rows <- vector("list", length(query_cells))
  for (i in seq_along(query_cells)) {
    index <- nn$nn.idx[i, ]
    distance <- nn$nn.dists[i, ]
    neighbor_weight <- 1 / pmax(distance, 1e-8)
    label_weight <- tapply(neighbor_weight, ref_cluster[index], sum)
    predicted_cluster <- names(label_weight)[order(-label_weight,
      names(label_weight), method = "radix")][[1L]]
    membership_score <- weighted.mean(ref_weight[index], neighbor_weight)
    l4_neighbors <- index[ref_weight[index] >=
      cfg$qualification$l4_gate$weight_positive_threshold & is.finite(ref_pt[index])]
    transferred_pt <- if (length(l4_neighbors) >=
      cfg$qualification$l4_gate$minimum_l4_neighbors_for_pseudotime) {
      local <- match(l4_neighbors, index)
      weighted.mean(ref_pt[l4_neighbors], neighbor_weight[local])
    } else NA_real_
    prediction[[i]] <- data.frame(candidate = candidate, query_cell_id = query_cells[[i]],
      query_type = query_type[[i]], true_cluster = truth_cluster[[i]],
      predicted_cluster = predicted_cluster, nearest_distance = min(distance),
      l4_membership_score = membership_score, true_l4_weight = truth_weight[[i]],
      true_native_l4 = truth_pt[[i]], transferred_native_l4 = transferred_pt,
      stringsAsFactors = FALSE)
    neighbor_rows[[i]] <- data.frame(candidate = candidate,
      query_cell_id = query_cells[[i]], rank = seq_len(k),
      reference_cell_id = reference_cells[index], distance = distance,
      reference_cluster = ref_cluster[index], reference_l4_weight = ref_weight[index],
      stringsAsFactors = FALSE)
  }
  prediction <- do.call(rbind, prediction)
  neighbors <- do.call(rbind, neighbor_rows)
  actual <- prediction$query_type == "held_out_single_cell"
  pseudo <- prediction$query_type == "bdev_pseudobulk"
  recalls <- tapply(prediction$predicted_cluster[actual] ==
    prediction$true_cluster[actual], prediction$true_cluster[actual], mean)
  latent_cor <- apply(query_latent[actual, , drop = FALSE], 2L, function(value)
    suppressWarnings(cor(value, log1p(query_depth[actual]), method = "spearman")))
  l4_truth <- prediction$true_l4_weight[actual] >=
    cfg$qualification$l4_gate$weight_positive_threshold
  threshold <- calibrate_threshold(l4_truth,
    prediction$l4_membership_score[actual],
    cfg$qualification$l4_gate$minimum_specificity)
  l4_call <- prediction$l4_membership_score[actual] >= threshold$selected$threshold[[1L]]
  pt_keep <- l4_truth & l4_call & is.finite(prediction$true_native_l4[actual]) &
    is.finite(prediction$transferred_native_l4[actual])
  metric_rows[[candidate]] <- data.frame(candidate = candidate,
    features = as.integer(fit_summary$n_features),
    reference_cells = as.integer(fit_summary$n_cells),
    held_out_cells = sum(actual), pseudobulks = sum(pseudo),
    single_cell_accuracy = mean(prediction$predicted_cluster[actual] ==
      prediction$true_cluster[actual]),
    macro_cluster_accuracy = mean(recalls),
    pseudobulk_accuracy = if (any(pseudo)) mean(prediction$predicted_cluster[pseudo] ==
      prediction$true_cluster[pseudo]) else NA_real_,
    max_abs_depth_latent_spearman = max(abs(latent_cor), na.rm = TRUE),
    neighbor_reference_cells_used = length(unique(neighbors$reference_cell_id)),
    maximum_reference_hub_fraction = max(table(neighbors$reference_cell_id)) / nrow(neighbors),
    l4_membership_auc = auc_rank(l4_truth, prediction$l4_membership_score[actual]),
    l4_threshold = threshold$selected$threshold[[1L]],
    l4_sensitivity = threshold$selected$sensitivity[[1L]],
    l4_specificity = threshold$selected$specificity[[1L]],
    l4_pseudotime_mae = if (any(pt_keep)) mean(abs(
      prediction$transferred_native_l4[actual][pt_keep] -
        prediction$true_native_l4[actual][pt_keep])) else NA_real_,
    l4_pseudotime_spearman = if (sum(pt_keep) >= 3L) suppressWarnings(cor(
      prediction$transferred_native_l4[actual][pt_keep],
      prediction$true_native_l4[actual][pt_keep], method = "spearman")) else NA_real_,
    latent_sha256 = as.character(fit_summary$latent_sha256), stringsAsFactors = FALSE)
  all_predictions[[candidate]] <- prediction
  all_neighbors[[candidate]] <- neighbors
  if (candidate == cfg$qualification$primary_candidate)
    hybrid_threshold <- list(selected = threshold$selected, audit = threshold$audit)
  log_event("Completed candidate", candidate, "macro accuracy",
    metric_rows[[candidate]]$macro_cluster_accuracy)
}

metrics <- do.call(rbind, metric_rows); rownames(metrics) <- NULL
predictions <- do.call(rbind, all_predictions); rownames(predictions) <- NULL
neighbors <- do.call(rbind, all_neighbors); rownames(neighbors) <- NULL
write_tsv(metrics, file.path(dirs[["tables"]], "05_candidate_metrics.tsv"))
write_tsv(predictions, file.path(dirs[["tables"]], "06_candidate_predictions.tsv"))
neighbor_path <- gzfile(file.path(dirs[["tables"]], "07_candidate_neighbors.tsv.gz"), "wt")
write.table(neighbors, neighbor_path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA"); close(neighbor_path)
write_tsv(hybrid_threshold$audit,
  file.path(dirs[["tables"]], "08_hybrid_l4_threshold_audit.tsv"))

hybrid <- metrics[metrics$candidate == cfg$qualification$primary_candidate, ]
full <- metrics[metrics$candidate == "full191433", ]
gate <- list(
  macro_accuracy = hybrid$macro_cluster_accuracy >=
    max(metrics$macro_cluster_accuracy) -
      cfg$qualification$hybrid_gate$maximum_macro_accuracy_loss,
  pseudobulk_accuracy = hybrid$pseudobulk_accuracy >=
    max(metrics$pseudobulk_accuracy) -
      cfg$qualification$hybrid_gate$maximum_pseudobulk_accuracy_loss,
  depth_correlation = hybrid$max_abs_depth_latent_spearman <=
    full$max_abs_depth_latent_spearman +
      cfg$qualification$hybrid_gate$maximum_depth_correlation_excess_over_full,
  l4_specificity = hybrid$l4_specificity >=
    cfg$qualification$l4_gate$minimum_specificity)
gate_passed <- all(unlist(gate, use.names = FALSE))
qualification <- list(
  schema_version = "bdev.full_peakvi_candidate_qualification.v1",
  status = if (gate_passed) "passed" else "failed",
  run_identity = cfg$run_identity,
  selected_candidate = if (gate_passed) cfg$qualification$primary_candidate else NULL,
  gate = gate,
  frozen_l4_rule = if (gate_passed) list(
    membership_source = "authenticated fitted Slingshot L4 weights",
    neighbor_k = cfg$qualification$neighbor_k,
    membership_threshold = hybrid$l4_threshold,
    reference_positive_weight = cfg$qualification$l4_gate$weight_positive_threshold,
    minimum_l4_neighbors_for_pseudotime =
      cfg$qualification$l4_gate$minimum_l4_neighbors_for_pseudotime,
    pseudotime_orientation = "native") else NULL,
  external_query_use = FALSE, real_immgen_use = FALSE,
  source_manifest_unchanged = identical(sha256(source_manifest), source_manifest_before))
write_json(qualification,
  file.path(dirs[["reports"]], "02_candidate_qualification_summary.json"))
stopifnot(qualification$source_manifest_unchanged)
if (!gate_passed)
  stop("The predeclared hybrid candidate gate failed; preserve results and do not freeze a full reference model")
log_event("Hybrid candidate gate passed; full-reference freeze is eligible")
