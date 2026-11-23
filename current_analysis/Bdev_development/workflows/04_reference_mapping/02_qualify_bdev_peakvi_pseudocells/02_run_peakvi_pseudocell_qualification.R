#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config"))
  stop("Usage: 02_run_peakvi_pseudocell_qualification.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbLatent)
  library(scwbRepresentation)
  library(scwbReferenceMap)
  library(Matrix)
  library(yaml)
  library(jsonlite)
  library(digest)
})

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
read_config <- function(path) {
  current <- yaml::read_yaml(path)
  base_name <- current$base_config
  current$base_config <- NULL
  if (is.null(base_name)) return(current)
  base <- read_config(normalizePath(file.path(dirname(path), base_name),
    mustWork = TRUE))
  utils::modifyList(base, current, keep.null = TRUE)
}
cfg <- read_config(cfg_path)
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_peakvi_pseudocell_qualification.log"),
  split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
write_tsv <- function(x, name) write.table(x, file.path(dirs[["tables"]], name),
  sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
package_identity <- function(package, expected) {
  observed <- as.character(packageVersion(package))
  stopifnot(identical(observed, as.character(expected)))
  observed
}
bind_rows <- function(rows) {
  fields <- unique(unlist(lapply(rows, names), use.names = FALSE))
  normalized <- lapply(rows, function(x) {
    for (field in setdiff(fields, names(x))) x[[field]] <- NA
    x[, fields, drop = FALSE]
  })
  out <- do.call(rbind, normalized); rownames(out) <- NULL; out
}
promoted_ref <- function(promoted, output, kind) {
  selected <- promoted$promoted$status == "ok" &
    promoted$promoted$output == output
  stopifnot(sum(selected) == 1L)
  wb_ref_artifact(promoted$promoted$artifact_id[selected][[1L]], kind)
}
metric_frame <- function(prediction) {
  types <- unique(as.character(prediction$query_type))
  do.call(rbind, lapply(types, function(type) {
    x <- prediction[prediction$query_type == type, , drop = FALSE]
    accepted <- x$support_call != "out_of_reference"
    complete <- accepted & is.finite(x$primary_continuous_value) &
      is.finite(x$true_native_l4)
    data.frame(query_type = type, n_query = nrow(x),
      n_supported = sum(x$support_call == "supported"),
      n_borderline = sum(x$support_call == "borderline"),
      n_out_of_reference = sum(x$support_call == "out_of_reference"),
      accepted_fraction = mean(accepted),
      l4_mae_accepted = if (any(complete)) mean(abs(
        x$primary_continuous_value[complete] - x$true_native_l4[complete])) else NA_real_,
      l4_rmse_accepted = if (any(complete)) sqrt(mean((
        x$primary_continuous_value[complete] - x$true_native_l4[complete])^2)) else NA_real_,
      l4_spearman_accepted = if (sum(complete) >= 3L) suppressWarnings(cor(
        x$primary_continuous_value[complete], x$true_native_l4[complete],
        method = "spearman")) else NA_real_, stringsAsFactors = FALSE)
  }))
}
one_column <- function(values, feature_ids, cell_id) {
  nz <- which(values != 0)
  Matrix::sparseMatrix(i = nz, j = rep.int(1L, length(nz)), x = values[nz],
    dims = c(length(feature_ids), 1L), dimnames = list(feature_ids, cell_id))
}

stopifnot(identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256),
  identical(cfg$inputs$authenticated_l4$orientation, "native"))
versions <- list(
  scworkbench = package_identity("scworkbench", cfg$runtime$scworkbench$version),
  scwbRepresentation = package_identity("scwbRepresentation",
    cfg$runtime$scwbRepresentation$version),
  scwbReferenceMap = package_identity("scwbReferenceMap",
    cfg$runtime$scwbReferenceMap$version),
  scwbLatent = package_identity("scwbLatent", cfg$runtime$scwbLatent$version))
for (package in c("scwbRepresentation", "scwbReferenceMap", "scwbLatent")) {
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  stopifnot(nzchar(manifest), identical(sha256(manifest),
    cfg$runtime[[package]]$manifest_sha256))
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
}
required_ops <- c("scwblatent::peakvi_fit", "scwblatent::peakvi_project",
  "scwblatent::model_bound_latent_union",
  "scwbrepresentation::knn_neighbors",
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
family <- readRDS(cfg$inputs$authenticated_l4$path)
selected <- as.character(family$key) == cfg$inputs$authenticated_l4$key &
  as.character(family$cell_id) %in% fixture_cells &
  is.finite(as.numeric(family$value))
l4 <- stats::setNames(as.numeric(family$value[selected]),
  as.character(family$cell_id[selected]))
stopifnot(length(l4) == length(fixture_cells), !anyDuplicated(names(l4)),
  all(fixture_cells %in% names(l4)))
fixture$native_l4 <- unname(l4[fixture_cells])
stopifnot(isTRUE(all.equal(fixture$pseudotime, fixture$native_l4,
  tolerance = 1e-12, check.attributes = FALSE)))

log_event("Reading exact Bdev ATAC count universe for", length(fixture_cells),
  "authenticated L4 fixture cells")
source <- wb_open_bundle(cfg$inputs$source_bundle$path, validate = FALSE,
  read_only = TRUE)
source_cells <- as.character(wb_cells(source, "all"))
stopifnot(length(source_cells) == cfg$inputs$source_bundle$expected_cells,
  all(fixture_cells %in% source_cells))
feature_ids <- wb_feature_ids(source, cfg$inputs$count_matrix)
stopifnot(length(feature_ids) == cfg$inputs$expected_features,
  !anyDuplicated(feature_ids))
metadata <- wb_meta(source, cells = fixture_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$sample_field,
  cfg$inputs$annotation_field))
counts <- wb_slice(source, cfg$inputs$count_matrix, features = feature_ids,
  cells = fixture_cells)
wb_close_bundle(source)
stopifnot(inherits(counts, "sparseMatrix"), identical(rownames(counts), feature_ids),
  identical(colnames(counts), fixture_cells), length(counts@x) > 0L,
  all(counts@x >= 0), all(counts@x == floor(counts@x)))

fixture$sample_id <- as.character(metadata[[cfg$inputs$sample_field]])
fixture$annotation <- as.character(metadata[[cfg$inputs$annotation_field]])
fixture$group_id <- sprintf("sample_%s_bin_%02d", fixture$sample_id,
  as.integer(fixture$bin))
group_sizes <- table(fixture$group_id)
eligible_groups <- names(group_sizes)[group_sizes >= cfg$split$minimum_group_cells]
stopifnot(length(eligible_groups) == cfg$split$expected_groups)
fixture$partition <- "reference"
for (group in eligible_groups) {
  idx <- which(fixture$group_id == group)
  idx <- idx[order(fixture$bundle_order[idx], method = "radix")]
  fixture$partition[idx[seq_along(idx) %% 4L == 0L]] <- "held_out"
}
reference_cells <- fixture$cell_id[fixture$partition == "reference"]
held_out_cells <- fixture$cell_id[fixture$partition == "held_out"]
stopifnot(length(reference_cells) + length(held_out_cells) == nrow(fixture),
  !length(intersect(reference_cells, held_out_cells)), length(held_out_cells) > 0L)
depths <- Matrix::colSums(counts)
training_depths <- depths[reference_cells]
depth_targets <- as.integer(round(stats::quantile(training_depths,
  probs = as.numeric(unlist(cfg$split$pseudocell_depth_quantiles)),
  names = FALSE, type = 7)))
stopifnot(length(depth_targets) == 3L, all(depth_targets > 0L))

group_rows <- split(which(fixture$partition == "held_out"),
  fixture$group_id[fixture$partition == "held_out"])
stopifnot(identical(sort(names(group_rows), method = "radix"),
  sort(eligible_groups, method = "radix")))
bulk_columns <- list(); pseudo_columns <- list(); design <- list(); cursor <- 0L
set.seed(as.integer(cfg$split$seed))
for (group in sort(names(group_rows), method = "radix")) {
  rows <- group_rows[[group]]
  cells <- fixture$cell_id[rows]
  aggregate <- Matrix::rowSums(counts[, cells, drop = FALSE])
  annotation_counts <- table(fixture$annotation[rows])
  truth_annotation <- names(annotation_counts)[order(-as.integer(annotation_counts),
    names(annotation_counts), method = "radix")][[1L]]
  truth_l4 <- stats::median(fixture$native_l4[rows])
  sample_id <- unique(fixture$sample_id[rows]); l4_bin <- unique(fixture$bin[rows])
  stopifnot(length(sample_id) == 1L, length(l4_bin) == 1L)
  bulk_id <- paste0("pseudobulk__", group)
  bulk_columns[[bulk_id]] <- one_column(aggregate, feature_ids, bulk_id)
  cursor <- cursor + 1L
  design[[cursor]] <- data.frame(query_cell_id = bulk_id,
    query_type = "raw_pseudobulk", group_id = group, sample_id = sample_id,
    l4_bin = as.integer(l4_bin), replicate = 1L, target_depth = NA_integer_,
    observed_depth = as.numeric(sum(aggregate)), constituent_count = length(cells),
    true_native_l4 = truth_l4, true_annotation = truth_annotation,
    stringsAsFactors = FALSE)
  nz <- which(aggregate > 0)
  probabilities <- aggregate[nz] / sum(aggregate[nz])
  for (replicate in seq_along(depth_targets)) {
    pseudo_id <- paste0("pseudocell_q", replicate, "__", group)
    sampled <- as.numeric(stats::rmultinom(1L, size = depth_targets[[replicate]],
      prob = probabilities))
    values <- numeric(length(feature_ids)); values[nz] <- sampled
    pseudo_columns[[pseudo_id]] <- one_column(values, feature_ids, pseudo_id)
    cursor <- cursor + 1L
    design[[cursor]] <- data.frame(query_cell_id = pseudo_id,
      query_type = "depth_matched_pseudocell", group_id = group,
      sample_id = sample_id, l4_bin = as.integer(l4_bin),
      replicate = replicate, target_depth = depth_targets[[replicate]],
      observed_depth = sum(sampled), constituent_count = length(cells),
      true_native_l4 = truth_l4, true_annotation = truth_annotation,
      stringsAsFactors = FALSE)
  }
}
synthetic_design <- bind_rows(design)
bulk_matrix <- do.call(cbind, bulk_columns)
pseudo_matrix <- do.call(cbind, pseudo_columns)
held_out_matrix <- counts[, held_out_cells, drop = FALSE]
actual_rows <- fixture[match(held_out_cells, fixture$cell_id), , drop = FALSE]
actual_design <- data.frame(query_cell_id = held_out_cells,
  query_type = "held_out_single_cell", group_id = actual_rows$group_id,
  sample_id = actual_rows$sample_id, l4_bin = as.integer(actual_rows$bin),
  replicate = 1L, target_depth = as.integer(depths[held_out_cells]),
  observed_depth = as.numeric(depths[held_out_cells]), constituent_count = 1L,
  true_native_l4 = actual_rows$native_l4,
  true_annotation = actual_rows$annotation, stringsAsFactors = FALSE)
query_design <- rbind(actual_design, synthetic_design)
query_matrix <- cbind(held_out_matrix, bulk_matrix, pseudo_matrix)
pseudobulk_group_count <- ncol(bulk_matrix)
query_matrix <- query_matrix[, query_design$query_cell_id, drop = FALSE]
stopifnot(identical(colnames(query_matrix), query_design$query_cell_id),
  all(query_design$observed_depth == Matrix::colSums(query_matrix)),
  all(synthetic_design$observed_depth[synthetic_design$query_type ==
    "depth_matched_pseudocell"] == synthetic_design$target_depth[
    synthetic_design$query_type == "depth_matched_pseudocell"]))
write_tsv(fixture, "01_reference_and_holdout_cells.tsv")
write_tsv(query_design, "02_synthetic_query_design.tsv")

log_event("Building fresh qualification bundle with", length(reference_cells),
  "Bdev-only reference cells and", ncol(query_matrix), "disjoint queries")
bundle_exists <- dir.exists(cfg$storage$qualification_bundle)
bundle <- if (bundle_exists) wb_open_bundle(cfg$storage$qualification_bundle,
  validate = FALSE, read_only = FALSE) else
  wb_init_bundle(cfg$storage$qualification_bundle, overwrite = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
all_bundle_cells <- c(reference_cells, query_design$query_cell_id)
reference_matrix <- counts[, reference_cells, drop = FALSE]
if (!bundle_exists) {
  wb_define_cells(bundle, all_bundle_cells)
  wb_attach_matrix(bundle, "Bdev_ATAC_reference_counts", reference_matrix,
    layer = "counts", feature_type = "peak", cell_ids = reference_cells,
    description = "Exact Bdev ATAC peak universe; Bdev-only PeakVI reference cells",
    source = "active_03c_bundle_ATAC_counts_fixture_slice",
    created_by = "scworkbench_Bdev Step 4 development qualification")
  wb_attach_matrix(bundle, "Bdev_ATAC_disjoint_queries", query_matrix,
    layer = "counts", feature_type = "peak", cell_ids = query_design$query_cell_id,
    description = "Held-out Bdev cells, raw pseudobulks, and depth-matched pseudo-cells",
    source = "active_03c_bundle_ATAC_counts_derived_validation_queries",
    created_by = "scworkbench_Bdev Step 4 development qualification")
  source_meta <- metadata[match(c(reference_cells, held_out_cells), metadata$cell_id),
    , drop = FALSE]
  synthetic_meta <- data.frame(cell_id = synthetic_design$query_cell_id,
    stringsAsFactors = FALSE, check.names = FALSE)
  synthetic_meta[[cfg$inputs$sample_field]] <- synthetic_design$sample_id
  synthetic_meta[[cfg$inputs$annotation_field]] <- synthetic_design$true_annotation
  all_meta <- rbind(source_meta, synthetic_meta)
  all_meta <- all_meta[match(all_bundle_cells, all_meta$cell_id), , drop = FALSE]
  stopifnot(identical(as.character(all_meta$cell_id), all_bundle_cells))
  wb_attach_metadata(bundle, all_meta, cell_ids = all_bundle_cells,
    source = "active_03c_metadata_plus_synthetic_query_truth",
    created_by = "scworkbench_Bdev Step 4 development qualification")
} else {
  stopifnot(identical(as.character(wb_cells(bundle, "all")), all_bundle_cells),
    all(c("Bdev_ATAC_reference_counts", "Bdev_ATAC_disjoint_queries") %in%
      wb_matrices(bundle)$matrix_id))
}
rm(counts, reference_matrix, query_matrix, held_out_matrix, bulk_matrix,
  pseudo_matrix, bulk_columns, pseudo_columns); invisible(gc())

log_event("Fitting Bdev-only PeakVI reference model")
artifacts <- wb_art_list(bundle)
retained_projection_artifact <- artifacts$source_op ==
  "scwblatent::peakvi_project" & artifacts$type == "obs_table" &
  artifacts$lifecycle_state == "current"
if (sum(retained_projection_artifact) == 1L) {
  projection_lineage <- wb_prov_artifact_lineage(bundle,
    artifacts$artifact_id[retained_projection_artifact][[1L]],
    direction = "upstream", max_depth = 5L)
  lineage_model_ids <- unique(na.omit(projection_lineage$nodes$artifact_id[
    projection_lineage$nodes$artifact_type == "model_ref" &
      projection_lineage$nodes$op_name == "scwblatent::peakvi_fit"]))
  lineage_summary_ids <- unique(na.omit(projection_lineage$nodes$artifact_id[
    projection_lineage$nodes$artifact_type == "obs_table" &
      projection_lineage$nodes$op_name == "scwblatent::peakvi_fit"]))
  stopifnot(length(lineage_model_ids) == 1L,
    length(lineage_summary_ids) == 1L)
  retained_fit_model <- artifacts$artifact_id == lineage_model_ids[[1L]]
  retained_fit_summary <- artifacts$artifact_id == lineage_summary_ids[[1L]]
} else {
  retained_fit_model <- artifacts$source_op == "scwblatent::peakvi_fit" &
    artifacts$type == "model_ref" & artifacts$lifecycle_state == "current"
  retained_fit_summary <- artifacts$source_op == "scwblatent::peakvi_fit" &
    artifacts$type == "obs_table" & artifacts$lifecycle_state == "current"
}
fit_retained <- sum(retained_fit_model) == 1L &&
  sum(retained_fit_summary) == 1L
if (fit_retained) {
  model_id <- artifacts$artifact_id[retained_fit_model][[1L]]
  fit_summary_id <- artifacts$artifact_id[retained_fit_summary][[1L]]
  model_value <- wb_art_load(bundle, model_id, expected_kind = "model_ref")
  reference_space_id <- as.character(model_value$model_space_id)[[1L]]
  model_ref <- wb_ref_artifact(model_id, "model")
  fit_summary_ref <- wb_ref_artifact(fit_summary_id, "obs_table")
  fit_summary <- wb_art_load(bundle, fit_summary_id,
    expected_kind = "obs_table")$table
  stopifnot(reference_space_id %in% wb_spaces(bundle)$space_id)
  log_event("Recovered retained PeakVI fit from bundle registry")
} else {
  stopifnot(sum(retained_fit_model) == 0L, sum(retained_fit_summary) == 0L)
  fit_workspace <- file.path(cfg$storage$sessions,
    paste0(cfg$run_identity, "_peakvi_fit"))
  fit_session <- if (dir.exists(fit_workspace))
    wb_session_resume(bundle, fit_workspace) else
    wb_session_create(bundle, cfg$storage$sessions,
      paste0(cfg$run_identity, "_peakvi_fit"))
  if (length(fit_session$runs)) {
    stopifnot(length(fit_session$runs) == 1L)
    fit <- list(session = fit_session, run = fit_session$runs[[1L]])
    log_event("Recovered completed staged PeakVI fit without refitting")
  } else
  fit_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_peakvi_fit_retry"))
  if (!exists("fit", inherits = FALSE))
    fit <- wb_session_run_op(fit_session, bundle, "scwblatent::peakvi_fit",
    inputs = list(matrix = "Bdev_ATAC_reference_counts",
      n_latent = as.integer(cfg$peakvi$n_latent),
      n_layers = as.integer(cfg$peakvi$n_layers),
      n_hidden = as.integer(cfg$peakvi$n_hidden),
      dropout = as.numeric(cfg$peakvi$dropout),
      max_epochs = as.integer(cfg$peakvi$max_epochs),
      lr = as.numeric(cfg$peakvi$learning_rate),
      early_stopping = isTRUE(cfg$peakvi$early_stopping),
      seed = as.integer(cfg$peakvi$seed),
      chunk_size = as.integer(cfg$peakvi$chunk_size),
      batch_key = NULL, use_gpu = FALSE, cell_set = "all", min_cells = 1L),
    notes = "Step 4 Bdev-only PeakVI fit; no ImmGen observations")
  reference_space_id <- as.character(fit$run$outputs$space$space_id %||%
    fit$run$outputs$space$allocated_id)[[1L]]
  stopifnot(nzchar(reference_space_id))
  fit_summary <- fit$run$outputs$fit_summary$table
  fit$session <- wb_session_resume(bundle, fit$session$durable$workspace)
  fit_promoted <- wb_session_promote_current(bundle, fit$session,
    outputs = c("model", "fit_summary"), name_prefix = "step4_peakvi_reference",
    description = "Bdev-only PeakVI reference model",
    created_by = "scworkbench_Bdev Step 4 development qualification",
    fail_fast = TRUE)
  model_ref <- promoted_ref(fit_promoted, "model", "model")
  fit_summary_ref <- promoted_ref(fit_promoted, "fit_summary", "obs_table")
  fit_promoted$session <- wb_session_discard(fit_promoted$session)
}
stopifnot(nrow(fit_summary) == 1L,
  as.integer(fit_summary$n_latent) == as.integer(cfg$peakvi$n_latent))

log_event("Projecting held-out and synthetic queries without refitting")
artifacts <- wb_art_list(bundle)
retained_projection <- artifacts$source_op == "scwblatent::peakvi_project" &
  artifacts$type == "obs_table" & artifacts$lifecycle_state == "current"
projection_space_retained <- "step4_peakvi_disjoint_queries" %in%
  wb_spaces(bundle)$space_id
if (sum(retained_projection) == 1L && projection_space_retained) {
  projection_summary_ref <- wb_ref_artifact(
    artifacts$artifact_id[retained_projection][[1L]], "obs_table")
  projection_summary <- wb_art_load(bundle,
    artifacts$artifact_id[retained_projection][[1L]],
    expected_kind = "obs_table")$table
  log_event("Recovered retained frozen-model projection from bundle registry")
} else {
  stopifnot(sum(retained_projection) == 0L, !projection_space_retained)
  project_workspace <- file.path(cfg$storage$sessions,
    paste0(cfg$run_identity, "_peakvi_project"))
  project_session <- if (dir.exists(project_workspace))
    wb_session_resume(bundle, project_workspace) else
    wb_session_create(bundle, cfg$storage$sessions,
      paste0(cfg$run_identity, "_peakvi_project"))
  stopifnot(length(project_session$runs) == 0L)
  project <- wb_session_run_op(project_session, bundle,
    "scwblatent::peakvi_project", inputs = list(model = model_ref,
      query_matrix = "Bdev_ATAC_disjoint_queries",
      space_id_out = "step4_peakvi_disjoint_queries",
      chunk_size = as.integer(cfg$peakvi$chunk_size), batch_key = NULL,
      use_gpu = FALSE,
      pickle_execution_trust = cfg$peakvi$pickle_execution_trust,
      cell_set = "all", min_cells = 1L),
    notes = "Step 4 held-out Bdev projection through frozen PeakVI")
  projection_summary <- project$run$outputs$projection_summary$table
  project$session <- wb_session_resume(bundle,
    project$session$durable$workspace)
  project_promoted <- wb_session_promote_current(bundle, project$session,
    outputs = "projection_summary", name_prefix = "step4_peakvi_projection",
    description = "Frozen-model projection of disjoint Bdev validation queries",
    created_by = "scworkbench_Bdev Step 4 development qualification",
    fail_fast = TRUE)
  projection_summary_ref <- promoted_ref(project_promoted,
    "projection_summary", "obs_table")
  project_promoted$session <- wb_session_discard(project_promoted$session)
}
stopifnot(identical(as.character(projection_summary$transformation),
  "transfer_without_refit"), nrow(projection_summary) == 1L,
  as.integer(projection_summary$n_latent) ==
    as.integer(cfg$peakvi$n_latent))

log_event("Creating authenticated joint latent space")
wb_close_bundle(bundle)
bundle <- wb_open_bundle(cfg$storage$qualification_bundle,
  validate = FALSE, read_only = FALSE)
joint_space_id <- "step4_peakvi_reference_and_queries"
if (!joint_space_id %in% wb_spaces(bundle)$space_id) {
  union_attempt_id <- paste0(cfg$run_identity,
    "_latent_union_registry_recovery_2")
  union_session <- wb_session_create(bundle, cfg$storage$sessions,
    union_attempt_id)
  union <- tryCatch(wb_session_run_op(union_session, bundle,
    "scwblatent::model_bound_latent_union", inputs = list(
      reference_space = reference_space_id,
      query_space = "step4_peakvi_disjoint_queries", model = model_ref,
      projection_summary = projection_summary_ref,
      space_id_out = joint_space_id,
      chunk_size = as.integer(cfg$mapping$query_chunk_size)),
    notes = "Exact union under one retained Bdev-only PeakVI model"),
    error = identity)
  if (inherits(union, "error")) {
    if (isTRUE(cfg$acceptance$require_authenticated_model_bound_union))
      stop(conditionMessage(union), call. = FALSE)
    stopifnot(grepl("Complete typed-space union source changed after authentication",
      conditionMessage(union), fixed = TRUE))
    reference_coordinates <- wb_space(bundle, reference_space_id,
      cells = reference_cells)
    query_coordinates <- wb_space(bundle, "step4_peakvi_disjoint_queries",
      cells = query_design$query_cell_id)
    stopifnot(ncol(reference_coordinates) == ncol(query_coordinates),
      !length(intersect(rownames(reference_coordinates),
        rownames(query_coordinates))))
    joint_coordinates <- rbind(reference_coordinates, query_coordinates)
    joint_frame <- data.frame(cell_id = rownames(joint_coordinates),
      as.data.frame(joint_coordinates, check.names = FALSE),
      check.names = FALSE, stringsAsFactors = FALSE)
    wb_attach_space(bundle, joint_space_id, joint_frame,
      cell_ids = rownames(joint_coordinates),
      source = paste0("exact row union of retained model-bound spaces; model=",
        model_ref$artifact_id),
      created_by = "scworkbench_Bdev Step 4 development qualification")
    observed_joint <- wb_space(bundle, joint_space_id,
      cells = rownames(joint_coordinates))
    stopifnot(identical(dim(observed_joint), dim(joint_coordinates)),
      identical(rownames(observed_joint), rownames(joint_coordinates)),
      isTRUE(all.equal(as.numeric(as.matrix(observed_joint)),
        as.numeric(as.matrix(joint_coordinates)), tolerance = 0,
        check.attributes = FALSE)))
    union_mode <- "exact_direct_attachment_after_model_bound_union_contract_rejection"
  } else {
    union$session <- wb_session_resume(bundle, union$session$durable$workspace)
    union_promoted <- wb_session_promote_current(bundle, union$session,
      outputs = "union_summary", name_prefix = "step4_peakvi_union",
      description = "Model-bound union of Bdev reference and validation queries",
      created_by = "scworkbench_Bdev Step 4 development qualification",
      fail_fast = TRUE)
    union_promoted$session <- wb_session_discard(union_promoted$session)
    union_mode <- "qualified_model_bound_latent_union"
  }
} else {
  retained_union <- wb_art_list(bundle)
  retained_union <- retained_union$source_op ==
    "scwblatent::model_bound_latent_union" &
    retained_union$type == "obs_table" &
    retained_union$lifecycle_state == "current"
  if (isTRUE(cfg$acceptance$require_authenticated_model_bound_union))
    stopifnot(sum(retained_union) == 1L)
  union_mode <- if (sum(retained_union) == 1L)
    "qualified_model_bound_latent_union_recovered" else
    "reused_retained_joint_space"
}
reference_coordinates_check <- wb_space(bundle, reference_space_id,
  cells = reference_cells)
query_coordinates_check <- wb_space(bundle, "step4_peakvi_disjoint_queries",
  cells = query_design$query_cell_id)
expected_joint_check <- rbind(reference_coordinates_check,
  query_coordinates_check)
observed_joint_check <- wb_space(bundle, joint_space_id,
  cells = rownames(expected_joint_check))
coordinate_columns <- setdiff(colnames(expected_joint_check), "cell_id")
expected_coordinate_count <- ncol(expected_joint_check) -
  as.integer(is.data.frame(expected_joint_check) &&
    "cell_id" %in% names(expected_joint_check))
observed_cell_ids <- if (is.data.frame(observed_joint_check) &&
    "cell_id" %in% names(observed_joint_check))
  as.character(observed_joint_check$cell_id) else
  rownames(observed_joint_check)
expected_cell_ids <- if (is.data.frame(expected_joint_check) &&
    "cell_id" %in% names(expected_joint_check))
  as.character(expected_joint_check$cell_id) else
  rownames(expected_joint_check)
stopifnot(identical(dim(observed_joint_check), dim(expected_joint_check)),
  identical(rownames(observed_joint_check), rownames(expected_joint_check)),
  identical(observed_cell_ids, expected_cell_ids),
  length(coordinate_columns) == expected_coordinate_count,
  length(coordinate_columns) == as.integer(cfg$peakvi$n_latent),
  isTRUE(all.equal(as.numeric(as.matrix(
    observed_joint_check[, coordinate_columns, drop = FALSE])),
    as.numeric(as.matrix(expected_joint_check[, coordinate_columns,
      drop = FALSE])), tolerance = 0,
    check.attributes = FALSE)))
wb_close_bundle(bundle)
bundle <- wb_open_bundle(cfg$storage$qualification_bundle,
  validate = FALSE, read_only = FALSE)

reference_annotation <- wb_prim_cell_vector(stats::setNames(
  fixture$annotation[match(reference_cells, fixture$cell_id)], reference_cells),
  value_type = "categorical")
reference_l4 <- wb_prim_cell_vector_family(
  rep(cfg$inputs$authenticated_l4$key, length(reference_cells)), reference_cells,
  unname(l4[reference_cells]), value_type = "continuous")
log_event("Calibrating reference support and native-L4 transfer")
calibration_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_calibrate"))
calibration <- wb_session_run_op(calibration_session, bundle,
  "scwbreferencemap::reference_map_calibrate", inputs = list(
    reference_matrix = "Bdev_ATAC_reference_counts",
    reference_space = reference_space_id,
    reference_annotation = reference_annotation,
    reference_continuous = reference_l4,
    holdout_group_field = cfg$inputs$sample_field,
    k_values = as.integer(unlist(cfg$mapping$k_values)),
    primary_k = as.integer(cfg$mapping$primary_k),
    support_quantiles = as.list(as.numeric(unlist(cfg$mapping$support_quantiles))),
    primary_continuous_key = cfg$inputs$authenticated_l4$key,
    space_chunk_size = as.integer(cfg$mapping$space_chunk_size),
    query_chunk_size = as.integer(cfg$mapping$query_chunk_size)),
  notes = "Step 4 PeakVI latent support calibration")
calibration$session <- wb_session_resume(bundle,
  calibration$session$durable$workspace)
calibration_promoted <- wb_session_promote_current(bundle, calibration$session,
  outputs = names(calibration$run$outputs), name_prefix = "step4_peakvi_calibration",
  description = "Bdev PeakVI support and transfer calibration",
  created_by = "scworkbench_Bdev Step 4 development qualification",
  fail_fast = TRUE)
calibration_ref <- promoted_ref(calibration_promoted, "calibration", "model")
thresholds <- calibration$run$outputs$thresholds$table
supported_threshold <- thresholds$threshold[thresholds$role == "supported_upper"][[1L]]
out_of_reference_threshold <- thresholds$threshold[
  thresholds$role == "borderline_upper"][[1L]]
calibration_promoted$session <- wb_session_discard(calibration_promoted$session)
wb_close_bundle(bundle)
bundle <- wb_open_bundle(cfg$storage$qualification_bundle,
  validate = FALSE, read_only = FALSE)

knn_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_knn"))
knn <- wb_session_run_op(knn_session, bundle,
  "scwbrepresentation::knn_neighbors", inputs = list(
    space = "step4_peakvi_reference_and_queries",
    query_cells = query_design$query_cell_id,
    candidate_cells = reference_cells, k = as.integer(cfg$mapping$neighbor_k),
    include_query = FALSE, backend = cfg$mapping$backend,
    space_chunk_size = as.integer(cfg$mapping$space_chunk_size),
    query_chunk_size = as.integer(cfg$mapping$query_chunk_size)),
  notes = "Exact query-to-reference neighbors in frozen PeakVI latent space")

transfer_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_transfer"))
transfer <- wb_session_run_op(transfer_session, bundle,
  "scwbreferencemap::reference_map_transfer", inputs = list(
    in_space = "step4_peakvi_disjoint_queries",
    reference_space = reference_space_id, neighbors = knn$run$outputs$neighbors,
    calibration = calibration_ref, support_output_key = "step4_peakvi_support",
    annotation_output_key = "step4_peakvi_mapped_cluster",
    continuous_output_keys = "mapped_native_L4",
    supported_threshold = as.numeric(supported_threshold),
    out_of_reference_threshold = as.numeric(out_of_reference_threshold),
    cell_ids = query_design$query_cell_id, min_cells = 1L,
    space_chunk_size = as.integer(cfg$mapping$space_chunk_size)),
  notes = "Step 4 frozen-PeakVI support and native-L4 transfer")
prediction <- transfer$run$outputs$per_query_summary$table
truth_index <- match(prediction$query_cell_id, query_design$query_cell_id)
stopifnot(!anyNA(truth_index), identical(as.character(prediction$query_cell_id),
  query_design$query_cell_id))
prediction <- cbind(prediction, query_design[truth_index, setdiff(
  names(query_design), "query_cell_id"), drop = FALSE])
prediction$absolute_l4_error <- abs(prediction$primary_continuous_value -
  prediction$true_native_l4)
prediction$annotation_correct <- !is.na(prediction$accepted_annotation) &
  prediction$accepted_annotation == prediction$true_annotation
metrics <- metric_frame(prediction)
write_tsv(prediction, "03_query_transfer_predictions.tsv")
write_tsv(metrics, "04_metrics_by_query_type.tsv")
write_tsv(calibration$run$outputs$calibration_metrics$table,
  "05_calibration_metrics.tsv")
write_tsv(thresholds, "06_support_thresholds.tsv")
write_tsv(transfer$run$outputs$grouped_summary$table,
  "07_support_grouped_summary.tsv")
knn$session <- wb_session_resume(bundle, knn$session$durable$workspace)
transfer$session <- wb_session_resume(bundle, transfer$session$durable$workspace)
knn$session <- wb_session_discard(knn$session)
transfer$session <- wb_session_discard(transfer$session)

metric <- function(type, field) metrics[metrics$query_type == type, field][[1L]]
performance_acceptance <- list(
  held_out_single_cell = isTRUE(metric("held_out_single_cell", "accepted_fraction") >=
      cfg$acceptance$held_out_single_cell$minimum_accepted_fraction) &&
    isTRUE(metric("held_out_single_cell", "l4_spearman_accepted") >=
      cfg$acceptance$held_out_single_cell$minimum_l4_spearman),
  depth_matched_pseudocell = isTRUE(metric("depth_matched_pseudocell",
      "accepted_fraction") >=
      cfg$acceptance$depth_matched_pseudocell$minimum_accepted_fraction) &&
    isTRUE(metric("depth_matched_pseudocell", "l4_spearman_accepted") >=
      cfg$acceptance$depth_matched_pseudocell$minimum_l4_spearman) &&
    isTRUE(metric("depth_matched_pseudocell", "l4_mae_accepted") <=
      cfg$acceptance$depth_matched_pseudocell$maximum_l4_mae),
  raw_pseudobulk = isTRUE(metric("raw_pseudobulk", "accepted_fraction") >=
      cfg$acceptance$raw_pseudobulk$minimum_accepted_fraction) &&
    isTRUE(metric("raw_pseudobulk", "l4_spearman_accepted") >=
      cfg$acceptance$raw_pseudobulk$minimum_l4_spearman))
source_manifest_after <- sha256(source_manifest)
structural_acceptance <- list(
  source_manifest_unchanged = identical(source_manifest_before,
    source_manifest_after),
  exact_feature_axis = nrow(fit_summary) == 1L &&
    as.integer(fit_summary$n_features) == cfg$inputs$expected_features,
  disjoint_reference_and_query = !length(intersect(reference_cells,
    query_design$query_cell_id)),
  projection_without_refit = identical(as.character(
    projection_summary$transformation), "transfer_without_refit"),
  requested_latent_dimensions = nrow(fit_summary) == 1L &&
    nrow(projection_summary) == 1L &&
    as.integer(fit_summary$n_latent) == as.integer(cfg$peakvi$n_latent) &&
    as.integer(projection_summary$n_latent) ==
      as.integer(cfg$peakvi$n_latent) &&
    length(coordinate_columns) == as.integer(cfg$peakvi$n_latent),
  authenticated_model_bound_union =
    !isTRUE(cfg$acceptance$require_authenticated_model_bound_union) ||
    union_mode %in% c("qualified_model_bound_latent_union",
      "qualified_model_bound_latent_union_recovered"),
  native_l4_used = identical(cfg$inputs$authenticated_l4$orientation, "native"),
  pseudocell_depths_exact = all(synthetic_design$observed_depth[
    synthetic_design$query_type == "depth_matched_pseudocell"] ==
    synthetic_design$target_depth[
    synthetic_design$query_type == "depth_matched_pseudocell"]),
  real_immgen_projected = FALSE, new_sif_built = FALSE,
  implementation_unlocked = FALSE)
stopifnot(all(unlist(structural_acceptance[c("source_manifest_unchanged",
  "exact_feature_axis", "disjoint_reference_and_query",
  "projection_without_refit", "requested_latent_dimensions",
  "authenticated_model_bound_union", "native_l4_used",
  "pseudocell_depths_exact")])))

finalized <- wb_finalize_recipient_bundle(bundle)
stopifnot(isTRUE(finalized$ready))
wb_close_bundle(bundle)
bundle <- NULL
bundle_manifest_sha <- sha256(file.path(cfg$storage$qualification_bundle,
  "manifest.json"))
bundle_files <- sort(list.files(cfg$storage$qualification_bundle,
  recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE),
  method = "radix")
bundle_files <- bundle_files[file.info(bundle_files)$isdir %in% FALSE]
bundle_inventory <- list(
  schema_version = "bdev.peakvi_pseudocell_bundle_inventory.v1",
  bundle = cfg$storage$qualification_bundle,
  manifest_sha256 = bundle_manifest_sha,
  file_count = length(bundle_files),
  files = lapply(bundle_files, function(path) list(
    path = substring(path, nchar(cfg$storage$qualification_bundle) + 2L),
    bytes = unname(file.info(path)$size), sha256 = sha256(path))))
write_json(bundle_inventory, file.path(dirs[["manifests"]],
  "03_bundle_inventory.json"))

gate_passed <- all(unlist(performance_acceptance))
summary <- list(
  schema_version = "bdev.peakvi_pseudocell_qualification_summary.v1",
  status = if (gate_passed) "passed" else "failed_performance_gate",
  run_identity = cfg$run_identity, scope = cfg$scope,
  runtime_versions = versions, reference_cells = length(reference_cells),
  held_out_single_cells = length(held_out_cells),
  pseudobulk_groups = pseudobulk_group_count,
  depth_matched_pseudocells = sum(query_design$query_type ==
    "depth_matched_pseudocell"),
  feature_count = length(feature_ids), depth_targets = as.list(depth_targets),
  peakvi_fit_summary = fit_summary, projection_summary = projection_summary,
  metrics_by_query_type = split(metrics, seq_len(nrow(metrics))),
  structural_acceptance = structural_acceptance,
  performance_acceptance = performance_acceptance,
  gate_passed = gate_passed,
  source_bundle_manifest_sha256 = source_manifest_after,
  qualification_bundle_manifest_sha256 = bundle_manifest_sha,
  joint_space_mode = union_mode,
  decision = if (gate_passed)
    "PeakVI pseudobulk/pseudo-cell gate passed; real ImmGen remains a later explicitly controlled step."
  else
    "PeakVI pseudobulk/pseudo-cell gate failed; do not project real ImmGen libraries.",
  boundaries = list(real_immgen_projected = FALSE, sif_built = FALSE,
    implementation_unlocked = FALSE))
write_json(summary, file.path(dirs[["reports"]],
  "01_peakvi_pseudocell_qualification_summary.json"))
write_json(list(
  schema_version = "bdev.peakvi_pseudocell_execution_receipt.v1",
  status = "complete", run_identity = cfg$run_identity,
  config_path = cfg_path, config_sha256 = sha256(cfg_path),
  source_manifest_sha256_before = source_manifest_before,
  source_manifest_sha256_after = source_manifest_after,
  qualification_bundle_manifest_sha256 = bundle_manifest_sha,
  gate_passed = gate_passed, archive_accessed = FALSE,
  real_immgen_projected = FALSE, sif_built = FALSE),
  file.path(dirs[["manifests"]], "02_execution_receipt.json"))
log_event("Step 4 PeakVI pseudocell qualification complete; gate_passed =",
  gate_passed)
