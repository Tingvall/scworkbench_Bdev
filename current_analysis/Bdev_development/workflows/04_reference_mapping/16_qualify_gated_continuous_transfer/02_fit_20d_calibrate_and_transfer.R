#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_index <- match("--config", args)
if (is.na(config_index) || config_index == length(args)) {
  stop("Usage: 02_fit_20d_calibrate_and_transfer.R --config CONFIG", call. = FALSE)
}
config_path <- normalizePath(args[[config_index + 1L]], mustWork = TRUE)

options(stringsAsFactors = FALSE, warn = 1, width = 180)
suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbLatent)
  library(scwbRepresentation)
  library(scwbReferenceMap)
  library(data.table)
  library(digest)
  library(jsonlite)
  library(yaml)
})

cfg <- yaml::read_yaml(config_path)
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
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
  } else {
    x <- as.matrix(x)
  }
  x
}

root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")), c("tables", "figures",
  "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_fit_calibrate_transfer.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(sha256(cfg$runtime$sif) == cfg$runtime$sif_sha256)
for (package in c("scwbLatent", "scwbRepresentation", "scwbReferenceMap")) {
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  stopifnot(nzchar(manifest), sha256(manifest) == cfg$runtime[[package]]$manifest_sha256)
  namespace <- switch(package, scwbLatent = "scwblatent",
    scwbRepresentation = "scwbrepresentation",
    scwbReferenceMap = "scwbreferencemap")
  try(wb_operation_pack_unload(namespace), silent = TRUE)
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
}
required_ops <- c("scwblatent::peakvi_fit", "scwblatent::peakvi_project",
  "scwblatent::model_bound_latent_union",
  "scwbrepresentation::knn_neighbors",
  "scwbreferencemap::reference_map_calibrate",
  "scwbreferencemap::reference_map_transfer")
stopifnot(all(required_ops %in% wb_list_ops()))

source_manifest <- file.path(cfg$inputs$qualification_parent_bundle$path,
  "manifest.json")
stopifnot(sha256(source_manifest) ==
  cfg$inputs$qualification_parent_bundle$manifest_sha256)
bundle <- wb_open_bundle(cfg$storage$qualification_bundle,
  validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)

all_cells <- as.character(wb_cells(bundle, "all"))
query_cells <- all_cells[grepl("^external_query_[0-9]+$", all_cells)]
reference_cells <- setdiff(all_cells, query_cells)
stopifnot(length(reference_cells) == cfg$acceptance$require_exact_reference_cells,
  length(query_cells) == cfg$acceptance$require_exact_query_cells)
matrix_table <- wb_matrices(bundle)
required_matrices <- unlist(cfg$inputs$qualification_parent_bundle[
  c("reference_matrix", "query_matrix", "authorizer_matrix")], use.names = FALSE)
stopifnot(all(required_matrices %in% matrix_table$matrix_id))
reference_matrix_row <- matrix_table[
  matrix_table$matrix_id == cfg$inputs$qualification_parent_bundle$reference_matrix, ]
query_matrix_row <- matrix_table[
  matrix_table$matrix_id == cfg$inputs$qualification_parent_bundle$query_matrix, ]
stopifnot(nrow(reference_matrix_row) == 1L, nrow(query_matrix_row) == 1L,
  reference_matrix_row$n_features[[1L]] == 20000L,
  query_matrix_row$n_features[[1L]] == 20000L,
  reference_matrix_row$n_cells[[1L]] == length(reference_cells),
  query_matrix_row$n_cells[[1L]] == length(query_cells))

metadata <- wb_meta(bundle, cells = reference_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$sample_field,
  cfg$inputs$l4_pseudotime_field, cfg$inputs$l4_weight_field))
stopifnot(identical(as.character(metadata$cell_id), reference_cells))
sample_id <- as.character(metadata[[cfg$inputs$sample_field]])
native_l4 <- as.numeric(metadata[[cfg$inputs$l4_pseudotime_field]])
l4_weight <- as.numeric(metadata[[cfg$inputs$l4_weight_field]])
stopifnot(!anyNA(sample_id), all(nzchar(sample_id)),
  length(unique(sample_id)) == 4L, all(is.finite(l4_weight)))

in_domain <- is.finite(native_l4) &
  l4_weight >= as.numeric(cfg$calibration$minimum_l4_weight)
domain <- ifelse(in_domain, cfg$calibration$in_domain_value,
  cfg$calibration$out_domain_value)
stopifnot(!anyNA(domain), all(nzchar(domain)), all(table(domain) > 0L))

trajectory_bin <- rep("out_domain", length(reference_cells))
in_index <- which(in_domain)
in_rank <- rank(native_l4[in_index], ties.method = "first")
n_bins <- as.integer(cfg$calibration$trajectory_bins)
in_bin <- pmin(n_bins, ceiling(in_rank / length(in_rank) * n_bins))
trajectory_bin[in_index] <- sprintf("native_l4_q%02d", in_bin)
stratum <- paste(sample_id, domain, trajectory_bin, sep = "|")
fold <- integer(length(reference_cells))
n_folds <- as.integer(cfg$calibration$folds)
for (level in sort(unique(stratum))) {
  index <- which(stratum == level)
  hash <- vapply(reference_cells[index], function(cell) digest::digest(
    paste(cell, cfg$calibration$fold_seed, sep = "|"), algo = "xxhash64",
    serialize = FALSE), character(1))
  ordered <- index[order(hash, reference_cells[index], method = "radix")]
  fold[ordered] <- (seq_along(ordered) - 1L) %% n_folds + 1L
}
fold_label <- sprintf("fold_%02d", fold)
fold_balance <- as.data.table(data.frame(fold = fold_label,
  sample_id = sample_id, l4_domain = domain,
  trajectory_stratum = trajectory_bin, stringsAsFactors = FALSE))[, .N,
  by = .(fold, sample_id, l4_domain, trajectory_stratum)]
fold_domain <- as.data.table(data.frame(fold = fold_label,
  sample_id = sample_id, l4_domain = domain, stringsAsFactors = FALSE))[, .N,
  by = .(fold, sample_id, l4_domain)]
expected_combinations <- n_folds * length(unique(sample_id)) *
  length(unique(domain))
stopifnot(length(unique(fold_label)) == n_folds,
  nrow(fold_domain) == expected_combinations, all(fold_domain$N > 0L))

qualification_metadata <- data.frame(cell_id = reference_cells,
  reference_map_balanced_fold_5 = fold_label,
  reference_map_l4_domain = domain,
  reference_map_l4_trajectory_stratum = trajectory_bin,
  stringsAsFactors = FALSE)
wb_attach_metadata(bundle, qualification_metadata, cell_ids = reference_cells,
  source = "deterministic sample-domain-native-L4 balanced qualification folds",
  created_by = "Bdev development reference mapping step 16")
fwrite(fold_balance, file.path(dirs[["tables"]],
  "01_reference_fold_stratum_counts.tsv.gz"), sep = "\t", quote = TRUE,
  na = "NA")
fwrite(fold_domain, file.path(dirs[["tables"]],
  "02_reference_fold_sample_domain_counts.tsv"), sep = "\t", quote = TRUE,
  na = "NA")

log_event("Fitting 20-dimensional HVP20k PeakVI reference")
fit_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_fit_20d"))
fit <- wb_session_run_op(fit_session, bundle, "scwblatent::peakvi_fit",
  inputs = list(
    matrix = cfg$inputs$qualification_parent_bundle$reference_matrix,
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
  notes = "Bdev-only 20-dimensional HVP20k gated-transfer qualification")
reference_space_id <- fit$run$outputs$space$space_id
if (is.null(reference_space_id)) reference_space_id <- fit$run$outputs$space$allocated_id
reference_space_id <- as.character(reference_space_id)[[1L]]
fit_summary <- fit$run$outputs$fit_summary$table
stopifnot(nrow(fit_summary) == 1L,
  as.integer(fit_summary$n_latent[[1L]]) == cfg$peakvi$n_latent)
fit$session <- wb_session_resume(bundle, fit$session$durable$workspace)
fit_promoted <- wb_session_promote_current(bundle, fit$session,
  outputs = c("model", "fit_summary"),
  name_prefix = "gated_transfer_hvp20k_20d_reference",
  description = "Bdev-only HVP20k 20-dimensional qualification reference",
  created_by = "Bdev development reference mapping step 16", fail_fast = TRUE)
model_artifact_id <- promoted_id(fit_promoted, "model")
model_ref <- wb_ref_artifact(model_artifact_id, "model")
fit_promoted$session <- wb_session_discard(fit_promoted$session)

query_space_id <- "gated_transfer_hvp20k_20d_query_fixture"
log_event("Projecting balanced 28-cell KO fixture without model refit")
project_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_project_fixture"))
project <- wb_session_run_op(project_session, bundle,
  "scwblatent::peakvi_project", inputs = list(model = model_ref,
    query_matrix = cfg$inputs$qualification_parent_bundle$query_matrix,
    space_id_out = query_space_id,
    chunk_size = as.integer(cfg$peakvi$chunk_size), batch_key = NULL,
    use_gpu = FALSE,
    pickle_execution_trust = cfg$peakvi$pickle_execution_trust,
    cell_set = "all", min_cells = 1L),
  notes = "Balanced KO fixture projection through frozen 20-dimensional model")
projection_summary <- project$run$outputs$projection_summary$table
stopifnot(nrow(projection_summary) == 1L,
  projection_summary$transformation[[1L]] == "transfer_without_refit",
  as.integer(projection_summary$n_latent[[1L]]) == cfg$peakvi$n_latent)
project$session <- wb_session_resume(bundle, project$session$durable$workspace)
project_promoted <- wb_session_promote_current(bundle, project$session,
  outputs = "projection_summary",
  name_prefix = "gated_transfer_hvp20k_20d_query_fixture",
  description = "Balanced KO fixture projected without reference refit",
  created_by = "Bdev development reference mapping step 16", fail_fast = TRUE)
projection_artifact_id <- promoted_id(project_promoted, "projection_summary")
projection_ref <- wb_ref_artifact(projection_artifact_id, "obs_table")
project_promoted$session <- wb_session_discard(project_promoted$session)

joint_space_id <- "gated_transfer_hvp20k_20d_reference_query_union"
union_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_latent_union"))
union <- wb_session_run_op(union_session, bundle,
  "scwblatent::model_bound_latent_union", inputs = list(
    reference_space = reference_space_id, query_space = query_space_id,
    model = model_ref, projection_summary = projection_ref,
    space_id_out = joint_space_id, chunk_size = 1024L),
  notes = "Authenticated 20-dimensional Bdev/KO qualification union")
union$session <- wb_session_resume(bundle, union$session$durable$workspace)
union_promoted <- wb_session_promote_current(bundle, union$session,
  outputs = "union_summary", name_prefix = "gated_transfer_hvp20k_20d_union",
  description = "Authenticated 20-dimensional Bdev/KO qualification union",
  created_by = "Bdev development reference mapping step 16", fail_fast = TRUE)
union_artifact_id <- promoted_id(union_promoted, "union_summary")
union_promoted$session <- wb_session_discard(union_promoted$session)

reference_annotation <- wb_prim_cell_vector(stats::setNames(domain,
  reference_cells), value_type = "categorical")
reference_l4 <- wb_prim_cell_vector_family(rep(cfg$calibration$continuous_key,
  length(reference_cells)), reference_cells, native_l4,
  value_type = "continuous")

log_event("Calibrating five balanced reference folds")
calibration_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_calibrate"))
calibration <- wb_session_run_op(calibration_session, bundle,
  "scwbreferencemap::reference_map_calibrate", inputs = list(
    reference_matrix = cfg$inputs$qualification_parent_bundle$authorizer_matrix,
    reference_space = reference_space_id,
    reference_annotation = reference_annotation,
    reference_continuous = reference_l4,
    holdout_group_field = cfg$calibration$holdout_group_field,
    k_values = as.integer(unlist(cfg$calibration$k_values)),
    primary_k = as.integer(cfg$calibration$primary_k),
    support_quantiles = as.list(as.numeric(unlist(
      cfg$calibration$support_quantiles))),
    primary_continuous_key = cfg$calibration$continuous_key,
    space_chunk_size = 4096L, query_chunk_size = 512L),
  notes = "Five balanced sample-domain-native-L4 reference folds")
heldout <- calibration$run$outputs$heldout_predictions$table
calibration_metrics <- calibration$run$outputs$calibration_metrics$table
thresholds <- calibration$run$outputs$thresholds$table
coverage <- calibration$run$outputs$coverage_limitations$table
stopifnot(nrow(heldout) == length(reference_cells),
  identical(sort(unique(as.character(heldout$holdout_category))),
    sort(unique(fold_label))), nrow(coverage) == n_folds,
  all(!nzchar(as.character(coverage$absent_annotation_categories))))
calibration$session <- wb_session_resume(bundle,
  calibration$session$durable$workspace)
calibration_promoted <- wb_session_promote_current(bundle,
  calibration$session, outputs = names(calibration$run$outputs),
  name_prefix = "gated_transfer_balanced_fold_calibration",
  description = "Five-fold balanced support and L4-domain calibration",
  created_by = "Bdev development reference mapping step 16", fail_fast = TRUE)
calibration_artifact_id <- promoted_id(calibration_promoted, "calibration")
calibration_ref <- wb_ref_artifact(calibration_artifact_id, "model")
calibration_promoted$session <- wb_session_discard(
  calibration_promoted$session)

log_event("Finding 30 neighbors and applying the generic continuous gate")
knn_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_knn"))
knn <- wb_session_run_op(knn_session, bundle,
  "scwbrepresentation::knn_neighbors", inputs = list(space = joint_space_id,
    query_cells = query_cells, candidate_cells = reference_cells,
    k = as.integer(cfg$calibration$primary_k), include_query = FALSE,
    space_chunk_size = 4096L),
  notes = "Balanced KO fixture to 20-dimensional HVP20k Bdev reference")

transfer_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_transfer"))
transfer <- wb_session_run_op(transfer_session, bundle,
  "scwbreferencemap::reference_map_transfer", inputs = list(
    in_space = query_space_id, reference_space = reference_space_id,
    neighbors = knn$run$outputs$neighbors, calibration = calibration_ref,
    support_output_key = "gated_transfer_support",
    annotation_output_key = "mapped_l4_domain",
    continuous_output_keys = cfg$calibration$continuous_output_key,
    continuous_acceptance = list(mode = "support_and_annotation",
      accepted_annotation_values = cfg$calibration$in_domain_value,
      minimum_annotation_confidence =
        as.numeric(cfg$calibration$minimum_annotation_confidence)),
    supported_threshold = as.numeric(thresholds$threshold[
      thresholds$role == "supported_upper"][[1L]]),
    out_of_reference_threshold = as.numeric(thresholds$threshold[
      thresholds$role == "borderline_upper"][[1L]]),
    cell_ids = query_cells, min_cells = 1L, space_chunk_size = 4096L),
  notes = "Generic support-and-annotation gate for native-L4 transfer")

prediction <- transfer$run$outputs$per_query_summary$table
categorical <- transfer$run$outputs$categorical_estimates$table
continuous <- transfer$run$outputs$continuous_estimates$table
acceptance <- transfer$run$outputs$continuous_acceptance$table
numeric_outputs <- transfer$run$outputs$numeric_outputs
stopifnot(nrow(prediction) == length(query_cells),
  identical(as.character(prediction$query_cell_id), query_cells),
  nrow(acceptance) == length(query_cells),
  identical(as.character(acceptance$query_cell_id), query_cells),
  !anyNA(acceptance$acceptance_reason),
  all(nzchar(as.character(acceptance$acceptance_reason))),
  all(acceptance$accepted == is.finite(acceptance$accepted_value)),
  all(is.na(acceptance$accepted_value[!acceptance$accepted])))

primary_categorical <- categorical[categorical$k ==
  as.integer(cfg$calibration$primary_k), ]
primary_continuous <- continuous[continuous$k ==
  as.integer(cfg$calibration$primary_k) & continuous$reference_key ==
  cfg$calibration$continuous_key, ]
stopifnot(nrow(primary_categorical) == length(query_cells),
  nrow(primary_continuous) == length(query_cells))
primary_categorical <- primary_categorical[match(query_cells,
  primary_categorical$query_cell_id), ]
primary_continuous <- primary_continuous[match(query_cells,
  primary_continuous$query_cell_id), ]
stopifnot(identical(as.character(primary_categorical$query_cell_id), query_cells),
  identical(as.character(primary_continuous$query_cell_id), query_cells))
expected_accept <- prediction$support_call != "out_of_reference" &
  primary_categorical$raw_prediction == cfg$calibration$in_domain_value &
  is.finite(primary_categorical$confidence) &
  primary_categorical$confidence >=
    as.numeric(cfg$calibration$minimum_annotation_confidence) &
  is.finite(primary_continuous$raw_estimate)
stopifnot(identical(as.logical(acceptance$accepted), expected_accept))

accepted_numeric_cells <- as.character(numeric_outputs$cell_id[
  numeric_outputs$key == cfg$calibration$continuous_output_key &
  is.finite(as.numeric(numeric_outputs$value))])
stopifnot(setequal(accepted_numeric_cells,
  query_cells[acceptance$accepted]))

slot_map <- fread(cfg$inputs$query_slot_map$path)
stopifnot(identical(as.character(slot_map$slot_id), query_cells))
query_result <- data.table(slot_id = query_cells,
  support_call = as.character(prediction$support_call),
  predicted_domain = as.character(primary_categorical$raw_prediction),
  domain_confidence = as.numeric(primary_categorical$confidence),
  raw_native_l4 = as.numeric(primary_continuous$raw_estimate),
  accepted_native_l4 = as.numeric(acceptance$accepted_value),
  accepted = as.logical(acceptance$accepted),
  acceptance_reason = as.character(acceptance$acceptance_reason))
query_result <- merge(slot_map, query_result, by = "slot_id", sort = FALSE)
query_result <- query_result[match(query_cells, slot_id)]

fwrite(heldout, file.path(dirs[["tables"]],
  "03_balanced_fold_heldout_predictions.tsv.gz"), sep = "\t", quote = TRUE,
  na = "NA")
fwrite(calibration_metrics, file.path(dirs[["tables"]],
  "04_calibration_metrics.tsv"), sep = "\t", quote = TRUE, na = "NA")
fwrite(thresholds, file.path(dirs[["tables"]],
  "05_support_thresholds.tsv"), sep = "\t", quote = TRUE, na = "NA")
fwrite(coverage, file.path(dirs[["tables"]],
  "06_calibration_coverage.tsv"), sep = "\t", quote = TRUE, na = "NA")
fwrite(query_result, file.path(dirs[["tables"]],
  "07_gated_fixture_mapping.tsv.gz"), sep = "\t", quote = TRUE, na = "NA")
fwrite(acceptance, file.path(dirs[["tables"]],
  "08_continuous_acceptance.tsv"), sep = "\t", quote = TRUE, na = "NA")

knn$session <- wb_session_resume(bundle, knn$session$durable$workspace)
transfer$session <- wb_session_resume(bundle, transfer$session$durable$workspace)
knn$session <- wb_session_discard(knn$session)
transfer$session <- wb_session_discard(transfer$session)

reason_counts <- as.list(table(factor(acceptance$acceptance_reason,
  levels = c("accepted", "out_of_reference", "annotation_not_accepted",
    "annotation_confidence_too_low", "continuous_value_unavailable"))))
report <- list(
  schema_version = "bdev.gated_continuous_transfer_qualification_report.v1",
  status = "passed",
  gate_passed = TRUE,
  run_identity = cfg$run_identity,
  reference_cells = length(reference_cells),
  query_cells = length(query_cells),
  selected_features = as.integer(reference_matrix_row$n_features[[1L]]),
  latent_dimensions = as.integer(cfg$peakvi$n_latent),
  model_artifact_id = model_artifact_id,
  projection_artifact_id = projection_artifact_id,
  union_artifact_id = union_artifact_id,
  calibration_artifact_id = calibration_artifact_id,
  fold_count = n_folds,
  fold_strategy = "deterministic_round_robin_within_sample_domain_global_native_l4_quintile",
  in_domain_reference_cells = sum(in_domain),
  out_domain_reference_cells = sum(!in_domain),
  minimum_l4_weight = as.numeric(cfg$calibration$minimum_l4_weight),
  primary_k = as.integer(cfg$calibration$primary_k),
  minimum_annotation_confidence =
    as.numeric(cfg$calibration$minimum_annotation_confidence),
  accepted_query_cells = sum(acceptance$accepted),
  acceptance_reason_counts = reason_counts,
  raw_estimate_rows = nrow(acceptance),
  source_bundle_manifest_sha256 = sha256(source_manifest),
  qualification_bundle_manifest_sha256 = sha256(file.path(
    cfg$storage$qualification_bundle, "manifest.json")),
  genotype_contrast_performed = FALSE,
  full_query_mapping_performed = FALSE,
  paper_analysis = FALSE,
  sif_built = FALSE)
write_json(report, file.path(dirs[["reports"]],
  "01_gated_continuous_transfer_qualification.json"))
write_json(list(
  schema_version = "bdev.gated_continuous_transfer_execution.v1",
  config_path = config_path,
  config_sha256 = sha256(config_path),
  scworkbench = as.character(packageVersion("scworkbench")),
  scwbLatent = as.character(packageVersion("scwbLatent")),
  scwbRepresentation = as.character(packageVersion("scwbRepresentation")),
  scwbReferenceMap = as.character(packageVersion("scwbReferenceMap")),
  scwbReferenceMap_manifest_sha256 = sha256(system.file("scworkbench",
    "operation-pack.json", package = "scwbReferenceMap")),
  source_manifest_unchanged = sha256(source_manifest) ==
    cfg$inputs$qualification_parent_bundle$manifest_sha256,
  source_archives_modified = FALSE), file.path(dirs[["manifests"]],
  "02_execution_receipt.json"))
log_event("Gated continuous-transfer qualification passed")
