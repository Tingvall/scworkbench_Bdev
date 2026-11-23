#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_index <- match("--config", args)
if (is.na(config_index) || config_index == length(args)) {
  stop("Usage: 03_resume_gated_transfer.R --config CONFIG", call. = FALSE)
}
config_path <- normalizePath(args[[config_index + 1L]], mustWork = TRUE)

options(stringsAsFactors = FALSE, warn = 1, width = 180)
suppressPackageStartupMessages({
  library(scworkbench)
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

root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "03_reports", "04_logs",
  "06_manifests")), c("tables", "reports", "logs", "manifests"))
sink(file.path(dirs[["logs"]], "02_resume_gated_transfer.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

for (package in c("scwbRepresentation", "scwbReferenceMap")) {
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  stopifnot(nzchar(manifest), sha256(manifest) == cfg$runtime[[package]]$manifest_sha256)
  namespace <- switch(package,
    scwbRepresentation = "scwbrepresentation",
    scwbReferenceMap = "scwbreferencemap")
  try(wb_operation_pack_unload(namespace), silent = TRUE)
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
}

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

spaces <- wb_spaces(bundle)
reference_candidates <- spaces$source == "operation_pack_v3_matrix_model" &
  spaces$dims == cfg$peakvi$n_latent
stopifnot(sum(reference_candidates) == 1L)
reference_space_id <- as.character(spaces$space_id[reference_candidates][[1L]])
query_space_id <- "gated_transfer_hvp20k_20d_query_fixture"
joint_space_id <- "gated_transfer_hvp20k_20d_reference_query_union"
stopifnot(all(c(query_space_id, joint_space_id) %in% spaces$space_id),
  all(spaces$dims[match(c(reference_space_id, query_space_id, joint_space_id),
    spaces$space_id)] == cfg$peakvi$n_latent))

artifacts <- wb_art_list(bundle)
current <- artifacts$lifecycle_state == "current"
take_one <- function(pattern, type, source_op) {
  take <- current & artifacts$type == type & artifacts$source_op == source_op &
    grepl(pattern, artifacts$name)
  stopifnot(sum(take) == 1L)
  as.character(artifacts$artifact_id[take][[1L]])
}
model_artifact_id <- take_one(
  "^gated_transfer_hvp20k_20d_reference .*\\(model\\)$", "model_ref",
  "scwblatent::peakvi_fit")
projection_artifact_id <- take_one(
  "^gated_transfer_hvp20k_20d_query_fixture .*\\(projection_summary\\)$",
  "obs_table", "scwblatent::peakvi_project")
union_artifact_id <- take_one(
  "^gated_transfer_hvp20k_20d_union .*\\(union_summary\\)$", "obs_table",
  "scwblatent::model_bound_latent_union")
calibration_artifact_id <- take_one(
  "^gated_transfer_balanced_fold_calibration .*\\(calibration\\)$",
  "model_ref", "scwbreferencemap::reference_map_calibrate")
threshold_artifact_id <- take_one(
  "^gated_transfer_balanced_fold_calibration .*\\(thresholds\\)$",
  "obs_table", "scwbreferencemap::reference_map_calibrate")
heldout_artifact_id <- take_one(
  "^gated_transfer_balanced_fold_calibration .*\\(heldout_predictions\\)$",
  "obs_table", "scwbreferencemap::reference_map_calibrate")
metrics_artifact_id <- take_one(
  "^gated_transfer_balanced_fold_calibration .*\\(calibration_metrics\\)$",
  "obs_table", "scwbreferencemap::reference_map_calibrate")
coverage_artifact_id <- take_one(
  "^gated_transfer_balanced_fold_calibration .*\\(coverage_limitations\\)$",
  "obs_table", "scwbreferencemap::reference_map_calibrate")

thresholds <- wb_art_load(bundle, threshold_artifact_id,
  expected_kind = "obs_table")$table
heldout <- wb_art_load(bundle, heldout_artifact_id,
  expected_kind = "obs_table")$table
calibration_metrics <- wb_art_load(bundle, metrics_artifact_id,
  expected_kind = "obs_table")$table
coverage <- wb_art_load(bundle, coverage_artifact_id,
  expected_kind = "obs_table")$table
supported_threshold <- as.numeric(thresholds$threshold[
  thresholds$role == "supported_upper"][[1L]])
out_threshold <- as.numeric(thresholds$threshold[
  thresholds$role == "borderline_upper"][[1L]])
stopifnot(nrow(heldout) == length(reference_cells),
  nrow(coverage) == cfg$calibration$folds,
  all(!nzchar(as.character(coverage$absent_annotation_categories))),
  is.finite(supported_threshold), is.finite(out_threshold),
  out_threshold >= supported_threshold)

cat(format(Sys.time(), tz = "UTC", usetz = TRUE),
  "Resuming exact neighbors and gated transfer with explicit thresholds\n")
knn_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_knn_resume1"))
knn <- wb_session_run_op(knn_session, bundle,
  "scwbrepresentation::knn_neighbors", inputs = list(space = joint_space_id,
    query_cells = query_cells, candidate_cells = reference_cells,
    k = as.integer(cfg$calibration$primary_k), include_query = FALSE,
    space_chunk_size = 4096L),
  notes = "Resumed balanced KO fixture exact neighbors")

transfer_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_transfer_resume1"))
transfer <- wb_session_run_op(transfer_session, bundle,
  "scwbreferencemap::reference_map_transfer", inputs = list(
    in_space = query_space_id, reference_space = reference_space_id,
    neighbors = knn$run$outputs$neighbors,
    calibration = wb_ref_artifact(calibration_artifact_id, "model"),
    support_output_key = "gated_transfer_support",
    annotation_output_key = "mapped_l4_domain",
    continuous_output_keys = cfg$calibration$continuous_output_key,
    continuous_acceptance = list(mode = "support_and_annotation",
      accepted_annotation_values = cfg$calibration$in_domain_value,
      minimum_annotation_confidence =
        as.numeric(cfg$calibration$minimum_annotation_confidence)),
    supported_threshold = supported_threshold,
    out_of_reference_threshold = out_threshold,
    cell_ids = query_cells, min_cells = 1L, space_chunk_size = 4096L),
  notes = "Resumed generic support-and-annotation native-L4 gate")

prediction <- transfer$run$outputs$per_query_summary$table
categorical <- transfer$run$outputs$categorical_estimates$table
continuous <- transfer$run$outputs$continuous_estimates$table
acceptance <- transfer$run$outputs$continuous_acceptance$table
numeric_outputs <- transfer$run$outputs$numeric_outputs
primary_categorical <- categorical[categorical$k ==
  as.integer(cfg$calibration$primary_k), ]
primary_continuous <- continuous[continuous$k ==
  as.integer(cfg$calibration$primary_k) & continuous$reference_key ==
  cfg$calibration$continuous_key, ]
primary_categorical <- primary_categorical[match(query_cells,
  primary_categorical$query_cell_id), ]
primary_continuous <- primary_continuous[match(query_cells,
  primary_continuous$query_cell_id), ]
stopifnot(nrow(prediction) == length(query_cells),
  identical(as.character(prediction$query_cell_id), query_cells),
  nrow(acceptance) == length(query_cells),
  identical(as.character(acceptance$query_cell_id), query_cells),
  identical(as.character(primary_categorical$query_cell_id), query_cells),
  identical(as.character(primary_continuous$query_cell_id), query_cells),
  !anyNA(acceptance$acceptance_reason),
  all(acceptance$accepted == is.finite(acceptance$accepted_value)),
  all(is.na(acceptance$accepted_value[!acceptance$accepted])))

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

metadata <- wb_meta(bundle, cells = reference_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$calibration$domain_annotation_field))
domain <- as.character(metadata[[cfg$calibration$domain_annotation_field]])
reason_counts <- as.list(table(factor(acceptance$acceptance_reason,
  levels = c("accepted", "out_of_reference", "annotation_not_accepted",
    "annotation_confidence_too_low", "continuous_value_unavailable"))))
report <- list(
  schema_version = "bdev.gated_continuous_transfer_qualification_report.v1",
  status = "passed", gate_passed = TRUE, run_identity = cfg$run_identity,
  reference_cells = length(reference_cells), query_cells = length(query_cells),
  selected_features = 20000L, latent_dimensions = as.integer(cfg$peakvi$n_latent),
  model_artifact_id = model_artifact_id,
  projection_artifact_id = projection_artifact_id,
  union_artifact_id = union_artifact_id,
  calibration_artifact_id = calibration_artifact_id,
  fold_count = as.integer(cfg$calibration$folds),
  fold_strategy = "deterministic_round_robin_within_sample_domain_global_native_l4_quintile",
  in_domain_reference_cells = sum(domain == cfg$calibration$in_domain_value),
  out_domain_reference_cells = sum(domain == cfg$calibration$out_domain_value),
  minimum_l4_weight = as.numeric(cfg$calibration$minimum_l4_weight),
  primary_k = as.integer(cfg$calibration$primary_k),
  supported_threshold = supported_threshold,
  out_of_reference_threshold = out_threshold,
  minimum_annotation_confidence =
    as.numeric(cfg$calibration$minimum_annotation_confidence),
  accepted_query_cells = sum(acceptance$accepted),
  acceptance_reason_counts = reason_counts,
  raw_estimate_rows = nrow(acceptance),
  source_bundle_manifest_sha256 = sha256(source_manifest),
  qualification_bundle_manifest_sha256 = sha256(file.path(
    cfg$storage$qualification_bundle, "manifest.json")),
  resumed_after_explicit_retained_threshold_gate = TRUE,
  genotype_contrast_performed = FALSE, full_query_mapping_performed = FALSE,
  paper_analysis = FALSE, sif_built = FALSE)
write_json(report, file.path(dirs[["reports"]],
  "01_gated_continuous_transfer_qualification.json"))
write_json(list(
  schema_version = "bdev.gated_continuous_transfer_execution.v1",
  config_path = config_path, config_sha256 = sha256(config_path),
  scworkbench = as.character(packageVersion("scworkbench")),
  scwbRepresentation = as.character(packageVersion("scwbRepresentation")),
  scwbReferenceMap = as.character(packageVersion("scwbReferenceMap")),
  scwbReferenceMap_manifest_sha256 = sha256(system.file("scworkbench",
    "operation-pack.json", package = "scwbReferenceMap")),
  source_manifest_unchanged = sha256(source_manifest) ==
    cfg$inputs$qualification_parent_bundle$manifest_sha256,
  source_archives_modified = FALSE), file.path(dirs[["manifests"]],
  "02_execution_receipt.json"))
cat("Resumed gated continuous-transfer qualification passed\n")

