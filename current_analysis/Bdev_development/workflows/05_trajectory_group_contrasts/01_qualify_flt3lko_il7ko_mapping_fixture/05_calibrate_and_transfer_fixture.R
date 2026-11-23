#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: 05_calibrate_and_transfer_fixture.R BUNDLE SLOT_MAP RESULTS_ROOT SESSIONS_DIR",
       call. = FALSE)
}

bundle_path <- normalizePath(args[[1L]], mustWork = TRUE)
slot_map_path <- normalizePath(args[[2L]], mustWork = TRUE)
results_root <- normalizePath(args[[3L]], mustWork = FALSE)
sessions_dir <- normalizePath(args[[4L]], mustWork = FALSE)

options(stringsAsFactors = FALSE, warn = 1, width = 180)
suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbLatent)
  library(scwbRepresentation)
  library(scwbReferenceMap)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(digest)
})

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
promoted_ref <- function(promoted, output, kind) {
  take <- promoted$promoted$status == "ok" & promoted$promoted$output == output
  stopifnot(sum(take) == 1L)
  wb_ref_artifact(as.character(promoted$promoted$artifact_id[take][[1L]]), kind)
}

dirs <- setNames(file.path(results_root, c("01_tables", "03_reports", "04_logs", "06_manifests")),
                 c("tables", "reports", "logs", "manifests"))
invisible(lapply(c(dirs, sessions = sessions_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "03_support_transfer.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

for (package in c("scwbLatent", "scwbRepresentation", "scwbReferenceMap")) {
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  stopifnot(nzchar(manifest))
  namespace <- switch(package, scwbLatent = "scwblatent",
                      scwbRepresentation = "scwbrepresentation",
                      scwbReferenceMap = "scwbreferencemap")
  try(wb_operation_pack_unload(namespace), silent = TRUE)
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
}
required_ops <- c("scwblatent::model_bound_latent_union",
  "scwbrepresentation::knn_neighbors",
  "scwbreferencemap::reference_map_calibrate",
  "scwbreferencemap::reference_map_transfer")
stopifnot(all(required_ops %in% wb_list_ops()))

bundle <- wb_open_bundle(bundle_path, validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
all_cells <- as.character(wb_cells(bundle, "all"))
query_cells <- all_cells[grepl("^external_query_[0-9]+$", all_cells)]
reference_cells <- setdiff(all_cells, query_cells)
stopifnot(length(query_cells) == 28L, length(reference_cells) == 26723L)
slot_map <- fread(slot_map_path)
stopifnot(identical(as.character(slot_map$slot_id), query_cells),
          all(table(slot_map$library_id) == 7L))

reference_space_id <- "peakvi_space_id_n0cp39pv1pyi"
query_space_id <- "frozen_hvp_depth20k_flt3lko_il7ko_fixture"
joint_space_id <- "hvp_depth20k_bdev_flt3lko_il7ko_fixture_union"
model_artifact_id <- "model_ref_20261117_051703_wi5ql81y"
projection_artifact_id <- "obs_table_20261117_181001_2fw4603c"
model_ref <- wb_ref_artifact(model_artifact_id, "model")
projection_ref <- wb_ref_artifact(projection_artifact_id, "obs_table")

if (!(joint_space_id %in% wb_spaces(bundle)$space_id)) {
  union_session <- wb_session_create(bundle, sessions_dir,
    "flt3lko_il7ko_mapping_fixture_20261117T175427Z_latent_union_retry")
  union <- wb_session_run_op(union_session, bundle,
    "scwblatent::model_bound_latent_union", inputs = list(
      reference_space = reference_space_id,
      query_space = query_space_id,
      model = model_ref,
      projection_summary = projection_ref,
      space_id_out = joint_space_id,
      chunk_size = 1024L),
    notes = "Authenticated frozen-model Bdev/KO bounded-fixture latent union")
  union$session <- wb_session_resume(bundle, union$session$durable$workspace)
  union_promoted <- wb_session_promote_current(bundle, union$session,
    outputs = "union_summary", name_prefix = "flt3lko_il7ko_fixture_union",
    description = "Authenticated Bdev/KO bounded-fixture latent union",
    created_by = "Bdev development capability 05 bounded mapping qualification",
    fail_fast = TRUE)
  union_promoted$session <- wb_session_discard(union_promoted$session)
}

metadata <- wb_meta(bundle, cells = reference_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "orig.ident", "native_l4"))
stopifnot(identical(as.character(metadata$cell_id), reference_cells),
          length(unique(as.character(metadata$orig.ident))) == 4L,
          sum(is.finite(as.numeric(metadata$native_l4))) > 10000L)
reference_l4 <- wb_prim_cell_vector_family(
  rep("L4", length(reference_cells)), reference_cells,
  as.numeric(metadata$native_l4), value_type = "continuous")

# The calibration operation uses its matrix input only to authenticate the
# reference cell axis and requested metadata field. A one-feature sparse matrix
# avoids materializing the unrelated 20k-peak count payload during that
# authorization step while retaining the exact 26,723-cell reference axis.
authorizer_matrix_id <- "hvp_depth20k_reference_mapping_authorizer"
if (!(authorizer_matrix_id %in% wb_matrices(bundle)$matrix_id)) {
  authorizer <- Matrix::sparseMatrix(
    i = rep.int(1L, length(reference_cells)),
    j = seq_along(reference_cells), x = 1,
    dims = c(1L, length(reference_cells)),
    dimnames = list("mapping_authorizer_peak", reference_cells))
  wb_attach_matrix(bundle, authorizer_matrix_id, authorizer,
    layer = "counts", feature_type = "peak", cell_ids = reference_cells,
    description = "One-feature sparse authorizer for Bdev mapping metadata and cell axis",
    source = "derived_from_frozen_hvp_depth20k_reference_cell_axis",
    created_by = "Bdev development capability 05 bounded mapping qualification")
}

artifacts <- wb_art_list(bundle)
calibration_take <- artifacts$source_op ==
  "scwbreferencemap::reference_map_calibrate" &
  artifacts$type == "model_ref" & artifacts$lifecycle_state == "current" &
  grepl("^flt3lko_il7ko_fixture_calibration - ", artifacts$name)
if (sum(calibration_take) == 1L) {
  calibration_id <- as.character(artifacts$artifact_id[calibration_take][[1L]])
  calibration_ref <- wb_ref_artifact(calibration_id, "model")
  threshold_take <- artifacts$source_op ==
    "scwbreferencemap::reference_map_calibrate" &
    artifacts$type == "obs_table" & artifacts$lifecycle_state == "current" &
    grepl("^flt3lko_il7ko_fixture_calibration - ", artifacts$name) &
    grepl("\\(thresholds\\)$", artifacts$name)
  metric_take <- artifacts$source_op ==
    "scwbreferencemap::reference_map_calibrate" &
    artifacts$type == "obs_table" & artifacts$lifecycle_state == "current" &
    grepl("^flt3lko_il7ko_fixture_calibration - ", artifacts$name) &
    grepl("\\(calibration_metrics\\)$", artifacts$name)
  stopifnot(sum(threshold_take) == 1L, sum(metric_take) == 1L)
  thresholds <- wb_art_load(bundle,
    as.character(artifacts$artifact_id[threshold_take][[1L]]),
    expected_kind = "obs_table")$table
  calibration_metrics <- wb_art_load(bundle,
    as.character(artifacts$artifact_id[metric_take][[1L]]),
    expected_kind = "obs_table")$table
} else {
  stopifnot(sum(calibration_take) == 0L)
  calibration_session <- wb_session_create(bundle, sessions_dir,
    "flt3lko_il7ko_mapping_fixture_20261117T175427Z_calibrate_retry1")
  calibration <- wb_session_run_op(calibration_session, bundle,
    "scwbreferencemap::reference_map_calibrate", inputs = list(
      reference_matrix = authorizer_matrix_id,
      reference_space = reference_space_id,
      reference_continuous = reference_l4,
      holdout_group_field = "orig.ident",
      k_values = as.integer(c(15L, 30L, 50L)),
      primary_k = 30L,
      support_quantiles = as.list(c(0.95, 0.99)),
      primary_continuous_key = "L4",
      space_chunk_size = 4096L,
      query_chunk_size = 512L),
    notes = "Bdev-only HVP-depth20k support and native-L4 calibration")
  calibration$session <- wb_session_resume(bundle,
    calibration$session$durable$workspace)
  calibration_promoted <- wb_session_promote_current(bundle,
    calibration$session, outputs = names(calibration$run$outputs),
    name_prefix = "flt3lko_il7ko_fixture_calibration",
    description = "Bdev-only HVP-depth20k mapping calibration",
    created_by = "Bdev development capability 05 bounded mapping qualification",
    fail_fast = TRUE)
  calibration_ref <- promoted_ref(calibration_promoted, "calibration", "model")
  thresholds <- calibration$run$outputs$thresholds$table
  calibration_metrics <- calibration$run$outputs$calibration_metrics$table
  calibration_promoted$session <- wb_session_discard(
    calibration_promoted$session)
}
supported_threshold <- thresholds$threshold[
  thresholds$role == "supported_upper"][[1L]]
out_threshold <- thresholds$threshold[
  thresholds$role == "borderline_upper"][[1L]]

knn_session <- wb_session_create(bundle, sessions_dir,
  "flt3lko_il7ko_mapping_fixture_20261117T175427Z_knn_retry2")
knn <- wb_session_run_op(knn_session, bundle,
  "scwbrepresentation::knn_neighbors", inputs = list(
    space = joint_space_id,
    query_cells = query_cells,
    candidate_cells = reference_cells,
    k = 50L,
    include_query = FALSE,
    space_chunk_size = 4096L),
  notes = "Bounded KO fixture to frozen HVP-depth20k Bdev reference neighbors")

transfer_session <- wb_session_create(bundle, sessions_dir,
  "flt3lko_il7ko_mapping_fixture_20261117T175427Z_transfer_retry1")
transfer <- wb_session_run_op(transfer_session, bundle,
  "scwbreferencemap::reference_map_transfer", inputs = list(
    in_space = query_space_id,
    reference_space = reference_space_id,
    neighbors = knn$run$outputs$neighbors,
    calibration = calibration_ref,
    support_output_key = "flt3lko_il7ko_fixture_support",
    annotation_output_key = "unused_no_categorical_transfer",
    continuous_output_keys = "mapped_native_L4",
    supported_threshold = as.numeric(supported_threshold),
    out_of_reference_threshold = as.numeric(out_threshold),
    cell_ids = query_cells,
    min_cells = 1L,
    space_chunk_size = 4096L),
  notes = "Support-aware native-L4 transfer for bounded KO fixture")

prediction <- transfer$run$outputs$per_query_summary$table
continuous <- transfer$run$outputs$continuous_estimates$table
primary <- continuous$k == 30L & continuous$reference_key == "L4"
primary_continuous <- continuous[primary, c("query_cell_id", "weighted_q25",
  "raw_estimate", "weighted_q75", "iqr", "cross_k_sensitivity")]
names(primary_continuous)[names(primary_continuous) == "raw_estimate"] <- "raw_native_l4"
result <- merge(slot_map, prediction, by.x = "slot_id", by.y = "query_cell_id", sort = FALSE)
result <- merge(result, primary_continuous, by.x = "slot_id", by.y = "query_cell_id", sort = FALSE)
result <- result[match(query_cells, result$slot_id), ]
result$accepted_native_l4 <- ifelse(result$support_call == "out_of_reference",
                                    NA_real_, result$raw_native_l4)
stopifnot(nrow(result) == 28L,
          all(is.na(result$accepted_native_l4[result$support_call == "out_of_reference"])))

fwrite(result, file.path(dirs[["tables"]], "07_support_aware_native_l4.tsv.gz"),
       sep = "\t", quote = TRUE, na = "NA")
fwrite(thresholds, file.path(dirs[["tables"]], "08_support_thresholds.tsv"),
       sep = "\t", quote = TRUE, na = "NA")
fwrite(calibration_metrics,
       file.path(dirs[["tables"]], "09_calibration_metrics.tsv"),
       sep = "\t", quote = TRUE, na = "NA")

knn$session <- wb_session_resume(bundle, knn$session$durable$workspace)
transfer$session <- wb_session_resume(bundle, transfer$session$durable$workspace)
knn$session <- wb_session_discard(knn$session)
transfer$session <- wb_session_discard(transfer$session)

support_counts <- as.list(table(factor(result$support_call,
  levels = c("supported", "borderline", "out_of_reference"))))
group_counts <- result[, .(
  cells = .N,
  supported = sum(support_call == "supported"),
  borderline = sum(support_call == "borderline"),
  out_of_reference = sum(support_call == "out_of_reference"),
  accepted_l4 = sum(is.finite(accepted_native_l4))
), by = .(knockout_target, genotype, library_id)]
fwrite(group_counts, file.path(dirs[["tables"]], "10_support_by_library.tsv"),
       sep = "\t", quote = TRUE, na = "NA")
report <- list(
  status = "passed",
  query_cells = nrow(result),
  reference_cells = length(reference_cells),
  support_counts = support_counts,
  supported_threshold = supported_threshold,
  out_of_reference_threshold = out_threshold,
  accepted_native_l4 = sum(is.finite(result$accepted_native_l4)),
  categorical_transfer = FALSE,
  trajectory_key = "L4",
  trajectory_orientation = "native",
  genotype_contrast_performed = FALSE,
  source_archives_modified = FALSE,
  paper_analysis = FALSE
)
write_json(report, file.path(dirs[["reports"]], "05_support_transfer_audit.json"))
write_json(list(
  scworkbench = as.character(packageVersion("scworkbench")),
  scwbLatent = as.character(packageVersion("scwbLatent")),
  scwbRepresentation = as.character(packageVersion("scwbRepresentation")),
  scwbReferenceMap = as.character(packageVersion("scwbReferenceMap")),
  bundle_manifest_sha256 = sha256(file.path(bundle_path, "manifest.json"))
), file.path(dirs[["manifests"]], "04_support_transfer_runtime.json"))
