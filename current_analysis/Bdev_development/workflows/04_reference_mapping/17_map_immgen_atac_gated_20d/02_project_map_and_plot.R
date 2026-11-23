#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_index <- match("--config", args)
if (is.na(config_index) || config_index == length(args)) {
  stop("Usage: 02_project_map_and_plot.R --config CONFIG", call. = FALSE)
}
config_path <- normalizePath(args[[config_index + 1L]], mustWork = TRUE)

options(stringsAsFactors = FALSE, warn = 1, width = 180)
suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbLatent)
  library(scwbRepresentation)
  library(scwbReferenceMap)
  library(Matrix)
  library(data.table)
  library(digest)
  library(jsonlite)
  library(yaml)
  library(slingshot)
  library(ggplot2)
  library(ggrepel)
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
  } else x <- as.matrix(x)
  x
}

root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")), c("tables", "figures",
  "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_project_map_and_plot.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

qualified_manifest <- file.path(cfg$inputs$qualified_bundle$path,
  "manifest.json")
bulk_manifest <- file.path(cfg$inputs$bulk_query_bundle$path, "manifest.json")
fixed_umap_manifest <- file.path(cfg$inputs$fixed_umap_bundle$path,
  "manifest.json")
stopifnot(
  sha256(cfg$runtime$sif) == cfg$runtime$sif_sha256,
  sha256(qualified_manifest) ==
    cfg$inputs$qualified_bundle$manifest_sha256,
  sha256(bulk_manifest) == cfg$inputs$bulk_query_bundle$manifest_sha256,
  sha256(fixed_umap_manifest) ==
    cfg$inputs$fixed_umap_bundle$manifest_sha256,
  sha256(cfg$inputs$library_manifest$path) ==
    cfg$inputs$library_manifest$sha256,
  sha256(cfg$inputs$qualification_report$path) ==
    cfg$inputs$qualification_report$sha256,
  sha256(cfg$inputs$authenticated_l4$model_path) ==
    cfg$inputs$authenticated_l4$model_sha256,
  sha256(cfg$inputs$authenticated_l4$pseudotime_path) ==
    cfg$inputs$authenticated_l4$pseudotime_sha256)

for (package in c("scwbLatent", "scwbRepresentation", "scwbReferenceMap")) {
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  stopifnot(nzchar(manifest),
    sha256(manifest) == cfg$runtime[[package]]$manifest_sha256)
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

libraries <- fread(cfg$inputs$library_manifest$path)
stopifnot(nrow(libraries) == cfg$scope$expected_query_libraries,
  !anyDuplicated(libraries$gsm), all(nzchar(libraries$population)),
  all(nzchar(libraries$lineage)))
slots <- sprintf("external_bulk20d_query_%03d", seq_len(nrow(libraries)))
slot_map <- data.table(slot_id = slots, libraries)
fwrite(slot_map, file.path(dirs[["tables"]], "01_bulk_query_slot_map.tsv"),
  sep = "\t", quote = TRUE, na = "NA")

bundle <- wb_open_bundle(cfg$storage$mapping_bundle,
  validate = FALSE, read_only = FALSE)
on.exit(if (!is.null(bundle)) try(wb_close_bundle(bundle), silent = TRUE),
  add = TRUE)
current_cells <- as.character(wb_cells(bundle, "all"))
bulk_pattern <- "^external_bulk20d_query_[0-9]+$"
source_cells <- current_cells[!grepl(bulk_pattern, current_cells)]
reference_cells <- source_cells[!grepl("^external_query_[0-9]+$", source_cells)]
stopifnot(length(reference_cells) == cfg$scope$expected_reference_cells)

reference_axis <- as.character(wb_feature_ids(bundle,
  cfg$inputs$qualified_bundle$reference_matrix))
stopifnot(length(reference_axis) == 20000L)
model <- wb_art_load(bundle, cfg$inputs$qualified_bundle$model_artifact_id,
  expected_kind = "model_ref")
model_payload <- file.path(model$model_tree_root, model$model_primary_path)
stopifnot(
  sha256(model_payload) == cfg$inputs$qualified_bundle$model_payload_sha256,
  as.integer(model$model_linkage$latent_dimension_count) == 20L,
  as.character(model$model_linkage$feature_order_identity$sha256) ==
    cfg$inputs$qualified_bundle$feature_order_sha256)
model_ref <- wb_ref_artifact(cfg$inputs$qualified_bundle$model_artifact_id,
  "model")

log_event("Reading exact HVP20k axis for", length(slots), "bulk libraries")
query_bundle <- wb_open_bundle(cfg$inputs$bulk_query_bundle$path,
  validate = FALSE, read_only = TRUE)
bulk_counts <- wb_slice(query_bundle,
  cfg$inputs$bulk_query_bundle$matrix, features = reference_axis,
  cells = as.character(libraries$gsm))
wb_close_bundle(query_bundle)
stopifnot(nrow(bulk_counts) == 20000L, ncol(bulk_counts) == length(slots),
  identical(as.character(rownames(bulk_counts)), reference_axis),
  identical(as.character(colnames(bulk_counts)), as.character(libraries$gsm)))
colnames(bulk_counts) <- slots

if (identical(current_cells, source_cells)) {
  query_meta <- data.frame(cell_id = slots,
    orig.ident = as.character(libraries$gsm),
    population = as.character(libraries$population),
    lineage = as.character(libraries$lineage),
    native_l4 = rep(NA_real_, length(slots)),
    l4_weight = rep(NA_real_, length(slots)),
    cell_role = rep("immgen_bulk_atac_gated20d_query", length(slots)),
    stringsAsFactors = FALSE)
  wb_ingest_apply(bundle,
    set_name = sub("^set:", "", cfg$mapping$query_cell_set),
    cell_ids = slots, metadata = query_meta,
    source = "checksummed_immgen_bulk_atac_exact_hvp20k_axis",
    created_by = "Bdev development reference mapping step 17",
    apply_mode = "atomic", emit_report_artifact = TRUE)
  current_cells <- as.character(wb_cells(bundle, "all"))
}
stopifnot(identical(current_cells, c(source_cells, slots)),
  identical(as.character(wb_cells(bundle, cfg$mapping$query_cell_set)), slots))

if (!(cfg$mapping$query_matrix_id %in% wb_matrices(bundle)$matrix_id)) {
  wb_attach_matrix(bundle, cfg$mapping$query_matrix_id,
    as(bulk_counts, "dgCMatrix"), layer = "counts", feature_type = "peak",
    cell_ids = slots,
    description = "ImmGen bulk ATAC on selected exact HVP-depth20k axis",
    source = "checksummed_existing_immgen_full_bdev_peak_counts",
    created_by = "Bdev development reference mapping step 17")
}
rm(bulk_counts); invisible(gc())

artifacts <- wb_art_list(bundle)
projection_take <- artifacts$source_op == "scwblatent::peakvi_project" &
  artifacts$type == "obs_table" & artifacts$lifecycle_state == "current" &
  grepl("^immgen_bulk_gated20d_projection - ", artifacts$name)
if (cfg$mapping$query_space_id %in% wb_spaces(bundle)$space_id &&
    sum(projection_take) == 1L) {
  projection_artifact_id <- as.character(
    artifacts$artifact_id[projection_take][[1L]])
  projection_ref <- wb_ref_artifact(projection_artifact_id, "obs_table")
  log_event("Recovered completed bulk projection without refit")
} else {
  stopifnot(!(cfg$mapping$query_space_id %in% wb_spaces(bundle)$space_id),
    sum(projection_take) == 0L)
  log_event("Projecting bulk libraries without model refit")
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
    notes = "ImmGen bulk projection through selected frozen 20D HVP20k model")
  projection_summary <- project$run$outputs$projection_summary$table
  stopifnot(nrow(projection_summary) == 1L,
    projection_summary$transformation[[1L]] == "transfer_without_refit",
    as.integer(projection_summary$n_latent[[1L]]) == 20L,
    as.integer(projection_summary$n_query_cells[[1L]]) == length(slots))
  project$session <- wb_session_resume(bundle, project$session$durable$workspace)
  project_promoted <- wb_session_promote_current(bundle, project$session,
    outputs = "projection_summary",
    name_prefix = "immgen_bulk_gated20d_projection",
    description = "ImmGen bulk projected through selected frozen 20D model",
    created_by = "Bdev development reference mapping step 17", fail_fast = TRUE)
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
      source = "selected_20d_reference_and_immgen_bulk_union",
      created_by = "Bdev development reference mapping step 17",
      description = "Selected reference and ImmGen bulk union", min_cells = 1L),
    notes = "Create authenticated ImmGen bulk/reference union")
  union_set$session <- wb_session_resume(bundle,
    union_set$session$durable$workspace)
  union_promoted <- wb_session_promote_current(bundle, union_set$session,
    outputs = c("manifest", "rule_table"),
    name_prefix = "immgen_bulk_gated20d_union_cell_set",
    description = "Selected reference and ImmGen bulk union",
    created_by = "Bdev development reference mapping step 17", fail_fast = TRUE)
  union_promoted$session <- wb_session_discard(union_promoted$session)
}
stopifnot(identical(as.character(wb_cells(bundle, cfg$mapping$union_cell_set)),
  c(reference_cells, slots)))

artifacts <- wb_art_list(bundle)
union_take <- artifacts$source_op == "scwblatent::model_bound_latent_union" &
  artifacts$type == "obs_table" & artifacts$lifecycle_state == "current" &
  grepl("^immgen_bulk_gated20d_union - ", artifacts$name)
if (cfg$mapping$joint_space_id %in% wb_spaces(bundle)$space_id &&
    sum(union_take) == 1L) {
  union_artifact_id <- as.character(artifacts$artifact_id[union_take][[1L]])
  log_event("Recovered authenticated reference/bulk latent union")
} else {
  stopifnot(!(cfg$mapping$joint_space_id %in% wb_spaces(bundle)$space_id),
    sum(union_take) == 0L)
  log_event("Authenticating reference/bulk latent union")
  union_session <- wb_session_create(bundle, cfg$storage$sessions,
    paste0(cfg$run_identity, "_latent_union"))
  union <- wb_session_run_op(union_session, bundle,
    "scwblatent::model_bound_latent_union", inputs = list(
      reference_space = cfg$inputs$qualified_bundle$reference_space,
      query_space = cfg$mapping$query_space_id, model = model_ref,
      projection_summary = projection_ref,
      space_id_out = cfg$mapping$joint_space_id, chunk_size = 1024L,
      cell_set = cfg$mapping$union_cell_set, min_cells = 1L),
    notes = "Authenticated selected 20D Bdev/ImmGen bulk latent union")
  union$session <- wb_session_resume(bundle, union$session$durable$workspace)
  union_promoted <- wb_session_promote_current(bundle, union$session,
    outputs = "union_summary", name_prefix = "immgen_bulk_gated20d_union",
    description = "Authenticated selected 20D Bdev/ImmGen bulk union",
    created_by = "Bdev development reference mapping step 17", fail_fast = TRUE)
  union_artifact_id <- promoted_id(union_promoted, "union_summary")
  union_promoted$session <- wb_session_discard(union_promoted$session)
}

artifacts <- wb_art_list(bundle)
calibration <- wb_art_load(bundle,
  cfg$inputs$qualified_bundle$calibration_artifact_id,
  expected_kind = "model_ref")
stopifnot(as.character(calibration$model_space_id) ==
  cfg$inputs$qualified_bundle$reference_space)
calibration_ref <- wb_ref_artifact(
  cfg$inputs$qualified_bundle$calibration_artifact_id, "model")
threshold_take <- artifacts$source_op ==
  "scwbreferencemap::reference_map_calibrate" &
  artifacts$type == "obs_table" & artifacts$lifecycle_state == "current" &
  startsWith(artifacts$name, "gated_transfer_balanced_fold_calibration - ") &
  grepl("\\(thresholds\\)$", artifacts$name)
stopifnot(sum(threshold_take) == 1L)
thresholds <- wb_art_load(bundle, artifacts$artifact_id[threshold_take][[1L]],
  expected_kind = "obs_table")$table
supported_threshold <- as.numeric(thresholds$threshold[
  thresholds$role == "supported_upper"][[1L]])
out_threshold <- as.numeric(thresholds$threshold[
  thresholds$role == "borderline_upper"][[1L]])
stopifnot(isTRUE(all.equal(supported_threshold,
  as.numeric(cfg$mapping$supported_threshold), tolerance = 1e-12)),
  isTRUE(all.equal(out_threshold,
    as.numeric(cfg$mapping$out_of_reference_threshold), tolerance = 1e-12)))

log_event("Applying exact 30-neighbor support/native-L4 gate")
map_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_map"))
knn <- wb_session_run_op(map_session, bundle,
  "scwbrepresentation::knn_neighbors", inputs = list(
    space = cfg$mapping$joint_space_id, query_cells = slots,
    candidate_cells = reference_cells, k = as.integer(cfg$mapping$primary_k),
    include_query = FALSE,
    space_chunk_size = as.integer(cfg$mapping$space_chunk_size)),
  notes = "ImmGen bulk to selected 20D Bdev reference neighbors")
transfer <- wb_session_run_op(knn$session, bundle,
  "scwbreferencemap::reference_map_transfer", inputs = list(
    in_space = cfg$mapping$query_space_id,
    reference_space = cfg$inputs$qualified_bundle$reference_space,
    neighbors = knn$run$outputs$neighbors,
    calibration = calibration_ref,
    support_output_key = "immgen_bulk_gated20d_support",
    annotation_output_key = "immgen_bulk_gated20d_l4_domain",
    continuous_output_keys = cfg$mapping$continuous_output_key,
    continuous_acceptance = list(mode = "support_and_annotation",
      accepted_annotation_values = cfg$mapping$accepted_annotation_value,
      minimum_annotation_confidence =
        as.numeric(cfg$mapping$minimum_annotation_confidence)),
    supported_threshold = supported_threshold,
    out_of_reference_threshold = out_threshold,
    cell_ids = slots, min_cells = 1L,
    space_chunk_size = as.integer(cfg$mapping$space_chunk_size)),
  notes = "Qualified support-and-annotation native-L4 gate for ImmGen bulk")

neighbors <- as.data.table(knn$run$outputs$neighbors$table)
prediction <- as.data.table(transfer$run$outputs$per_query_summary$table)
categorical <- as.data.table(transfer$run$outputs$categorical_estimates$table)
continuous <- as.data.table(transfer$run$outputs$continuous_estimates$table)
acceptance <- as.data.table(transfer$run$outputs$continuous_acceptance$table)
transfer$session <- wb_session_discard(transfer$session)
neighbors <- neighbors[order(match(query_cell_id, slots), neighbor_rank)]
prediction <- prediction[match(slots, query_cell_id)]
acceptance <- acceptance[match(slots, query_cell_id)]
primary_categorical <- categorical[k == as.integer(cfg$mapping$primary_k)]
primary_categorical <- primary_categorical[match(slots, query_cell_id)]
primary_continuous <- continuous[k == as.integer(cfg$mapping$primary_k) &
  reference_key == "L4"]
primary_continuous <- primary_continuous[match(slots, query_cell_id)]
stopifnot(nrow(neighbors) == length(slots) * cfg$mapping$primary_k,
  identical(as.character(prediction$query_cell_id), slots),
  identical(as.character(acceptance$query_cell_id), slots),
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

fixed_bundle <- wb_open_bundle(cfg$inputs$fixed_umap_bundle$path,
  validate = FALSE, read_only = TRUE)
fixed <- as_space_matrix(wb_space(fixed_bundle,
  cfg$inputs$fixed_umap_bundle$space, cells = reference_cells))
wb_close_bundle(fixed_bundle)
stopifnot(identical(rownames(fixed), reference_cells), ncol(fixed) == 2L)
colnames(fixed) <- c("UMAP_1", "UMAP_2")
neighbors[, `:=`(UMAP_1 = fixed[cell_id, 1L], UMAP_2 = fixed[cell_id, 2L])]
placement <- neighbors[, {
  scale <- max(distance)
  weight <- if (scale == 0) as.numeric(distance == 0) else
    exp(-((distance / scale)^2))
  .(UMAP_1 = weighted.mean(UMAP_1, weight),
    UMAP_2 = weighted.mean(UMAP_2, weight),
    effective_neighbors = 1 / sum((weight / sum(weight))^2))
}, by = query_cell_id]
placement <- placement[match(slots, query_cell_id)]
placement <- cbind(result, placement[, .(UMAP_1, UMAP_2,
  effective_neighbors)])
placement[, detail := paste0(population, " · ", gsm)]

fwrite(result, file.path(dirs[["tables"]],
  "02_immgen_bulk_gated_20d_mapping.tsv"), sep = "\t", quote = TRUE,
  na = "NA")
fwrite(neighbors, file.path(dirs[["tables"]],
  "03_immgen_bulk_gated_20d_neighbors.tsv.gz"), sep = "\t", quote = TRUE,
  na = "NA")
fwrite(placement, file.path(dirs[["tables"]],
  "04_immgen_bulk_fixed_umap_placements.tsv"), sep = "\t", quote = TRUE,
  na = "NA")
fwrite(as.data.table(thresholds), file.path(dirs[["tables"]],
  "05_qualified_support_thresholds.tsv"), sep = "\t", quote = TRUE,
  na = "NA")

reference_plot <- data.frame(cell_id = reference_cells,
  UMAP_1 = fixed[, 1L], UMAP_2 = fixed[, 2L])
trajectory_model <- readRDS(cfg$inputs$authenticated_l4$model_path)
curve_object <- slingCurves(trajectory_model)[[
  cfg$inputs$authenticated_l4$curve_name]]
curve_xy <- curve_object$s[curve_object$ord, , drop = FALSE]
curve <- data.frame(UMAP_1 = curve_xy[, 1L], UMAP_2 = curve_xy[, 2L])
accepted <- placement[accepted == TRUE]
rejected <- placement[accepted == FALSE]
pad <- as.numeric(cfg$visualization$fixed_coordinate_padding_fraction)
x_range <- range(reference_plot$UMAP_1)
y_range <- range(reference_plot$UMAP_2)
x_limits <- x_range + c(-1, 1) * diff(x_range) * pad
y_limits <- y_range + c(-1, 1) * diff(y_range) * pad
caption <- sprintf(
  "Frozen 20D HVP20k projection; %d/%d libraries pass the complete native-L4 gate. Rejected libraries are dark grey.",
  nrow(accepted), nrow(placement))
p <- ggplot() +
  geom_point(data = reference_plot, aes(UMAP_1, UMAP_2),
    color = cfg$visualization$reference_color, size = 0.24, alpha = 0.42) +
  geom_path(data = curve, aes(UMAP_1, UMAP_2), color = "white",
    linewidth = 2.3) +
  geom_path(data = curve, aes(UMAP_1, UMAP_2), color = "#171717",
    linewidth = 0.9) +
  geom_point(data = rejected, aes(UMAP_1, UMAP_2, shape = lineage),
    fill = cfg$visualization$rejected_color, color = "#111111", size = 4.2,
    stroke = 0.65) +
  geom_point(data = accepted,
    aes(UMAP_1, UMAP_2, fill = accepted_native_l4, shape = lineage),
    color = "#111111", size = 4.4, stroke = 0.7) +
  geom_label_repel(data = placement,
    aes(UMAP_1, UMAP_2, label = detail), seed = cfg$visualization$label_seed,
    max.overlaps = Inf, max.time = 30, size = 2.55, label.size = 0.15,
    fill = adjustcolor("white", 0.9), color = "#222222") +
  scale_fill_gradientn(colors = unlist(cfg$visualization$l4_palette),
    name = "Accepted native L4") +
  scale_shape_manual(values = c(B = 21, "Stem&Prog" = 24), name = "Lineage") +
  coord_equal(xlim = x_limits, ylim = y_limits, clip = "off") +
  labs(title = "ImmGen bulk ATAC on the selected Bdev reference",
    subtitle = "Exact frozen 20D model and qualified support/native-L4 gate",
    x = "ATAC UMAP 1", y = "ATAC UMAP 2", caption = caption) +
  theme_classic(base_size = 11.5) +
  theme(plot.title = element_text(face = "bold", size = 16),
    legend.position = "bottom", plot.margin = margin(10, 55, 10, 10))
ggsave(file.path(dirs[["figures"]],
  "01_immgen_bulk_gated_20d_fixed_atac_umap.png"), p,
  width = 15, height = 10.5, dpi = 260, bg = "white")
ggsave(file.path(dirs[["figures"]],
  "01_immgen_bulk_gated_20d_fixed_atac_umap.pdf"), p,
  width = 15, height = 10.5, device = cairo_pdf, bg = "white")

support_counts <- as.list(table(factor(result$support_call,
  levels = c("supported", "borderline", "out_of_reference"))))
reason_counts <- as.list(table(factor(result$acceptance_reason,
  levels = c("accepted", "out_of_reference", "annotation_not_accepted",
    "annotation_confidence_too_low", "continuous_value_unavailable"))))
stopifnot(sha256(model_payload) ==
    cfg$inputs$qualified_bundle$model_payload_sha256,
  sha256(qualified_manifest) ==
    cfg$inputs$qualified_bundle$manifest_sha256,
  sha256(bulk_manifest) == cfg$inputs$bulk_query_bundle$manifest_sha256,
  sha256(fixed_umap_manifest) ==
    cfg$inputs$fixed_umap_bundle$manifest_sha256)
wb_finalize_recipient_bundle(bundle)
wb_close_bundle(bundle)
bundle <- NULL

report <- list(
  schema_version = "bdev.immgen_bulk_gated_20d_report.v1",
  status = "passed", gate_passed = TRUE, run_identity = cfg$run_identity,
  reference_cells = length(reference_cells), query_libraries = nrow(result),
  selected_features = 20000L, latent_dimensions = 20L,
  model_artifact_id = cfg$inputs$qualified_bundle$model_artifact_id,
  model_payload_sha256 = sha256(model_payload),
  projection_artifact_id = projection_artifact_id,
  union_artifact_id = union_artifact_id,
  calibration_artifact_id = cfg$inputs$qualified_bundle$calibration_artifact_id,
  primary_k = as.integer(cfg$mapping$primary_k),
  supported_threshold = supported_threshold,
  out_of_reference_threshold = out_threshold,
  minimum_annotation_confidence =
    as.numeric(cfg$mapping$minimum_annotation_confidence),
  support_counts = support_counts,
  acceptance_reason_counts = reason_counts,
  accepted_query_libraries = sum(result$accepted),
  withheld_query_libraries = sum(!result$accepted),
  raw_estimate_rows = nrow(result),
  fixed_umap_placements = nrow(placement),
  qualified_source_manifest_sha256 = sha256(qualified_manifest),
  bulk_source_manifest_sha256 = sha256(bulk_manifest),
  fixed_umap_source_manifest_sha256 = sha256(fixed_umap_manifest),
  mapping_bundle_manifest_sha256 = sha256(file.path(
    cfg$storage$mapping_bundle, "manifest.json")),
  model_refit = FALSE, umap_refit = FALSE, biological_inference = FALSE,
  paper_analysis = FALSE, sif_built = FALSE)
write_json(report, file.path(dirs[["reports"]],
  "01_immgen_bulk_gated_20d_mapping.json"))
write_json(list(
  schema_version = "bdev.immgen_bulk_gated_20d_execution.v1",
  config_path = config_path, config_sha256 = sha256(config_path),
  scworkbench = as.character(packageVersion("scworkbench")),
  scwbLatent = as.character(packageVersion("scwbLatent")),
  scwbRepresentation = as.character(packageVersion("scwbRepresentation")),
  scwbReferenceMap = as.character(packageVersion("scwbReferenceMap")),
  qualified_source_manifest_unchanged =
    sha256(qualified_manifest) == cfg$inputs$qualified_bundle$manifest_sha256,
  bulk_source_manifest_unchanged =
    sha256(bulk_manifest) == cfg$inputs$bulk_query_bundle$manifest_sha256,
  fixed_umap_source_manifest_unchanged =
    sha256(fixed_umap_manifest) == cfg$inputs$fixed_umap_bundle$manifest_sha256,
  model_payload_unchanged = sha256(model_payload) ==
    cfg$inputs$qualified_bundle$model_payload_sha256,
  source_archives_modified = FALSE), file.path(dirs[["manifests"]],
  "02_execution_receipt.json"))
log_event("ImmGen bulk gated 20D mapping and fixed-UMAP figure passed")
