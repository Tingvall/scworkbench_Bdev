#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 05_map_immgen_trajectory.R --config <config.yaml>")

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

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
read_config <- function(path) {
  current <- yaml::read_yaml(path)
  base_name <- current$base_config
  current$base_config <- NULL
  if (is.null(base_name)) return(current)
  utils::modifyList(read_config(normalizePath(file.path(dirname(path), base_name))),
    current, keep.null = TRUE)
}
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
promoted_ref <- function(promoted, output, kind) {
  selected <- promoted$promoted$status == "ok" & promoted$promoted$output == output
  stopifnot(sum(selected) == 1L)
  wb_ref_artifact(promoted$promoted$artifact_id[selected][[1L]], kind)
}
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- read_config(cfg_path)
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "05_immgen_trajectory_mapping.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256),
  identical(cfg$inputs$authenticated_l4$orientation, "native"),
  isFALSE(cfg$scope$categorical_cluster_mapping))
for (package in c("scwbRepresentation", "scwbReferenceMap", "scwbLatent")) {
  stopifnot(as.character(packageVersion(package)) == cfg$runtime[[package]]$version)
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  stopifnot(identical(sha256(manifest), cfg$runtime[[package]]$manifest_sha256))
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
}
required_ops <- c("scwblatent::peakvi_fit", "scwblatent::peakvi_project",
  "scwblatent::model_bound_latent_union", "scwbrepresentation::knn_neighbors",
  "scwbreferencemap::reference_map_calibrate",
  "scwbreferencemap::reference_map_transfer")
stopifnot(all(required_ops %in% wb_list_ops()))

library_manifest <- read.delim(file.path(dirs[["tables"]],
  "01_immgen_library_run_manifest.tsv"), check.names = FALSE, quote = "")
qc <- read.delim(file.path(dirs[["tables"]], "03_immgen_preprocessing_qc.tsv"),
  check.names = FALSE, quote = "")
stopifnot(nrow(library_manifest) == cfg$inputs$selected_libraries$expected_libraries,
  nrow(qc) == nrow(library_manifest),
  setequal(library_manifest$gsm, qc$gsm), all(qc$status == "complete"),
  all(qc$fastq_md5_valid))

fixture <- read.delim(cfg$inputs$fixture$path, check.names = FALSE,
  quote = "\"", na.strings = "NA")
family <- readRDS(cfg$inputs$authenticated_l4$path)
selected <- as.character(family$key) == cfg$inputs$authenticated_l4$key &
  as.character(family$cell_id) %in% fixture$cell_id &
  is.finite(as.numeric(family$value))
l4 <- stats::setNames(as.numeric(family$value[selected]),
  as.character(family$cell_id[selected]))
stopifnot(length(l4) == nrow(fixture), !anyDuplicated(names(l4)))
fixture$native_l4 <- unname(l4[fixture$cell_id])
fixture$group_id <- sprintf("sample_%s_bin_%02d", fixture$sample_id,
  as.integer(fixture$bin))
group_sizes <- table(fixture$group_id)
eligible <- names(group_sizes)[group_sizes >= cfg$split$minimum_group_cells]
fixture$partition <- "reference"
for (group in eligible) {
  index <- which(fixture$group_id == group)
  index <- index[order(fixture$bundle_order[index], method = "radix")]
  fixture$partition[index[seq_along(index) %% 4L == 0L]] <- "held_out"
}
reference_cells <- as.character(fixture$cell_id[fixture$partition == "reference"])
stopifnot(length(reference_cells) == 1382L)

source_manifest <- file.path(cfg$inputs$source_bundle$path, "manifest.json")
source_manifest_before <- sha256(source_manifest)
stopifnot(identical(source_manifest_before,
  cfg$inputs$source_bundle$manifest_sha256))
source <- wb_open_bundle(cfg$inputs$source_bundle$path, validate = FALSE,
  read_only = TRUE)
feature_ids <- wb_feature_ids(source, cfg$inputs$count_matrix)
reference_counts <- wb_slice(source, cfg$inputs$count_matrix,
  features = feature_ids, cells = reference_cells)
reference_meta <- wb_meta(source, cells = reference_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$sample_field))
wb_close_bundle(source)
stopifnot(length(feature_ids) == cfg$inputs$expected_features,
  identical(rownames(reference_counts), feature_ids),
  identical(colnames(reference_counts), reference_cells))

log_event("Reading", nrow(library_manifest), "ImmGen peak-count profiles")
query_cells <- as.character(library_manifest$gsm)
query_i <- list(); query_j <- list(); query_x <- list()
for (j in seq_along(query_cells)) {
  path <- file.path(cfg$storage$peak_counts,
    paste0(query_cells[[j]], ".bdev_peak_counts.tsv.gz"))
  values <- read.delim(gzfile(path), check.names = FALSE, quote = "")
  stopifnot(identical(as.character(values$feature_id), feature_ids),
    all(values$count >= 0), all(values$count == floor(values$count)))
  nonzero <- which(values$count != 0)
  query_i[[j]] <- nonzero
  query_j[[j]] <- rep.int(j, length(nonzero))
  query_x[[j]] <- values$count[nonzero]
}
query_counts <- Matrix::sparseMatrix(i = unlist(query_i, use.names = FALSE),
  j = unlist(query_j, use.names = FALSE), x = unlist(query_x, use.names = FALSE),
  dims = c(length(feature_ids), length(query_cells)),
  dimnames = list(feature_ids, query_cells))
stopifnot(all(Matrix::colSums(query_counts) > 0))

log_event("Creating derived bundle with", length(reference_cells),
  "Bdev references and", length(query_cells), "ImmGen libraries")
bundle <- wb_init_bundle(cfg$storage$mapping_bundle, overwrite = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
all_cells <- c(reference_cells, query_cells)
wb_define_cells(bundle, all_cells)
wb_attach_matrix(bundle, "Bdev_ATAC_reference_counts", reference_counts,
  layer = "counts", feature_type = "peak", cell_ids = reference_cells,
  description = "Qualified Step 4 Bdev PeakVI reference counts",
  source = cfg$inputs$source_bundle$path,
  created_by = "scworkbench_Bdev Step 4 ImmGen trajectory mapping")
wb_attach_matrix(bundle, "ImmGen_ATAC_Bdev_peak_counts", query_counts,
  layer = "counts", feature_type = "peak", cell_ids = query_cells,
  description = "ImmGen bulk ATAC libraries reprocessed on exact Bdev peak axis",
  source = "GSE100738/SRP110978 raw paired FASTQ",
  created_by = "scworkbench_Bdev Step 4 ImmGen trajectory mapping")
query_meta <- data.frame(cell_id = query_cells, stringsAsFactors = FALSE)
query_meta[[cfg$inputs$sample_field]] <- query_cells
metadata <- rbind(reference_meta, query_meta)
metadata <- metadata[match(all_cells, metadata$cell_id), , drop = FALSE]
wb_attach_metadata(bundle, metadata, cell_ids = all_cells,
  source = "active Bdev metadata plus frozen ImmGen library manifest",
  created_by = "scworkbench_Bdev Step 4 ImmGen trajectory mapping")
rm(reference_counts, query_counts, query_i, query_j, query_x); invisible(gc())

log_event("Fitting deterministic qualified 10D Bdev PeakVI reference")
fit_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_peakvi_fit"))
fit <- wb_session_run_op(fit_session, bundle, "scwblatent::peakvi_fit",
  inputs = list(matrix = "Bdev_ATAC_reference_counts",
    n_latent = as.integer(cfg$peakvi$n_latent),
    n_layers = as.integer(cfg$peakvi$n_layers),
    n_hidden = as.integer(cfg$peakvi$n_hidden),
    dropout = as.numeric(cfg$peakvi$dropout),
    max_epochs = as.integer(cfg$peakvi$max_epochs),
    lr = as.numeric(cfg$peakvi$learning_rate),
    early_stopping = isTRUE(cfg$peakvi$early_stopping),
    seed = as.integer(cfg$peakvi$seed), chunk_size = as.integer(cfg$peakvi$chunk_size),
    batch_key = NULL, use_gpu = FALSE, cell_set = "all", min_cells = 1L),
  notes = "Step 4 deterministic Bdev-only 10D PeakVI reference fit")
reference_space_id <- as.character(fit$run$outputs$space$space_id %||%
  fit$run$outputs$space$allocated_id)[[1L]]
fit_summary <- fit$run$outputs$fit_summary$table
stopifnot(nrow(fit_summary) == 1L,
  fit_summary$latent_sha256[[1L]] == cfg$peakvi$expected_reference_latent_sha256)
fit$session <- wb_session_resume(bundle, fit$session$durable$workspace)
fit_promoted <- wb_session_promote_current(bundle, fit$session,
  outputs = c("model", "fit_summary"), name_prefix = "step4_immgen_peakvi_reference",
  description = "Qualified deterministic Bdev-only 10D PeakVI reference",
  created_by = "scworkbench_Bdev Step 4 ImmGen trajectory mapping", fail_fast = TRUE)
model_ref <- promoted_ref(fit_promoted, "model", "model")
fit_promoted$session <- wb_session_discard(fit_promoted$session)

log_event("Projecting ImmGen libraries without refitting")
project_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_peakvi_project"))
project <- wb_session_run_op(project_session, bundle, "scwblatent::peakvi_project",
  inputs = list(model = model_ref, query_matrix = "ImmGen_ATAC_Bdev_peak_counts",
    space_id_out = "step4_immgen_peakvi_queries",
    chunk_size = as.integer(cfg$peakvi$chunk_size), batch_key = NULL,
    use_gpu = FALSE, pickle_execution_trust = cfg$peakvi$pickle_execution_trust,
    cell_set = "all", min_cells = 1L),
  notes = "Frozen-model projection of 28 ImmGen ATAC libraries")
projection_summary <- project$run$outputs$projection_summary$table
stopifnot(projection_summary$transformation[[1L]] == "transfer_without_refit",
  projection_summary$n_latent[[1L]] == cfg$peakvi$n_latent)
project$session <- wb_session_resume(bundle, project$session$durable$workspace)
project_promoted <- wb_session_promote_current(bundle, project$session,
  outputs = "projection_summary", name_prefix = "step4_immgen_peakvi_projection",
  description = "Frozen PeakVI projection of ImmGen ATAC libraries",
  created_by = "scworkbench_Bdev Step 4 ImmGen trajectory mapping", fail_fast = TRUE)
projection_ref <- promoted_ref(project_promoted, "projection_summary", "obs_table")
project_promoted$session <- wb_session_discard(project_promoted$session)

log_event("Authenticating model-bound reference/query latent union")
union_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_latent_union"))
union <- wb_session_run_op(union_session, bundle,
  "scwblatent::model_bound_latent_union", inputs = list(
    reference_space = reference_space_id,
    query_space = "step4_immgen_peakvi_queries", model = model_ref,
    projection_summary = projection_ref,
    space_id_out = "step4_immgen_reference_and_queries",
    chunk_size = as.integer(cfg$mapping$query_chunk_size)),
  notes = "Authenticated model-bound Bdev and ImmGen latent union")
union$session <- wb_session_resume(bundle, union$session$durable$workspace)
union_promoted <- wb_session_promote_current(bundle, union$session,
  outputs = "union_summary", name_prefix = "step4_immgen_peakvi_union",
  description = "Authenticated Bdev/ImmGen model-bound latent union",
  created_by = "scworkbench_Bdev Step 4 ImmGen trajectory mapping", fail_fast = TRUE)
union_promoted$session <- wb_session_discard(union_promoted$session)

reference_l4 <- wb_prim_cell_vector_family(
  rep(cfg$inputs$authenticated_l4$key, length(reference_cells)), reference_cells,
  unname(l4[reference_cells]), value_type = "continuous")
log_event("Calibrating support and native L4 only; no cluster labels")
calibration_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_calibrate"))
calibration <- wb_session_run_op(calibration_session, bundle,
  "scwbreferencemap::reference_map_calibrate", inputs = list(
    reference_matrix = "Bdev_ATAC_reference_counts",
    reference_space = reference_space_id,
    reference_continuous = reference_l4,
    holdout_group_field = cfg$inputs$sample_field,
    k_values = as.integer(unlist(cfg$mapping$k_values)),
    primary_k = as.integer(cfg$mapping$primary_k),
    support_quantiles = as.list(as.numeric(unlist(cfg$mapping$support_quantiles))),
    primary_continuous_key = cfg$inputs$authenticated_l4$key,
    space_chunk_size = as.integer(cfg$mapping$space_chunk_size),
    query_chunk_size = as.integer(cfg$mapping$query_chunk_size)),
  notes = "Trajectory-only Bdev PeakVI support and native L4 calibration")
calibration$session <- wb_session_resume(bundle,
  calibration$session$durable$workspace)
calibration_promoted <- wb_session_promote_current(bundle, calibration$session,
  outputs = names(calibration$run$outputs), name_prefix = "step4_immgen_calibration",
  description = "Trajectory-only Bdev PeakVI support and L4 calibration",
  created_by = "scworkbench_Bdev Step 4 ImmGen trajectory mapping", fail_fast = TRUE)
calibration_ref <- promoted_ref(calibration_promoted, "calibration", "model")
thresholds <- calibration$run$outputs$thresholds$table
supported_threshold <- thresholds$threshold[thresholds$role == "supported_upper"][[1L]]
out_threshold <- thresholds$threshold[thresholds$role == "borderline_upper"][[1L]]
calibration_promoted$session <- wb_session_discard(calibration_promoted$session)

knn_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_knn"))
knn <- wb_session_run_op(knn_session, bundle,
  "scwbrepresentation::knn_neighbors", inputs = list(
    space = "step4_immgen_reference_and_queries", query_cells = query_cells,
    candidate_cells = reference_cells, k = as.integer(cfg$mapping$neighbor_k),
    include_query = FALSE, backend = cfg$mapping$backend,
    space_chunk_size = as.integer(cfg$mapping$space_chunk_size),
    query_chunk_size = as.integer(cfg$mapping$query_chunk_size)),
  notes = "ImmGen-to-Bdev neighbors in authenticated 10D PeakVI space")
transfer_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_transfer"))
transfer <- wb_session_run_op(transfer_session, bundle,
  "scwbreferencemap::reference_map_transfer", inputs = list(
    in_space = "step4_immgen_peakvi_queries",
    reference_space = reference_space_id, neighbors = knn$run$outputs$neighbors,
    calibration = calibration_ref, support_output_key = "immgen_bdev_support",
    annotation_output_key = "unused_no_cluster_transfer",
    continuous_output_keys = "mapped_native_L4",
    supported_threshold = as.numeric(supported_threshold),
    out_of_reference_threshold = as.numeric(out_threshold),
    cell_ids = query_cells, min_cells = 1L,
    space_chunk_size = as.integer(cfg$mapping$space_chunk_size)),
  notes = "Trajectory-only native L4 transfer; categorical mapping omitted")

prediction <- transfer$run$outputs$per_query_summary$table
continuous <- transfer$run$outputs$continuous_estimates$table
primary <- continuous$k == as.integer(cfg$mapping$primary_k) &
  continuous$reference_key == cfg$inputs$authenticated_l4$key
primary_continuous <- continuous[primary, c("query_cell_id", "weighted_q25",
  "raw_estimate", "weighted_q75", "iqr", "cross_k_sensitivity")]
names(primary_continuous)[names(primary_continuous) == "raw_estimate"] <-
  "raw_native_l4"
library_results <- merge(library_manifest, prediction,
  by.x = "gsm", by.y = "query_cell_id", sort = FALSE)
library_results <- merge(library_results, primary_continuous,
  by.x = "gsm", by.y = "query_cell_id", sort = FALSE)
library_results <- library_results[match(query_cells, library_results$gsm), ]
library_results$accepted_native_l4 <- ifelse(
  library_results$support_call == "out_of_reference", NA_real_,
  library_results$raw_native_l4)
stopifnot(all(is.na(library_results$accepted_annotation)))

population_names <- unique(library_manifest$population)
population_rows <- lapply(population_names, function(population) {
  x <- library_results[library_results$population == population, , drop = FALSE]
  accepted <- is.finite(x$accepted_native_l4)
  data.frame(population = population, n_libraries = nrow(x),
    n_supported = sum(x$support_call == "supported"),
    n_borderline = sum(x$support_call == "borderline"),
    n_out_of_reference = sum(x$support_call == "out_of_reference"),
    accepted_fraction = mean(x$support_call != "out_of_reference"),
    median_native_l4 = if (any(accepted)) median(x$accepted_native_l4[accepted]) else NA_real_,
    min_library_q25 = if (any(accepted)) min(x$weighted_q25[accepted]) else NA_real_,
    max_library_q75 = if (any(accepted)) max(x$weighted_q75[accepted]) else NA_real_,
    median_cross_k_sensitivity = if (any(accepted))
      median(x$cross_k_sensitivity[accepted]) else NA_real_,
    stringsAsFactors = FALSE)
})
population_results <- do.call(rbind, population_rows)
write_tsv(library_results, file.path(dirs[["tables"]],
  "04_immgen_library_l4_mapping.tsv"))
write_tsv(population_results, file.path(dirs[["tables"]],
  "05_immgen_population_l4_summary.tsv"))
write_tsv(continuous, file.path(dirs[["tables"]],
  "06_immgen_l4_cross_k_estimates.tsv"))
write_tsv(calibration$run$outputs$calibration_metrics$table,
  file.path(dirs[["tables"]], "07_trajectory_calibration_metrics.tsv"))
write_tsv(thresholds, file.path(dirs[["tables"]],
  "08_support_thresholds.tsv"))
write_tsv(transfer$run$outputs$grouped_summary$table,
  file.path(dirs[["tables"]], "09_support_grouped_summary.tsv"))

plot_l4 <- function(device, path, width, height, res = NULL) {
  if (is.null(res)) device(path, width = width, height = height) else
    device(path, width = width, height = height, units = "in", res = res)
  on.exit(dev.off(), add = TRUE)
  x <- library_results[order(library_results$population,
    library_results$accepted_native_l4, na.last = TRUE), ]
  colors <- c(supported = "#2166AC", borderline = "#FDAE61",
    out_of_reference = "#B2182B")
  y <- rev(seq_len(nrow(x)))
  par(mar = c(4.5, 10, 2, 1))
  plot(range(c(x$weighted_q25, x$weighted_q75), na.rm = TRUE), range(y),
    type = "n", yaxt = "n", xlab = "Transferred native L4 pseudotime",
    ylab = "", main = "ImmGen ATAC trajectory mapping")
  segments(x$weighted_q25, y, x$weighted_q75, y,
    col = colors[x$support_call], lwd = 2)
  points(x$raw_native_l4, y, pch = 19, col = colors[x$support_call])
  axis(2, at = y, labels = paste(x$population, x$gsm, sep = "  "),
    las = 2, cex.axis = 0.65)
  legend("bottomright", legend = names(colors), col = colors, pch = 19,
    bty = "n", cex = 0.8)
}
plot_l4(grDevices::pdf, file.path(dirs[["figures"]],
  "01_immgen_library_native_l4.pdf"), 9, 9)
plot_l4(grDevices::png, file.path(dirs[["figures"]],
  "01_immgen_library_native_l4.png"), 9, 9, 180)

plot_support <- function(device, path, width, height, res = NULL) {
  if (is.null(res)) device(path, width = width, height = height) else
    device(path, width = width, height = height, units = "in", res = res)
  on.exit(dev.off(), add = TRUE)
  colors <- c(supported = "#2166AC", borderline = "#FDAE61",
    out_of_reference = "#B2182B")
  y <- seq_len(nrow(library_results))
  par(mar = c(4.5, 10, 2, 1))
  plot(library_results$distance_statistic, y, pch = 19,
    col = colors[library_results$support_call], yaxt = "n",
    xlab = "Reference-support distance statistic", ylab = "",
    main = "ImmGen support against the Bdev reference")
  abline(v = supported_threshold, lty = 2, col = colors[["supported"]])
  abline(v = out_threshold, lty = 2, col = colors[["out_of_reference"]])
  axis(2, at = y, labels = paste(library_results$population,
    library_results$gsm, sep = "  "), las = 2, cex.axis = 0.65)
}
plot_support(grDevices::pdf, file.path(dirs[["figures"]],
  "02_immgen_reference_support.pdf"), 9, 9)
plot_support(grDevices::png, file.path(dirs[["figures"]],
  "02_immgen_reference_support.png"), 9, 9, 180)

knn$session <- wb_session_resume(bundle, knn$session$durable$workspace)
transfer$session <- wb_session_resume(bundle, transfer$session$durable$workspace)
knn$session <- wb_session_discard(knn$session)
transfer$session <- wb_session_discard(transfer$session)
source_manifest_after <- sha256(source_manifest)
stopifnot(identical(source_manifest_before, source_manifest_after))
finalized <- wb_finalize_recipient_bundle(bundle)
stopifnot(isTRUE(finalized$ready))
wb_close_bundle(bundle); bundle <- NULL

support_counts <- as.list(table(factor(library_results$support_call,
  levels = c("supported", "borderline", "out_of_reference"))))
summary <- list(schema_version = "bdev.immgen_atac_trajectory_mapping_summary.v1",
  status = "complete", run_identity = cfg$run_identity,
  query_libraries = nrow(library_results),
  query_populations = nrow(population_results),
  bdev_reference_cells = length(reference_cells), feature_count = length(feature_ids),
  latent_dimensions = cfg$peakvi$n_latent,
  reference_latent_sha256 = fit_summary$latent_sha256[[1L]],
  projection_transformation = projection_summary$transformation[[1L]],
  support_counts = support_counts,
  accepted_libraries = sum(library_results$support_call != "out_of_reference"),
  categorical_cluster_transfer = FALSE,
  trajectory = list(key = cfg$inputs$authenticated_l4$key,
    orientation = cfg$inputs$authenticated_l4$orientation,
    primary_k = cfg$mapping$primary_k),
  interpretation = "L4 estimates are reportable only for supported or borderline libraries; out-of-reference libraries retain no accepted L4 value.",
  source_bundle_manifest_unchanged = TRUE,
  archive_modified = FALSE, paper_analysis = FALSE,
  implementation_unlocked = FALSE, sif_built = FALSE)
write_json(summary, file.path(dirs[["reports"]],
  "01_immgen_atac_trajectory_mapping_summary.json"))
write_json(list(schema_version = "bdev.immgen_atac_trajectory_execution.v1",
  status = "complete", run_identity = cfg$run_identity,
  config_path = cfg_path, config_sha256 = sha256(cfg_path),
  source_manifest_sha256_before = source_manifest_before,
  source_manifest_sha256_after = source_manifest_after,
  mapping_bundle = cfg$storage$mapping_bundle,
  mapping_bundle_manifest_sha256 = sha256(file.path(cfg$storage$mapping_bundle,
    "manifest.json")), categorical_cluster_transfer = FALSE,
  real_immgen_projected = TRUE),
  file.path(dirs[["manifests"]], "05_execution_receipt.json"))
log_event("ImmGen trajectory-only mapping complete; support:",
  paste(names(support_counts), unlist(support_counts), collapse = ", "))
