#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 01_prepare_attribution.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(yaml); library(jsonlite); library(digest); library(Matrix)
})
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
read_config <- function(path) {
  current <- yaml::read_yaml(path); base <- current$base_config; current$base_config <- NULL
  if (is.null(base)) return(current)
  utils::modifyList(read_config(normalizePath(file.path(dirname(path), base))), current, keep.null = TRUE)
}
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
write_tsv_gz <- function(x, path) { con <- gzfile(path, "wt"); on.exit(close(con)); write.table(x, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA") }
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
promoted_ref <- function(promoted, output, kind) {
  take <- promoted$promoted$status == "ok" & promoted$promoted$output == output
  stopifnot(sum(take) == 1L); wb_ref_artifact(promoted$promoted$artifact_id[take][[1L]], kind)
}

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE); cfg <- read_config(cfg_path)
.libPaths(c(cfg$runtime$run_library, .libPaths()))
suppressPackageStartupMessages({ library(scworkbench); library(scwbLatent) })
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports", "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_prepare_attribution.log"), split = TRUE); on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256),
  identical(sha256(file.path(cfg$inputs$source_mapping_bundle$path, "manifest.json")), cfg$inputs$source_mapping_bundle$manifest_sha256),
  identical(sha256(file.path(cfg$inputs$source_bundle$path, "manifest.json")), cfg$inputs$source_bundle$manifest_sha256),
  identical(sha256(cfg$inputs$provisional_rna_labels$path), cfg$inputs$provisional_rna_labels$sha256),
  identical(sha256(cfg$inputs$label_bootstrap_stability$path), cfg$inputs$label_bootstrap_stability$sha256),
  identical(sha256(cfg$inputs$neighbor_driver_ledger$path), cfg$inputs$neighbor_driver_ledger$sha256))
for (package in "scwbLatent") {
  stopifnot(as.character(packageVersion(package)) == cfg$runtime[[package]]$version)
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  stopifnot(identical(sha256(manifest), cfg$runtime[[package]]$manifest_sha256))
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
}

source_map_before <- sha256(file.path(cfg$inputs$source_mapping_bundle$path, "manifest.json"))
source_bdev_before <- sha256(file.path(cfg$inputs$source_bundle$path, "manifest.json"))
labels <- read.delim(cfg$inputs$provisional_rna_labels$path, check.names = FALSE, quote = "\"", na.strings = "NA")
stability <- read.delim(gzfile(cfg$inputs$label_bootstrap_stability$path), check.names = FALSE, quote = "\"", na.strings = "NA")
labels <- merge(labels, stability[, c("cell_id", "same_as_baseline_fraction", "assigned_fraction", "stable_at_threshold")], by = "cell_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(labels) == cfg$acceptance$expected_label_fixture_cells, !anyDuplicated(labels$cell_id))

mapping <- wb_open_bundle(cfg$inputs$source_mapping_bundle$path, validate = FALSE, read_only = TRUE)
reference_space <- wb_space(mapping, cfg$inputs$source_mapping_bundle$reference_space)
reference_cells <- rownames(reference_space)
feature_ids <- wb_feature_ids(mapping, cfg$inputs$source_mapping_bundle$reference_matrix)
reference_counts <- wb_slice(mapping, cfg$inputs$source_mapping_bundle$reference_matrix, features = feature_ids, cells = reference_cells)
library_manifest <- read.delim(cfg$inputs$selected_libraries$path, check.names = FALSE, quote = "\"", na.strings = "NA")
bulk_cells <- as.character(library_manifest$gsm)
bulk_counts <- wb_slice(mapping, cfg$inputs$source_mapping_bundle$bulk_matrix, features = feature_ids, cells = bulk_cells)
wb_close_bundle(mapping)
stopifnot(length(reference_cells) == cfg$acceptance$expected_reference_cells, length(bulk_cells) == cfg$acceptance$expected_immgen_libraries)

bdev <- wb_open_bundle(cfg$inputs$source_bundle$path, validate = FALSE, read_only = TRUE)
stopifnot(identical(wb_feature_ids(bdev, cfg$inputs$count_matrix), feature_ids))
label_counts <- wb_slice(bdev, cfg$inputs$count_matrix, features = feature_ids, cells = labels$cell_id)
tech_fields <- unlist(cfg$diagnostic$technical_fields, use.names = FALSE)
all_bdev_cells <- unique(c(reference_cells, labels$cell_id))
bdev_meta <- wb_meta(bdev, cells = all_bdev_cells, cols = c(list(wb_ref_metadata_field("core:cell_id")), as.list(tech_fields)))
wb_close_bundle(bdev)

new_label_cells <- setdiff(labels$cell_id, reference_cells)
label_query_counts <- label_counts[, new_label_cells, drop = FALSE]
query_counts <- cbind(bulk_counts, label_query_counts)
query_cells <- colnames(query_counts)
all_cells <- c(reference_cells, query_cells)
if (dir.exists(cfg$storage$diagnostic_bundle)) stop("Refusing to overwrite diagnostic bundle")
bundle <- wb_init_bundle(cfg$storage$diagnostic_bundle, overwrite = FALSE); on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
wb_define_cells(bundle, all_cells)
wb_attach_matrix(bundle, "Bdev_ATAC_reference_counts", reference_counts, layer = "counts", feature_type = "peak", cell_ids = reference_cells,
  description = "Qualified Bdev-only PeakVI reference counts", source = cfg$inputs$source_mapping_bundle$path, created_by = "Bdev position-driver diagnostic")
wb_attach_matrix(bundle, "Position_driver_queries", query_counts, layer = "counts", feature_type = "peak", cell_ids = query_cells,
  description = "ImmGen bulk libraries plus non-reference provisional RNA-labelled Bdev diagnostic cells", source = "pinned Bdev and ImmGen inputs", created_by = "Bdev position-driver diagnostic")
meta <- bdev_meta[match(all_cells, bdev_meta$cell_id), , drop = FALSE]
meta$cell_id[seq_along(all_cells)] <- all_cells
meta$orig.ident[match(bulk_cells, all_cells)] <- bulk_cells
wb_attach_metadata(bundle, meta, cell_ids = all_cells, source = "pinned Bdev metadata; ImmGen library IDs", created_by = "Bdev position-driver diagnostic")

fit_session <- wb_session_create(bundle, cfg$storage$sessions, paste0(cfg$run_identity, "_fit"))
fit <- wb_session_run_op(fit_session, bundle, "scwblatent::peakvi_fit", inputs = list(
  matrix = "Bdev_ATAC_reference_counts", n_latent = as.integer(cfg$peakvi$n_latent), n_layers = as.integer(cfg$peakvi$n_layers),
  n_hidden = as.integer(cfg$peakvi$n_hidden), dropout = as.numeric(cfg$peakvi$dropout), max_epochs = as.integer(cfg$peakvi$max_epochs),
  lr = as.numeric(cfg$peakvi$learning_rate), early_stopping = isTRUE(cfg$peakvi$early_stopping), seed = as.integer(cfg$peakvi$seed),
  chunk_size = as.integer(cfg$peakvi$chunk_size), batch_key = NULL, use_gpu = FALSE, cell_set = "all", min_cells = 1L),
  notes = "Deterministic Bdev-only PeakVI refit for isolated position-driver diagnostic")
reference_space_id <- as.character(fit$run$outputs$space$space_id %||% fit$run$outputs$space$allocated_id)[[1L]]
stopifnot(fit$run$outputs$fit_summary$table$latent_sha256[[1L]] == cfg$peakvi$expected_reference_latent_sha256)
fit$session <- wb_session_resume(bundle, fit$session$durable$workspace)
prom <- wb_session_promote_current(bundle, fit$session, outputs = c("model", "fit_summary"), name_prefix = "position_driver_reference",
  description = "Qualified deterministic Bdev-only PeakVI reference", created_by = "Bdev position-driver diagnostic", fail_fast = TRUE)
model_ref <- promoted_ref(prom, "model", "model"); prom$session <- wb_session_discard(prom$session)

project_session <- wb_session_create(bundle, cfg$storage$sessions, paste0(cfg$run_identity, "_project"))
project <- wb_session_run_op(project_session, bundle, "scwblatent::peakvi_project", inputs = list(model = model_ref,
  query_matrix = "Position_driver_queries", space_id_out = "position_driver_queries", chunk_size = as.integer(cfg$peakvi$chunk_size),
  batch_key = NULL, use_gpu = FALSE, pickle_execution_trust = cfg$peakvi$pickle_execution_trust, cell_set = "all", min_cells = 1L),
  notes = "Projection without refit of ImmGen and labelled Bdev diagnostic cells")
stopifnot(project$run$outputs$projection_summary$table$transformation[[1L]] == "transfer_without_refit")
project$session <- wb_session_resume(bundle, project$session$durable$workspace)
pp <- wb_session_promote_current(bundle, project$session, outputs = "projection_summary", name_prefix = "position_driver_projection",
  description = "Frozen-model position-driver projection", created_by = "Bdev position-driver diagnostic", fail_fast = TRUE)
projection_ref <- promoted_ref(pp, "projection_summary", "obs_table"); pp$session <- wb_session_discard(pp$session)

union_session <- wb_session_create(bundle, cfg$storage$sessions, paste0(cfg$run_identity, "_union"))
union <- wb_session_run_op(union_session, bundle, "scwblatent::model_bound_latent_union", inputs = list(reference_space = reference_space_id,
  query_space = "position_driver_queries", model = model_ref, projection_summary = projection_ref, space_id_out = "position_driver_union",
  chunk_size = as.integer(cfg$mapping$query_chunk_size), cell_set = "all", min_cells = 10L), notes = "Authenticated diagnostic latent union")
union$session <- wb_session_resume(bundle, union$session$durable$workspace)
up <- wb_session_promote_current(bundle, union$session, outputs = "union_summary", name_prefix = "position_driver_union",
  description = "Authenticated diagnostic model-bound latent union", created_by = "Bdev position-driver diagnostic", fail_fast = TRUE)
up$session <- wb_session_discard(up$session)
coords <- as.data.frame(wb_space(bundle, "position_driver_union")); coords$cell_id <- rownames(coords); rownames(coords) <- NULL
write_tsv(coords[, c("cell_id", setdiff(names(coords), "cell_id"))], file.path(dirs[["tables"]], "01_position_driver_latent_coordinates.tsv"))

label_meta <- merge(labels, bdev_meta, by = "cell_id", all.x = TRUE, sort = FALSE)
label_meta$projection_role <- ifelse(label_meta$cell_id %in% reference_cells, "mapping_reference", "diagnostic_projection")
write_tsv_gz(label_meta, file.path(dirs[["tables"]], "02_provisional_label_diagnostic_cells.tsv.gz"))
ref_labels <- labels[match(reference_cells, labels$cell_id), c("cell_id", "label", "score_margin", "same_as_baseline_fraction", "stable_at_threshold")]
ref_export <- merge(data.frame(cell_id = reference_cells), bdev_meta, by = "cell_id", all.x = TRUE, sort = FALSE)
ref_export <- merge(ref_export, ref_labels, by = "cell_id", all.x = TRUE, sort = FALSE)
write_tsv_gz(ref_export, file.path(dirs[["tables"]], "03_reference_cell_technical_metadata.tsv.gz"))

# Empirical peak candidates: peaks where each bulk library and its observed
# neighbor centroid deviate in the same direction from the expected RNA-label centroid.
ledger <- read.delim(cfg$inputs$neighbor_driver_ledger$path, check.names = FALSE, quote = "\"", na.strings = "NA")
scale <- as.numeric(cfg$diagnostic$normalization_scale)
norm_cols <- function(x) { totals <- Matrix::colSums(x); log1p(t(t(x) / totals) * scale) }
ref_norm <- norm_cols(reference_counts); bulk_norm <- norm_cols(bulk_counts); label_norm <- norm_cols(label_counts)
ref_binary <- reference_counts > 0; label_binary <- label_counts > 0; bulk_binary <- bulk_counts > 0
reference_prevalence <- Matrix::rowMeans(ref_binary); bulk_prevalence <- Matrix::rowMeans(bulk_binary)
top_n <- as.integer(cfg$diagnostic$top_peaks_per_library); peak_rows <- vector("list", length(bulk_cells))
for (i in seq_along(bulk_cells)) {
  gsm <- bulk_cells[[i]]; population <- ledger$population[match(gsm, ledger$query_cell_id)]
  target_cells <- labels$cell_id[labels$label == population]
  if (!length(target_cells)) next
  x <- ledger[ledger$query_cell_id == gsm, , drop = FALSE]; w <- x$normalized_weight / sum(x$normalized_weight)
  expected <- if (length(target_cells)) Matrix::rowMeans(label_norm[, target_cells, drop = FALSE]) else rep(NA_real_, length(feature_ids))
  expected_prev <- if (length(target_cells)) Matrix::rowMeans(label_binary[, target_cells, drop = FALSE]) else rep(NA_real_, length(feature_ids))
  observed <- as.numeric(ref_norm[, x$cell_id, drop = FALSE] %*% w)
  observed_prev <- as.numeric(ref_binary[, x$cell_id, drop = FALSE] %*% w)
  query <- as.numeric(bulk_norm[, gsm]); alignment <- (query - expected) * (observed - expected)
  alignment[!is.finite(alignment)] <- -Inf
  rank <- head(order(alignment, decreasing = TRUE, method = "radix"), top_n)
  parts <- do.call(rbind, strsplit(feature_ids[rank], "[:-]"))
  peak_rows[[i]] <- data.frame(gsm = gsm, population = population, rank = seq_along(rank), feature_id = feature_ids[rank],
    chromosome = parts[, 1L], start = suppressWarnings(as.numeric(parts[, 2L])), end = suppressWarnings(as.numeric(parts[, 3L])),
    query_log_cp10k = query[rank], expected_label_log_cp10k = expected[rank], observed_neighbor_log_cp10k = observed[rank],
    expected_label_prevalence = expected_prev[rank], observed_neighbor_prevalence = observed_prev[rank],
    reference_prevalence = reference_prevalence[rank], immgen_bulk_prevalence = bulk_prevalence[rank], alignment_score = alignment[rank],
    stringsAsFactors = FALSE)
}
peak_candidates <- do.call(rbind, peak_rows); peak_candidates$width <- peak_candidates$end - peak_candidates$start
write_tsv_gz(peak_candidates, file.path(dirs[["tables"]], "04_peak_driver_candidates.tsv.gz"))
target_summary <- as.data.frame(table(labels$label), stringsAsFactors = FALSE); names(target_summary) <- c("population", "diagnostic_cells")
stable_flag <- tolower(as.character(labels$stable_at_threshold)) == "true"
stable_counts <- as.data.frame(table(factor(labels$label[stable_flag], levels = unique(labels$label))), stringsAsFactors = FALSE); names(stable_counts) <- c("population", "stable_cells")
target_summary <- merge(target_summary, stable_counts, by = "population", all.x = TRUE); target_summary$stable_cells[is.na(target_summary$stable_cells)] <- 0L
write_tsv(target_summary, file.path(dirs[["tables"]], "05_label_target_coverage.tsv"))
query_matrix_qc <- data.frame(gsm = bulk_cells, matrix_total = as.numeric(Matrix::colSums(bulk_counts)),
  nonzero_peaks = as.numeric(Matrix::colSums(bulk_counts > 0)), stringsAsFactors = FALSE)
query_matrix_qc$nonzero_fraction <- query_matrix_qc$nonzero_peaks / nrow(bulk_counts)
write_tsv(query_matrix_qc, file.path(dirs[["tables"]], "06_immgen_matrix_qc.tsv"))
axis_audit <- data.frame(feature_id = feature_ids, reference_prevalence = reference_prevalence,
  immgen_bulk_prevalence = bulk_prevalence, prevalence_gap = bulk_prevalence - reference_prevalence,
  stringsAsFactors = FALSE)
write_tsv_gz(axis_audit, file.path(dirs[["tables"]], "15_peak_axis_prevalence_audit.tsv.gz"))

source_map_after <- sha256(file.path(cfg$inputs$source_mapping_bundle$path, "manifest.json")); source_bdev_after <- sha256(file.path(cfg$inputs$source_bundle$path, "manifest.json"))
stopifnot(identical(source_map_before, source_map_after), identical(source_bdev_before, source_bdev_after))
summary <- list(schema_version = "bdev.immgen_position_driver_preparation.v1", status = "prepared", run_identity = cfg$run_identity,
  reference_cells = length(reference_cells), immgen_libraries = length(bulk_cells), label_fixture_cells = nrow(labels),
  label_cells_projected = length(new_label_cells), label_cells_already_in_reference = sum(labels$cell_id %in% reference_cells),
  candidate_peaks = nrow(peak_candidates), labels_accepted = FALSE, source_manifests_unchanged = TRUE)
write_json(summary, file.path(dirs[["reports"]], "01_preparation_summary.json"))
write_json(c(summary, list(config_path = cfg_path, config_sha256 = sha256(cfg_path))), file.path(dirs[["manifests"]], "01_preparation_receipt.json"))
wb_close_bundle(bundle)
cat("Prepared position-driver diagnostic with", length(new_label_cells), "new labelled-cell projections\n")
