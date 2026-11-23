#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste("Usage: 04_project_fixture.R BUNDLE MATRIX_DIR FIXTURE_TSV",
             "RESULTS_ROOT SESSIONS_DIR MODEL_PAYLOAD"), call. = FALSE)
}

bundle_path <- normalizePath(args[[1L]], mustWork = TRUE)
matrix_path <- normalizePath(args[[2L]], mustWork = TRUE)
fixture_path <- normalizePath(args[[3L]], mustWork = TRUE)
results_root <- normalizePath(args[[4L]], mustWork = FALSE)
sessions_dir <- normalizePath(args[[5L]], mustWork = FALSE)
model_payload_path <- normalizePath(args[[6L]], mustWork = TRUE)

options(stringsAsFactors = FALSE, warn = 1, width = 180)
suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbLatent)
  library(BPCells)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(digest)
})

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
as_space_matrix <- function(x) {
  if (is.data.frame(x)) {
    ids <- as.character(x$cell_id)
    x <- as.matrix(x[, setdiff(names(x), "cell_id"), drop = FALSE])
    rownames(x) <- ids
  } else x <- as.matrix(x)
  x
}
promoted_id <- function(x, output) {
  take <- x$promoted$status == "ok" & x$promoted$output == output
  stopifnot(sum(take) == 1L)
  as.character(x$promoted$artifact_id[take][[1L]])
}

dirs <- setNames(file.path(results_root, c("01_tables", "03_reports", "04_logs", "06_manifests")),
                 c("tables", "reports", "logs", "manifests"))
invisible(lapply(c(dirs, sessions = sessions_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "02_projection.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

manifest <- system.file("scworkbench", "operation-pack.json", package = "scwbLatent")
stopifnot(nzchar(manifest))
try(wb_operation_pack_unload("scwblatent"), silent = TRUE)
wb_operation_pack_activate(manifest, lib.loc = .libPaths())
stopifnot("scwblatent::peakvi_project" %in% wb_list_ops())

bundle <- wb_open_bundle(bundle_path, validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
all_cells <- as.character(wb_cells(bundle, "all"))
slots <- all_cells[grepl("^external_query_[0-9]+$", all_cells)]
stopifnot(length(slots) == 28L)

fixture <- fread(fixture_path)
fixture <- fixture[order(library_id, cell_id)]
projection_cells <- fixture[, head(.SD, 7L), by = library_id]
stopifnot(nrow(projection_cells) == 28L, all(table(projection_cells$library_id) == 7L))
projection_cells[, slot_id := slots]
setcolorder(projection_cells, c("slot_id", setdiff(names(projection_cells), "slot_id")))
fwrite(projection_cells, file.path(dirs[["tables"]], "05_projection_slot_map.tsv"),
       sep = "\t", quote = TRUE, na = "NA")

matrix <- BPCells::open_matrix_dir(matrix_path)
feature_ids <- as.character(rownames(matrix))
stopifnot(nrow(matrix) == 20000L, ncol(matrix) == 400L,
          all(projection_cells$cell_id %in% colnames(matrix)))
query <- matrix[, projection_cells$cell_id]
query <- as(query, "dgCMatrix")
colnames(query) <- slots
stopifnot(identical(as.character(rownames(query)), feature_ids),
          identical(as.character(colnames(query)), slots))

matrix_id <- "hvp_depth20k_flt3lko_il7ko_fixture_counts"
stopifnot(!(matrix_id %in% wb_matrices(bundle)$matrix_id))
wb_attach_matrix(bundle, matrix_id, query, layer = "counts", feature_type = "peak",
  cell_ids = slots,
  description = "Balanced 28-cell FLT3LKO/IL7KO bounded projection fixture",
  source = "active_development_exact_liftover_quantification_with_33_feature_structural_zero_mask",
  created_by = "Bdev development capability 05 bounded mapping qualification")

model_artifact_id <- "model_ref_20261117_051703_wi5ql81y"
model_payload_before <- sha256(model_payload_path)
stopifnot(identical(model_payload_before,
  "024d11e3efb4822c31cd75c5371d6320860026f46dadd9ff80d680712df134a8"))
model <- wb_art_load(bundle, model_artifact_id, expected_kind = "model_ref")
stopifnot(as.integer(model$model_linkage$latent_dimension_count) == 10L)

space_id <- "frozen_hvp_depth20k_flt3lko_il7ko_fixture"
session <- wb_session_create(bundle, sessions_dir,
  "flt3lko_il7ko_mapping_fixture_20261117T175427Z_peakvi_project")
projected <- wb_session_run_op(session, bundle, "scwblatent::peakvi_project",
  inputs = list(
    model = wb_ref_artifact(model_artifact_id, "model"),
    query_matrix = matrix_id,
    space_id_out = space_id,
    chunk_size = 32L,
    batch_key = NULL,
    use_gpu = FALSE,
    pickle_execution_trust = "operator_verified_source_v1",
    cell_set = "all",
    min_cells = 1L
  ),
  notes = "Bounded balanced FLT3LKO/IL7KO projection through frozen HVP-depth20k model without refit")
summary <- projected$run$outputs$projection_summary$table
stopifnot(nrow(summary) == 1L, summary$transformation[[1L]] == "transfer_without_refit",
          as.integer(summary$n_latent[[1L]]) == 10L)
projected$session <- wb_session_resume(bundle, projected$session$durable$workspace)
promoted <- wb_session_promote_current(bundle, projected$session,
  outputs = "projection_summary", name_prefix = "flt3lko_il7ko_fixture_projection",
  description = "Bounded balanced KO projection without reference refit",
  created_by = "Bdev development capability 05 bounded mapping qualification",
  fail_fast = TRUE)
projection_artifact_id <- promoted_id(promoted, "projection_summary")
promoted$session <- wb_session_discard(promoted$session)

coordinates <- as_space_matrix(wb_space(bundle, space_id, cells = slots))
stopifnot(identical(rownames(coordinates), slots), ncol(coordinates) == 10L)
coordinate_table <- data.frame(slot_id = slots, coordinates,
                               check.names = FALSE, stringsAsFactors = FALSE)
coordinate_table <- merge(projection_cells[, .(slot_id, cell_id, library_id,
  knockout_target, genotype, tissue)], coordinate_table, by = "slot_id",
  sort = FALSE)
coordinate_table <- coordinate_table[match(slots, coordinate_table$slot_id), ]
fwrite(coordinate_table, file.path(dirs[["tables"]], "06_projected_coordinates.tsv.gz"),
       sep = "\t", quote = TRUE, na = "NA")

model_payload_after <- sha256(model_payload_path)
stopifnot(identical(model_payload_after, model_payload_before))
report <- list(
  status = "passed",
  projected_cells = nrow(coordinate_table),
  cells_per_library = as.list(table(coordinate_table$library_id)),
  latent_dimensions = ncol(coordinates),
  transformation = "transfer_without_refit",
  model_artifact_id = model_artifact_id,
  model_payload_sha256_before = model_payload_before,
  model_payload_sha256_after = model_payload_after,
  projection_artifact_id = projection_artifact_id,
  source_archives_modified = FALSE,
  paper_analysis = FALSE
)
write_json(report, file.path(dirs[["reports"]], "04_projection_audit.json"))
write_json(list(scworkbench = as.character(packageVersion("scworkbench")),
                scwbLatent = as.character(packageVersion("scwbLatent")),
                scwbLatent_manifest_sha256 = sha256(manifest)),
           file.path(dirs[["manifests"]], "03_projection_runtime.json"))
