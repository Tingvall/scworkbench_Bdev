#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 02_quantify_full_axis.R --config <config.yaml>", call. = FALSE)

suppressPackageStartupMessages({
  library(BPCells)
  library(Signac)
  library(GenomicRanges)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(yaml)
  library(digest)
})

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
root <- cfg$storage$results_root
matrix_dir <- cfg$storage$quantified_matrix
if (dir.exists(matrix_dir) || file.exists(matrix_dir))
  stop("Refusing to overwrite quantified matrix: ", matrix_dir, call. = FALSE)
stopifnot(
  sha256(cfg$inputs$qc_metadata$path) == cfg$inputs$qc_metadata$sha256,
  sha256(cfg$inputs$fragments$path) ==
    "006e5937205e2522f5955e8baf672c5ddde4fcd1384e8126a94e516c907114ef",
  sha256(cfg$inputs$selected_axis$path) == cfg$inputs$selected_axis$sha256,
  sha256(cfg$inputs$reciprocal_liftover_map$path) ==
    cfg$inputs$reciprocal_liftover_map$sha256,
  file.exists(paste0(cfg$inputs$fragments$path, ".tbi")))

sink(file.path(root, "04_logs", "01_quantify_full_axis.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)
cells <- as.data.table(read.delim(gzfile(file.path(root, "01_tables",
  "01_all_qc_selected_cells.tsv.gz")), check.names = FALSE, quote = "\"",
  na.strings = "NA"))
axis <- as.data.table(read.delim(gzfile(file.path(root, "01_tables",
  "02_hvp_depth20k_axis.tsv.gz")), check.names = FALSE, quote = "\"",
  na.strings = "NA"))
mapping <- as.data.table(read.delim(gzfile(
  cfg$inputs$reciprocal_liftover_map$path), check.names = FALSE,
  quote = "\"", na.strings = "NA"))
stopifnot(
  nrow(cells) == cfg$scope$expected_query_cells,
  !anyNA(cells$cell_id), !anyDuplicated(cells$cell_id),
  identical(sort(unique(cells$library_id)),
    sort(c("FLT3LKO_BM_WT", "FLT3LKO_BM_KO", "IL7KO_BM_WT", "IL7KO_BM_KO"))),
  nrow(axis) == 20000L, identical(axis$feature_order, seq_len(20000L)),
  nrow(mapping) == 20000L,
  identical(mapping$feature_id, axis$feature_id),
  sum(mapping$status == "accepted") ==
    cfg$inputs$reciprocal_liftover_map$expected_accepted_features)

accepted <- mapping[status == "accepted"]
features <- GRanges(
  seqnames = accepted$mm39_seqname,
  ranges = IRanges(start = accepted$mm39_start_1based,
                  end = accepted$mm39_end_1based))
names(features) <- accepted$feature_id
query_cells <- as.character(cells$cell_id)
log_event("Creating read-only fragment view for", length(query_cells), "cells")
fragment_object <- Signac::CreateFragmentObject(
  path = cfg$inputs$fragments$path,
  cells = setNames(query_cells, query_cells),
  validate.fragments = FALSE, verbose = FALSE)
log_event("Quantifying", length(features), "reciprocally lifted peaks")
counts <- Signac::FeatureMatrix(
  fragments = fragment_object, features = features, cells = query_cells,
  process_n = as.integer(cfg$mapping$process_n), keep_all_features = FALSE,
  verbose = TRUE)
stopifnot(setequal(colnames(counts), query_cells))
counts <- counts[, query_cells, drop = FALSE]
stopifnot(nrow(counts) == nrow(accepted))
rownames(counts) <- accepted$feature_id

triplet <- summary(counts)
full_counts <- Matrix::sparseMatrix(
  i = accepted$feature_order[triplet$i], j = triplet$j, x = triplet$x,
  dims = c(nrow(axis), length(query_cells)),
  dimnames = list(axis$feature_id, query_cells))
counts <- as(full_counts, "dgCMatrix")
rm(full_counts, triplet); invisible(gc())
stopifnot(
  identical(rownames(counts), as.character(axis$feature_id)),
  identical(colnames(counts), query_cells),
  all(Matrix::rowSums(counts[mapping$status != "accepted", , drop = FALSE]) == 0))

stage <- paste0(matrix_dir, ".tmp-", Sys.getpid())
dir.create(dirname(matrix_dir), recursive = TRUE, showWarnings = FALSE)
on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)
log_event("Writing BPCells uint32 matrix")
BPCells::write_matrix_dir(BPCells::convert_matrix_type(counts, "uint32_t"), stage)
persisted <- BPCells::open_matrix_dir(stage)
stopifnot(identical(dim(persisted), c(20000L, as.integer(cfg$scope$expected_query_cells))),
  identical(as.character(rownames(persisted)), as.character(axis$feature_id)),
  identical(as.character(colnames(persisted)), as.character(query_cells)))
if (!file.rename(stage, matrix_dir)) stop("Could not publish matrix", call. = FALSE)

cell_totals <- as.numeric(Matrix::colSums(counts))
cell_detected <- as.integer(Matrix::colSums(counts != 0))
summary_table <- cells[, .(cell_id, library_id, knockout_target, genotype, tissue)]
summary_table[, total_hvp_depth20k_fragments := cell_totals]
summary_table[, detected_hvp_depth20k_peaks := cell_detected]
fwrite(summary_table, file.path(root, "01_tables",
  "03_full_quantification_summary.tsv.gz"), sep = "\t", quote = TRUE, na = "NA")

report <- list(
  schema_version = "bdev.flt3lko_il7ko_full_quantification.v1",
  status = "passed", run_identity = cfg$run_identity,
  dimensions = unname(dim(counts)), nonzero_entries = length(counts@x),
  exact_lifted_features = nrow(accepted),
  structural_zero_features = sum(mapping$status != "accepted"),
  cells_with_zero_counts = sum(cell_totals == 0),
  cell_total_summary = as.list(summary(cell_totals)),
  cell_detected_feature_summary = as.list(summary(cell_detected)),
  exact_model_feature_order = TRUE, complete_query_axis = TRUE,
  model_refit = FALSE, umap_refit = FALSE, paper_analysis = FALSE,
  sif_built = FALSE)
write_json(report, file.path(root, "03_reports", "02_full_quantification_audit.json"))
capture.output(sessionInfo(), file = file.path(root, "06_manifests",
  "02_quantification_session_info.txt"))
log_event("Published", matrix_dir)
