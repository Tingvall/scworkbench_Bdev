#!/usr/bin/env Rscript

# Prefilter the outside-bundle Parse Seurat object by total raw RNA counts.
#
# Usage:
#   Rscript --vanilla 02_prefilter_min500counts.R INPUT_RDS OUTPUT_RDS
#
# The input is never modified. Cells are retained exactly when the column sum
# of the RNA counts layer is at least 500; no sample or additional QC filter is
# applied. The script refuses to overwrite its RDS or validation sidecars.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: 02_prefilter_min500counts.R INPUT_RDS OUTPUT_RDS",
    call. = FALSE
  )
}

input_rds <- normalizePath(args[[1L]], mustWork = TRUE)
output_rds <- normalizePath(args[[2L]], mustWork = FALSE)
if (identical(input_rds, output_rds)) {
  stop("Input and output RDS paths must differ", call. = FALSE)
}
if (file.exists(output_rds)) {
  stop("Output RDS already exists; refusing to overwrite: ", output_rds, call. = FALSE)
}
if (!dir.exists(dirname(output_rds))) {
  stop("Output directory does not exist: ", dirname(output_rds), call. = FALSE)
}

output_stem <- tools::file_path_sans_ext(output_rds)
sidecar_paths <- c(
  validation = paste0(output_stem, "_validation.json"),
  sample_summary = paste0(output_stem, "_sample_summary.tsv"),
  environment = paste0(output_stem, "_environment.tsv"),
  session_info = paste0(output_stem, "_session_info.txt")
)
existing_sidecars <- sidecar_paths[file.exists(sidecar_paths)]
if (length(existing_sidecars)) {
  stop(
    "Validation sidecar(s) already exist; refusing to overwrite: ",
    paste(existing_sidecars, collapse = ", "), call. = FALSE
  )
}

options(width = 160, warn = 1)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(digest)
  library(jsonlite)
})

min_counts <- 500
thread_limit <- suppressWarnings(as.integer(Sys.getenv("PARSE_THREADS", unset = "4")))
if (length(thread_limit) != 1L || is.na(thread_limit) || thread_limit < 1L) {
  stop("PARSE_THREADS must be a positive integer", call. = FALSE)
}
memory_limit_gib <- suppressWarnings(as.numeric(Sys.getenv("PARSE_MAX_MEMORY_GIB", unset = NA_character_)))

required_versions <- c(
  Seurat = "5.0.0",
  SeuratObject = "5.0.0",
  Matrix = "1.5.0",
  digest = "0.6.0",
  jsonlite = "1.8.0"
)
for (package_name in names(required_versions)) {
  if (packageVersion(package_name) < numeric_version(required_versions[[package_name]])) {
    stop(
      package_name, " >= ", required_versions[[package_name]], " is required; found ",
      as.character(packageVersion(package_name)), call. = FALSE
    )
  }
}

log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256")
}

write_rds_atomic <- function(object, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  saveRDS(object, temporary_path, version = 3, compress = "gzip")
  if (!file.rename(temporary_path, path)) {
    stop("Could not rename temporary RDS into place: ", path, call. = FALSE)
  }
  invisible(path)
}

write_text_atomic <- function(lines, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  writeLines(lines, temporary_path, useBytes = TRUE)
  if (!file.rename(temporary_path, path)) {
    stop("Could not rename temporary text file into place: ", path, call. = FALSE)
  }
  invisible(path)
}

write_tsv_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  utils::write.table(
    value, file = temporary_path, sep = "\t", quote = TRUE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
  if (!file.rename(temporary_path, path)) {
    stop("Could not rename temporary TSV into place: ", path, call. = FALSE)
  }
  invisible(path)
}

write_json_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  jsonlite::write_json(
    value, temporary_path, pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null", digits = NA
  )
  if (!file.rename(temporary_path, path)) {
    stop("Could not rename temporary JSON into place: ", path, call. = FALSE)
  }
  invisible(path)
}

summarize_counts <- function(x) {
  quantiles <- stats::quantile(
    x, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
    names = FALSE, type = 7
  )
  stats::setNames(
    as.list(as.numeric(quantiles)),
    c("min", "p01", "p05", "p25", "median", "p75", "p95", "p99", "max")
  )
}

input_stat_before <- file.info(input_rds)[, c("size", "mtime"), drop = FALSE]
log_event("Checksumming the unfiltered input RDS")
input_sha256 <- sha256_file(input_rds)

log_event("Reading the unfiltered Seurat RDS:", input_rds)
parse_seurat <- readRDS(input_rds)
if (!inherits(parse_seurat, "Seurat")) stop("Input RDS is not a Seurat object", call. = FALSE)
if (!"RNA" %in% names(parse_seurat@assays)) stop("Input Seurat object has no RNA assay", call. = FALSE)

input_cells <- colnames(parse_seurat)
input_features <- rownames(parse_seurat)
if (anyNA(input_cells) || any(!nzchar(input_cells)) || anyDuplicated(input_cells)) {
  stop("Input cell identifiers are missing, empty, or duplicated", call. = FALSE)
}
if (anyNA(input_features) || any(!nzchar(input_features)) || anyDuplicated(input_features)) {
  stop("Input feature identifiers are missing, empty, or duplicated", call. = FALSE)
}
input_metadata <- parse_seurat@meta.data
if (!identical(rownames(input_metadata), input_cells)) {
  stop("Input Seurat metadata is not aligned to the cell order", call. = FALSE)
}
if (!"sample_id" %in% names(input_metadata)) {
  stop("Input Seurat metadata lacks sample_id", call. = FALSE)
}

log_event("Calculating total raw RNA counts per cell from the counts layer")
counts <- SeuratObject::LayerData(parse_seurat[["RNA"]], layer = "counts")
if (!inherits(counts, "sparseMatrix")) stop("RNA counts layer is not a sparse Matrix", call. = FALSE)
if (!identical(dim(counts), c(length(input_features), length(input_cells))) ||
    !identical(rownames(counts), input_features) ||
    !identical(colnames(counts), input_cells)) {
  stop("RNA counts layer dimensions or identifiers do not align to the Seurat object", call. = FALSE)
}
raw_rna_counts <- Matrix::colSums(counts)
if (length(raw_rna_counts) != length(input_cells) || anyNA(raw_rna_counts) ||
    any(!is.finite(raw_rna_counts)) || any(raw_rna_counts < 0)) {
  stop("Calculated raw RNA count totals are invalid", call. = FALSE)
}

ncount_metadata_matches <- FALSE
if ("nCount_RNA" %in% names(input_metadata)) {
  ncount_metadata <- as.numeric(input_metadata[["nCount_RNA"]])
  ncount_metadata_matches <- length(ncount_metadata) == length(raw_rna_counts) &&
    !anyNA(ncount_metadata) && all(ncount_metadata == as.numeric(raw_rna_counts))
  if (!ncount_metadata_matches) {
    stop("nCount_RNA does not exactly match the raw RNA counts-layer column sums", call. = FALSE)
  }
}

retain <- as.numeric(raw_rna_counts) >= min_counts
if (!any(retain)) stop("No cells meet the minimum-count threshold", call. = FALSE)
retained_cells <- input_cells[retain]
cells_before <- length(input_cells)
cells_retained <- length(retained_cells)
cells_removed <- cells_before - cells_retained
log_event(
  "Retaining", cells_retained, "of", cells_before,
  "cells with raw RNA counts >=", min_counts, "; removing", cells_removed
)

sample_ids <- as.character(input_metadata$sample_id)
if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
  stop("sample_id contains missing or empty values", call. = FALSE)
}
sample_levels <- unique(sample_ids)
sample_index <- match(sample_ids, sample_levels)
sample_summary <- data.frame(
  sample_id = sample_levels,
  cells_before = tabulate(sample_index, nbins = length(sample_levels)),
  cells_retained = tabulate(sample_index[retain], nbins = length(sample_levels)),
  stringsAsFactors = FALSE
)
sample_summary$cells_removed <- sample_summary$cells_before - sample_summary$cells_retained
if (sum(sample_summary$cells_before) != cells_before ||
    sum(sample_summary$cells_retained) != cells_retained ||
    sum(sample_summary$cells_removed) != cells_removed) {
  stop("Sample-level cell counts do not reconcile to the object totals", call. = FALSE)
}

input_feature_metadata <- parse_seurat[["RNA"]][[]]
log_event("Subsetting the Seurat object; no feature or additional cell filter is applied")
filtered_seurat <- parse_seurat[, retained_cells]
if (ncol(filtered_seurat) != cells_retained || nrow(filtered_seurat) != length(input_features)) {
  stop("Filtered Seurat object dimensions are incorrect", call. = FALSE)
}
if (!identical(colnames(filtered_seurat), retained_cells)) {
  stop("Filtering changed the retained cell order", call. = FALSE)
}
if (!identical(rownames(filtered_seurat), input_features)) {
  stop("Filtering changed the feature identifiers or feature order", call. = FALSE)
}
if (!identical(filtered_seurat@meta.data, input_metadata[retain, , drop = FALSE])) {
  stop("Filtering changed retained cell metadata values or order", call. = FALSE)
}
if (!identical(filtered_seurat[["RNA"]][[]], input_feature_metadata)) {
  stop("Filtering changed RNA feature metadata values or order", call. = FALSE)
}

filtered_counts <- SeuratObject::LayerData(filtered_seurat[["RNA"]], layer = "counts")
filtered_raw_rna_counts <- Matrix::colSums(filtered_counts)
if (!identical(as.numeric(filtered_raw_rna_counts), as.numeric(raw_rna_counts[retain]))) {
  stop("Filtered raw RNA counts differ from the selected input columns", call. = FALSE)
}
if (any(filtered_raw_rna_counts < min_counts)) {
  stop("Filtered object contains a cell below the minimum-count threshold", call. = FALSE)
}
if (any(as.numeric(raw_rna_counts[!retain]) >= min_counts)) {
  stop("A removed cell meets the minimum-count threshold", call. = FALSE)
}

package_names <- names(required_versions)
package_versions <- vapply(package_names, function(x) as.character(packageVersion(x)), character(1L))
completed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
filtered_seurat@misc$prefiltering <- list(
  completed_utc = completed_utc,
  method = "RNA counts-layer column sum",
  threshold = list(metric = "total_raw_RNA_counts", operator = ">=", value = min_counts),
  additional_cell_filters = NULL,
  feature_filter = NULL,
  input = list(path = input_rds, sha256 = input_sha256, cells = cells_before, features = length(input_features)),
  result = list(cells_retained = cells_retained, cells_removed = cells_removed, features_retained = length(input_features)),
  environment = list(
    image = Sys.getenv("SCWORKBENCH_IMAGE", unset = NA_character_),
    image_sha256 = Sys.getenv("SCWORKBENCH_IMAGE_SHA256", unset = NA_character_),
    R = as.character(getRversion()),
    packages = as.list(package_versions),
    thread_limit = thread_limit,
    virtual_memory_limit_gib = if (is.na(memory_limit_gib)) NULL else memory_limit_gib
  )
)

rm(counts, parse_seurat, input_metadata, input_feature_metadata, retained_cells)
invisible(gc())

log_event("Writing the prefiltered Seurat RDS atomically:", output_rds)
write_rds_atomic(filtered_seurat, output_rds)
log_event("Checksumming the prefiltered Seurat RDS")
output_sha256 <- sha256_file(output_rds)

input_stat_after <- file.info(input_rds)[, c("size", "mtime"), drop = FALSE]
input_unchanged <- identical(as.numeric(input_stat_before$size), as.numeric(input_stat_after$size)) &&
  identical(as.numeric(input_stat_before$mtime), as.numeric(input_stat_after$mtime))
if (!input_unchanged) stop("Input RDS size or modification time changed during filtering", call. = FALSE)

environment_manifest <- data.frame(
  component = c("R", package_names, "thread_limit", "virtual_memory_limit_gib"),
  version = c(
    as.character(getRversion()), unname(package_versions), as.character(thread_limit),
    if (is.na(memory_limit_gib)) "not_reported" else as.character(memory_limit_gib)
  ),
  stringsAsFactors = FALSE
)
validation <- list(
  passed = TRUE,
  completed_utc = completed_utc,
  filter = list(
    assay = "RNA",
    layer = "counts",
    metric = "column sum of raw RNA counts",
    operator = ">=",
    threshold = min_counts,
    additional_cell_filters = NULL,
    feature_filter = NULL
  ),
  input = list(
    path = input_rds,
    size_bytes = as.numeric(input_stat_before$size),
    sha256 = input_sha256,
    unchanged = input_unchanged,
    cells = cells_before,
    features = length(input_features)
  ),
  output = list(
    path = output_rds,
    size_bytes = as.numeric(file.info(output_rds)$size),
    sha256 = output_sha256,
    cells = cells_retained,
    cells_removed = cells_removed,
    features = nrow(filtered_seurat),
    samples = nrow(sample_summary)
  ),
  checks = list(
    nCount_RNA_matches_counts_layer = ncount_metadata_matches,
    all_retained_cells_meet_threshold = TRUE,
    all_removed_cells_are_below_threshold = TRUE,
    cell_order_is_input_subset_order = TRUE,
    feature_order_unchanged = TRUE,
    cell_metadata_unchanged = TRUE,
    feature_metadata_unchanged = TRUE,
    sample_counts_reconcile = TRUE
  ),
  raw_RNA_count_distribution_before = summarize_counts(raw_rna_counts),
  raw_RNA_count_distribution_retained = summarize_counts(filtered_raw_rna_counts),
  environment = list(
    R = as.character(getRversion()),
    packages = as.list(package_versions),
    image = Sys.getenv("SCWORKBENCH_IMAGE", unset = NA_character_),
    image_sha256 = Sys.getenv("SCWORKBENCH_IMAGE_SHA256", unset = NA_character_),
    thread_limit = thread_limit,
    virtual_memory_limit_gib = if (is.na(memory_limit_gib)) NULL else memory_limit_gib
  )
)

write_tsv_atomic(sample_summary, sidecar_paths[["sample_summary"]])
write_tsv_atomic(environment_manifest, sidecar_paths[["environment"]])
write_text_atomic(capture.output(sessionInfo()), sidecar_paths[["session_info"]])
write_json_atomic(validation, sidecar_paths[["validation"]])
log_event("Complete:", output_rds)
