#!/usr/bin/env Rscript

# Prepare the unfiltered Parse RNA matrix outside a scWorkbench bundle.
#
# Usage:
#   Rscript --vanilla 01_readparsebio_unfiltered.R ALL_SAMPLE_DIR OUTPUT_DIR
#
# ALL_SAMPLE_DIR may be the Parse `all-sample` directory or its
# `DGE_unfiltered` child. The script never reads `DGE_filtered` and refuses to
# overwrite a non-empty output directory.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: 01_readparsebio_unfiltered.R ALL_SAMPLE_DIR OUTPUT_DIR",
    call. = FALSE
  )
}

input_arg <- normalizePath(args[[1L]], mustWork = TRUE)
input_dir <- if (basename(input_arg) == "DGE_unfiltered") {
  input_arg
} else {
  file.path(input_arg, "DGE_unfiltered")
}
input_dir <- normalizePath(input_dir, mustWork = TRUE)
if (basename(input_dir) != "DGE_unfiltered") {
  stop("The selected Parse input must be DGE_unfiltered", call. = FALSE)
}

output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty; refusing to overwrite: ", output_dir, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

options(width = 160, warn = 1)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(digest)
  library(jsonlite)
})

thread_limit <- suppressWarnings(as.integer(Sys.getenv("PARSE_THREADS", unset = "4")))
if (length(thread_limit) != 1L || is.na(thread_limit) || thread_limit < 1L) {
  stop("PARSE_THREADS must be a positive integer", call. = FALSE)
}
data.table::setDTthreads(thread_limit)
memory_limit_gib <- suppressWarnings(as.numeric(Sys.getenv("PARSE_MAX_MEMORY_GIB", unset = NA_character_)))

required_versions <- c(
  Seurat = "5.0.0",
  SeuratObject = "5.0.0",
  Matrix = "1.5.0",
  data.table = "1.14.0",
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

write_tsv <- function(value, path) {
  utils::write.table(
    value, file = path, sep = "\t", quote = TRUE, row.names = FALSE,
    col.names = TRUE, na = "NA"
  )
}

write_json <- function(value, path) {
  jsonlite::write_json(
    value, path, pretty = TRUE, auto_unbox = TRUE, null = "null",
    na = "null", digits = NA
  )
}

write_rds_atomic <- function(object, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  saveRDS(object, temporary_path, version = 3, compress = "gzip")
  if (!file.rename(temporary_path, path)) {
    stop("Could not rename temporary RDS into place", call. = FALSE)
  }
  invisible(path)
}

read_gzip_csv <- function(path) {
  data.table::fread(
    cmd = paste("gzip -cd --", shQuote(path)),
    sep = ",", header = TRUE, showProgress = TRUE
  )
}

decompress_gzip <- function(source, destination, chunk_bytes = 8L * 1024L * 1024L) {
  input_connection <- gzfile(source, open = "rb")
  output_connection <- file(destination, open = "wb")
  on.exit({
    close(input_connection)
    close(output_connection)
  }, add = TRUE)
  repeat {
    chunk <- readBin(input_connection, what = "raw", n = chunk_bytes)
    if (!length(chunk)) break
    writeBin(chunk, output_connection)
  }
  invisible(destination)
}

read_parsebio_compressed_delivery <- function(input_dir, output_dir, source_files) {
  # Seurat 5.5.1 detects a compressed Matrix Market file but looks for the
  # metadata using the literal names cell_metadata.csv and all_genes.csv.
  # Build those names temporarily without altering or copying the count matrix.
  staging_dir <- tempfile(pattern = ".readparsebio-input-", tmpdir = output_dir)
  dir.create(staging_dir)
  on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)

  staged_matrix <- file.path(staging_dir, "count_matrix.mtx.gz")
  if (!file.symlink(source_files[["count_matrix"]], staged_matrix)) {
    stop("Could not create the temporary read-only matrix symlink", call. = FALSE)
  }
  decompress_gzip(source_files[["cell_metadata"]], file.path(staging_dir, "cell_metadata.csv"))
  decompress_gzip(source_files[["feature_metadata"]], file.path(staging_dir, "all_genes.csv"))
  Seurat::ReadParseBio(staging_dir)
}

read_matrix_market_dimensions <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(line)) stop("Matrix Market dimensions line was not found: ", path)
    if (!startsWith(line, "%")) break
  }
  fields <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
  if (length(fields) != 3L || anyNA(suppressWarnings(as.numeric(fields)))) {
    stop("Invalid Matrix Market dimensions line: ", line, call. = FALSE)
  }
  stats::setNames(as.numeric(fields), c("source_cells", "source_features", "nonzero_entries"))
}

source_files <- c(
  count_matrix = file.path(input_dir, "count_matrix.mtx.gz"),
  cell_metadata = file.path(input_dir, "cell_metadata.csv.gz"),
  feature_metadata = file.path(input_dir, "all_genes.csv.gz")
)
missing_files <- source_files[!file.exists(source_files)]
if (length(missing_files)) {
  stop("Missing required unfiltered Parse input(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
}

mtx_candidates <- list.files(input_dir, pattern = "\\.mtx(?:\\.gz)?$", full.names = TRUE)
if (length(mtx_candidates) != 1L || normalizePath(mtx_candidates) != normalizePath(source_files[["count_matrix"]])) {
  stop("DGE_unfiltered must contain exactly one Matrix Market file: count_matrix.mtx.gz", call. = FALSE)
}

source_stat_before <- file.info(source_files)[, c("size", "mtime"), drop = FALSE]
log_event("Checksumming the three DGE_unfiltered source files")
source_inventory <- data.frame(
  role = names(source_files),
  path = unname(source_files),
  size_bytes = as.numeric(source_stat_before$size),
  modified_utc = format(source_stat_before$mtime, tz = "UTC", usetz = TRUE),
  sha256 = vapply(source_files, sha256_file, character(1L)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(source_inventory, file.path(output_dir, "source_inventory.tsv"))

matrix_market <- read_matrix_market_dimensions(source_files[["count_matrix"]])
log_event(
  "Declared unfiltered matrix:", matrix_market[["source_cells"]], "cells x",
  matrix_market[["source_features"]], "features with",
  matrix_market[["nonzero_entries"]], "stored entries"
)

log_event("Reading feature metadata")
feature_metadata <- read_gzip_csv(source_files[["feature_metadata"]])
required_feature_columns <- c("gene_id", "gene_name", "genome")
if (!all(required_feature_columns %in% names(feature_metadata))) {
  stop(
    "Feature metadata lacks required column(s): ",
    paste(setdiff(required_feature_columns, names(feature_metadata)), collapse = ", "),
    call. = FALSE
  )
}
if (nrow(feature_metadata) != matrix_market[["source_features"]]) {
  stop("Feature metadata row count does not match the Matrix Market header", call. = FALSE)
}

# Match the feature naming performed internally by Seurat::ReadParseBio():
# column 2 (gene_name), missing values replaced from column 1, then make.unique().
expected_feature_names <- as.character(feature_metadata$gene_name)
missing_feature_names <- is.na(expected_feature_names)
expected_feature_names[missing_feature_names] <- as.character(feature_metadata$gene_id[missing_feature_names])
expected_feature_names <- make.unique(expected_feature_names)

log_event("Running Seurat::ReadParseBio on DGE_unfiltered; no cell or feature filtering is applied")
counts <- read_parsebio_compressed_delivery(input_dir, output_dir, source_files)
if (!inherits(counts, "sparseMatrix")) stop("ReadParseBio did not return a sparse matrix", call. = FALSE)
if (!identical(as.numeric(dim(counts)), c(matrix_market[["source_features"]], matrix_market[["source_cells"]]))) {
  stop("ReadParseBio matrix dimensions do not match the Matrix Market header", call. = FALSE)
}
if (!identical(rownames(counts), expected_feature_names)) {
  stop("ReadParseBio feature names do not match all_genes.csv in source order", call. = FALSE)
}
if (anyDuplicated(rownames(counts))) stop("ReadParseBio returned duplicate feature names", call. = FALSE)
if (length(counts@x) != matrix_market[["nonzero_entries"]]) {
  stop("ReadParseBio stored-entry count does not match the Matrix Market header", call. = FALSE)
}

log_event("Reading and validating complete unfiltered cell metadata")
cell_metadata <- read_gzip_csv(source_files[["cell_metadata"]])
required_cell_columns <- c(
  "bc_wells", "sample", "species", "gene_count", "tscp_count",
  "mread_count", "bc1_well", "bc2_well", "bc3_well"
)
if (!all(required_cell_columns %in% names(cell_metadata))) {
  stop(
    "Cell metadata lacks required column(s): ",
    paste(setdiff(required_cell_columns, names(cell_metadata)), collapse = ", "),
    call. = FALSE
  )
}
if (nrow(cell_metadata) != ncol(counts)) {
  stop("Cell metadata row count does not match the ReadParseBio matrix", call. = FALSE)
}
cell_ids <- as.character(cell_metadata$bc_wells)
if (anyNA(cell_ids) || any(!nzchar(cell_ids))) stop("Cell identifiers contain missing or empty values", call. = FALSE)
if (anyDuplicated(cell_ids)) stop("Cell identifiers are not unique", call. = FALSE)
if (!identical(cell_ids, colnames(counts))) {
  stop("Cell metadata and ReadParseBio matrix cell orders differ", call. = FALSE)
}

# Preserve the delivered fields and add explicit fields required for later
# bundle ingestion. Sample spelling is retained exactly; only the terminal
# replicate suffix is parsed.
sample_values <- as.character(cell_metadata$sample)
replicate_pattern <- "_(?:rep|repl)([0-9]+)$"
replicate_values <- sub(paste0("^.*", replicate_pattern), "\\1", sample_values, perl = TRUE)
condition_values <- sub(replicate_pattern, "", sample_values, perl = TRUE)
unparsed_samples <- !grepl(replicate_pattern, sample_values, perl = TRUE)
if (any(unparsed_samples)) {
  stop(
    "Could not parse condition/replicate suffix for sample(s): ",
    paste(unique(sample_values[unparsed_samples]), collapse = ", "),
    call. = FALSE
  )
}
cell_metadata[, `:=`(
  orig.ident = sample_values,
  sample_id = sample_values,
  condition = condition_values,
  replicate = replicate_values,
  dataset_id = "parse_cytokines",
  source_matrix_set = "DGE_unfiltered"
)]

sample_summary <- cell_metadata[, .(n_cells = .N), by = .(sample_id, condition, replicate)]
data.table::setorder(sample_summary, condition, replicate, sample_id)
write_tsv(as.data.frame(sample_summary), file.path(output_dir, "sample_summary.tsv"))

metadata_df <- as.data.frame(cell_metadata, stringsAsFactors = FALSE)
rownames(metadata_df) <- cell_ids
rm(cell_metadata)
invisible(gc())

log_event("Creating the unfiltered Seurat object")
parse_seurat <- SeuratObject::CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  project = "parse_cytokines_unfiltered",
  meta.data = metadata_df,
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_"
)
if (ncol(parse_seurat) != matrix_market[["source_cells"]] ||
    nrow(parse_seurat) != matrix_market[["source_features"]]) {
  stop("CreateSeuratObject changed the unfiltered matrix dimensions", call. = FALSE)
}
if (!identical(colnames(parse_seurat), cell_ids)) {
  stop("CreateSeuratObject changed cell identifiers or cell order", call. = FALSE)
}
if (!identical(rownames(parse_seurat), expected_feature_names)) {
  stop("CreateSeuratObject changed feature identifiers or feature order", call. = FALSE)
}
seurat_metadata <- parse_seurat[[]]
if (!identical(rownames(seurat_metadata), cell_ids)) {
  stop("Seurat metadata cell order differs from the source cell order", call. = FALSE)
}
if (!all(names(metadata_df) %in% names(seurat_metadata))) {
  stop("CreateSeuratObject did not retain every source and derived metadata column", call. = FALSE)
}
metadata_columns_match <- vapply(
  names(metadata_df),
  function(column_name) identical(metadata_df[[column_name]], seurat_metadata[[column_name]]),
  logical(1L)
)
if (!all(metadata_columns_match)) {
  stop(
    "Metadata values changed during Seurat construction for column(s): ",
    paste(names(metadata_columns_match)[!metadata_columns_match], collapse = ", "),
    call. = FALSE
  )
}
rm(seurat_metadata)

feature_df <- data.frame(
  gene_id = as.character(feature_metadata$gene_id),
  gene_name_original = as.character(feature_metadata$gene_name),
  genome = as.character(feature_metadata$genome),
  stringsAsFactors = FALSE,
  row.names = expected_feature_names,
  check.names = FALSE
)
parse_seurat[["RNA"]] <- SeuratObject::AddMetaData(parse_seurat[["RNA"]], metadata = feature_df)
retained_feature_metadata <- parse_seurat[["RNA"]][[]]
if (!identical(rownames(retained_feature_metadata), expected_feature_names) ||
    !identical(retained_feature_metadata$gene_id, feature_df$gene_id) ||
    !identical(retained_feature_metadata$gene_name_original, feature_df$gene_name_original) ||
    !identical(retained_feature_metadata$genome, feature_df$genome)) {
  stop("Seurat assay feature metadata does not match all_genes.csv", call. = FALSE)
}
rm(retained_feature_metadata)

package_names <- names(required_versions)
package_versions <- vapply(package_names, function(x) as.character(packageVersion(x)), character(1L))
environment_manifest <- data.frame(
  component = c("R", package_names, "thread_limit", "virtual_memory_limit_gib"),
  version = c(
    as.character(getRversion()), unname(package_versions), as.character(thread_limit),
    if (is.na(memory_limit_gib)) "not_reported" else as.character(memory_limit_gib)
  ),
  stringsAsFactors = FALSE
)
write_tsv(environment_manifest, file.path(output_dir, "environment.tsv"))
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))

parse_seurat@misc$source_provenance <- list(
  dataset_id = "parse_cytokines",
  source_selection = "DGE_unfiltered",
  source_directory = input_dir,
  inputs = split(source_inventory, seq_len(nrow(source_inventory))),
  reader = list(
    function_name = "Seurat::ReadParseBio",
    seurat_version = unname(package_versions[["Seurat"]]),
    delivery_adapter = "temporary decompression of .csv.gz metadata to the literal filenames required by ReadParseBio; count matrix linked read-only",
    feature_column = "gene_name",
    matrix_transposed_by_reader = TRUE
  ),
  filtering = list(cell_filter = NULL, feature_filter = NULL),
  environment = list(
    image = Sys.getenv("SCWORKBENCH_IMAGE", unset = NA_character_),
    image_sha256 = Sys.getenv("SCWORKBENCH_IMAGE_SHA256", unset = NA_character_),
    thread_limit = thread_limit,
    virtual_memory_limit_gib = if (is.na(memory_limit_gib)) NULL else memory_limit_gib
  )
)

rm(counts, metadata_df, feature_metadata, feature_df)
invisible(gc())

output_rds <- file.path(output_dir, "parse_cytokines_unfiltered_seurat.rds")
log_event("Writing Seurat RDS atomically:", output_rds)
write_rds_atomic(parse_seurat, output_rds)

log_event("Checksumming the Seurat RDS")
output_sha256 <- sha256_file(output_rds)
source_stat_after <- file.info(source_files)[, c("size", "mtime"), drop = FALSE]
source_unchanged <- identical(as.numeric(source_stat_before$size), as.numeric(source_stat_after$size)) &&
  identical(as.numeric(source_stat_before$mtime), as.numeric(source_stat_after$mtime))
if (!source_unchanged) stop("One or more source file sizes or modification times changed during the run", call. = FALSE)

validation <- list(
  passed = TRUE,
  source_selection = list(
    requested = "DGE_unfiltered",
    resolved_directory = input_dir,
    filtered_input_read = FALSE,
    source_files_unchanged = source_unchanged
  ),
  matrix = list(
    class = class(SeuratObject::LayerData(parse_seurat[["RNA"]], layer = "counts")),
    features = nrow(parse_seurat),
    cells = ncol(parse_seurat),
    stored_entries = matrix_market[["nonzero_entries"]],
    feature_order_matches_source = TRUE,
    cell_order_matches_source = TRUE,
    duplicate_feature_ids = 0L,
    duplicate_cell_ids = 0L
  ),
  metadata = list(
    delivered_columns = required_cell_columns,
    added_columns = c("orig.ident", "sample_id", "condition", "replicate", "dataset_id", "source_matrix_set"),
    sample_count = nrow(sample_summary),
    cells_retained = ncol(parse_seurat),
    cells_removed = 0L,
    features_removed = 0L
  ),
  output = list(
    path = output_rds,
    size_bytes = as.numeric(file.info(output_rds)$size),
    sha256 = output_sha256
  ),
  environment = list(
    R = as.character(getRversion()),
    packages = as.list(package_versions),
    image = Sys.getenv("SCWORKBENCH_IMAGE", unset = NA_character_),
    image_sha256 = Sys.getenv("SCWORKBENCH_IMAGE_SHA256", unset = NA_character_),
    thread_limit = thread_limit,
    virtual_memory_limit_gib = if (is.na(memory_limit_gib)) NULL else memory_limit_gib
  )
)
write_json(validation, file.path(output_dir, "validation.json"))
write_json(
  list(
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    status = "complete",
    validation_passed = TRUE,
    output_rds = output_rds,
    output_sha256 = output_sha256
  ),
  file.path(output_dir, "completion.json")
)
log_event("Complete:", output_rds)
