#!/usr/bin/env Rscript

# Count existing GRCm39 KO fragments over the lifted lookup intervals for the
# canonical Bdev mm10 peaks. Output rows retain the original Bdev feature IDs.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "Usage: 07_quantify_bdev_peaks_in_ko.R ACCEPTED_MAPPING_TSV KO_FRAGMENTS QUERY_CELL_IDS OUTPUT_DIR PROCESS_N",
    call. = FALSE
  )
}

mapping_path <- normalizePath(args[[1L]], mustWork = TRUE)
fragments_path <- normalizePath(args[[2L]], mustWork = TRUE)
cells_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
process_n <- as.integer(args[[5L]])
if (is.na(process_n) || process_n < 1L) stop("PROCESS_N must be positive", call. = FALSE)
fragments_index_path <- paste0(fragments_path, ".tbi")
if (!file.exists(fragments_index_path)) stop("Missing fragment index: ", fragments_index_path)
if (file.exists(output_dir) || dir.exists(output_dir)) {
  stop("Output directory already exists; refusing to overwrite: ", output_dir, call. = FALSE)
}

options(width = 160, warn = 1)
suppressPackageStartupMessages({
  library(Signac)
  library(GenomicRanges)
  library(Matrix)
  library(BPCells)
  library(data.table)
  library(jsonlite)
  library(scworkbench)
})

log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}

write_json <- function(value, path) {
  jsonlite::write_json(value, path, pretty = TRUE, auto_unbox = TRUE,
                       null = "null", na = "null", digits = NA)
}

main <- function() {
  mapping_connection <- gzfile(mapping_path, open = "rt")
  mapping <- utils::read.delim(
    mapping_connection, header = TRUE, sep = "\t", quote = "\"",
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA"
  )
  close(mapping_connection)
  required <- c("bdev_row", "bdev_feature_id", "mm39_seqname",
                "mm39_start_1based", "mm39_end_1based", "status")
  stopifnot(all(required %in% names(mapping)), nrow(mapping) > 0L,
            all(mapping$status == "accepted"),
            !anyDuplicated(mapping$bdev_feature_id))
  query_cells <- readLines(cells_path, warn = FALSE)
  stopifnot(length(query_cells) == 23252L, !anyNA(query_cells),
            !anyDuplicated(query_cells), all(nzchar(query_cells)))

  lookup <- GRanges(
    seqnames = mapping$mm39_seqname,
    ranges = IRanges(start = mapping$mm39_start_1based,
                     end = mapping$mm39_end_1based)
  )
  names(lookup) <- mapping$bdev_feature_id
  stopifnot(length(lookup) == nrow(mapping), all(width(lookup) > 0L))

  parent <- dirname(output_dir)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  stage_dir <- file.path(parent, paste0(".", basename(output_dir), ".tmp-", Sys.getpid()))
  if (dir.exists(stage_dir)) unlink(stage_dir, recursive = TRUE)
  dir.create(stage_dir, recursive = FALSE)
  on.exit(if (dir.exists(stage_dir)) unlink(stage_dir, recursive = TRUE), add = TRUE)

  log_event("Creating read-only KO fragment object for", length(query_cells), "accepted query cells")
  fragment_object <- Signac::CreateFragmentObject(
    path = fragments_path,
    cells = setNames(query_cells, query_cells),
    validate.fragments = FALSE,
    verbose = FALSE
  )
  log_event("Quantifying", length(lookup), "Bdev-defined peaks in KO GRCm39 fragments")
  counts <- Signac::FeatureMatrix(
    fragments = fragment_object,
    features = lookup,
    cells = query_cells,
    process_n = process_n,
    # Signac 1.17.1 pads every chunk to the full feature universe before
    # row-binding when keep_all_features=TRUE, after which duplicate row-name
    # selection retains only the first chunk. All accepted lookup intervals
    # use standard chromosomes present in this fragment file, so FALSE is both
    # correct and required for complete chunk assembly.
    keep_all_features = FALSE,
    verbose = TRUE
  )
  if (!setequal(colnames(counts), query_cells)) {
    stop("Quantified cells do not match the requested query cells", call. = FALSE)
  }
  counts <- counts[, query_cells, drop = FALSE]
  stopifnot(nrow(counts) == length(lookup), ncol(counts) == length(query_cells))
  rownames(counts) <- mapping$bdev_feature_id
  stopifnot(all(as.character(rownames(counts)) == as.character(mapping$bdev_feature_id)),
            all(as.character(colnames(counts)) == as.character(query_cells)),
            all(counts@x >= 0), all(counts@x == floor(counts@x)))

  feature_nonzero <- Matrix::rowSums(counts != 0)
  cell_totals <- Matrix::colSums(counts)
  cell_detected <- Matrix::colSums(counts != 0)
  nonzero_features <- sum(feature_nonzero > 0)
  if (nonzero_features < 0.5 * nrow(counts) ||
      max(which(feature_nonzero > 0), 0L) <= process_n) {
    stop(
      "Feature coverage gate failed; counting may be incomplete across chunks: ",
      nonzero_features, " of ", nrow(counts), " features have signal",
      call. = FALSE
    )
  }

  matrix_dir <- file.path(stage_dir, "ATAC_bdev_mm10_peaks_ko_counts")
  log_event("Writing quantified integer matrix as BPCells")
  typed_counts <- BPCells::convert_matrix_type(counts, "uint32_t")
  BPCells::write_matrix_dir(typed_counts, matrix_dir)
  persisted <- BPCells::open_matrix_dir(matrix_dir)
  stopifnot(identical(dim(persisted), dim(counts)),
            all(as.character(rownames(persisted)) == as.character(rownames(counts))),
            all(as.character(colnames(persisted)) == as.character(colnames(counts))))

  cell_summary <- data.table(
    cell_id = query_cells,
    total_bdev_peak_fragments = as.numeric(cell_totals),
    detected_bdev_peaks = as.integer(cell_detected)
  )
  fwrite(cell_summary, file.path(stage_dir, "cell_quantification_summary.tsv.gz"),
         sep = "\t", quote = TRUE, na = "NA")
  writeLines(query_cells, file.path(stage_dir, "query_cell_ids.txt"))
  file.copy(mapping_path, file.path(stage_dir, "accepted_bdev_peak_mapping.tsv.gz"),
            overwrite = FALSE)

  completion <- list(
    status = "complete",
    stage = "quantify_bdev_peaks_in_ko",
    canonical_feature_identity = "Bdev peaks_counts mm10 feature IDs",
    lookup_coordinate_system = "mm39/GRCm39",
    count_method = "Signac::FeatureMatrix fragment overlap counting",
    matrix_id_proposed = "ATAC_bdev_mm10_peaks",
    matrix_dir = file.path(output_dir, "ATAC_bdev_mm10_peaks_ko_counts"),
    dimensions = unname(dim(counts)),
    nonzero_entries = length(counts@x),
    nonzero_features = nonzero_features,
    zero_features = sum(feature_nonzero == 0),
    cells_with_zero_counts = sum(cell_totals == 0),
    cell_total_summary = as.list(summary(as.numeric(cell_totals))),
    cell_detected_feature_summary = as.list(summary(as.numeric(cell_detected))),
    inputs = list(
      accepted_mapping = list(path = mapping_path,
                              sha256 = scworkbench::wb_source_sha256(mapping_path)),
      fragments = list(path = fragments_path,
                       sha256 = scworkbench::wb_source_sha256(fragments_path)),
      fragments_index = list(path = fragments_index_path,
                             sha256 = scworkbench::wb_source_sha256(fragments_index_path)),
      query_cells = list(path = cells_path,
                         sha256 = scworkbench::wb_source_sha256(cells_path))
    ),
    bdev_parent_modified = FALSE,
    ko_parent_modified = FALSE
  )
  write_json(completion, file.path(stage_dir, "completion.json"))
  capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

  if (!file.rename(stage_dir, output_dir)) stop("Could not publish output directory", call. = FALSE)
  log_event("Quantification complete:", nrow(counts), "peaks x", ncol(counts), "cells")
}

main()
