#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    paste(
      "Usage: 08_validate_bdev_peak_quantification.R QUANT_DIR QC_SELECTION_TSV",
      "BDEV_PEAKS_BPCELLS BDEV_BUNDLE KO_BUNDLE OUTPUT_DIR EXPECTED_BDEV_FINGERPRINT"
    ),
    call. = FALSE
  )
}

quant_dir <- normalizePath(args[[1L]], mustWork = TRUE)
qc_path <- normalizePath(args[[2L]], mustWork = TRUE)
bdev_matrix_dir <- normalizePath(args[[3L]], mustWork = TRUE)
bdev_bundle_dir <- normalizePath(args[[4L]], mustWork = TRUE)
ko_bundle_dir <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
expected_bdev_fingerprint <- args[[7L]]
if (file.exists(output_dir) || dir.exists(output_dir)) {
  stop("Output directory already exists; refusing to overwrite: ", output_dir, call. = FALSE)
}

options(width = 160, warn = 1)
suppressPackageStartupMessages({
  library(BPCells)
  library(jsonlite)
  library(scworkbench)
})

read_gz_tsv <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))
  utils::read.delim(
    connection, header = TRUE, sep = "\t", quote = "\"",
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA"
  )
}

fingerprint <- function(path) {
  bundle <- scworkbench::wb_open_bundle(path, read_only = TRUE, validate = FALSE)
  on.exit(scworkbench::wb_close_bundle(bundle))
  scworkbench::wb_manifest(bundle)$current_state$fingerprint
}

main <- function() {
  completion <- jsonlite::read_json(file.path(quant_dir, "completion.json"), simplifyVector = TRUE)
  mapping <- read_gz_tsv(file.path(quant_dir, "accepted_bdev_peak_mapping.tsv.gz"))
  cell_summary <- read_gz_tsv(file.path(quant_dir, "cell_quantification_summary.tsv.gz"))
  query_cells <- readLines(file.path(quant_dir, "query_cell_ids.txt"), warn = FALSE)
  qc <- read_gz_tsv(qc_path)
  matrix <- BPCells::open_matrix_dir(file.path(quant_dir, "ATAC_bdev_mm10_peaks_ko_counts"))
  bdev <- BPCells::open_matrix_dir(bdev_matrix_dir)

  qc_match <- match(query_cells, qc$cell_id)
  bdev_match <- match(mapping$bdev_feature_id, rownames(bdev))
  checks <- list(
    completion_complete = identical(completion$status, "complete"),
    dimensions_expected = identical(as.integer(dim(matrix)), c(97788L, 23252L)),
    matrix_stored_as_uint32 = identical(BPCells::matrix_type(matrix), "uint32_t"),
    feature_ids_exact = identical(as.character(rownames(matrix)),
                                  as.character(mapping$bdev_feature_id)),
    feature_order_matches_bdev_rows = !anyNA(bdev_match) &&
      identical(as.integer(bdev_match), as.integer(mapping$bdev_row)),
    query_cell_ids_exact = identical(as.character(colnames(matrix)),
                                     as.character(query_cells)),
    cell_summary_ids_exact = identical(as.character(cell_summary$cell_id),
                                       as.character(query_cells)),
    qc_cells_exact = !anyNA(qc_match) && all(qc$qc_use_downstream[qc_match]),
    no_zero_features = identical(as.integer(completion$zero_features), 0L),
    no_zero_cells = identical(as.integer(completion$cells_with_zero_counts), 0L),
    nonnegative_cell_counts = all(cell_summary$total_bdev_peak_fragments >= 0),
    bdev_parent_fingerprint_unchanged = identical(
      fingerprint(bdev_bundle_dir), expected_bdev_fingerprint
    ),
    ko_parent_fingerprint_unchanged = identical(
      fingerprint(ko_bundle_dir),
      "c102c6b3c00a022ea108889924cc2094e6cebfefb534a76aa730dd982b73723b"
    )
  )
  stopifnot(all(unlist(checks)))

  source_fragments <- qc$atac_fragments[qc_match]
  quantified <- cell_summary$total_bdev_peak_fragments
  capture_fraction <- quantified / source_fragments
  source_checks <- list(
    pearson_count_vs_total_fragments = unname(cor(quantified, source_fragments)),
    spearman_count_vs_total_fragments = unname(cor(quantified, source_fragments,
                                                   method = "spearman")),
    capture_fraction_summary = as.list(summary(capture_fraction)),
    quantified_counts_do_not_exceed_total_fragments = all(capture_fraction <= 1),
    cells_exceeding_total_fragments = sum(capture_fraction > 1)
  )
  stopifnot(source_checks$spearman_count_vs_total_fragments > 0.95,
            source_checks$quantified_counts_do_not_exceed_total_fragments)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  report <- list(
    status = "passed",
    stage = "validate_bdev_peak_quantification",
    checks = checks,
    source_checks = source_checks,
    matrix = list(
      path = file.path(quant_dir, "ATAC_bdev_mm10_peaks_ko_counts"),
      dimensions = unname(dim(matrix)),
      storage_type = BPCells::matrix_type(matrix),
      canonical_features = "Bdev peaks_counts mm10 feature IDs",
      query_cells = length(query_cells)
    ),
    parent_bundles_modified = FALSE
  )
  jsonlite::write_json(
    report, file.path(output_dir, "validation.json"), pretty = TRUE,
    auto_unbox = TRUE, null = "null", na = "null", digits = NA
  )
  capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
}

main()

