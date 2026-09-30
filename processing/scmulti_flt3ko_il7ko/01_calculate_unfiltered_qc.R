#!/usr/bin/env Rscript

# Calculate factual RNA/ATAC QC metrics for every cell in a Cell Ranger ARC
# aggregate. This script does not select, filter, or remove cells.
# Usage: Rscript --vanilla 01_calculate_unfiltered_qc.R \
#   CELLRANGER_OUTS PER_BARCODE_METRICS_DIR REFERENCE_DIR OUTPUT_DIR

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: 01_calculate_unfiltered_qc.R CELLRANGER_OUTS PER_BARCODE_METRICS_DIR REFERENCE_DIR OUTPUT_DIR")
outs_dir <- normalizePath(args[[1L]], mustWork = TRUE)
per_barcode_dir <- normalizePath(args[[2L]], mustWork = TRUE)
reference_dir <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- args[[4L]]
options(width = 180, warn = 1)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
completion_path <- file.path(output_dir, "completion.json")
if (file.exists(completion_path)) stop("Refusing to overwrite completed QC output: ", output_dir)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(jsonlite)
  library(Matrix)
  library(rhdf5)
  library(SeuratObject)
  library(Signac)
})
log_event <- function(...) { cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n"); flush.console() }

required_files <- c(
  matrix_h5 = file.path(outs_dir, "filtered_feature_bc_matrix.h5"),
  aggr_csv = file.path(outs_dir, "aggr.csv"),
  fragments = file.path(outs_dir, "atac_fragments.tsv.gz"),
  fragments_index = file.path(outs_dir, "atac_fragments.tsv.gz.tbi"),
  peaks = file.path(outs_dir, "atac_peaks.bed"),
  summary = file.path(outs_dir, "summary.csv"),
  web_summary = file.path(outs_dir, "web_summary.html"),
  tss = file.path(reference_dir, "regions", "tss.bed"),
  reference_json = file.path(reference_dir, "reference.json")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing required inputs: ", paste(missing_files, collapse = ", "))

sample_def <- data.table(
  aggregate_suffix = 1:4,
  library_id = c("FLT3LKO_BM_WT", "FLT3LKO_BM_KO", "IL7KO_BM_WT", "IL7KO_BM_KO"),
  knockout_target = c("Flt3l", "Flt3l", "Il7", "Il7"),
  genotype = c("WT", "KO", "WT", "KO"),
  tissue = "BM"
)
sample_def[, per_barcode_path := file.path(per_barcode_dir, paste0(library_id, "_per_barcode_metrics.csv"))]
if (any(!file.exists(sample_def$per_barcode_path))) stop("Missing per-barcode metrics")

aggr <- fread(required_files[["aggr_csv"]])
stopifnot(
  identical(as.character(aggr$library_id), sample_def$library_id),
  identical(as.character(basename(aggr$per_barcode_metrics)), basename(sample_def$per_barcode_path))
)
reference_metadata <- fromJSON(required_files[["reference_json"]], simplifyVector = TRUE)
stopifnot(identical(as.character(reference_metadata$genomes), "GRCm39"))

checkpoint_path <- file.path(output_dir, "pre_fragment_checkpoint.rds")
if (file.exists(checkpoint_path)) {
  log_event("Resuming from validated pre-fragment checkpoint")
  checkpoint <- readRDS(checkpoint_path)
  cell_qc <- checkpoint$cell_qc
  h5_shape <- checkpoint$h5_shape
  barcodes <- checkpoint$barcodes
  feature_type <- checkpoint$feature_type
  n_cells <- checkpoint$n_cells
  suffix_counts <- checkpoint$suffix_counts
} else {
log_event("Reading aggregate HDF5 structure")
h5_shape <- as.integer(h5read(required_files[["matrix_h5"]], "matrix/shape"))
barcodes <- as.character(h5read(required_files[["matrix_h5"]], "matrix/barcodes"))
feature_type <- as.character(h5read(required_files[["matrix_h5"]], "matrix/features/feature_type"))
feature_name <- as.character(h5read(required_files[["matrix_h5"]], "matrix/features/name"))
feature_interval <- as.character(h5read(required_files[["matrix_h5"]], "matrix/features/interval"))
indptr <- as.numeric(h5read(required_files[["matrix_h5"]], "matrix/indptr"))
stopifnot(
  length(barcodes) == h5_shape[[2L]], length(feature_type) == h5_shape[[1L]],
  length(indptr) == length(barcodes) + 1L, !anyDuplicated(barcodes), all(diff(indptr) >= 0)
)

barcode_suffix <- as.integer(sub("^.*-", "", barcodes))
cell_qc <- data.table(
  cell_id = barcodes,
  original_barcode = paste0(sub("-[0-9]+$", "", barcodes), "-1"),
  aggregate_suffix = barcode_suffix
)
cell_qc <- sample_def[cell_qc, on = "aggregate_suffix"]
stopifnot(nrow(cell_qc) == length(barcodes), !anyNA(cell_qc$library_id))
suffix_counts <- cell_qc[, .N, by = .(aggregate_suffix, library_id)]
setorder(suffix_counts, aggregate_suffix)
stopifnot(identical(as.integer(suffix_counts$N), c(6753L, 8477L, 8521L, 4790L)))

log_event("Calculating count-matrix metrics in HDF5 chunks")
is_gene <- feature_type == "Gene Expression"
is_peak <- feature_type == "Peaks"
is_mito_gene <- is_gene & (grepl("^mt-", feature_name, ignore.case = TRUE) | grepl("^chrM[:_]", feature_interval))
is_ribosomal_gene <- is_gene & grepl("^Rp[sl][0-9A-Za-z-]*$", feature_name)
n_cells <- length(barcodes)
matrix_metrics <- data.table(
  cell_id = barcodes, matrix_rna_umis = numeric(n_cells), matrix_detected_genes = integer(n_cells),
  matrix_mito_umis = numeric(n_cells), matrix_ribosomal_umis = numeric(n_cells),
  matrix_atac_peak_counts = numeric(n_cells), matrix_detected_peaks = integer(n_cells)
)
sum_by_group <- function(values, groups, n_groups) {
  out <- numeric(n_groups)
  if (!length(values)) return(out)
  aggregated <- rowsum(values, groups, reorder = FALSE)
  out[as.integer(rownames(aggregated))] <- aggregated[, 1L]
  out
}
chunk_size <- 256L
for (cell_start in seq.int(1L, n_cells, by = chunk_size)) {
  cell_end <- min(n_cells, cell_start + chunk_size - 1L)
  n_chunk <- cell_end - cell_start + 1L
  ptr <- indptr[cell_start:(cell_end + 1L)]
  entry_start <- ptr[[1L]] + 1
  entry_end <- ptr[[length(ptr)]]
  if (entry_end >= entry_start) {
    entry_index <- seq.int(entry_start, entry_end)
    rows <- as.integer(h5read(required_files[["matrix_h5"]], "matrix/indices", index = list(entry_index))) + 1L
    values <- as.numeric(h5read(required_files[["matrix_h5"]], "matrix/data", index = list(entry_index)))
    groups <- rep.int(seq_len(n_chunk), diff(ptr))
    gene_entries <- is_gene[rows]; peak_entries <- is_peak[rows]
    mito_entries <- is_mito_gene[rows]; ribosomal_entries <- is_ribosomal_gene[rows]
    destination <- cell_start:cell_end
    matrix_metrics$matrix_rna_umis[destination] <- sum_by_group(values[gene_entries], groups[gene_entries], n_chunk)
    matrix_metrics$matrix_detected_genes[destination] <- tabulate(groups[gene_entries], nbins = n_chunk)
    matrix_metrics$matrix_mito_umis[destination] <- sum_by_group(values[mito_entries], groups[mito_entries], n_chunk)
    matrix_metrics$matrix_ribosomal_umis[destination] <- sum_by_group(values[ribosomal_entries], groups[ribosomal_entries], n_chunk)
    matrix_metrics$matrix_atac_peak_counts[destination] <- sum_by_group(values[peak_entries], groups[peak_entries], n_chunk)
    matrix_metrics$matrix_detected_peaks[destination] <- tabulate(groups[peak_entries], nbins = n_chunk)
  }
  if (cell_end %% 2048L < chunk_size || cell_end == n_cells) log_event("Processed matrix cells", cell_end, "of", n_cells)
}
matrix_metrics[, matrix_mito_percent := fifelse(matrix_rna_umis > 0, 100 * matrix_mito_umis / matrix_rna_umis, NA_real_)]
matrix_metrics[, matrix_ribosomal_percent := fifelse(matrix_rna_umis > 0, 100 * matrix_ribosomal_umis / matrix_rna_umis, NA_real_)]
cell_qc <- matrix_metrics[cell_qc, on = "cell_id"]
stopifnot(
  length(cell_qc$cell_id) == length(barcodes),
  all(as.character(cell_qc$cell_id) == as.character(barcodes))
)

source_metric_columns <- c(
  "barcode", "gex_barcode", "atac_barcode", "is_cell", "excluded_reason",
  "gex_raw_reads", "gex_mapped_reads", "gex_umis_count", "gex_genes_count",
  "atac_raw_reads", "atac_unmapped_reads", "atac_lowmapq", "atac_dup_reads",
  "atac_chimeric_reads", "atac_mitochondrial_reads", "atac_fragments",
  "atac_TSS_fragments", "atac_peak_region_fragments", "atac_peak_region_cutsites"
)
log_event("Reading original per-barcode Cell Ranger metrics")
source_metrics <- rbindlist(lapply(seq_len(nrow(sample_def)), function(i) {
  definition <- sample_def[i]
  metrics <- fread(definition$per_barcode_path, select = source_metric_columns, showProgress = FALSE)
  wanted <- cell_qc[library_id == definition$library_id, .(cell_id, original_barcode)]
  matched <- match(wanted$original_barcode, metrics$barcode)
  if (anyNA(matched)) stop("Unmatched barcodes for ", definition$library_id, ": ", sum(is.na(matched)))
  metrics <- metrics[matched]
  stopifnot(identical(as.character(metrics$barcode), wanted$original_barcode))
  metrics[, cell_id := wanted$cell_id]
  metrics[, library_id_check := definition$library_id]
  metrics
}), use.names = TRUE)
stopifnot(nrow(source_metrics) == n_cells, !anyDuplicated(source_metrics$cell_id), all(source_metrics$is_cell == 1L))
cell_qc <- source_metrics[cell_qc, on = "cell_id"]
stopifnot(identical(cell_qc$cell_id, barcodes), identical(cell_qc$library_id, cell_qc$library_id_check))
cell_qc[, library_id_check := NULL]
cell_qc[, atac_frip := fifelse(atac_fragments > 0, atac_peak_region_fragments / atac_fragments, NA_real_)]
cell_qc[, atac_tss_fragment_fraction := fifelse(atac_fragments > 0, atac_TSS_fragments / atac_fragments, NA_real_)]
cell_qc[, atac_mitochondrial_read_fraction := fifelse(atac_raw_reads > 0, atac_mitochondrial_reads / atac_raw_reads, NA_real_)]
saveRDS(
  list(
    cell_qc = cell_qc, h5_shape = h5_shape, barcodes = barcodes,
    feature_type = feature_type, n_cells = n_cells, suffix_counts = suffix_counts
  ),
  checkpoint_path
)
}

log_event("Constructing transient fragment-backed Signac object")
tss_table <- fread(required_files[["tss"]], col.names = c("chrom", "start", "end", "gene_id", "score", "strand"), showProgress = FALSE)
tss_table <- unique(tss_table, by = c("chrom", "start", "end", "strand"))
tss_positions <- makeGRangesFromDataFrame(
  as.data.frame(tss_table), keep.extra.columns = TRUE, ignore.strand = FALSE,
  seqnames.field = "chrom", start.field = "start", end.field = "end",
  strand.field = "strand", starts.in.df.are.0based = TRUE
)
tss_positions <- tss_positions[as.character(seqnames(tss_positions)) != "chrM"]
fragment_object <- CreateFragmentObject(
  path = required_files[["fragments"]], cells = setNames(barcodes, barcodes),
  validate.fragments = FALSE, verbose = FALSE
)
dummy_counts <- sparseMatrix(
  i = rep.int(1:2, n_cells), j = rep(seq_len(n_cells), each = 2L),
  x = rep.int(1, 2L * n_cells), dims = c(2L, n_cells),
  dimnames = list(c("chr1-1-2", "chr1-4-5"), barcodes)
)
chromatin_assay <- CreateChromatinAssay(
  counts = dummy_counts, sep = c("-", "-"), fragments = fragment_object,
  min.cells = 0, min.features = 0
)
qc_object <- CreateSeuratObject(counts = chromatin_assay, assay = "ATAC")
tss_checkpoint_path <- file.path(output_dir, "tss_checkpoint.rds")
if (file.exists(tss_checkpoint_path)) {
  log_event("Resuming from completed TSS-enrichment checkpoint")
  qc_object <- readRDS(tss_checkpoint_path)
} else {
  log_event("Calculating TSS enrichment against GRCm39-2024-A TSS positions")
  qc_object <- suppressWarnings(TSSEnrichment(
    object = qc_object, tss.positions = tss_positions, fast = TRUE, assay = "ATAC",
    process_n = 2000L, verbose = TRUE, region_extension = 1000L
  ))
  saveRDS(qc_object, tss_checkpoint_path)
}
log_event("Calculating nucleosome signal using the standard Signac sampling target")
nucleosome_n <- n_cells * 5000L
qc_object <- suppressWarnings(NucleosomeSignal(object = qc_object, assay = "ATAC", n = nucleosome_n, verbose = TRUE))
signac_metrics <- as.data.table(qc_object[[]], keep.rownames = "cell_id")[
  , .(cell_id, TSS.enrichment, TSS.percentile, nucleosome_signal, nucleosome_percentile)
]
setnames(signac_metrics, c("TSS.enrichment", "TSS.percentile"), c("tss_enrichment", "tss_percentile"))
signac_match <- match(cell_qc$cell_id, signac_metrics$cell_id)
stopifnot(!anyNA(signac_match), !anyDuplicated(signac_metrics$cell_id))
cell_qc[, `:=`(
  tss_enrichment = signac_metrics$tss_enrichment[signac_match],
  tss_percentile = signac_metrics$tss_percentile[signac_match],
  nucleosome_signal = signac_metrics$nucleosome_signal[signac_match],
  nucleosome_percentile = signac_metrics$nucleosome_percentile[signac_match]
)]
stopifnot(
  length(cell_qc$cell_id) == length(barcodes),
  all(as.character(cell_qc$cell_id) == as.character(barcodes))
)

metric_labels <- c(
  gex_umis_count = "RNA UMIs (Cell Ranger)", gex_genes_count = "Detected genes (Cell Ranger)",
  matrix_mito_percent = "Mitochondrial RNA (%)", matrix_ribosomal_percent = "Ribosomal RNA (%)",
  atac_fragments = "ATAC fragments (Cell Ranger)", atac_frip = "ATAC FRiP",
  atac_tss_fragment_fraction = "ATAC TSS fragment fraction", tss_enrichment = "TSS enrichment (Signac)",
  nucleosome_signal = "Nucleosome signal (Signac)",
  atac_mitochondrial_read_fraction = "ATAC mitochondrial read fraction",
  matrix_detected_peaks = "Detected aggregate peaks"
)
summary_long <- melt(
  cell_qc, id.vars = c("library_id", "knockout_target", "genotype", "tissue"),
  measure.vars = names(metric_labels), variable.name = "metric", value.name = "value"
)[, .(
  cells = .N, missing = sum(is.na(value)), minimum = min(value, na.rm = TRUE),
  q1 = quantile(value, 0.25, na.rm = TRUE, names = FALSE), median = median(value, na.rm = TRUE),
  q3 = quantile(value, 0.75, na.rm = TRUE, names = FALSE), maximum = max(value, na.rm = TRUE)
), by = .(library_id, knockout_target, genotype, tissue, metric)]
summary_long[, metric_label := unname(metric_labels[metric])]
setcolorder(summary_long, c("library_id", "knockout_target", "genotype", "tissue", "metric", "metric_label", "cells", "missing", "minimum", "q1", "median", "q3", "maximum"))

source_comparison <- list(
  gex_umis_exact = all(cell_qc$matrix_rna_umis == cell_qc$gex_umis_count),
  gex_umis_mismatch_cells = sum(cell_qc$matrix_rna_umis != cell_qc$gex_umis_count),
  detected_genes_exact = all(cell_qc$matrix_detected_genes == cell_qc$gex_genes_count),
  detected_genes_mismatch_cells = sum(cell_qc$matrix_detected_genes != cell_qc$gex_genes_count),
  peak_counts_equal_peak_region_cutsites = all(cell_qc$matrix_atac_peak_counts == cell_qc$atac_peak_region_cutsites),
  peak_count_mismatch_cells = sum(cell_qc$matrix_atac_peak_counts != cell_qc$atac_peak_region_cutsites)
)

fwrite(cell_qc, file.path(output_dir, "cell_qc_metrics.tsv.gz"), sep = "\t", compress = "gzip", na = "NA")
fwrite(summary_long, file.path(output_dir, "sample_qc_summary.tsv"), sep = "\t", na = "NA")
fwrite(
  sample_def[, setdiff(names(sample_def), "per_barcode_path"), with = FALSE],
  file.path(output_dir, "sample_manifest.tsv"), sep = "\t", na = "NA"
)
inventory_paths <- c(required_files, setNames(sample_def$per_barcode_path, paste0("per_barcode_", sample_def$library_id)))
inventory_info <- file.info(inventory_paths)
log_event("Calculating SHA-256 checksums for retained source inputs")
inventory_sha256 <- vapply(unname(inventory_paths), function(path) {
  result <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) stop("sha256sum failed for ", path)
  strsplit(result[[1L]], "[[:space:]]+")[[1L]][[1L]]
}, character(1L))
inventory <- data.table(
  role = names(inventory_paths), path = unname(inventory_paths), bytes = as.numeric(inventory_info$size),
  modified = format(inventory_info$mtime, "%Y-%m-%dT%H:%M:%S%z"), sha256 = inventory_sha256
)
fwrite(inventory, file.path(output_dir, "input_inventory.tsv"), sep = "\t")

validation <- list(
  status = "ok", filtering_performed = FALSE, thresholding_performed = FALSE,
  input_cells = n_cells, reported_cells = nrow(cell_qc), unique_cell_ids = uniqueN(cell_qc$cell_id),
  all_cells_are_cellranger_cells = all(cell_qc$is_cell == 1L),
  sample_counts = split(suffix_counts$N, suffix_counts$library_id), matrix_shape = unname(h5_shape),
  feature_counts = as.list(table(feature_type)), tss_positions = length(tss_positions),
  nucleosome_fragment_line_target = nucleosome_n, source_matrix_comparison = source_comparison
)
write_json(validation, file.path(output_dir, "validation.json"), pretty = TRUE, auto_unbox = TRUE, null = "null", digits = 16)
context <- list(
  analysis = "Unfiltered per-library QC for FLT3LKO and IL7KO bone-marrow multiome",
  technical_bundle_id = "scmulti_flt3ko_il7ko", corrected_biological_identity = "scmulti_flt3lko_il7ko",
  knockout_targets = c("Flt3l", "Il7"), cellranger_pipeline = "cellranger-arc aggr",
  cellranger_version = "2.2.0", normalization = "none", genome = "GRCm39",
  reference_version = as.character(reference_metadata$version),
  reference_gtf_hash = as.character(reference_metadata[["gtf_hash.gz"]]),
  filtering_performed = FALSE, thresholds_selected = FALSE,
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_json(context, file.path(output_dir, "qc_context.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
unlink(checkpoint_path)
unlink(tss_checkpoint_path)
write_json(
  list(
    status = "complete", completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    filtering_performed = FALSE, input_cells = n_cells, output_cells = nrow(cell_qc),
    outputs = c("cell_qc_metrics.tsv.gz", "sample_qc_summary.tsv", "sample_manifest.tsv", "input_inventory.tsv", "validation.json", "qc_context.json", "session_info.txt")
  ), completion_path, pretty = TRUE, auto_unbox = TRUE
)
log_event("Completed unfiltered QC for", nrow(cell_qc), "cells at", output_dir)
