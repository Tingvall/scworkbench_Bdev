#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 160)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(
  "Usage: 03_create_source_bundle.R FEATURECOUNTS_TSV CONSENSUS_BED SAMPLE_MANIFEST BUNDLE_DIR REPORT_DIR",
  call. = FALSE
)

count_path <- normalizePath(args[[1L]], mustWork = TRUE)
peaks_path <- normalizePath(args[[2L]], mustWork = TRUE)
manifest_path <- normalizePath(args[[3L]], mustWork = TRUE)
bundle_dir <- normalizePath(args[[4L]], mustWork = FALSE)
report_dir <- normalizePath(args[[5L]], mustWork = FALSE)
if (file.exists(bundle_dir) || dir.exists(bundle_dir))
  stop("Bundle destination already exists: ", bundle_dir, call. = FALSE)
if (file.exists(report_dir) || dir.exists(report_dir))
  stop("Report destination already exists: ", report_dir, call. = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(jsonlite)
  library(Matrix)
  library(digest)
})

write_json_atomic <- function(value, path) {
  tmp <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  write_json(value, tmp, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  if (!file.rename(tmp, path)) stop("Could not publish ", path, call. = FALSE)
}

manifest <- fread(manifest_path, data.table = FALSE)
stopifnot(
  nrow(manifest) == 6L, !anyDuplicated(manifest$sample_id),
  sum(manifest$condition == "control") == 3L,
  sum(manifest$condition == "TCF3_depleted") == 3L
)
manifest$observation_unit <- "bulk_ATAC_library"
manifest$cell_id <- manifest$sample_id

fc <- fread(count_path, skip = "Geneid")
annotation_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
stopifnot(all(annotation_columns %in% names(fc)))
count_columns <- setdiff(names(fc), annotation_columns)
count_samples <- sub("\\.filtered\\.bam$", "", basename(count_columns))
stopifnot(setequal(count_samples, manifest$sample_id), !anyDuplicated(fc$Geneid))
count_columns <- count_columns[match(manifest$sample_id, count_samples)]
counts <- as.matrix(fc[, ..count_columns])
storage.mode(counts) <- "integer"
rownames(counts) <- fc$Geneid
colnames(counts) <- manifest$sample_id
stopifnot(!anyNA(counts), all(counts >= 0L))
counts <- methods::as(counts, "dgCMatrix")

source_inputs <- list(
  wb_source_input(
    path = count_path, source_id = "consensus_peak_featurecounts",
    accession = basename(count_path), role = "raw_count_matrix", format = "tsv",
    scope = "external_file", sha256 = digest(count_path, algo = "sha256", file = TRUE)
  ),
  wb_source_input(
    path = peaks_path, source_id = "consensus_peak_coordinates",
    accession = basename(peaks_path), role = "feature_coordinates", format = "bed",
    scope = "external_file", sha256 = digest(peaks_path, algo = "sha256", file = TRUE)
  ),
  wb_source_input(
    path = manifest_path, source_id = "gse259355_sample_manifest",
    accession = basename(manifest_path), role = "sample_metadata", format = "tsv",
    scope = "external_file", sha256 = digest(manifest_path, algo = "sha256", file = TRUE)
  )
)
source_identity <- list(
  dataset_name = "GSE259355 pro-B TCF3 acute-depletion ATAC-seq",
  technical_bundle_id = "gse259355_tcf3_prob",
  biological_label = "bone marrow pro-B cells",
  source_type = "bulk paired-end ATAC-seq consensus-peak fragment counts",
  external_preprocessing_boundary = paste(
    "SRA paired-end reads aligned to UCSC mm10; proper-pair MAPQ>=30 alignments retained;",
    "duplicates, non-primary chromosomes, chrM, and Boyle mm10 blacklist regions removed;",
    "MACS3 BAMPE sample peaks merged and counted with featureCounts."
  ),
  genome = list(organism = "Mus musculus", build = "GRCm38", reference = "UCSC mm10"),
  inputs = source_inputs,
  observation_semantics = list(
    columns_are = "bulk_ATAC_libraries",
    columns_are_not = "single_cells",
    libraries = ncol(counts),
    comparison = "three acute TCF3/E2A-depleted libraries versus three untreated controls"
  )
)
wb_source_preflight(source_identity)

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
metadata <- manifest[, c("cell_id", setdiff(names(manifest), c("cell_id", "sample_id")), "sample_id")]
matrices <- list(
  ATAC_consensus_peak_counts = list(
    mat = counts, layer = "counts", feature_type = "peak",
    description = paste(
      "Raw fragment counts in the merged GSE259355 consensus-peak universe.",
      "Columns are six independent bulk ATAC libraries, not single cells; coordinates use UCSC mm10."
    )
  )
)

bundle <- wb_init_bundle(
  bundle_dir, bundle_id = "gse259355_tcf3_prob",
  bundle_name = "GSE259355 pro-B TCF3 acute-depletion bulk ATAC", overwrite = FALSE
)
on.exit(if (!is.null(bundle)) try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
ingest <- wb_ingest_apply(
  bundle, set_name = "gse259355_tcf3_prob", cell_ids = manifest$sample_id,
  metadata = metadata, matrices = matrices,
  source = "GSE259355 SRA libraries processed by the repository bulk ATAC workflow",
  source_identity = source_identity,
  created_by = Sys.getenv("SCWORKBENCH_CREATED_BY", unset = "jtg"),
  apply_mode = "atomic", emit_report_artifact = TRUE
)
write_json_atomic(attr(ingest, "ingest_report"), file.path(report_dir, "ingest_report.json"))
wb_matrix_set_genome(bundle, "ATAC_consensus_peak_counts", "Mus musculus", "GRCm38", "UCSC mm10")
wb_close_bundle(bundle)
bundle <- NULL

bundle <- wb_open_bundle(bundle_dir)
stopifnot(
  identical(wb_cells(bundle, "all"), manifest$sample_id),
  identical(wb_cells(bundle, "set:gse259355_tcf3_prob"), manifest$sample_id)
)
imported <- wb_matrix(bundle, "ATAC_consensus_peak_counts")
stopifnot(
  identical(dim(imported), dim(counts)), identical(rownames(imported), rownames(counts)),
  identical(colnames(imported), colnames(counts))
)
imported_values <- as.matrix(wb_matrix_read_features(
  bundle, "ATAC_consensus_peak_counts", feature_ids = rownames(counts),
  cell_ids = colnames(counts), missing = "NA"
))
stopifnot(identical(as.numeric(imported_values), as.numeric(counts)))
imported_metadata <- wb_meta(bundle, cells = manifest$sample_id)
stopifnot(identical(as.character(imported_metadata$sample_id), manifest$sample_id))
write_json_atomic(wb_prov_source_inputs(bundle), file.path(report_dir, "source_provenance.json"))
write_json_atomic(
  wb_entity_certificate(bundle, "matrix", "ATAC_consensus_peak_counts"),
  file.path(report_dir, "matrix_certificate.json")
)
readiness <- wb_finalize_recipient_bundle(bundle)
stopifnot(isTRUE(readiness$ready))
write_json_atomic(readiness, file.path(report_dir, "recipient_readiness.json"))
wb_close_bundle(bundle)
bundle <- NULL

validation <- wb_bundle_validation_report(bundle_dir, level = "recipient", progress = TRUE)
stopifnot(isTRUE(validation$ok))
write_json_atomic(validation, file.path(report_dir, "validation.json"))
recipient <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = TRUE)
bundle_manifest <- wb_manifest(recipient)
stopifnot(identical(wb_cells(recipient, "all"), manifest$sample_id))
wb_close_bundle(recipient)
capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
write_json_atomic(list(
  status = "complete",
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  bundle_id = "gse259355_tcf3_prob", bundle_dir = bundle_dir,
  observation_unit = "bulk_ATAC_library", libraries = ncol(counts),
  consensus_peaks = nrow(counts), matrix_id = "ATAC_consensus_peak_counts",
  bundle_state_fingerprint = bundle_manifest$current_state$fingerprint,
  recipient_validation_fingerprint = validation$post_validation_fingerprint,
  validation_status = validation$status
), file.path(report_dir, "completion.json"))
