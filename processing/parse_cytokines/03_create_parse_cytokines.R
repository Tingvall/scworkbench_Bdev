#!/usr/bin/env Rscript

# Create and validate the parse_cytokines scWorkbench bundle from the accepted
# pre-filtered Seurat RDS.
#
# Usage:
#   Rscript --vanilla 03_create_parse_cytokines.R INPUT_RDS BUNDLE_DIR REPORT_DIR

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 03_create_parse_cytokines.R INPUT_RDS BUNDLE_DIR REPORT_DIR", call. = FALSE)
}

input_rds <- normalizePath(args[[1L]], mustWork = TRUE)
bundle_dir <- normalizePath(args[[2L]], mustWork = FALSE)
report_dir <- normalizePath(args[[3L]], mustWork = FALSE)
if (file.exists(bundle_dir) || dir.exists(bundle_dir)) {
  stop("Bundle destination already exists; refusing to overwrite: ", bundle_dir, call. = FALSE)
}
if (file.exists(report_dir) || dir.exists(report_dir)) {
  stop("Report destination already exists; refusing to overwrite: ", report_dir, call. = FALSE)
}
if (!identical(basename(bundle_dir), "bundle") ||
    !identical(basename(dirname(bundle_dir)), "parse_cytokines")) {
  stop("Expected bundle destination .../parse_cytokines/bundle", call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1L]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "..", "..", "src", "parse_cytokines_validation.R"))

options(width = 160, warn = 1)
suppressPackageStartupMessages({
  library(scworkbench)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(jsonlite)
})

expected_sha256 <- "1e4da7149190e6c93f71bc28fc47d1f58a20e22ee9cbec6894f0824394e9bc40"
expected_dimensions <- c(features = 78298L, cells = 251599L)
required_metadata <- c("sample_id", "replicate", "condition")
created_by <- Sys.getenv("SCWORKBENCH_CREATED_BY", unset = "jtg")

log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}

write_json_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  jsonlite::write_json(
    value, temporary_path, pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null", digits = NA
  )
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish JSON report: ", path, call. = FALSE)
  }
  invisible(path)
}

write_tsv_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  utils::write.table(
    value, temporary_path, sep = "\t", quote = TRUE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish TSV report: ", path, call. = FALSE)
  }
  invisible(path)
}

main <- function() {
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(report_dir)) stop("Could not create report directory: ", report_dir, call. = FALSE)

  input_stat_before <- file.info(input_rds)[, c("size", "mtime", "ctime"), drop = FALSE]
  log_event("Verifying accepted pre-filtered RDS checksum")
  input_sha256 <- scworkbench::wb_source_sha256(input_rds)
  stopifnot(identical(input_sha256, expected_sha256))

  source_input <- scworkbench::wb_source_input(
    path = input_rds,
    source_id = "parse_prefiltered_seurat_rds",
    accession = basename(input_rds),
    role = "prefiltered_source_object",
    format = "rds",
    scope = "external_file",
    sha256 = input_sha256
  )
  source_identity <- list(
    dataset_name = "parse_cytokines",
    source_type = "pre-filtered Seurat RDS",
    source_object = "Seurat object returned by readRDS()",
    extraction = "SeuratObject::LayerData(object[['RNA']], layer = 'counts')",
    inputs = list(source_input),
    prefilter = list(
      metric = "column sum of raw RNA counts",
      operator = ">=",
      threshold = 500L,
      additional_cell_filters = NULL,
      feature_filter = NULL
    )
  )
  scworkbench::wb_source_preflight(source_identity)

  log_event("Reading accepted pre-filtered Seurat RDS through its read-only mount")
  parse_seurat <- readRDS(input_rds)
  stopifnot(inherits(parse_seurat, "Seurat"), "RNA" %in% SeuratObject::Assays(parse_seurat))
  cells <- colnames(parse_seurat)
  features <- rownames(parse_seurat)
  source_metadata <- parse_seurat[[]]
  counts <- SeuratObject::LayerData(parse_seurat[["RNA"]], layer = "counts")
  stopifnot(
    inherits(counts, "sparseMatrix"),
    identical(dim(counts), unname(expected_dimensions)),
    identical(colnames(counts), cells),
    identical(rownames(counts), features),
    identical(rownames(source_metadata), cells),
    !anyNA(cells), all(nzchar(cells)), !anyDuplicated(cells),
    !anyNA(features), all(nzchar(features)), !anyDuplicated(features),
    all(required_metadata %in% names(source_metadata))
  )
  for (column in required_metadata) {
    values <- as.character(source_metadata[[column]])
    if (anyNA(values) || any(!nzchar(values))) {
      stop("Required metadata contains missing or empty values: ", column, call. = FALSE)
    }
  }

  metadata <- data.frame(
    cell_id = cells,
    sample_id = as.character(source_metadata$sample_id),
    replicate = as.character(source_metadata$replicate),
    condition = as.character(source_metadata$condition),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  sample_summary <- stats::aggregate(
    rep.int(1L, nrow(metadata)),
    by = metadata[c("sample_id", "condition", "replicate")],
    FUN = sum
  )
  names(sample_summary)[ncol(sample_summary)] <- "n_cells"
  sample_summary <- sample_summary[order(sample_summary$condition, sample_summary$replicate, sample_summary$sample_id), ]
  rownames(sample_summary) <- NULL
  stopifnot(nrow(sample_summary) == 96L, !anyDuplicated(sample_summary$sample_id), sum(sample_summary$n_cells) == nrow(metadata))

  source_identity$matrix_mapping <- list(
    RNA_counts = list(
      matrix_id = "RNA_counts",
      source_assay = "RNA",
      source_layer = "counts",
      layer = "counts",
      feature_type = "gene",
      dimensions = dim(counts)
    )
  )
  source_identity$metadata_columns <- required_metadata
  source_identity$cell_universe <- list(
    cells = length(cells), features = length(features), samples = nrow(sample_summary),
    selection = "accepted pre-filtered RDS; total raw RNA counts >= 500"
  )
  write_json_atomic(source_identity, file.path(report_dir, "source_identity.json"))
  write_json_atomic(list(
    bundle_id = "parse_cytokines",
    source_sha256 = input_sha256,
    matrices = source_identity$matrix_mapping,
    metadata_columns = required_metadata,
    cells = length(cells),
    features = length(features),
    samples = nrow(sample_summary)
  ), file.path(report_dir, "import_plan.json"))
  write_tsv_atomic(sample_summary, file.path(report_dir, "sample_summary.tsv"))

  matrices <- list(
    RNA_counts = list(
      mat = counts,
      layer = "counts",
      feature_type = "gene",
      description = paste0(
        "Raw RNA counts from accepted pre-filtered Parse Seurat RDS; all source features retained; ",
        "cells have total raw RNA counts >= 500; source SHA-256 ", input_sha256
      )
    )
  )

  log_event("Creating parse_cytokines bundle with raw RNA counts and required metadata")
  bundle <- scworkbench::wb_init_bundle(
    bundle_dir,
    bundle_id = "parse_cytokines",
    bundle_name = "parse_cytokines",
    overwrite = FALSE
  )
  on.exit(if (!is.null(bundle)) try(scworkbench::wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  source_tag <- paste0(basename(input_rds), ";sha256=", input_sha256)
  ingest <- scworkbench::wb_ingest_apply(
    bundle,
    set_name = "parse_cytokines",
    cell_ids = cells,
    metadata = metadata,
    matrices = matrices,
    source = source_tag,
    source_identity = source_identity,
    created_by = created_by,
    apply_mode = "atomic",
    emit_report_artifact = TRUE
  )
  ingest_report <- attr(ingest, "ingest_report")
  ingest_report_artifact_id <- attr(ingest, "ingest_report_artifact_id")
  write_json_atomic(ingest_report, file.path(report_dir, "ingest_report.json"))

  scworkbench::wb_close_bundle(bundle)
  bundle <- NULL
  bundle <- scworkbench::wb_open_bundle(bundle_dir)
  log_event("Comparing imported identifiers, metadata, and distributed raw-count values to the RDS")
  source_comparison <- validate_parse_cytokines_import(
    bundle, parse_seurat, required_metadata = required_metadata
  )
  write_json_atomic(source_comparison, file.path(report_dir, "source_comparison.json"))

  source_records <- scworkbench::wb_prov_source_inputs(bundle)
  stopifnot(
    "parse_prefiltered_seurat_rds" %in% source_records$source_id,
    input_sha256 %in% source_records$checksum_value
  )
  write_json_atomic(source_records, file.path(report_dir, "source_provenance.json"))

  metadata_fields <- scworkbench::wb_metadata_fields(bundle, include_core = FALSE)
  certificates <- list(
    RNA_counts = scworkbench::wb_entity_certificate(bundle, "matrix", "RNA_counts")
  )
  for (i in seq_len(nrow(metadata_fields))) {
    certificate_name <- paste0("metadata:", metadata_fields$source_name[[i]])
    certificates[[certificate_name]] <- scworkbench::wb_entity_certificate(
      bundle, "metadata_field", metadata_fields$field_id[[i]]
    )
  }
  if (length(ingest_report_artifact_id) == 1L && !is.na(ingest_report_artifact_id) &&
      nzchar(ingest_report_artifact_id)) {
    certificates$ingest_report <- scworkbench::wb_entity_certificate(
      bundle, "artifact", ingest_report_artifact_id
    )
  }
  write_json_atomic(certificates, file.path(report_dir, "entity_certificates.json"))
  write_json_atomic(scworkbench::wb_manifest(bundle), file.path(report_dir, "manifest_before_finalization.json"))
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))

  input_stat_after <- file.info(input_rds)[, c("size", "mtime", "ctime"), drop = FALSE]
  stopifnot(identical(input_stat_before, input_stat_after))

  log_event("Finalizing bundle for recipient read-only validation")
  readiness <- scworkbench::wb_finalize_recipient_bundle(bundle)
  stopifnot(isTRUE(readiness$ready))
  write_json_atomic(readiness, file.path(report_dir, "recipient_readiness.json"))
  scworkbench::wb_close_bundle(bundle)
  bundle <- NULL

  validation <- scworkbench::wb_bundle_validation_report(
    bundle_dir, level = "recipient", progress = TRUE
  )
  write_json_atomic(validation, file.path(report_dir, "validation.json"))
  stopifnot(isTRUE(validation$ok))

  recipient <- scworkbench::wb_open_bundle(bundle_dir, read_only = TRUE, validate = TRUE)
  manifest <- scworkbench::wb_manifest(recipient)
  scworkbench::wb_close_bundle(recipient)
  write_json_atomic(manifest, file.path(report_dir, "manifest.json"))

  completion <- list(
    status = "complete",
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    bundle_id = "parse_cytokines",
    bundle_dir = bundle_dir,
    source_rds = input_rds,
    source_sha256 = input_sha256,
    source_stat_unchanged = TRUE,
    matrix = list(id = "RNA_counts", layer = "counts", feature_type = "gene"),
    metadata_columns = required_metadata,
    cells = length(cells),
    features = length(features),
    samples = nrow(sample_summary),
    ingest_report_artifact_id = ingest_report_artifact_id,
    bundle_state_fingerprint = manifest$current_state$fingerprint,
    recipient_validation_fingerprint = validation$post_validation_fingerprint,
    validation_status = validation$status,
    source_comparison = source_comparison
  )
  write_json_atomic(completion, file.path(report_dir, "completion.json"))
  log_event("Complete:", bundle_dir)
  invisible(validation)
}

main()
