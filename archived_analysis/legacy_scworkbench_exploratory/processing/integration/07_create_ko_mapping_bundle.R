#!/usr/bin/env Rscript

# Create an isolated working integration bundle containing one shared,
# Bdev-anchored ATAC peak matrix for accepted Bdev reference cells and accepted
# FLT3LKO/IL7KO query cells. Parent bundles are opened read-only.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    paste(
      "Usage: 07_create_ko_mapping_bundle.R BDEV_BUNDLE KO_BUNDLE PEAK_MAPPING",
      "KO_BPCELLS OUTPUT_BUNDLE REPORT_DIR CREATED_BY"
    ), call. = FALSE
  )
}

bdev_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
ko_bundle_dir <- normalizePath(args[[2L]], mustWork = TRUE)
mapping_path <- normalizePath(args[[3L]], mustWork = TRUE)
ko_matrix_dir <- normalizePath(args[[4L]], mustWork = TRUE)
output_bundle_dir <- normalizePath(args[[5L]], mustWork = FALSE)
report_dir <- normalizePath(args[[6L]], mustWork = FALSE)
created_by <- args[[7L]]
if (file.exists(output_bundle_dir) || dir.exists(output_bundle_dir)) {
  stop("Output bundle already exists; refusing to overwrite: ", output_bundle_dir, call. = FALSE)
}
if (file.exists(report_dir) || dir.exists(report_dir)) {
  stop("Report directory already exists; refusing to overwrite: ", report_dir, call. = FALSE)
}

options(width = 160, warn = 1)
suppressPackageStartupMessages({
  library(scworkbench)
  library(BPCells)
  library(data.table)
  library(jsonlite)
  library(Matrix)
})

expected_bdev_fingerprint <- "1027f6390c6b853acc2bf40d3fae9abb59393f4256cb39285685eb2b2973f7ec"
expected_ko_fingerprint <- "c102c6b3c00a022ea108889924cc2094e6cebfefb534a76aa730dd982b73723b"

log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}

write_json <- function(value, name) {
  jsonlite::write_json(
    value, file.path(report_dir, name), pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null", digits = NA
  )
}

read_gz_tsv <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))
  utils::read.delim(
    connection, header = TRUE, sep = "\t", quote = "\"",
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA"
  )
}

manifest_fingerprint <- function(bundle) wb_manifest(bundle)$current_state$fingerprint

distributed_indices <- function(n, target) {
  unique(as.integer(round(seq.int(1L, n, length.out = min(as.integer(target), n)))))
}

main <- function() {
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    if (!file.exists(file.path(report_dir, "creation.json"))) {
      unlink(report_dir, recursive = TRUE)
    }
  }, add = TRUE)

  bdev <- wb_open_bundle(bdev_bundle_dir, read_only = TRUE, validate = FALSE)
  on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
  ko <- wb_open_bundle(ko_bundle_dir, read_only = TRUE, validate = FALSE)
  on.exit(try(wb_close_bundle(ko), silent = TRUE), add = TRUE)
  stopifnot(
    identical(manifest_fingerprint(bdev), expected_bdev_fingerprint),
    identical(manifest_fingerprint(ko), expected_ko_fingerprint)
  )

  mapping <- read_gz_tsv(mapping_path)
  stopifnot(nrow(mapping) == 97788L, all(mapping$status == "accepted"),
            !anyDuplicated(mapping$bdev_feature_id),
            identical(as.integer(mapping$bdev_row), sort(as.integer(mapping$bdev_row))))

  bdev_source_cells <- wb_cells(bdev, "all")
  ko_source_cells <- colnames(BPCells::open_matrix_dir(ko_matrix_dir))
  ko_accepted_cells <- wb_cells(ko, "set:qc_use_downstream")
  stopifnot(length(bdev_source_cells) == 26723L,
            length(ko_source_cells) == 23252L,
            setequal(as.character(ko_source_cells), as.character(ko_accepted_cells)))

  bdev_cells <- paste0("bdev::", bdev_source_cells)
  ko_cells <- paste0("ko::", ko_source_cells)
  all_cells <- c(bdev_cells, ko_cells)
  stopifnot(!anyDuplicated(all_cells))

  bdev_md <- wb_meta(bdev, cells = bdev_source_cells)
  ko_md <- wb_meta(ko, cells = ko_source_cells)
  stopifnot(identical(as.character(bdev_md$cell_id), bdev_source_cells),
            identical(as.character(ko_md$cell_id), ko_source_cells))
  metadata <- data.frame(
    cell_id = all_cells,
    source_bundle = c(rep("bdev_reference", length(bdev_cells)),
                      rep("scmulti_flt3ko_il7ko", length(ko_cells))),
    source_cell_id = c(bdev_source_cells, ko_source_cells),
    mapping_role = c(rep("reference", length(bdev_cells)),
                     rep("query", length(ko_cells))),
    source_genome = c(rep("mm10/GRCm38", length(bdev_cells)),
                      rep("mm39/GRCm39", length(ko_cells))),
    sample_id = c(as.character(bdev_md$orig.ident), as.character(ko_md$library_id)),
    genotype = c(rep(NA_character_, length(bdev_cells)), as.character(ko_md$genotype)),
    knockout_target = c(rep(NA_character_, length(bdev_cells)),
                        as.character(ko_md$knockout_target)),
    tissue = c(rep(NA_character_, length(bdev_cells)), as.character(ko_md$tissue)),
    qc_use_downstream = c(rep(TRUE, length(bdev_cells)),
                          as.logical(ko_md$qc_use_downstream)),
    reference_pt = c(as.numeric(bdev_md$pt), rep(NA_real_, length(ko_cells))),
    reference_phase = c(as.character(bdev_md$Phase), rep(NA_character_, length(ko_cells))),
    reference_peak_cluster = c(as.character(bdev_md[["peaks_snn_res.1.2"]]),
                               rep(NA_character_, length(ko_cells))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  stopifnot(identical(metadata$cell_id, all_cells),
            identical(as.numeric(metadata$reference_pt[seq_along(bdev_cells)]),
                      as.numeric(bdev_md$pt)),
            identical(is.na(metadata$reference_pt[seq_along(bdev_cells)]),
                      is.na(bdev_md$pt)))

  log_event("Opening and aligning Bdev and KO count matrices")
  bdev_counts <- wb_matrix(bdev, "peaks_counts")
  ko_counts <- BPCells::open_matrix_dir(ko_matrix_dir)
  stopifnot(nrow(bdev_counts) == 97846L,
            all(as.character(rownames(bdev_counts)[mapping$bdev_row]) ==
                  as.character(mapping$bdev_feature_id)),
            all(as.character(rownames(ko_counts)) == as.character(mapping$bdev_feature_id)),
            all(as.character(colnames(ko_counts)) == as.character(ko_source_cells)))
  bdev_selected <- bdev_counts[mapping$bdev_row, , drop = FALSE]
  bdev_selected <- BPCells::convert_matrix_type(bdev_selected, "uint32_t")
  colnames(bdev_selected) <- bdev_cells
  colnames(ko_counts) <- ko_cells
  common_counts <- cbind(bdev_selected, ko_counts)
  stopifnot(identical(dim(common_counts), c(97788L, 49975L)),
            all(as.character(rownames(common_counts)) == as.character(mapping$bdev_feature_id)),
            all(as.character(colnames(common_counts)) == as.character(all_cells)),
            identical(BPCells::matrix_type(common_counts), "uint32_t"))
  log_event("Materializing the combined sparse matrix for the scWorkbench ingest boundary")
  common_sparse <- methods::as(common_counts, "dgCMatrix")
  stopifnot(identical(dim(common_sparse), dim(common_counts)),
            all(common_sparse@x >= 0), all(common_sparse@x == floor(common_sparse@x)))

  quant_completion <- file.path(dirname(ko_matrix_dir), "completion.json")
  quant_validation <- file.path(
    dirname(dirname(dirname(ko_matrix_dir))), "reports",
    "bdev_peak_quantification_v1_20260930", "validation.json"
  )
  stopifnot(file.exists(quant_completion), file.exists(quant_validation))
  source_identity <- list(
    dataset_name = "bdev_integrated_ko_mapping_v1",
    phase = "bdev_trained_peakvi_projection",
    reference_parent = list(bundle_id = "bdev_reference",
                            fingerprint = expected_bdev_fingerprint),
    query_parent = list(bundle_id = "scmulti_flt3ko_il7ko",
                        fingerprint = expected_ko_fingerprint),
    canonical_features = list(
      source = "bdev_reference peaks_counts mm10",
      retained = nrow(mapping),
      query_lookup = "unique exact-reciprocal mm39 liftOver intervals",
      canonical_ids = "original Bdev mm10 peak IDs"
    ),
    cell_identity = list(
      reference_prefix = "bdev::", query_prefix = "ko::",
      reference_cells = length(bdev_cells), query_cells = length(ko_cells)
    ),
    inputs = list(
      wb_source_input(path = mapping_path, source_id = "accepted_bdev_peak_mapping",
                      accession = basename(mapping_path), role = "feature_mapping",
                      format = "tsv_gz", scope = "external_file",
                      sha256 = wb_source_sha256(mapping_path)),
      wb_source_input(path = quant_completion, source_id = "ko_peak_quantification",
                      accession = basename(dirname(ko_matrix_dir)), role = "matrix_creation_receipt",
                      format = "json", scope = "external_file",
                      sha256 = wb_source_sha256(quant_completion)),
      wb_source_input(path = quant_validation, source_id = "ko_peak_quantification_validation",
                      accession = "bdev_peak_quantification_v1_20260930", role = "validation",
                      format = "json", scope = "external_file",
                      sha256 = wb_source_sha256(quant_validation))
    )
  )
  wb_source_preflight(source_identity)
  write_json(source_identity, "source_identity.json")

  log_event("Creating isolated Bdev/KO mapping bundle with", length(all_cells), "cells")
  bundle <- wb_init_bundle(
    output_bundle_dir,
    bundle_id = "bdev_integrated_ko_mapping_v1",
    bundle_name = "Bdev-trained PeakVI mapping of FLT3LKO/IL7KO query",
    overwrite = FALSE
  )
  on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  ingest <- wb_ingest_apply(
    bundle,
    set_name = "bdev_ko_mapping",
    cell_ids = all_cells,
    metadata = metadata,
    matrices = list(
      ATAC_bdev_mm10_peaks = list(
        mat = common_sparse,
        layer = "counts",
        feature_type = "peak",
        description = paste(
          "Bdev mm10 peaks_counts reference universe; KO GRCm39 fragments quantified",
          "over exact reciprocal lookup intervals; canonical rows are Bdev mm10 peak IDs"
        )
      )
    ),
    source = "accepted bdev_reference and scmulti_flt3ko_il7ko parents",
    source_identity = source_identity,
    created_by = created_by,
    apply_mode = "atomic",
    emit_report_artifact = TRUE
  )
  write_json(attr(ingest, "ingest_report"), "ingest_report.json")

  wb_cell_set_derive(
    bundle, cell_set = "set:bdev_reference", cells = bdev_cells,
    role = "analysis_universe", status = "accepted",
    source = "accepted bdev_reference parent",
    created_by = created_by,
    description = "Bdev cells used as the PeakVI reference",
    overwrite = FALSE, promote_manifest = TRUE
  )
  wb_cell_set_derive(
    bundle, cell_set = "set:ko_query", cells = ko_cells,
    role = "analysis_universe", status = "accepted",
    source = "accepted scmulti_flt3ko_il7ko QC universe",
    created_by = created_by,
    description = "FLT3LKO/IL7KO cells accepted for Bdev projection",
    overwrite = FALSE, promote_manifest = TRUE
  )

  stopifnot(identical(wb_cells(bundle, "all"), all_cells),
            identical(wb_cells(bundle, "set:bdev_reference"), bdev_cells),
            identical(wb_cells(bundle, "set:ko_query"), ko_cells),
            all(as.character(wb_feature_ids(bundle, "ATAC_bdev_mm10_peaks")) ==
                  as.character(mapping$bdev_feature_id)))
  feature_idx <- distributed_indices(nrow(mapping), 64L)
  ref_idx <- distributed_indices(length(bdev_cells), 24L)
  query_idx <- distributed_indices(length(ko_cells), 24L)
  imported_ref <- wb_matrix_read_features(
    bundle, "ATAC_bdev_mm10_peaks", mapping$bdev_feature_id[feature_idx],
    bdev_cells[ref_idx], missing = "NA"
  )
  source_ref <- bdev_counts[mapping$bdev_row[feature_idx], ref_idx, drop = FALSE]
  imported_query <- wb_matrix_read_features(
    bundle, "ATAC_bdev_mm10_peaks", mapping$bdev_feature_id[feature_idx],
    ko_cells[query_idx], missing = "NA"
  )
  source_query <- ko_counts[feature_idx, query_idx, drop = FALSE]
  stopifnot(identical(as.numeric(as.matrix(imported_ref)),
                      as.numeric(as.matrix(source_ref))),
            identical(as.numeric(as.matrix(imported_query)),
                      as.numeric(as.matrix(source_query))))

  manifest <- wb_manifest(bundle)
  persisted_matrix <- wb_matrix(bundle, "ATAC_bdev_mm10_peaks")
  creation <- list(
    status = "complete",
    bundle_id = manifest$bundle_id,
    bundle_dir = output_bundle_dir,
    bundle_state_fingerprint = manifest$current_state$fingerprint,
    parents = list(bdev_reference = expected_bdev_fingerprint,
                   scmulti_flt3ko_il7ko = expected_ko_fingerprint),
    matrix = list(
      id = "ATAC_bdev_mm10_peaks",
      dimensions = unname(dim(persisted_matrix)),
      input_storage_type = BPCells::matrix_type(common_counts),
      persisted_storage_type = BPCells::matrix_type(persisted_matrix),
      values_are_nonnegative_integers = TRUE
    ),
    reference_cells = length(bdev_cells),
    query_cells = length(ko_cells),
    exact_reference_pt_and_missingness = TRUE,
    exact_ordered_feature_identity = TRUE,
    distributed_source_count_checks = TRUE,
    parent_bundles_modified = FALSE
  )
  write_json(creation, "creation.json")
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
  log_event("Mapping bundle creation complete")
}

main()
