#!/usr/bin/env Rscript

# Create the unfiltered scmulti_flt3ko_il7ko source bundle from the aggregated
# Cell Ranger ARC matrix. All cells are retained; QC is represented by metadata
# and a separate accepted analysis-universe cell set.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "Usage: 05_create_source_bundle.R OUTS_DIR QC_SELECTION_TSV BUNDLE_DIR REPORT_DIR INPUT_INVENTORY",
    call. = FALSE
  )
}

outs_dir <- normalizePath(args[[1L]], mustWork = TRUE)
qc_selection_path <- normalizePath(args[[2L]], mustWork = TRUE)
bundle_dir <- normalizePath(args[[3L]], mustWork = FALSE)
report_dir <- normalizePath(args[[4L]], mustWork = FALSE)
inventory_path <- normalizePath(args[[5L]], mustWork = TRUE)

if (file.exists(bundle_dir) || dir.exists(bundle_dir)) {
  stop("Bundle destination already exists; refusing to overwrite: ", bundle_dir, call. = FALSE)
}
if (file.exists(report_dir) || dir.exists(report_dir)) {
  stop("Report destination already exists; refusing to overwrite: ", report_dir, call. = FALSE)
}

options(width = 160, warn = 1)
suppressPackageStartupMessages({
  library(scworkbench)
  library(BPCells)
  library(rhdf5)
  library(data.table)
  library(jsonlite)
  library(Matrix)
})

created_by <- Sys.getenv("SCWORKBENCH_CREATED_BY", unset = "jtg")
matrix_h5 <- file.path(outs_dir, "filtered_feature_bc_matrix.h5")
fragments_path <- file.path(outs_dir, "atac_fragments.tsv.gz")
fragments_index_path <- paste0(fragments_path, ".tbi")
peaks_path <- file.path(outs_dir, "atac_peaks.bed")
aggr_path <- file.path(outs_dir, "aggr.csv")
reference_json <- "/home/jtg/resources/refs/refdata-cellranger-arc-GRCm39-2024-A/reference.json"
downstream_ids_path <- file.path(dirname(qc_selection_path), "downstream_cell_ids.txt")
qc_validation_path <- file.path(dirname(qc_selection_path), "validation.json")

required_paths <- c(
  matrix_h5, fragments_path, fragments_index_path, peaks_path, aggr_path,
  reference_json, qc_selection_path, downstream_ids_path, qc_validation_path,
  inventory_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) stop("Missing required inputs: ", paste(missing_paths, collapse = ", "), call. = FALSE)

log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}

write_json_atomic <- function(value, path) {
  tmp <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  jsonlite::write_json(value, tmp, pretty = TRUE, auto_unbox = TRUE,
                       null = "null", na = "null", digits = NA)
  if (!file.rename(tmp, path)) stop("Could not publish ", path, call. = FALSE)
  invisible(path)
}

write_tsv_atomic <- function(value, path) {
  compression_suffix <- if (endsWith(path, ".gz")) ".gz" else ""
  tmp <- file.path(
    dirname(path),
    paste0(".", basename(path), ".tmp-", Sys.getpid(), compression_suffix)
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  data.table::fwrite(value, tmp, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(tmp, path)) stop("Could not publish ", path, call. = FALSE)
  invisible(path)
}

distributed_indices <- function(n, target) {
  unique(as.integer(round(seq.int(1L, n, length.out = min(as.integer(target), n)))))
}

values_equal <- function(source, imported) {
  if (is.logical(source)) source <- as.integer(source)
  if (is.logical(imported)) imported <- as.integer(imported)
  if (is.numeric(source) && is.numeric(imported)) {
    return(identical(as.numeric(source), as.numeric(imported)))
  }
  identical(as.character(source), as.character(imported))
}

source_input_from_inventory <- function(inventory, role, source_id, format) {
  row <- inventory[inventory$role == role, , drop = FALSE]
  if (nrow(row) != 1L) stop("Inventory role must occur exactly once: ", role, call. = FALSE)
  wb_source_input(
    path = row$path[[1L]], source_id = source_id,
    accession = basename(row$path[[1L]]), role = role, format = format,
    scope = "external_file", sha256 = row$sha256[[1L]]
  )
}

main <- function() {
  inventory <- data.table::fread(inventory_path, data.table = FALSE)
  stopifnot(all(c("role", "path", "sha256") %in% names(inventory)))

  source_inputs <- list(
    source_input_from_inventory(inventory, "matrix_h5", "arc_filtered_feature_bc_matrix", "10x_hdf5"),
    source_input_from_inventory(inventory, "fragments", "arc_atac_fragments", "tsv_gz"),
    source_input_from_inventory(inventory, "fragments_index", "arc_atac_fragments_index", "tbi"),
    source_input_from_inventory(inventory, "peaks", "arc_atac_peaks", "bed"),
    source_input_from_inventory(inventory, "aggr_csv", "arc_aggregation_manifest", "csv"),
    source_input_from_inventory(inventory, "reference_json", "arc_reference_manifest", "json"),
    wb_source_input(
      path = qc_selection_path, source_id = "qc_cell_selection", accession = basename(qc_selection_path),
      role = "cell_metadata", format = "tsv_gz", scope = "external_file"
    ),
    wb_source_input(
      path = downstream_ids_path, source_id = "qc_downstream_cell_ids", accession = basename(downstream_ids_path),
      role = "cell_selection", format = "txt", scope = "external_file"
    ),
    wb_source_input(
      path = qc_validation_path, source_id = "qc_selection_validation", accession = basename(qc_validation_path),
      role = "validation", format = "json", scope = "external_file"
    ),
    wb_source_input(
      path = inventory_path, source_id = "arc_input_inventory", accession = basename(inventory_path),
      role = "input_inventory", format = "tsv", scope = "external_file"
    )
  )

  source_identity <- list(
    dataset_name = "scmulti_flt3lko_il7ko",
    technical_bundle_id = "scmulti_flt3ko_il7ko",
    biological_label = "FLT3LKO_IL7KO_BM",
    source_type = "aggregated Cell Ranger ARC 2.2.0 output",
    external_preprocessing_boundary = paste(
      "Cell Ranger ARC mixed RNA/ATAC HDF5 opened with BPCells; feature types split explicitly;",
      "QC metrics and user-approved selection were computed before bundle creation."
    ),
    genome = list(organism = "Mus musculus", build = "GRCm39", reference = "refdata-cellranger-arc-GRCm39-2024-A"),
    inputs = source_inputs,
    qc_policy = list(
      physical_filtering_performed = FALSE,
      all_cells_retained = TRUE,
      downstream_set = "set:qc_use_downstream",
      all_must_pass = TRUE,
      strict_inequalities = TRUE,
      thresholds = list(
        rna_umi = "> 1000", detected_genes = "> 500", mitochondrial_percent = "< 20",
        atac_fragments = "> 1000", tss_enrichment = "> 2", nucleosome_signal = "< 2"
      )
    )
  )

  log_event("Running strict source preflight")
  wb_source_preflight(source_identity)

  log_event("Opening Gene Expression and Peaks matrices as disk-backed BPCells matrices")
  rna <- BPCells::open_matrix_10x_hdf5(matrix_h5, feature_type = "Gene Expression")
  atac <- BPCells::open_matrix_10x_hdf5(matrix_h5, feature_type = "Peaks")
  h5_feature_type <- rhdf5::h5read(matrix_h5, "matrix/features/feature_type")
  h5_feature_id <- rhdf5::h5read(matrix_h5, "matrix/features/id")
  h5_feature_name <- rhdf5::h5read(matrix_h5, "matrix/features/name")

  gene_mask <- h5_feature_type == "Gene Expression"
  peak_mask <- h5_feature_type == "Peaks"
  gene_symbols <- h5_feature_name[gene_mask]
  bundle_gene_names <- make.unique(gene_symbols)
  bundle_peak_names <- h5_feature_name[peak_mask]
  rownames(rna) <- bundle_gene_names
  rownames(atac) <- bundle_peak_names

  cells <- colnames(rna)
  stopifnot(
    identical(colnames(atac), cells), length(cells) == 28541L,
    nrow(rna) == 33696L, nrow(atac) == 200063L,
    !anyDuplicated(cells), !anyDuplicated(rownames(rna)), !anyDuplicated(rownames(atac))
  )

  qc_connection <- gzfile(qc_selection_path, open = "rt")
  qc <- utils::read.delim(
    qc_connection, header = TRUE, sep = "\t", quote = "\"",
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA"
  )
  close(qc_connection)
  stopifnot("cell_id" %in% names(qc), !anyDuplicated(qc$cell_id), setequal(qc$cell_id, cells))
  qc <- qc[match(cells, qc$cell_id), , drop = FALSE]
  stopifnot(identical(qc$cell_id, cells), !anyNA(qc$qc_use_downstream))
  if (anyNA(qc$qc_failure_reasons)) qc$qc_failure_reasons[is.na(qc$qc_failure_reasons)] <- ""
  downstream_cells <- cells[as.logical(qc$qc_use_downstream)]
  expected_downstream <- readLines(downstream_ids_path, warn = FALSE)
  stopifnot(identical(downstream_cells, expected_downstream), length(downstream_cells) == 23252L)

  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(report_dir)) stop("Could not create report directory: ", report_dir, call. = FALSE)

  gene_mapping <- data.frame(
    h5_row = which(gene_mask), ensembl_id = h5_feature_id[gene_mask],
    source_symbol = gene_symbols, bundle_feature_id = bundle_gene_names,
    symbol_was_disambiguated = bundle_gene_names != gene_symbols,
    stringsAsFactors = FALSE
  )
  peak_mapping <- data.frame(
    h5_row = which(peak_mask), source_feature_id = h5_feature_id[peak_mask],
    bundle_feature_id = bundle_peak_names, stringsAsFactors = FALSE
  )
  write_tsv_atomic(gene_mapping, file.path(report_dir, "rna_feature_mapping.tsv.gz"))
  write_tsv_atomic(peak_mapping, file.path(report_dir, "atac_feature_mapping.tsv.gz"))

  source_identity$matrix_mapping <- list(
    RNA_counts = list(
      source_feature_type = "Gene Expression", layer = "counts", feature_type = "gene",
      dimensions = dim(rna), row_ids = "make.unique(Cell Ranger feature name)",
      original_gene_ids = "rna_feature_mapping.tsv.gz"
    ),
    ATAC_aggregate_peaks = list(
      source_feature_type = "Peaks", layer = "counts", feature_type = "peak",
      dimensions = dim(atac), row_ids = "Cell Ranger peak feature name",
      original_peak_ids = "atac_feature_mapping.tsv.gz"
    )
  )
  source_identity$metadata_columns <- setdiff(names(qc), "cell_id")
  source_identity$cell_universe <- list(
    retained_cells = length(cells), downstream_cells = length(downstream_cells),
    retained_not_selected = length(cells) - length(downstream_cells), samples = length(unique(qc$library_id))
  )
  write_json_atomic(source_identity, file.path(report_dir, "source_identity.json"))
  write_json_atomic(list(
    bundle_id = "scmulti_flt3ko_il7ko", imported_cell_set = "set:scmulti_flt3ko_il7ko",
    matrices = source_identity$matrix_mapping, metadata_columns = source_identity$metadata_columns,
    retained_cells = length(cells), downstream_cell_set = "set:qc_use_downstream",
    downstream_cells = length(downstream_cells), physical_filtering_performed = FALSE
  ), file.path(report_dir, "import_plan.json"))

  matrices <- list(
    RNA_counts = list(
      mat = rna, layer = "counts", feature_type = "gene",
      description = paste(
        "Raw Gene Expression counts from aggregated Cell Ranger ARC HDF5; all cells and genes retained;",
        "gene symbols made unique with base::make.unique; original Ensembl IDs recorded in the import report."
      )
    ),
    ATAC_aggregate_peaks = list(
      mat = atac, layer = "counts", feature_type = "peak",
      description = paste(
        "Raw Peaks counts from aggregated Cell Ranger ARC HDF5; all cells and the Cell Ranger aggregate",
        "peak universe retained; GRCm39 coordinates."
      )
    )
  )

  # The public atomic ingest adapter accepts Matrix payloads, not BPCells S4
  # wrappers. Materialize exact sparse matrices only at this governed handoff.
  log_event("Materializing exact sparse Matrix payloads for governed ingestion")
  rna_sparse <- methods::as(rna, "dgCMatrix")
  atac_sparse <- methods::as(atac, "dgCMatrix")
  stopifnot(identical(dimnames(rna_sparse), dimnames(rna)))
  stopifnot(identical(dimnames(atac_sparse), dimnames(atac)))
  matrices$RNA_counts$mat <- rna_sparse
  matrices$ATAC_aggregate_peaks$mat <- atac_sparse
  rm(rna_sparse, atac_sparse)
  invisible(gc())

  log_event("Creating bundle and ingesting", length(cells), "cells with raw RNA and ATAC matrices")
  bundle <- wb_init_bundle(
    bundle_dir, bundle_id = "scmulti_flt3ko_il7ko",
    bundle_name = "scmulti_flt3ko_il7ko", overwrite = FALSE
  )
  on.exit(if (!is.null(bundle)) try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  ingest <- wb_ingest_apply(
    bundle, set_name = "scmulti_flt3ko_il7ko", cell_ids = cells,
    metadata = qc, matrices = matrices,
    source = "Cell Ranger ARC aggregate plus user-approved QC annotations",
    source_identity = source_identity, created_by = created_by,
    apply_mode = "atomic", emit_report_artifact = TRUE
  )
  write_json_atomic(attr(ingest, "ingest_report"), file.path(report_dir, "ingest_report.json"))

  wb_matrix_set_genome(bundle, "RNA_counts", "Mus musculus", "GRCm39",
                       "10x Genomics refdata-cellranger-arc-GRCm39-2024-A")
  wb_matrix_set_genome(bundle, "ATAC_aggregate_peaks", "Mus musculus", "GRCm39",
                       "10x Genomics refdata-cellranger-arc-GRCm39-2024-A")

  log_event("Deriving accepted downstream QC analysis universe without removing cells")
  set_result <- wb_cell_set_derive(
    bundle, cell_set = "set:qc_use_downstream", cells = downstream_cells,
    role = "analysis_universe", status = "accepted",
    source = "user_approved_multiome_qc_thresholds", created_by = created_by,
    description = paste(
      "Cells passing all strict thresholds: RNA UMI >1000; genes >500; mitochondrial percent <20;",
      "ATAC fragments >1000; TSS enrichment >2; nucleosome signal <2. All source cells remain retained."
    ),
    overwrite = FALSE, promote_manifest = TRUE
  )
  write_json_atomic(set_result$promoted, file.path(report_dir, "downstream_cell_set_promotion.json"))

  wb_close_bundle(bundle)
  bundle <- NULL
  bundle <- wb_open_bundle(bundle_dir)

  log_event("Validating identifiers, metadata, cell sets, and distributed raw matrix values")
  stopifnot(
    identical(wb_cells(bundle, cell_set = "all"), cells),
    identical(wb_cells(bundle, cell_set = "set:scmulti_flt3ko_il7ko"), cells),
    identical(wb_cells(bundle, cell_set = "set:qc_use_downstream"), downstream_cells)
  )
  matrix_checks <- list()
  for (matrix_id in names(matrices)) {
    source_matrix <- matrices[[matrix_id]]$mat
    imported <- wb_matrix(bundle, matrix_id)
    stopifnot(
      identical(dim(imported), dim(source_matrix)),
      identical(rownames(imported), rownames(source_matrix)),
      identical(colnames(imported), colnames(source_matrix))
    )
    feature_idx <- distributed_indices(nrow(source_matrix), 64L)
    cell_idx <- distributed_indices(ncol(source_matrix), 32L)
    by_feature <- as.matrix(wb_matrix_read_features(
      bundle, matrix_id, rownames(source_matrix)[feature_idx], cells, missing = "NA"
    ))
    expected_by_feature <- as.matrix(source_matrix[feature_idx, , drop = FALSE])
    by_cell <- as.matrix(wb_matrix_read_features(
      bundle, matrix_id, rownames(source_matrix), cells[cell_idx], missing = "NA"
    ))
    expected_by_cell <- as.matrix(source_matrix[, cell_idx, drop = FALSE])
    stopifnot(
      identical(as.numeric(by_feature), as.numeric(expected_by_feature)),
      identical(as.numeric(by_cell), as.numeric(expected_by_cell))
    )
    matrix_checks[[matrix_id]] <- list(
      dimensions = dim(source_matrix), ordered_feature_ids_equal = TRUE,
      ordered_cell_ids_equal = TRUE, sampled_features_across_all_cells = length(feature_idx),
      sampled_cells_across_all_features = length(cell_idx), sampled_values_exact = TRUE
    )
  }

  imported_metadata <- wb_meta(bundle, cells = cells)
  metadata_checks <- vapply(setdiff(names(qc), "cell_id"), function(column) {
    column %in% names(imported_metadata) && values_equal(qc[[column]], imported_metadata[[column]])
  }, logical(1L))
  stopifnot(all(metadata_checks))
  recipe <- wb_get_cell_set_recipe(bundle, "set:qc_use_downstream", parse = FALSE)
  stopifnot(
    nrow(recipe) == 1L, identical(as.character(recipe$role[[1L]]), "analysis_universe"),
    identical(as.character(recipe$status[[1L]]), "accepted")
  )
  write_tsv_atomic(recipe, file.path(report_dir, "downstream_cell_set_recipe.tsv"))

  source_comparison <- list(
    retained_cells = length(cells), downstream_cells = length(downstream_cells),
    retained_not_selected = length(cells) - length(downstream_cells),
    physical_filtering_performed = FALSE, exact_cell_partition = TRUE,
    matrices = matrix_checks, metadata_columns_exact = as.list(metadata_checks)
  )
  write_json_atomic(source_comparison, file.path(report_dir, "source_comparison.json"))
  write_json_atomic(wb_prov_source_inputs(bundle), file.path(report_dir, "source_provenance.json"))
  write_json_atomic(list(
    RNA_counts = wb_entity_certificate(bundle, "matrix", "RNA_counts"),
    ATAC_aggregate_peaks = wb_entity_certificate(bundle, "matrix", "ATAC_aggregate_peaks"),
    qc_use_downstream = wb_entity_certificate(bundle, "cell_set", "set:qc_use_downstream")
  ), file.path(report_dir, "entity_certificates.json"))
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))

  log_event("Finalizing bundle for recipient validation")
  readiness <- wb_finalize_recipient_bundle(bundle)
  stopifnot(isTRUE(readiness$ready))
  write_json_atomic(readiness, file.path(report_dir, "recipient_readiness.json"))
  wb_close_bundle(bundle)
  bundle <- NULL

  validation <- wb_bundle_validation_report(bundle_dir, level = "recipient", progress = TRUE)
  write_json_atomic(validation, file.path(report_dir, "validation.json"))
  stopifnot(isTRUE(validation$ok))

  recipient <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = TRUE)
  manifest <- wb_manifest(recipient)
  stopifnot(
    identical(wb_cells(recipient, "all"), cells),
    identical(wb_cells(recipient, "set:qc_use_downstream"), downstream_cells)
  )
  wb_close_bundle(recipient)
  write_json_atomic(manifest, file.path(report_dir, "manifest.json"))

  completion <- list(
    status = "complete", completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    bundle_id = "scmulti_flt3ko_il7ko", biological_dataset = "scmulti_flt3lko_il7ko",
    bundle_dir = bundle_dir, retained_cells = length(cells), samples = length(unique(qc$library_id)),
    downstream_cell_set = "set:qc_use_downstream", downstream_cells = length(downstream_cells),
    retained_not_selected = length(cells) - length(downstream_cells),
    physical_filtering_performed = FALSE,
    matrices = lapply(matrix_checks, function(x) x$dimensions),
    metadata_columns = setdiff(names(qc), "cell_id"),
    bundle_state_fingerprint = manifest$current_state$fingerprint,
    recipient_validation_fingerprint = validation$post_validation_fingerprint,
    validation_status = validation$status
  )
  write_json_atomic(completion, file.path(report_dir, "completion.json"))
  log_event("Complete:", bundle_dir)
}

main()
