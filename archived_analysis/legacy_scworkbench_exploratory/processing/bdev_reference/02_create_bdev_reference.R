# Create the reference bundle through the installed scWorkbench public API.
# Usage: Rscript --vanilla processing/bdev_reference/02_create_bdev_reference.R SOURCE BUNDLE REPORTS
# Bind the source read-only and the bundle's parent directory writable.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 02_create_bdev_reference.R SOURCE BUNDLE REPORTS")
source_path <- normalizePath(args[[1L]], mustWork = TRUE)
bundle_dir <- args[[2L]]
report_dir <- args[[3L]]
if (file.exists(bundle_dir) || dir.exists(bundle_dir)) stop("Bundle destination already exists; refusing to overwrite")
if (basename(bundle_dir) != "bdev_reference") stop("Expected bundle name bdev_reference")
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1L]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "..", "..", "src", "reference_validation.R"))
options(width = 160, warn = 1)
suppressPackageStartupMessages({
  library(scworkbench)
  library(Seurat)
  library(Signac)
  library(jsonlite)
})
log_event <- function(...) { cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n"); flush.console() }
write_report <- function(value, name) {
  jsonlite::write_json(value, file.path(report_dir, name), pretty = TRUE,
                       auto_unbox = TRUE, null = "null", na = "null", digits = NA)
}

main <- function() {
  expected_sha <- "c880abb2cd41a80a6a1af30f8065a71baee12101cb9bab3d7118d9719020b8f9"
  log_event("Verifying original source checksum")
  source_sha <- wb_source_sha256(source_path)
  stopifnot(identical(source_sha, expected_sha))
  source_before <- file.info(source_path)[, c("size", "mtime", "ctime")]
  original_source <- wb_source_input(
    path = source_path, source_id = "bdev_rdata", accession = "Bdev.Rdata",
    role = "source_object", format = "rdata", scope = "external_file", sha256 = source_sha
  )
  source_identity <- list(
    dataset_name = "bdev_reference", source_type = "Seurat in RData workspace",
    source_object = "Bdev", extraction = "load(source, envir = isolated_environment); get('Bdev')",
    inputs = list(original_source)
  )
  wb_source_preflight(source_identity)
  log_event("Loading Bdev source through a read-only connection")
  workspace <- new.env(parent = baseenv())
  connection <- gzfile(source_path, "rb")
  loaded <- load(connection, envir = workspace)
  close(connection)
  stopifnot(identical(loaded, "Bdev"))
  Bdev <- get("Bdev", envir = workspace, inherits = FALSE)
  stopifnot(inherits(Bdev, "Seurat"), ncol(Bdev) == 26723L)
  cells <- colnames(Bdev)
  metadata_source <- Bdev[[]]
  stopifnot(identical(rownames(metadata_source), cells), ncol(metadata_source) == 68L)
  stopifnot(is.numeric(metadata_source$pt), sum(is.na(metadata_source$pt)) == 14081L)
  stopifnot(all(is.finite(metadata_source$pt[!is.na(metadata_source$pt)])))

  matrices <- list()
  mapping <- list()
  for (assay in SeuratObject::Assays(Bdev)) {
    feature_type <- if (assay %in% c("RNA", "SCT")) "gene" else if (assay %in% c("ATAC", "peaks")) "peak" else "feature"
    for (layer in SeuratObject::Layers(Bdev[[assay]], search = NA)) {
      mat <- SeuratObject::LayerData(Bdev[[assay]], layer = layer)
      stopifnot(identical(colnames(mat), cells), !anyDuplicated(rownames(mat)))
      matrix_id <- paste(assay, gsub("\\.", "_", layer), sep = "_")
      matrices[[matrix_id]] <- list(
        mat = mat, layer = layer, feature_type = feature_type,
        description = paste0("Source Bdev.Rdata object Bdev: assay ", assay, ", layer ", layer,
                             "; retained source values and feature universe; SHA-256 ", source_sha)
      )
      mapping[[matrix_id]] <- list(matrix_id = matrix_id, assay = assay, layer = layer,
                                  feature_type = feature_type, dimensions = dim(mat))
    }
  }
  stopifnot(length(matrices) == 9L)
  reductions <- SeuratObject::Reductions(Bdev)
  stopifnot(setequal(reductions, c("pca", "lsi", "umap.atac", "umap.rna", "umap.wnn")))
  spaces <- setNames(lapply(reductions, function(nm) {
    emb <- SeuratObject::Embeddings(Bdev[[nm]])
    stopifnot(identical(rownames(emb), cells), all(is.finite(emb)))
    emb
  }), reductions)
  metadata <- data.frame(cell_id = cells, metadata_source, check.names = FALSE, stringsAsFactors = FALSE)
  for (nm in names(metadata)) if (is.factor(metadata[[nm]])) metadata[[nm]] <- as.character(metadata[[nm]])
  source_identity$matrix_mapping <- mapping
  source_identity$reductions <- lapply(reductions, function(nm) list(
    name = nm, source_location = paste0("Bdev@reductions[[\"", nm, "\"]]@cell.embeddings"),
    columns = colnames(spaces[[nm]]), assay = SeuratObject::DefaultAssay(Bdev[[nm]])
  ))
  source_identity$metadata_columns <- names(metadata_source)
  source_identity$reference_pseudotime <- list(source_field = "pt", selected_by = "user", missing_policy = "preserve")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  write_report(source_identity, "source_identity.json")
  write_report(list(matrices = mapping, reductions = source_identity$reductions,
                    metadata_columns = names(metadata_source), n_cells = length(cells)), "import_plan.json")
  log_event("Creating bundle and importing", length(matrices), "matrices,", length(spaces), "spaces, and 68 metadata columns")
  b <- wb_init_bundle(bundle_dir, bundle_id = "bdev_reference", bundle_name = "bdev_reference", overwrite = FALSE)
  on.exit(try(wb_close_bundle(b), silent = TRUE), add = TRUE)
  source_tag <- paste0("Bdev.Rdata::Bdev;sha256=", source_sha)
  ingest <- wb_ingest_apply(
    b, set_name = "bdev_reference", cell_ids = cells, metadata = metadata,
    spaces = spaces, matrices = matrices, source = source_tag,
    source_identity = source_identity, created_by = "jtg", apply_mode = "atomic",
    emit_report_artifact = TRUE
  )
  write_report(attr(ingest, "ingest_report"), "ingest_report.json")
  wb_close_bundle(b)
  b <- wb_open_bundle(bundle_dir)
  log_event("Verifying source IDs, all metadata, all coordinates, and distributed matrix values")
  source_comparison <- validate_reference_import(b, Bdev, mapping, reductions)
  write_report(source_comparison, "source_comparison.json")
  source_records <- wb_prov_source_inputs(b)
  write_report(source_records, "source_provenance.json")
  stopifnot("bdev_rdata" %in% source_records$source_id)
  stopifnot(source_sha %in% source_records$checksum_value)
  log_event("Materializing the user-selected pt pseudotime through a registered operation")
  interface <- wb_op_interface("metadata_cell_vector", bundle = b)
  write_report(interface, "pseudotime_interface.json")
  result <- wb_session_run_op(
    wb_session_new(), b, "metadata_cell_vector",
    inputs = list(source_column = "pt", value_type = "continuous", missing_policy = "drop",
                  vector_name = "bdev_reference_pseudotime", cell_set = "all", min_cells = 1L),
    notes = "User selected source metadata pt as reference pseudotime. Retain its original direction and values. Missing values remain in bundle metadata; the typed vector covers defined pseudotime only. No fitted trajectory model is inferred."
  )
  promotion <- wb_session_promote_current(
    b, result$session, outputs = "vector", name_prefix = "bdev_reference_pseudotime", created_by = "jtg", fail_fast = TRUE
  )
  promoted <- promotion$promoted
  artifact_id <- as.character(promoted$artifact_id[promoted$output == "vector"])
  stopifnot(length(artifact_id) == 1L, nzchar(artifact_id))
  finish_reference_import(
    b, metadata_source, mapping, reductions, source_path, source_sha, source_before,
    report_dir, artifact_id, promotion = promotion, run_id = result$run$run_id
  )
}
main()
