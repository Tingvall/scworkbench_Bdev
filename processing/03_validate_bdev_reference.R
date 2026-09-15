# Validate/finalize an existing reference without reimporting or refitting it.
# Usage: Rscript --vanilla processing/03_validate_bdev_reference.R SOURCE BUNDLE REPORTS
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 03_validate_bdev_reference.R SOURCE BUNDLE REPORTS")
source_path <- normalizePath(args[[1L]], mustWork = TRUE)
bundle_dir <- normalizePath(args[[2L]], mustWork = TRUE)
report_dir <- normalizePath(args[[3L]], mustWork = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1L]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "..", "src", "reference_validation.R"))
suppressPackageStartupMessages({
  library(scworkbench)
  library(Seurat)
  library(Signac)
  library(jsonlite)
})

main <- function() {
  expected_sha <- "c880abb2cd41a80a6a1af30f8065a71baee12101cb9bab3d7118d9719020b8f9"
  cat(format(Sys.time()), "Checking original source and existing bundle\n")
  source_sha <- wb_source_sha256(source_path)
  stopifnot(identical(source_sha, expected_sha))
  source_before <- file.info(source_path)[, c("size", "mtime", "ctime")]
  workspace <- new.env(parent = baseenv())
  connection <- gzfile(source_path, "rb")
  loaded <- load(connection, envir = workspace)
  close(connection)
  stopifnot(identical(loaded, "Bdev"))
  Bdev <- get("Bdev", envir = workspace, inherits = FALSE)
  plan <- jsonlite::read_json(file.path(report_dir, "import_plan.json"), simplifyVector = FALSE)
  mapping <- plan$matrices
  reductions <- SeuratObject::Reductions(Bdev)
  b <- wb_open_bundle(bundle_dir)
  on.exit(try(wb_close_bundle(b), silent = TRUE), add = TRUE)
  stopifnot(identical(wb_manifest(b)$bundle_id, "bdev_reference"))
  source_records <- wb_prov_source_inputs(b)
  stopifnot("bdev_rdata" %in% source_records$source_id,
            source_sha %in% source_records$checksum_value)
  cat(format(Sys.time()), "Verifying source values and retained pseudotime\n")
  comparison <- validate_reference_import(b, Bdev, mapping, reductions)
  jsonlite::write_json(comparison, file.path(report_dir, "source_comparison.json"),
                       pretty = TRUE, auto_unbox = TRUE, digits = NA)
  vectors <- wb_list_cell_vectors(b)
  vectors <- vectors[vectors$op_name == "metadata_cell_vector" &
                     grepl("bdev_reference_pseudotime", vectors$name, fixed = TRUE), , drop = FALSE]
  stopifnot(nrow(vectors) == 1L)
  finish_reference_import(
    b, Bdev[[]], mapping, reductions, source_path, source_sha, source_before,
    report_dir, vectors$artifact_id[[1L]]
  )
}
main()
