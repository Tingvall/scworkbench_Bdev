#!/usr/bin/env Rscript

# Usage: Rscript --vanilla 02_render_unfiltered_qc_report.R QC_DIR TEMPLATE OUTPUT_FILE
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 02_render_unfiltered_qc_report.R QC_DIR TEMPLATE OUTPUT_FILE")
qc_dir <- normalizePath(args[[1L]], mustWork = TRUE)
template <- normalizePath(args[[2L]], mustWork = TRUE)
output_file <- args[[3L]]
if (!file.exists(file.path(qc_dir, "completion.json"))) stop("QC calculations are incomplete: ", qc_dir)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
rmarkdown::render(
  input = template,
  output_format = "github_document",
  output_file = basename(output_file),
  output_dir = dirname(output_file),
  params = list(qc_dir = qc_dir),
  envir = new.env(parent = globalenv()),
  clean = TRUE,
  quiet = FALSE
)
