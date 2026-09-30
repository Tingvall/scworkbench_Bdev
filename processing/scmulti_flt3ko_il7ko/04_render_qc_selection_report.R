#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 04_render_qc_selection_report.R SELECTION_DIR TEMPLATE OUTPUT_FILE")
selection_dir <- normalizePath(args[[1L]], mustWork = TRUE)
template <- normalizePath(args[[2L]], mustWork = TRUE)
output_file <- args[[3L]]
if (!file.exists(file.path(selection_dir, "completion.json"))) stop("Selection is incomplete")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
rmarkdown::render(
  template, output_format = "github_document", output_file = basename(output_file),
  output_dir = dirname(output_file), params = list(selection_dir = selection_dir),
  envir = new.env(parent = globalenv()), clean = TRUE, quiet = FALSE
)
