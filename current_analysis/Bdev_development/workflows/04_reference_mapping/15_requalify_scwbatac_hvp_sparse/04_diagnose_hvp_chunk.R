#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1]] != "--config") {
  stop("usage: 04_diagnose_hvp_chunk.R --config <yaml>")
}
suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbATAC)
  library(Matrix)
  library(yaml)
})

stage <- function(name, expression) {
  cat("STAGE_START", name, "\n")
  value <- withCallingHandlers(force(expression), warning = function(warning) {
    cat("STAGE_WARNING", name, conditionMessage(warning), "\n")
    cat("STAGE_WARNING_CALL", paste(deparse(conditionCall(warning)),
      collapse = " "), "\n")
    stop(warning)
  })
  cat("STAGE_PASS", name, "class=", paste(class(value), collapse = ","),
    "object_bytes=", as.numeric(object.size(value)), "\n")
  value
}

config <- yaml::read_yaml(normalizePath(args[[2]], mustWork = TRUE))
bundle <- wb_open_bundle(config$inputs$source_bundle$path,
  validate = FALSE, read_only = TRUE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
features <- as.character(wb_feature_ids(bundle,
  config$inputs$source_bundle$count_matrix))
cells <- head(as.character(wb_cells(bundle, "all")),
  as.integer(config$feature_selection$chunk_size))
values <- stage("slice", wb_slice(bundle,
  config$inputs$source_bundle$count_matrix,
  features = features, cells = cells))
binary <- stage("binarize", values > 0)
depth <- stage("column_sums", Matrix::colSums(binary))
detected <- stage("row_sums", Matrix::rowSums(binary))
log_depth <- stage("log_depth", log1p(as.numeric(depth)))
sparse_log_depth <- stage("sparse_log_depth",
  Matrix::Matrix(log_depth, ncol = 1L, sparse = TRUE))
product <- stage("logical_sparse_product", binary %*% sparse_log_depth)
numeric_product <- stage("numeric_product", as.numeric(product))
numeric_binary <- stage("numeric_binary", methods::as(binary, "dMatrix"))
safe_product <- stage("numeric_sparse_product",
  numeric_binary %*% sparse_log_depth)
safe_numeric_product <- stage("safe_numeric_product", as.numeric(safe_product))
stopifnot(isTRUE(all.equal(numeric_product, safe_numeric_product,
  tolerance = 0)))
cat("DIAGNOSTIC_COMPLETE\n")
