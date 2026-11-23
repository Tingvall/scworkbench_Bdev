#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 03_export_bdev_peak_axis.R --config <config.yaml>")
suppressPackageStartupMessages({
  library(scworkbench)
  library(yaml)
  library(digest)
  library(jsonlite)
})

merge_config <- function(path) {
  current <- yaml::read_yaml(path)
  base_name <- current$base_config
  current$base_config <- NULL
  if (is.null(base_name)) return(current)
  utils::modifyList(merge_config(normalizePath(file.path(dirname(path), base_name))),
    current, keep.null = TRUE)
}
cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- merge_config(cfg_path)
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
source_manifest <- file.path(cfg$inputs$source_bundle$path, "manifest.json")
stopifnot(identical(sha256(source_manifest),
  cfg$inputs$source_bundle$manifest_sha256))

bundle <- wb_open_bundle(cfg$inputs$source_bundle$path, validate = FALSE,
  read_only = TRUE)
features <- wb_feature_ids(bundle, cfg$inputs$count_matrix)
wb_close_bundle(bundle)
stopifnot(length(features) == cfg$inputs$expected_features,
  !anyDuplicated(features))
matches <- regexec("^(.*)-([0-9]+)-([0-9]+)$", features)
parts <- regmatches(features, matches)
stopifnot(all(lengths(parts) == 4L))
peaks <- data.frame(
  chrom = vapply(parts, `[[`, character(1), 2L),
  start = as.integer(vapply(parts, `[[`, character(1), 3L)),
  end = as.integer(vapply(parts, `[[`, character(1), 4L)),
  feature_id = features,
  feature_order = seq_along(features), stringsAsFactors = FALSE)
stopifnot(all(peaks$start >= 0L), all(peaks$end > peaks$start))

tables <- file.path(cfg$storage$results_root, "01_tables")
manifests <- file.path(cfg$storage$results_root, "06_manifests")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(manifests, recursive = TRUE, showWarnings = FALSE)
bed <- file.path(tables, "02_bdev_exact_peak_axis.bed")
write.table(peaks[, c("chrom", "start", "end", "feature_id")], bed,
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
receipt <- list(schema_version = "bdev.exact_peak_axis_export.v1",
  status = "complete", run_identity = cfg$run_identity,
  source_bundle = cfg$inputs$source_bundle$path,
  source_bundle_manifest_sha256 = sha256(source_manifest),
  matrix_id = cfg$inputs$count_matrix, feature_count = length(features),
  peak_bed = bed, peak_bed_sha256 = sha256(bed),
  first_feature = features[[1L]], last_feature = tail(features, 1L),
  archive_modified = FALSE)
jsonlite::write_json(receipt, file.path(manifests, "03_peak_axis_export.json"),
  pretty = TRUE, auto_unbox = TRUE)
cat(jsonlite::toJSON(receipt, pretty = TRUE, auto_unbox = TRUE), "\n")
