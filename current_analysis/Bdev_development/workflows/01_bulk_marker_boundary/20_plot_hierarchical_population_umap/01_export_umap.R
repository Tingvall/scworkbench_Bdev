#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config")) {
  stop("Usage: 01_export_umap.R --config <config.yaml>")
}
suppressPackageStartupMessages({
  library(scworkbench)
  library(yaml)
  library(jsonlite)
  library(digest)
})

config_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(config_path)
root <- cfg$storage$results_root
tables <- file.path(root, "01_tables")
manifests <- file.path(root, "06_manifests")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(manifests, recursive = TRUE, showWarnings = FALSE)
stopifnot(digest::digest(file = file.path(cfg$bundle$path, "manifest.json"),
  algo = "sha256", serialize = FALSE) == cfg$bundle$manifest_sha256)

decisions <- read.delim(cfg$inputs$decisions$path, check.names = FALSE)
cells <- unique(as.character(decisions$cell_id))
stopifnot(length(cells) == 1746L)
bundle <- wb_open_bundle(cfg$bundle$path, validate = FALSE, read_only = TRUE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$current_state$fingerprint, cfg$bundle$fingerprint))
coordinates <- wb_space(bundle, cfg$bundle$space_id, cells = cells)
stopifnot(nrow(coordinates) == length(cells), identical(as.character(coordinates$cell_id), cells))
names(coordinates)[names(coordinates) == "x"] <- "UMAP_1"
names(coordinates)[names(coordinates) == "y"] <- "UMAP_2"
coordinate_connection <- gzfile(file.path(tables, "01_bdev_umap_atac_coordinates.tsv.gz"), "wt")
write.table(coordinates, coordinate_connection, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
close(coordinate_connection)
receipt <- list(
  schema_version = "bdev.hierarchical_population_umap_export.v1",
  status = "complete_read_only",
  run_identity = cfg$run_identity,
  bundle = cfg$bundle$path,
  bundle_fingerprint = manifest$current_state$fingerprint,
  space_id = cfg$bundle$space_id,
  cells = nrow(coordinates),
  bundle_mutated = FALSE)
jsonlite::write_json(receipt, file.path(manifests, "01_coordinate_export_receipt.json"),
  pretty = TRUE, auto_unbox = TRUE)
