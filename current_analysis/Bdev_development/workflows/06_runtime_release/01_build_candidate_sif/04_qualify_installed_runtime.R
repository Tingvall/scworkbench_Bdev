#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--output")
  stop("usage: 04_qualify_installed_runtime.R --output <json>")
suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
  library(digest)
})

expected_versions <- c(scworkbench = "0.1.2", scwbStatistics = "0.1.0",
  scwbGenePrograms = "0.1.0", scwbRNA = "0.1.0",
  scwbRepresentation = "0.1.1", scwbTrajectory = "0.1.8",
  scwbLatent = "0.1.3", scwbReferenceMap = "0.1.0", scwbATAC = "0.1.1")
expected_manifests <- c(
  scwbStatistics = "14b4610595e255f199c58748e387981b294825a145a3ffffa2e0c5e309c175ad",
  scwbGenePrograms = "29e071596984f0e4930e072db58c76386a058fb40be55a70ffb8a541921feaee",
  scwbRNA = "067db0625181775721a91f2188bb9f2ca72443fc8aa06d0b453669035ebad410",
  scwbRepresentation = "c7716855de8613da874fc022874e150b7dd3e231a66e26095fa4a42ecf461a65",
  scwbTrajectory = "526a945e656cd7e58e99153bd3923c547fbdc040b12c1ccc09992c45b87ca352",
  scwbLatent = "5075f5b31a23b6427af97e255560b617d06b299b4ff8c826f489ee57c4984f23",
  scwbATAC = "3a61bc277357939c02a825246c960b61e2ee2a22f5a7ceff7cff7d67287443d3",
  scwbReferenceMap = "b190a85b0a743af769e6625b6d94fa44f32cee477479bbda19281ee096cb0747")
activation_order <- c("scwbStatistics", "scwbGenePrograms", "scwbRNA",
  "scwbRepresentation", "scwbTrajectory", "scwbLatent", "scwbATAC",
  "scwbReferenceMap")

observed_versions <- vapply(names(expected_versions), function(package)
  as.character(packageVersion(package)), character(1))
stopifnot(identical(observed_versions, expected_versions))
stopifnot(as.character(packageVersion("fgsea")) == "1.34.2",
  as.character(packageVersion("BiocParallel")) == "1.42.2")
observed_manifests <- vapply(activation_order, function(package) {
  path <- system.file("scworkbench", "operation-pack.json", package = package)
  stopifnot(nzchar(path))
  digest::digest(path, algo = "sha256", file = TRUE)
}, character(1))
stopifnot(identical(observed_manifests, expected_manifests[activation_order]))
stopifnot(nrow(wb_operation_packs()) == 0L)
for (package in activation_order) {
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  wb_operation_pack_activate(manifest)
}
packs <- wb_operation_packs()
operations <- wb_list_ops()
stopifnot(nrow(packs) == 8L, length(unique(packs$namespace)) == 8L,
  !anyDuplicated(operations), "scwbatac::hvp_select" %in% operations,
  any(grepl("^scwbreferencemap::", operations)))
report <- list(
  schema_version = "bdev.candidate_sif_installed_runtime.v1",
  status = "passed",
  package_versions = as.list(observed_versions),
  runtime_dependencies = list(fgsea = "1.34.2", BiocParallel = "1.42.2"),
  operation_pack_manifest_sha256 = as.list(observed_manifests),
  active_packs = split(packs, seq_len(nrow(packs))),
  operation_count = length(operations),
  operations = as.list(sort(operations)),
  no_operation_collisions = !anyDuplicated(operations),
  gate_passed = TRUE)
jsonlite::write_json(report, args[[2L]], pretty = TRUE, auto_unbox = TRUE,
  null = "null", na = "null")
