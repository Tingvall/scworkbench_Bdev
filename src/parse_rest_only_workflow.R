# Shared runtime, paths and reproducibility for the separate rest-only review.
options(stringsAsFactors = FALSE, warn = 1, width = 180)
suppressPackageStartupMessages({library(scworkbench); library(BPCells); library(jsonlite); library(yaml)})
source("src/bundle_config.R")
source("src/bundle_validation.R")
registry <- load_bundle_registry()
parent_cfg <- read_project_yaml("config/parse_cytokines.yml")
cfg <- read_project_yaml("config/parse_cytokines_rest_only.yml")
stopifnot(identical(as.character(packageVersion("scworkbench")), registry$bundles$parse_cytokines$runtime$scworkbench_version))
root <- cfg$storage_root
report <- file.path(root, "reports", cfg$run_id)
work <- file.path(root, "work", cfg$run_id)
discovery_dir <- file.path(work, "discovery")
session_root <- file.path(root, "sessions", cfg$run_id)
dir.create(report, recursive = TRUE, showWarnings = FALSE)
dir.create(work, recursive = TRUE, showWarnings = FALSE)
set.seed(cfg$analysis$random_seed)
log_event <- function(...) {cat(format(Sys.time(), tz = "UTC"), ..., "\n"); flush.console()}
write_json <- function(x, name, directory = report) jsonlite::write_json(x, file.path(directory, name),
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
write_tsv <- function(x, name, directory = report) write.table(x, file.path(directory, name),
  sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
read_tsv <- function(path) read.delim(path, check.names = FALSE, na.strings = c("NA", ""))
parent_check <- function() {
  x <- validate_registered_bundle(registry, "parse_cytokines", run_recipient_validation = FALSE)
  path <- registry$bundles$parse_cytokines$validation_report
  stopifnot(isTRUE(jsonlite::read_json(path, simplifyVector = TRUE)$ok))
  x
}
auc <- function(values, positive) {
  stopifnot(length(values) == length(positive), all(is.finite(values)), !anyNA(positive))
  n1 <- sum(positive); n0 <- sum(!positive)
  stopifnot(n1 > 0L, n0 > 0L)
  (sum(rank(values, ties.method = "average")[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
jaccard <- function(sets) {
  pairs <- combn(seq_along(sets), 2)
  mean(apply(pairs, 2, function(pair) length(intersect(sets[[pair[1]]], sets[[pair[2]]])) /
    length(union(sets[[pair[1]]], sets[[pair[2]]]))))
}
slug <- function(x) gsub("(^_+|_+$)", "", tolower(gsub("[^A-Za-z0-9]+", "_", x)))
save_provenance <- function(stage, inputs, extra = list()) {
  paths <- unique(normalizePath(c(inputs, "config/parse_cytokines_rest_only.yml", "config/parse_cytokines.yml",
    "config/bundles.yml", "src/parse_rest_only_workflow.R", "src/bundle_config.R", "src/bundle_validation.R"), mustWork = TRUE))
  write_tsv(data.frame(path = paths, sha256 = vapply(paths, wb_source_sha256, character(1))), paste0(stage, "_inputs.tsv"))
  write_json(c(list(status = "complete", stage = stage, source = parent_check(),
    bundle = cfg$bundle_dir, dataset_id = cfg$bundle_id, namespace = cfg$namespace,
    runtime = registry$bundles$parse_cytokines$runtime, parameters = cfg$analysis,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), promoted = FALSE), extra), paste0(stage, ".json"))
  capture.output(sessionInfo(), file = file.path(report, paste0(stage, "_session_info.txt")))
}
write_sets <- function(members, ids, stem) {
  table <- members[members$signature_id %in% ids, c("signature_id", "feature_id")]
  names(table)[1] <- "feature_set"
  table$set_name <- table$feature_set
  table$description <- "Parse target-versus-other-single-stimuli signature; no treatment-versus-control requirement; pending review."
  write_tsv(table[, c("feature_set", "set_name", "feature_id", "description")], paste0(stem, ".tsv"))
  lines <- vapply(ids, function(id) paste(c(id, "Parse rest-only, pending review", table$feature_id[table$feature_set == id]), collapse = "\t"), character(1))
  writeLines(lines, file.path(report, paste0(stem, ".gmt")))
}
