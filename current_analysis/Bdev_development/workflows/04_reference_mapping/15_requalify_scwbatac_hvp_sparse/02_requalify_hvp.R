#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1]] != "--config") {
  stop("usage: 02_requalify_hvp.R --config <yaml>")
}
suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbATAC)
  library(yaml)
  library(jsonlite)
  library(digest)
})

sha <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
scientific_state <- function(bundle) {
  as.character(scworkbench:::.wb_session_source_identity(bundle)$scientific_state_fingerprint)
}
without_publication_time <- function(value) {
  if (!is.null(value$publication)) value$publication$published_at <- NULL
  value
}
write_tsv <- function(value, path) {
  connection <- if (grepl("[.]gz$", path)) gzfile(path, "wt") else path
  on.exit(if (inherits(connection, "connection")) close(connection), add = TRUE)
  write.table(value, connection, sep = "\t", quote = TRUE,
    row.names = FALSE, na = "NA")
}
write_json <- function(value, path) jsonlite::write_json(value, path,
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
max_difference <- function(left, right) {
  difference <- abs(as.numeric(left) - as.numeric(right))
  if (!length(difference) || all(is.na(difference))) return(NA_real_)
  max(difference, na.rm = TRUE)
}

config_path <- normalizePath(args[[2]], mustWork = TRUE)
config <- yaml::read_yaml(config_path)
stopifnot(config$scope$stage == "Bdev_development",
  config$scope$external_query_use == "prohibited",
  config$scope$knockout_query_use == "prohibited")
results <- config$storage$results_root
directories <- setNames(file.path(results, c("01_tables", "02_figures",
  "03_reports", "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(c(directories, sessions = config$storage$sessions),
  dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(directories[["logs"]], "01_hvp_sparse_requalification.log"),
  split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(sha(config$runtime$sif) == config$runtime$sif_sha256,
  as.character(packageVersion("scworkbench")) ==
    config$runtime$scworkbench$version,
  as.character(packageVersion("scwbATAC")) ==
    config$runtime$scwbATAC$version)
manifest <- system.file("scworkbench", "operation-pack.json",
  package = "scwbATAC")
stopifnot(sha(manifest) == config$runtime$scwbATAC$manifest_sha256)
wb_operation_pack_activate(manifest, lib.loc = .libPaths())
stopifnot(config$feature_selection$operation %in% wb_list_ops())

source_manifest <- file.path(config$inputs$source_bundle$path, "manifest.json")
source_manifest_before <- sha(source_manifest)
stopifnot(source_manifest_before ==
  config$inputs$source_bundle$manifest_sha256)
source_manifest_document_before <- jsonlite::read_json(source_manifest,
  simplifyVector = FALSE)
original_source <- config$inputs$source_bundle$original_source
original_manifest <- if (is.null(original_source)) NULL else
  file.path(original_source$path, "manifest.json")
original_manifest_before <- if (is.null(original_manifest)) NA_character_ else
  sha(original_manifest)
if (!is.null(original_source)) stopifnot(
  original_manifest_before == original_source$manifest_sha256)
bundle <- wb_open_bundle(config$inputs$source_bundle$path,
  validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
scientific_state_before <- scientific_state(bundle)
stopifnot(scientific_state_before ==
  config$inputs$source_bundle$scientific_state_fingerprint)
cells <- as.character(wb_cells(bundle, "all"))
features <- as.character(wb_feature_ids(bundle,
  config$inputs$source_bundle$count_matrix))
stopifnot(length(cells) == config$inputs$source_bundle$expected_cells,
  length(features) == config$inputs$source_bundle$expected_features)

new_ledgers <- list()
new_summaries <- list()
elapsed <- list()
session_evidence <- list()
for (mode in names(config$feature_selection$modes)) {
  specification <- config$feature_selection$modes[[mode]]
  log_event("Starting official", config$feature_selection$operation, mode)
  session <- wb_session_create(bundle, config$storage$sessions,
    paste0(config$run_identity, "_", mode))
  started <- proc.time()[["elapsed"]]
  result <- wb_session_run_op(session, bundle,
    config$feature_selection$operation, inputs = list(
      matrix = config$inputs$source_bundle$count_matrix,
      n_features = as.integer(config$feature_selection$n_features),
      depth_adjustment = isTRUE(specification$depth_adjustment),
      min_detection_pct =
        as.numeric(config$feature_selection$min_detection_pct),
      max_detection_pct =
        as.numeric(config$feature_selection$max_detection_pct),
      chunk_size = as.integer(config$feature_selection$chunk_size),
      cell_set = "all", min_cells = 2L),
    notes = paste("Sparse-safe full-Bdev HVP requalification", mode))
  elapsed[[mode]] <- proc.time()[["elapsed"]] - started
  ledger <- result$run$outputs$selected_features$table
  summary <- result$run$outputs$fit_summary$table
  stopifnot(nrow(ledger) == length(features),
    sum(ledger$selected) == config$feature_selection$n_features,
    all(ledger$method == specification$expected_method))
  new_ledgers[[mode]] <- ledger
  new_summaries[[mode]] <- summary
  write_tsv(ledger, file.path(directories[["tables"]],
    paste0(if (mode == "hvp_depth20k") "01" else "03", "_", mode,
      "_scores.tsv.gz")))
  write_tsv(summary, file.path(directories[["tables"]],
    paste0(if (mode == "hvp_depth20k") "02" else "04", "_", mode,
      "_summary.tsv")))
  session_evidence[[mode]] <- list(
    durable_workspace = result$session$durable$workspace,
    output_names = names(result$run$outputs))
  log_event("Completed", mode, "elapsed_seconds =", elapsed[[mode]])
  gc()
}

prior <- list(
  hvp_depth20k = read.delim(gzfile(
    config$inputs$prior_run$adjusted_ledger), check.names = FALSE),
  hvp_raw20k = read.delim(gzfile(
    config$inputs$prior_run$raw_ledger), check.names = FALSE))
numeric_columns <- c("mean_accessibility", "detected_cells", "detection_pct",
  "raw_variance", "residual_variance", "depth_slope",
  "depth_r_squared", "score")
comparisons <- list()
for (mode in names(new_ledgers)) {
  current <- new_ledgers[[mode]]
  earlier <- prior[[mode]]
  stopifnot(setequal(current$feature_id, earlier$feature_id))
  earlier <- earlier[match(current$feature_id, earlier$feature_id), , drop = FALSE]
  selected_current <- current$feature_id[current$selected]
  selected_earlier <- earlier$feature_id[earlier$selected]
  differences <- setNames(lapply(numeric_columns, function(column) {
    max_difference(current[[column]], earlier[[column]])
  }), numeric_columns)
  comparisons[[mode]] <- list(
    complete_feature_axis_identical = identical(current$feature_id,
      earlier$feature_id),
    selected_axis_identical = identical(selected_current, selected_earlier),
    selected_set_identical = setequal(selected_current, selected_earlier),
    rank_identical = identical(as.integer(current$rank),
      as.integer(earlier$rank)),
    method_identical = identical(as.character(current$method),
      as.character(earlier$method)),
    depth_adjustment_status_identical = identical(
      as.character(current$depth_adjustment_status),
      as.character(earlier$depth_adjustment_status)),
    maximum_absolute_differences = differences)
}

tolerance <- as.numeric(config$acceptance$numerical_tolerance)
gate_passed <- all(vapply(comparisons, function(comparison) {
  differences <- unlist(comparison$maximum_absolute_differences,
    use.names = FALSE)
  isTRUE(comparison$complete_feature_axis_identical) &&
    isTRUE(comparison$selected_axis_identical) &&
    isTRUE(comparison$rank_identical) &&
    isTRUE(comparison$method_identical) &&
    isTRUE(comparison$depth_adjustment_status_identical) &&
    all(is.na(differences) | differences <= tolerance)
}, logical(1)))
scientific_state_after <- scientific_state(bundle)
wb_close_bundle(bundle)
bundle <- NULL
source_manifest_after <- sha(source_manifest)
source_manifest_document_after <- jsonlite::read_json(source_manifest,
  simplifyVector = FALSE)
source_manifest_unchanged <- identical(source_manifest_before,
  source_manifest_after)
source_manifest_timestamp_only <- identical(
  without_publication_time(source_manifest_document_before),
  without_publication_time(source_manifest_document_after))
original_manifest_after <- if (is.null(original_manifest)) NA_character_ else
  sha(original_manifest)
original_source_unchanged <- is.null(original_source) || identical(
  original_manifest_before, original_manifest_after)
scientific_state_unchanged <- identical(scientific_state_before,
  scientific_state_after)
gate_passed <- gate_passed && scientific_state_unchanged &&
  original_source_unchanged &&
  (source_manifest_unchanged || (isTRUE(config$acceptance$
    permit_qualification_parent_publication_timestamp_refresh) &&
    source_manifest_timestamp_only))

summary <- list(
  schema_version = "bdev.scwbatac_hvp_sparse_requalification_summary.v1",
  status = if (gate_passed) "passed" else "failed",
  run_identity = config$run_identity,
  operation = config$feature_selection$operation,
  scworkbench_commit = config$runtime$scworkbench$source_commit,
  scwbATAC_commit = config$runtime$scwbATAC$source_commit,
  cells = length(cells),
  features = length(features),
  selected_features_per_mode = lapply(new_ledgers,
    function(ledger) sum(ledger$selected)),
  elapsed_seconds = elapsed,
  comparisons = comparisons,
  retained_session_evidence = session_evidence,
  source_manifest_before = source_manifest_before,
  source_manifest_after = source_manifest_after,
  source_manifest_unchanged = source_manifest_unchanged,
  source_manifest_timestamp_only_change = source_manifest_timestamp_only,
  source_scientific_state_before = scientific_state_before,
  source_scientific_state_after = scientific_state_after,
  source_scientific_state_unchanged = scientific_state_unchanged,
  original_source_manifest_before = original_manifest_before,
  original_source_manifest_after = original_manifest_after,
  original_source_manifest_unchanged = original_source_unchanged,
  implementation_unlocked = FALSE,
  sif_built = FALSE,
  gate_passed = gate_passed)
write_json(summary, file.path(directories[["reports"]],
  "01_hvp_sparse_requalification_summary.json"))
write_json(list(
  schema_version = "bdev.scwbatac_hvp_sparse_execution_receipt.v1",
  status = if (gate_passed) "complete_passed" else "complete_failed",
  config_path = config_path,
  config_sha256 = sha(config_path),
  sif_sha256 = sha(config$runtime$sif),
  operation_pack_manifest_sha256 = sha(manifest),
  source_manifest_sha256 = source_manifest_after,
  gate_passed = gate_passed), file.path(directories[["manifests"]],
  "02_execution_receipt.json"))
stopifnot(gate_passed)
log_event("Sparse-safe full-Bdev HVP requalification passed")
