#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config")) {
  stop("Usage: 02_score_bdev_fixture.R --config <config.yaml>")
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
intermediate <- cfg$storage$intermediate_root
dirs <- setNames(file.path(root, c("01_tables", "03_reports", "04_logs", "06_manifests")),
  c("tables", "reports", "logs", "manifests"))
session_root <- file.path(intermediate, "05_objects", "01_sessions")
invisible(lapply(c(dirs, sessions = session_root), dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_bdev_scoring.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
write_tsv <- function(x, path) {
  connection <- if (grepl("[.]gz$", path)) gzfile(path, "wt") else path
  on.exit(if (inherits(connection, "connection")) close(connection), add = TRUE)
  write.table(x, connection, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  output <- lapply(fields, function(x) unique(x[seq.int(3L, length(x))]))
  names(output) <- vapply(fields, `[[`, character(1), 1L)
  stopifnot(!anyDuplicated(names(output)))
  output
}
run_ref <- function(run, output) wb_ref_run(run$run_id, output)
input_path <- function(name) cfg$inputs[[name]]$path
fixture <- read.delim(input_path("fixture_cells"), check.names = FALSE)
broad <- read.delim(input_path("bdev_broad_decisions"), check.names = FALSE)
thresholds <- read.delim(file.path(dirs[["tables"]], "02_population_margin_calibration.tsv"), check.names = FALSE)
population_map <- read.delim(input_path("population_summary"), check.names = FALSE)
sets <- read_gmt(input_path("population_gmt"))
stopifnot(nrow(fixture) == as.integer(cfg$bdev_test$expected_cells),
  !anyDuplicated(fixture$cell_id), !anyDuplicated(broad$cell_id),
  identical(sort(fixture$cell_id), sort(broad$cell_id)), length(sets) == 64L,
  identical(sort(names(sets)), sort(population_map$geneset_id)))

set_to_target <- setNames(population_map$target, population_map$geneset_id)
set_to_lineage <- setNames(population_map$lineage, population_map$geneset_id)
threshold_by_target <- setNames(thresholds$calibrated_minimum_margin, thresholds$target)
available_by_target <- setNames(as.logical(thresholds$threshold_available), thresholds$target)

parent_hierarchy <- jsonlite::read_json(input_path("parent_hierarchy_manifest"), simplifyVector = TRUE)
source_bundle <- parent_hierarchy$parent_bundle$path
execution_bundle <- file.path(intermediate, "05_objects", "03_fresh_session_bundle")
stopifnot(dir.exists(execution_bundle))
log_event("Loading operation pack and opening fresh isolated development session bundle")
wb_operation_pack_load("scwbGenePrograms")
pack_record <- wb_operation_packs()
stopifnot(nrow(pack_record) == 1L,
  identical(as.character(pack_record$manifest_sha256[[1L]]),
    cfg$bdev_test$operation_pack_manifest_sha256))
log_event("Operation pack loaded; opening bundle")
bundle <- wb_open_bundle(execution_bundle, validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
log_event("Bundle opened; reading manifest")
manifest_before <- wb_manifest(bundle)
stopifnot(identical(manifest_before$current_state$fingerprint, cfg$bdev_test$bundle_fingerprint))

log_event("Manifest verified; reading fixture metadata")
fixture_md <- wb_meta(bundle, cells = fixture$cell_id)
stopifnot(identical(as.character(fixture_md$cell_id), fixture$cell_id))

execution_session_id <- paste0(cfg$bdev_test$session_id, "_fresh1")
log_event("Fixture metadata read; creating session", execution_session_id)
session <- wb_session_create(bundle, workspace_root = session_root,
  session_id = execution_session_id)
log_event("Scoring 64 global-local programs without import or promotion")
score_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::ucell",
  inputs = list(
    matrix = cfg$bdev_test$matrix_id,
    cell_ids = fixture$cell_id,
    min_cells = 10L,
    features = sets,
    name = "IMMGEN_GLOBAL_LOCAL_HIERARCHY",
    max_rank = as.integer(cfg$bdev_test$ucell_max_rank),
    ncores = 1L,
    parallel_backend = "serial",
    chunk_size = 500L,
    missing_features_policy = "drop"),
  notes = "Diagnostic global-local population scoring for hierarchical calibration")
session <- score_result$session
calibration_result <- wb_session_run_op(session, bundle,
  "scwbgeneprograms::cell_vector_family_calibrate", inputs = list(
    signal = run_ref(score_result$run, "signature_scores"),
    cell_ids = fixture$cell_id,
    min_cells = 10L,
    method = "percentile",
    ties_method = "average",
    na_rm = TRUE),
  notes = "Diagnostic percentile calibration for global-local programs")

score_audit <- score_result$run$outputs$signature_summary$table
score_audit$detected_fraction <- score_audit$n_features_matched / score_audit$n_features_input
stopifnot(all(score_audit$n_features_matched >= 5L), all(score_audit$detected_fraction >= 0.5))
write_tsv(score_audit, file.path(dirs[["tables"]], "03_bdev_score_overlap_audit.tsv"))

primitive <- calibration_result$run$outputs$calibrated_signal
calibrated <- data.frame(cell_id = as.character(primitive$cell_id),
  geneset_id = as.character(primitive$key), calibrated_score = as.numeric(primitive$value),
  stringsAsFactors = FALSE)
calibrated$target <- unname(set_to_target[calibrated$geneset_id])
calibrated$lineage <- unname(set_to_lineage[calibrated$geneset_id])
stopifnot(!anyNA(calibrated$target), !anyNA(calibrated$lineage))
write_tsv(calibrated, file.path(dirs[["tables"]], "04_bdev_calibrated_population_scores.tsv.gz"))

score_split <- split(calibrated, calibrated$cell_id)
broad_by_cell <- split(broad, broad$cell_id)
policies <- unlist(cfg$bdev_test$policies)
decision_rows <- list()
row_index <- 0L
for (cell_id in fixture$cell_id) {
  broad_row <- broad_by_cell[[cell_id]]
  states <- if (broad_row$eligible_broad_state_count[[1L]] == 0L) character() else
    strsplit(broad_row$eligible_broad_states[[1L]], ";", fixed = TRUE)[[1L]]
  rows <- score_split[[cell_id]]
  rows <- rows[rows$lineage %in% states, , drop = FALSE]
  if (nrow(rows)) rows <- rows[order(-rows$calibrated_score, rows$geneset_id), , drop = FALSE]
  for (policy in policies) {
    row_index <- row_index + 1L
    if (!nrow(rows)) {
      winner <- label <- NA_character_; winner_lineage <- NA_character_
      top <- runner <- margin <- threshold <- NA_real_; reason <- "broad_unassigned"
    } else {
      winner <- rows$target[[1L]]; winner_lineage <- rows$lineage[[1L]]
      top <- rows$calibrated_score[[1L]]
      runner <- if (nrow(rows) > 1L) rows$calibrated_score[[2L]] else NA_real_
      margin <- if (is.finite(runner)) top - runner else Inf
      threshold <- if (policy == "fixed_margin_0.1") 0.1 else if (policy == "population_calibrated" &&
        isTRUE(available_by_target[[winner]])) threshold_by_target[[winner]] else 0
      if (top < 0.8) reason <- "min_score"
      else if (is.finite(margin) && margin <= 1e-8) reason <- "tie"
      else if (policy == "population_calibrated" && !isTRUE(available_by_target[[winner]])) reason <- "no_population_calibration"
      else if (margin + 1e-12 < threshold) reason <- "min_margin"
      else reason <- "none"
      label <- if (reason == "none") winner else NA_character_
    }
    decision_rows[[row_index]] <- data.frame(
      policy = policy, cell_id = cell_id,
      broad_top_state = broad_row$broad_top_state[[1L]],
      eligible_broad_states = broad_row$eligible_broad_states[[1L]],
      eligible_broad_state_count = broad_row$eligible_broad_state_count[[1L]],
      winner_population = winner, winner_lineage = winner_lineage,
      fine_label = label, top_score = top, runner_up_score = runner,
      score_margin = margin, required_margin = threshold,
      decision = if (reason == "none") "fine_label" else if (length(states)) "broad_only" else "unassigned",
      abstain_reason = reason, stringsAsFactors = FALSE)
  }
}
decisions <- do.call(rbind, decision_rows)
decisions <- merge(decisions,
  fixture_md[, c("cell_id", "peaks_snn_res.1.2", "pt", "Phase", "orig.ident"), drop = FALSE],
  by = "cell_id", all.x = TRUE, sort = FALSE)
write_tsv(decisions, file.path(dirs[["tables"]], "05_bdev_hierarchical_policy_decisions.tsv.gz"))

summary_rows <- do.call(rbind, lapply(policies, function(policy) {
  rows <- decisions[decisions$policy == policy, , drop = FALSE]
  data.frame(policy = policy, fixture_cells = nrow(rows),
    fine_labels = sum(rows$decision == "fine_label"),
    broad_only = sum(rows$decision == "broad_only"),
    fully_unassigned = sum(rows$decision == "unassigned"),
    b_lineage_fine_labels = sum(rows$decision == "fine_label" & rows$winner_lineage == "B"),
    stem_prog_fine_labels = sum(rows$decision == "fine_label" & rows$winner_lineage == "Stem/Prog"),
    stringsAsFactors = FALSE)
}))
write_tsv(summary_rows, file.path(dirs[["tables"]], "06_bdev_policy_summary.tsv"))
label_summary <- as.data.frame(table(decisions$policy, decisions$fine_label, useNA = "no"), stringsAsFactors = FALSE)
names(label_summary) <- c("policy", "fine_label", "cells")
label_summary <- label_summary[label_summary$cells > 0L, , drop = FALSE]
write_tsv(label_summary, file.path(dirs[["tables"]], "07_bdev_fine_label_summary.tsv"))

manifest_after <- wb_manifest(bundle)
stopifnot(identical(manifest_after$current_state$fingerprint, cfg$bdev_test$bundle_fingerprint))
reference_summary <- jsonlite::read_json(file.path(dirs[["reports"]], "01_reference_calibration_summary.json"), simplifyVector = TRUE)
summary <- list(
  schema_version = "bdev.immgen_hierarchical_population_calibration_summary.v1",
  status = "complete_diagnostic_only",
  run_identity = cfg$run_identity,
  reference = reference_summary,
  bdev_policy_summary = summary_rows,
  exact_pt_preserved = TRUE,
  exact_pt_evidence = "inherited_from_checksum_verified_parent_hierarchy_manifest_without_bundle_mutation",
  execution_bundle = execution_bundle,
  session_bundle_registry_unchanged = TRUE,
  session_bundle_fingerprint = manifest_after$current_state$fingerprint,
  accepted_labels_emitted = FALSE,
  feature_sets_registered = FALSE,
  artifacts_promoted = FALSE)
write_json(summary, file.path(dirs[["reports"]], "02_combined_summary.json"))
report <- c(
  "# Hierarchical population-specific calibration test",
  "",
  paste0("- Held-out true lineage retained: ", reference_summary$true_lineage_eligible, "/", reference_summary$heldout_libraries),
  paste0("- Held-out no-margin local precision: ", round(100 * reference_summary$no_margin_precision, 1), "%"),
  paste0("- Populations with an exploratory threshold: ", reference_summary$populations_with_calibrated_threshold, "/64"),
  "",
  paste0("- Bdev no-margin fine calls: ", summary_rows$fine_labels[summary_rows$policy == "no_margin"]),
  paste0("- Bdev fixed-0.10 fine calls: ", summary_rows$fine_labels[summary_rows$policy == "fixed_margin_0.1"]),
  paste0("- Bdev population-calibrated fine calls: ", summary_rows$fine_labels[summary_rows$policy == "population_calibrated"]),
  "",
  "Population thresholds are optimistic resubstitution estimates with small support. All calls are diagnostic; broad-only retention is the intended fallback.")
writeLines(report, file.path(dirs[["reports"]], "03_hierarchical_calibration_report.md"))
write_json(list(
  schema_version = "bdev.immgen_hierarchical_population_session_receipt.v1",
  status = "complete_without_promotion",
  session_id = execution_session_id,
  score_run_id = score_result$run$run_id,
  calibration_run_id = calibration_result$run$run_id,
  execution_bundle = execution_bundle,
  bundle_fingerprint = manifest_after$current_state$fingerprint,
  registry_unchanged = TRUE,
  scworkbench_artifacts_promoted = FALSE),
  file.path(dirs[["manifests"]], "02_session_receipt.json"))
log_event("Completed diagnostic hierarchy calibration test")
