#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config")) {
  stop("Usage: 03_score_bdev_fixture.R --config <config.yaml>")
}

suppressPackageStartupMessages({
  library(scworkbench)
  library(yaml)
  library(jsonlite)
})

config_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(config_path)
root <- cfg$storage$results_root
intermediate <- cfg$storage$intermediate_root
dirs <- setNames(file.path(root, c("01_tables", "03_reports", "04_logs", "06_manifests")),
  c("tables", "reports", "logs", "manifests"))
session_root <- file.path(intermediate, "05_objects", "02_sessions")
invisible(lapply(c(dirs, sessions = session_root), dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_bdev_scoring.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

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
  output
}
run_ref <- function(run, output) wb_ref_run(run$run_id, output)
count_summary <- function(data, fields) {
  output <- as.data.frame(table(data[fields], useNA = "ifany"), stringsAsFactors = FALSE)
  names(output)[ncol(output)] <- "cells"
  output <- output[output$cells > 0L, , drop = FALSE]
  totals <- aggregate(output$cells, output[fields[1L]], sum)
  names(totals)[2L] <- "sample_total"
  output <- merge(output, totals, by = fields[1L], sort = FALSE)
  output$fraction <- output$cells / output$sample_total
  output
}

fixture <- read.delim(cfg$inputs$fixture_cells$path, check.names = FALSE)
stopifnot(nrow(fixture) == as.integer(cfg$bdev_test$expected_cells), !anyDuplicated(fixture$cell_id))
gmt_path <- file.path(dirs[["tables"]], "02_immgen_stable_broad_sets.gmt")
sets <- read_gmt(gmt_path)
stopifnot(length(sets) == 9L, all(lengths(sets) >= 5L))
lineage_by_set <- setNames(sub("^IMMGEN_STABLE_BROAD_", "", names(sets)), names(sets))
lineage_by_set <- gsub("STEM_PROG", "Stem&Prog", lineage_by_set, fixed = TRUE)
lineage_by_set <- gsub("MONOCYTE", "monocyte", lineage_by_set, fixed = TRUE)
lineage_by_set <- gsub("MYELOID", "myeloid", lineage_by_set, fixed = TRUE)
lineage_by_set <- gsub("STROMA", "stroma", lineage_by_set, fixed = TRUE)
lineage_by_set <- gsub("ABT", "abT", lineage_by_set, fixed = TRUE)
lineage_by_set <- gsub("GDT", "gdT", lineage_by_set, fixed = TRUE)

bundle_path <- file.path(intermediate, "05_objects", "01_session_bundle")
wb_operation_pack_load("scwbGenePrograms")
pack_record <- wb_operation_packs()
stopifnot(nrow(pack_record) == 1L,
  identical(as.character(pack_record$manifest_sha256[[1L]]), cfg$bdev_test$operation_pack_manifest_sha256))
bundle <- wb_open_bundle(bundle_path, validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest_before <- wb_manifest(bundle)
stopifnot(identical(manifest_before$current_state$fingerprint, cfg$inputs$source_bundle$fingerprint))
fixture_md <- wb_meta(bundle, cells = fixture$cell_id)
stopifnot(identical(as.character(fixture_md$cell_id), fixture$cell_id))

session <- wb_session_create(bundle, workspace_root = session_root,
  session_id = paste0(cfg$bdev_test$session_id, "_fresh1"))
score_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::ucell",
  inputs = list(
    matrix = cfg$bdev_test$matrix_id,
    cell_ids = fixture$cell_id,
    min_cells = 10L,
    features = sets,
    name = "IMMGEN_STABLE_BROAD",
    max_rank = as.integer(cfg$bdev_test$ucell_max_rank),
    ncores = 1L,
    parallel_backend = "serial",
    chunk_size = 500L,
    missing_features_policy = "drop"),
  notes = "Diagnostic ImmGen-only stable broad-lineage scoring")
session <- score_result$session
calibration_result <- wb_session_run_op(session, bundle,
  "scwbgeneprograms::cell_vector_family_calibrate", inputs = list(
    signal = run_ref(score_result$run, "signature_scores"),
    cell_ids = fixture$cell_id,
    min_cells = 10L,
    method = "percentile",
    ties_method = "average",
    na_rm = TRUE),
  notes = "Diagnostic percentile calibration of stable broad programs")

score_audit <- score_result$run$outputs$signature_summary$table
score_audit$detected_fraction <- score_audit$n_features_matched / score_audit$n_features_input
stopifnot(all(score_audit$n_features_matched >= 5L), all(score_audit$detected_fraction >= 0.5))
write_tsv(score_audit, file.path(dirs[["tables"]], "06_bdev_score_overlap_audit.tsv"))

primitive <- calibration_result$run$outputs$calibrated_signal
scores <- data.frame(cell_id = as.character(primitive$cell_id),
  geneset_id = as.character(primitive$key), calibrated_score = as.numeric(primitive$value),
  stringsAsFactors = FALSE)
scores$lineage <- unname(lineage_by_set[scores$geneset_id])
stopifnot(!anyNA(scores$lineage))
write_tsv(scores, file.path(dirs[["tables"]], "07_bdev_broad_scores.tsv.gz"))

score_minimum <- as.numeric(cfg$bdev_test$score_minimum)
fixed_margin <- as.numeric(cfg$bdev_test$fixed_margin)
score_split <- split(scores, scores$cell_id)
decision_rows <- vector("list", length(fixture$cell_id))
for (index in seq_along(fixture$cell_id)) {
  cell_id <- fixture$cell_id[[index]]
  rows <- score_split[[cell_id]]
  rows <- rows[order(-rows$calibrated_score, rows$lineage), , drop = FALSE]
  top <- rows$calibrated_score[[1L]]
  runner <- rows$calibrated_score[[2L]]
  margin <- top - runner
  candidates <- rows$lineage[
    rows$calibrated_score >= score_minimum & top - rows$calibrated_score <= fixed_margin + 1e-12
  ]
  status <- if (top < score_minimum) "unassigned_low_score" else if (margin < fixed_margin) "borderline" else "accepted"
  decision_rows[[index]] <- data.frame(
    cell_id = cell_id,
    top_lineage = rows$lineage[[1L]],
    broad_label = if (status == "accepted") rows$lineage[[1L]] else NA_character_,
    candidate_broad_states = if (length(candidates)) paste(candidates, collapse = ";") else NA_character_,
    candidate_state_count = length(candidates),
    top_score = top,
    runner_up_score = runner,
    score_margin = margin,
    status = status,
    stringsAsFactors = FALSE)
}
decisions <- do.call(rbind, decision_rows)
metadata_fields <- c("cell_id", cfg$bdev_test$sample_field, cfg$bdev_test$pseudotime_field)
decisions <- merge(decisions, fixture_md[, metadata_fields, drop = FALSE], by = "cell_id", all.x = TRUE, sort = FALSE)
names(decisions)[names(decisions) == cfg$bdev_test$sample_field] <- cfg$bdev_test$reported_sample_field
stopifnot(!anyNA(decisions$sampleId), length(unique(decisions$sampleId)) == 4L)
write_tsv(decisions, file.path(dirs[["tables"]], "08_bdev_broad_decisions.tsv.gz"))

label_rows <- decisions
label_rows$reported_label <- ifelse(label_rows$status == "accepted", label_rows$broad_label,
  ifelse(label_rows$status == "borderline", "borderline", "unassigned"))
label_summary <- count_summary(label_rows, c("sampleId", "reported_label"))
status_summary <- count_summary(decisions, c("sampleId", "status"))
top_summary <- count_summary(decisions, c("sampleId", "top_lineage"))
write_tsv(label_summary, file.path(dirs[["tables"]], "09_bdev_sampleid_label_summary.tsv"))
write_tsv(status_summary, file.path(dirs[["tables"]], "10_bdev_sampleid_status_summary.tsv"))
write_tsv(top_summary, file.path(dirs[["tables"]], "11_bdev_sampleid_forced_top_summary.tsv"))

overall <- as.data.frame(table(decisions$status), stringsAsFactors = FALSE)
names(overall) <- c("status", "cells")
sample_ids <- sort(unique(as.character(decisions$sampleId)))
sample_report <- lapply(sample_ids, function(sample_id) {
  rows <- decisions[as.character(decisions$sampleId) == sample_id, , drop = FALSE]
  list(sampleId = sample_id, cells = nrow(rows), accepted = sum(rows$status == "accepted"),
    borderline = sum(rows$status == "borderline"), unassigned = sum(rows$status == "unassigned_low_score"))
})
manifest_after <- wb_manifest(bundle)
stopifnot(identical(manifest_after$current_state$fingerprint, cfg$inputs$source_bundle$fingerprint))
summary <- list(
  schema_version = "bdev.stable_immgen_broad_bdev_summary.v1",
  status = "complete_diagnostic_only",
  run_identity = cfg$run_identity,
  fixture_cells = nrow(decisions),
  sample_field_source = cfg$bdev_test$sample_field,
  reported_sample_field = cfg$bdev_test$reported_sample_field,
  sample_ids = sample_report,
  overall_status = overall,
  exact_pt_preserved = TRUE,
  pt_defined = sum(!is.na(decisions$pt)),
  pt_missing = sum(is.na(decisions$pt)),
  bundle_fingerprint = manifest_after$current_state$fingerprint,
  bundle_registry_unchanged = TRUE,
  Bdev_used_for_gene_selection = FALSE,
  accepted_labels_emitted = FALSE,
  artifacts_promoted = FALSE)
write_json(summary, file.path(dirs[["reports"]], "02_bdev_summary.json"))
write_json(list(
  schema_version = "bdev.stable_immgen_broad_session_receipt.v1",
  status = "complete_without_promotion",
  session_id = paste0(cfg$bdev_test$session_id, "_fresh1"),
  score_run_id = score_result$run$run_id,
  calibration_run_id = calibration_result$run$run_id,
  execution_bundle = bundle_path,
  bundle_fingerprint = manifest_after$current_state$fingerprint,
  registry_unchanged = TRUE,
  scworkbench_artifacts_promoted = FALSE),
  file.path(dirs[["manifests"]], "03_session_receipt.json"))
