#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config")) stop("Usage: 01_reduce_margin001.R --config <config.yaml>")
suppressPackageStartupMessages({ library(scworkbench); library(yaml); library(jsonlite); library(digest) })

cfg <- yaml::read_yaml(normalizePath(args[[2L]], mustWork = TRUE))
root <- cfg$storage$results_root
intermediate <- cfg$storage$intermediate_root
dirs <- setNames(file.path(root, c("01_tables", "03_reports", "04_logs", "06_manifests")), c("tables", "reports", "logs", "manifests"))
session_root <- file.path(intermediate, "05_objects", "01_sessions")
invisible(lapply(c(dirs, sessions = session_root), dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_margin001_reduction.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

stopifnot(identical(sha256(cfg$inputs$calibrated_signal$path), cfg$inputs$calibrated_signal$sha256),
  identical(sha256(file.path(cfg$inputs$execution_bundle$path, "manifest.json")), cfg$inputs$execution_bundle$manifest_sha256))
fixture <- read.delim(cfg$inputs$fixture_cells$path, check.names = FALSE)
stopifnot(nrow(fixture) == as.integer(cfg$fixture$expected_cells), !anyDuplicated(fixture$cell_id))
signal <- readRDS(cfg$inputs$calibrated_signal$path)
stopifnot(inherits(signal, "wb_primitive"), identical(signal$kind, "cell_vector_family"),
  length(unique(signal$cell_id)) == nrow(fixture), length(unique(signal$key)) == 47L,
  length(signal$value) == nrow(fixture) * 47L)

wb_operation_pack_load("scwbGenePrograms")
pack_record <- wb_operation_packs()
stopifnot(nrow(pack_record) == 1L, identical(as.character(pack_record$manifest_sha256[[1L]]), cfg$runtime$operation_pack_manifest_sha256))
write_json(wb_op_interface("scwbgeneprograms::cell_vector_family_reduce"), file.path(dirs[["manifests"]], "interface_cell_vector_family_reduce.json"))

bundle <- wb_open_bundle(cfg$inputs$execution_bundle$path, validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest_before <- wb_manifest(bundle)
stopifnot(identical(manifest_before$current_state$fingerprint, cfg$inputs$execution_bundle$fingerprint))
fixture_md <- wb_meta(bundle, cells = fixture$cell_id)
stopifnot(identical(as.character(fixture_md$cell_id), fixture$cell_id))

session <- wb_session_create(bundle, workspace_root = session_root, session_id = "margin001_reduce_v1")
result <- wb_session_run_op(session, bundle, "scwbgeneprograms::cell_vector_family_reduce", inputs = list(
  signal = signal, cell_ids = fixture$cell_id, min_cells = 10L,
  rule = cfg$reduction$rule, min_score = as.numeric(cfg$reduction$minimum_calibrated_score),
  min_margin = as.numeric(cfg$reduction$minimum_calibrated_margin), unknown_label = cfg$reduction$unknown_label,
  ties = cfg$reduction$ties, tie_tol = as.numeric(cfg$reduction$tie_tolerance)),
  notes = "Exploratory native reduction of the pinned calibrated narrow47 scores at margin 0.01.")

decision <- result$run$outputs$decision_table$table
decision$population <- as.character(decision$label)
metadata <- fixture_md[, c("cell_id", cfg$fixture$sample_field, cfg$fixture$pseudotime_field), drop = FALSE]
names(metadata)[names(metadata) == cfg$fixture$sample_field] <- cfg$fixture$reported_sample_field
decision <- merge(decision, metadata, by = "cell_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(decision) == nrow(fixture), !anyNA(decision[[cfg$fixture$reported_sample_field]]))
write_tsv(decision, file.path(dirs[["tables"]], "01_cell_label_decisions_margin001.tsv"))

counts <- as.data.frame(table(decision[[cfg$fixture$reported_sample_field]], decision$population), stringsAsFactors = FALSE)
names(counts) <- c(cfg$fixture$reported_sample_field, "population", "cells")
counts <- counts[counts$cells > 0L, , drop = FALSE]
totals <- aggregate(counts$cells, counts[cfg$fixture$reported_sample_field], sum)
names(totals)[2L] <- "sample_total"
counts <- merge(counts, totals, by = cfg$fixture$reported_sample_field, sort = FALSE)
counts$fraction <- counts$cells / counts$sample_total
write_tsv(counts, file.path(dirs[["tables"]], "02_labels_by_sampleid_margin001.tsv"))

sample_ids <- sort(unique(as.character(decision[[cfg$fixture$reported_sample_field]])))
sample_report <- lapply(sample_ids, function(sample_id) {
  rows <- decision[as.character(decision[[cfg$fixture$reported_sample_field]]) == sample_id, , drop = FALSE]
  assigned <- rows[rows$population != cfg$reduction$unknown_label, , drop = FALSE]
  label_counts <- sort(table(assigned$population), decreasing = TRUE)
  top <- head(label_counts, 3L)
  list(sampleId = sample_id, cells = nrow(rows), assigned = nrow(assigned), unassigned = nrow(rows) - nrow(assigned),
    top3 = lapply(seq_along(top), function(i) list(label = names(top)[[i]], cells = as.integer(top[[i]]),
      fraction_all = as.numeric(top[[i]]) / nrow(rows), fraction_assigned = as.numeric(top[[i]]) / max(1L, nrow(assigned)))))
})
manifest_after <- wb_manifest(bundle)
stopifnot(identical(manifest_after$current_state$fingerprint, manifest_before$current_state$fingerprint))
summary <- list(schema_version = "bdev.user_narrow47_margin001_summary.v1", status = "complete_diagnostic_only",
  run_identity = cfg$run_identity, parent_run = cfg$inputs$parent_run, cells = nrow(decision),
  minimum_calibrated_score = cfg$reduction$minimum_calibrated_score,
  minimum_calibrated_margin = cfg$reduction$minimum_calibrated_margin,
  assigned = sum(decision$population != cfg$reduction$unknown_label),
  unassigned = sum(decision$population == cfg$reduction$unknown_label), sample_ids = sample_report,
  broad_or_lineage_gate_used = FALSE, native_scworkbench_reduction = TRUE,
  bundle_fingerprint_unchanged = TRUE, score_artifacts_promoted = FALSE, accepted_labels_emitted = FALSE)
write_json(summary, file.path(dirs[["reports"]], "01_margin001_summary.json"))
write_json(list(schema_version = "bdev.user_narrow47_margin001_session_receipt.v1", status = "complete_without_promotion",
  reduction_run_id = result$run$run_id, session_id = "margin001_reduce_v1", artifacts_promoted = FALSE),
  file.path(dirs[["manifests"]], "01_session_receipt.json"))

