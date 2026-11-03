#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config")) stop("Usage: 02_score_bdev_raw.R --config <config.yaml>")
suppressPackageStartupMessages({ library(scworkbench); library(yaml); library(jsonlite); library(digest) })

cfg <- yaml::read_yaml(normalizePath(args[[2L]], mustWork = TRUE))
root <- cfg$storage$results_root
intermediate <- cfg$storage$intermediate_root
dirs <- setNames(file.path(root, c("01_tables", "03_reports", "04_logs", "06_manifests")), c("tables", "reports", "logs", "manifests"))
session_root <- file.path(intermediate, "05_objects", "01_sessions")
invisible(lapply(c(dirs, sessions = session_root), dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_bdev_raw_scoring.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
run_ref <- function(run, output) wb_ref_run(run$run_id, output)
read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  out <- lapply(fields, function(x) unique(x[seq.int(3L, length(x))])); names(out) <- vapply(fields, `[[`, character(1), 1L); out
}

fixture <- read.delim(cfg$inputs$fixture_cells$path, check.names = FALSE)
stopifnot(nrow(fixture) == as.integer(cfg$fixture$expected_cells), !anyDuplicated(fixture$cell_id))
sets <- read_gmt(file.path(dirs[["tables"]], "03_simple_mean_contrast_sets.gmt"))
stopifnot(length(sets) == 64L, all(lengths(sets) == 20L), !anyDuplicated(names(sets)))
wb_operation_pack_load("scwbGenePrograms")
pack_record <- wb_operation_packs()
stopifnot(nrow(pack_record) == 1L, identical(as.character(pack_record$manifest_sha256[[1L]]), cfg$runtime$operation_pack_manifest_sha256))

bundle <- wb_open_bundle(cfg$inputs$execution_bundle$path, validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest_before <- wb_manifest(bundle)
stopifnot(identical(manifest_before$current_state$fingerprint, cfg$inputs$execution_bundle$fingerprint))
fixture_md <- wb_meta(bundle, cells = fixture$cell_id)
stopifnot(identical(as.character(fixture_md$cell_id), fixture$cell_id))
metadata <- fixture_md[, c("cell_id", cfg$fixture$sample_field, cfg$fixture$pseudotime_field), drop = FALSE]
names(metadata)[names(metadata) == cfg$fixture$sample_field] <- cfg$fixture$reported_sample_field

session <- wb_session_create(bundle, workspace_root = session_root, session_id = "simple_mean_contrast_bdev_v1")
score_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::ucell", inputs = list(
  matrix = cfg$scoring$matrix_id, cell_ids = fixture$cell_id, min_cells = 10L,
  features = sets, name = "SIMPLE_MEAN_CONTRAST", max_rank = as.integer(cfg$scoring$ucell_max_rank),
  ncores = as.integer(cfg$scoring$ncores), parallel_backend = cfg$scoring$parallel_backend,
  chunk_size = as.integer(cfg$scoring$chunk_size), missing_features_policy = cfg$scoring$missing_features_policy),
  notes = "Raw UCell scoring of simple ImmGen population-balanced mean-contrast top20 sets.")
reduce_result <- wb_session_run_op(score_result$session, bundle, "scwbgeneprograms::cell_vector_family_reduce", inputs = list(
  signal = run_ref(score_result$run, "signature_scores"), cell_ids = fixture$cell_id, min_cells = 10L,
  rule = cfg$scoring$reduction_rule, min_score = as.numeric(cfg$scoring$minimum_raw_score),
  min_margin = as.numeric(cfg$scoring$minimum_raw_margin), unknown_label = cfg$scoring$unknown_label,
  ties = cfg$scoring$ties, tie_tol = as.numeric(cfg$scoring$tie_tolerance)),
  notes = "Direct raw-UCell reduction without percentile calibration or lineage gates.")

write_tsv(score_result$run$outputs$signature_summary$table, file.path(dirs[["tables"]], "05_ucell_signature_summary.tsv"))
decision <- merge(reduce_result$run$outputs$decision_table$table, metadata, by = "cell_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(decision) == nrow(fixture), !anyNA(decision[[cfg$fixture$reported_sample_field]]))
write_tsv(decision, file.path(dirs[["tables"]], "06_bdev_raw_decisions.tsv"))

sample_ids <- sort(unique(as.character(decision[[cfg$fixture$reported_sample_field]])))
sample_report <- lapply(sample_ids, function(sample_id) {
  rows <- decision[as.character(decision[[cfg$fixture$reported_sample_field]]) == sample_id, , drop = FALSE]
  assigned <- rows[rows$label != cfg$scoring$unknown_label, , drop = FALSE]
  counts <- sort(table(assigned$label), decreasing = TRUE); top <- head(counts, 3L)
  list(sampleId = sample_id, cells = nrow(rows), assigned = nrow(assigned), unassigned = nrow(rows) - nrow(assigned),
    top3 = lapply(seq_along(top), function(i) list(label = names(top)[[i]], cells = as.integer(top[[i]]),
      fraction_all = as.numeric(top[[i]]) / nrow(rows), fraction_assigned = as.numeric(top[[i]]) / max(1L, nrow(assigned)))))
})
manifest_after <- wb_manifest(bundle)
stopifnot(identical(manifest_after$current_state$fingerprint, manifest_before$current_state$fingerprint))
summary <- list(schema_version = "bdev.simple_mean_contrast_bdev_summary.v1", status = "complete_diagnostic_only",
  run_identity = cfg$run_identity, cells = nrow(decision), sample_ids = sample_report,
  assigned = sum(decision$label != cfg$scoring$unknown_label), unassigned = sum(decision$label == cfg$scoring$unknown_label),
  minimum_raw_margin = cfg$scoring$minimum_raw_margin, percentile_calibration_used = FALSE,
  broad_or_lineage_gate_used = FALSE, bdev_used_for_gene_selection = FALSE,
  bundle_fingerprint_unchanged = TRUE, artifacts_promoted = FALSE, accepted_labels_emitted = FALSE)
write_json(summary, file.path(dirs[["reports"]], "02_bdev_summary.json"))
write_json(list(schema_version = "bdev.simple_mean_contrast_session.v1", status = "complete_without_promotion",
  session_id = "simple_mean_contrast_bdev_v1", score_run_id = score_result$run$run_id,
  reduction_run_id = reduce_result$run$run_id, artifacts_promoted = FALSE),
  file.path(dirs[["manifests"]], "01_session_receipt.json"))

