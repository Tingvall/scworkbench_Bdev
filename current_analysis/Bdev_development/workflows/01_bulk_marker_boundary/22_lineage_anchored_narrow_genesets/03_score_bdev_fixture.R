#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config")) stop("Usage: 03_score_bdev_fixture.R --config <config.yaml>")
suppressPackageStartupMessages({ library(scworkbench); library(yaml); library(jsonlite) })

config_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(config_path)
root <- cfg$storage$results_root
intermediate <- cfg$storage$intermediate_root
dirs <- setNames(file.path(root, c("01_tables", "03_reports", "04_logs", "06_manifests")), c("tables", "reports", "logs", "manifests"))
session_root <- file.path(intermediate, "05_objects", "02_sessions")
invisible(lapply(c(dirs, sessions = session_root), dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_bdev_scoring.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

write_tsv <- function(x, path) {
  connection <- if (grepl("[.]gz$", path)) gzfile(path, "wt") else path
  on.exit(if (inherits(connection, "connection")) close(connection), add = TRUE)
  write.table(x, connection, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
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
module_map <- read.delim(file.path(dirs[["tables"]], "01_module_map.tsv"), check.names = FALSE)
sets <- read_gmt(file.path(dirs[["tables"]], "03_score_collection.gmt"))
stopifnot(nrow(fixture) == as.integer(cfg$bdev_test$expected_cells), !anyDuplicated(fixture$cell_id),
  nrow(module_map) == 128L, length(sets) == 128L, identical(sort(names(sets)), sort(module_map$geneset_id)))

bundle_path <- file.path(intermediate, "05_objects", "01_session_bundle")
wb_operation_pack_load("scwbGenePrograms")
pack_record <- wb_operation_packs()
stopifnot(nrow(pack_record) == 1L, identical(as.character(pack_record$manifest_sha256[[1L]]), cfg$bdev_test$operation_pack_manifest_sha256))
bundle <- wb_open_bundle(bundle_path, validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest_before <- wb_manifest(bundle)
stopifnot(identical(manifest_before$current_state$fingerprint, cfg$inputs$source_bundle$fingerprint))
fixture_md <- wb_meta(bundle, cells = fixture$cell_id)
stopifnot(identical(as.character(fixture_md$cell_id), fixture$cell_id))

session_id <- paste0(cfg$bdev_test$session_id, "_fresh1")
session <- wb_session_create(bundle, workspace_root = session_root, session_id = session_id)
score_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::ucell", inputs = list(
  matrix = cfg$bdev_test$matrix_id, cell_ids = fixture$cell_id, min_cells = 10L,
  features = sets, name = "IMMGEN_LINEAGE_ANCHORED_NARROW",
  max_rank = as.integer(cfg$bdev_test$ucell_max_rank), ncores = 1L,
  parallel_backend = "serial", chunk_size = 500L, missing_features_policy = "drop"),
  notes = "Diagnostic scoring of 64 ImmGen anchor and local module pairs")
session <- score_result$session
calibration_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::cell_vector_family_calibrate", inputs = list(
  signal = run_ref(score_result$run, "signature_scores"), cell_ids = fixture$cell_id,
  min_cells = 10L, method = "percentile", ties_method = "average", na_rm = TRUE),
  notes = "Diagnostic percentile calibration of narrow modules")

score_audit <- score_result$run$outputs$signature_summary$table
score_audit$detected_fraction <- score_audit$n_features_matched / score_audit$n_features_input
stopifnot(all(score_audit$n_features_matched >= 5L), all(score_audit$detected_fraction >= 0.5))
write_tsv(score_audit, file.path(dirs[["tables"]], "08_bdev_score_overlap_audit.tsv"))
primitive <- calibration_result$run$outputs$calibrated_signal
scores <- data.frame(cell_id = as.character(primitive$cell_id), geneset_id = as.character(primitive$key), calibrated_score = as.numeric(primitive$value), stringsAsFactors = FALSE)
scores <- merge(scores, module_map, by = "geneset_id", all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(scores$target), !anyNA(scores$module))
write_tsv(scores, file.path(dirs[["tables"]], "09_bdev_module_scores.tsv.gz"))

populations <- unique(module_map$target)
lineage_by_target <- setNames(module_map$lineage[match(populations, module_map$target)], populations)
anchor_weight <- as.numeric(cfg$composite$anchor_weight)
local_weight <- as.numeric(cfg$composite$local_weight)
conflict_weight <- as.numeric(cfg$composite$conflict_weight)
minimum_anchor <- as.numeric(cfg$acceptance$minimum_anchor_score)
minimum_contrast <- as.numeric(cfg$acceptance$minimum_anchor_contrast)
minimum_margin <- as.numeric(cfg$acceptance$minimum_composite_margin)
score_split <- split(scores, scores$cell_id)
decision_rows <- vector("list", length(fixture$cell_id))
for (index in seq_along(fixture$cell_id)) {
  cell_id <- fixture$cell_id[[index]]
  rows <- score_split[[cell_id]]
  anchor <- setNames(rows$calibrated_score[rows$module == "anchor"], rows$target[rows$module == "anchor"])
  local <- setNames(rows$calibrated_score[rows$module == "local"], rows$target[rows$module == "local"])
  stopifnot(identical(sort(names(anchor)), sort(populations)), identical(sort(names(local)), sort(populations)))
  conflict <- vapply(populations, function(target) max(anchor[lineage_by_target[names(anchor)] != lineage_by_target[[target]]]), numeric(1))
  names(conflict) <- populations
  composite <- anchor_weight * anchor[populations] + local_weight * local[populations] - conflict_weight * conflict[populations]
  ordered <- order(-composite, populations)
  winner <- populations[ordered[[1L]]]
  runner <- populations[ordered[[2L]]]
  margin <- composite[[winner]] - composite[[runner]]
  contrast <- anchor[[winner]] - conflict[[winner]]
  status <- if (anchor[[winner]] < minimum_anchor) "unknown_low_anchor" else if (contrast < minimum_contrast) "unknown_conflicting_anchor" else if (margin < minimum_margin) "borderline_population_margin" else "accepted"
  alternatives <- populations[composite >= composite[[winner]] - minimum_margin - 1e-12]
  decision_rows[[index]] <- data.frame(
    cell_id = cell_id, winner_population = winner, winner_lineage = lineage_by_target[[winner]],
    fine_label = if (status == "accepted") winner else NA_character_,
    alternative_populations = paste(alternatives, collapse = ";"), alternative_count = length(alternatives),
    anchor_score = anchor[[winner]], local_score = local[[winner]], conflict_score = conflict[[winner]],
    anchor_contrast = contrast, top_composite_score = composite[[winner]],
    runner_up_population = runner, runner_up_composite_score = composite[[runner]], composite_margin = margin,
    status = status, stringsAsFactors = FALSE)
}
decisions <- do.call(rbind, decision_rows)
metadata_fields <- c("cell_id", cfg$bdev_test$sample_field, cfg$bdev_test$pseudotime_field)
decisions <- merge(decisions, fixture_md[, metadata_fields, drop = FALSE], by = "cell_id", all.x = TRUE, sort = FALSE)
names(decisions)[names(decisions) == cfg$bdev_test$sample_field] <- cfg$bdev_test$reported_sample_field
stopifnot(!anyNA(decisions$sampleId), length(unique(decisions$sampleId)) == 4L)
write_tsv(decisions, file.path(dirs[["tables"]], "10_bdev_narrow_decisions.tsv.gz"))

reported <- decisions
reported$reported_label <- ifelse(reported$status == "accepted", reported$fine_label,
  ifelse(reported$status == "borderline_population_margin", "borderline", "unknown"))
label_summary <- count_summary(reported, c("sampleId", "reported_label"))
status_summary <- count_summary(decisions, c("sampleId", "status"))
accepted <- decisions[decisions$status == "accepted", , drop = FALSE]
lineage_summary <- if (nrow(accepted)) count_summary(accepted, c("sampleId", "winner_lineage")) else data.frame()
winner_summary <- count_summary(decisions, c("sampleId", "winner_lineage"))
write_tsv(label_summary, file.path(dirs[["tables"]], "11_bdev_sampleid_label_summary.tsv"))
write_tsv(status_summary, file.path(dirs[["tables"]], "12_bdev_sampleid_status_summary.tsv"))
write_tsv(lineage_summary, file.path(dirs[["tables"]], "13_bdev_sampleid_accepted_lineage_summary.tsv"))
write_tsv(winner_summary, file.path(dirs[["tables"]], "14_bdev_sampleid_forced_winner_lineage_summary.tsv"))

sample_ids <- sort(unique(as.character(decisions$sampleId)))
sample_report <- lapply(sample_ids, function(sample_id) {
  rows <- decisions[as.character(decisions$sampleId) == sample_id, , drop = FALSE]
  list(sampleId = sample_id, cells = nrow(rows), accepted = sum(rows$status == "accepted"),
    borderline = sum(rows$status == "borderline_population_margin"),
    unknown = sum(grepl("^unknown", rows$status)),
    accepted_B = sum(rows$status == "accepted" & rows$winner_lineage == "B"))
})
overall_status <- as.data.frame(table(decisions$status), stringsAsFactors = FALSE)
names(overall_status) <- c("status", "cells")
manifest_after <- wb_manifest(bundle)
stopifnot(identical(manifest_after$current_state$fingerprint, cfg$inputs$source_bundle$fingerprint))
summary <- list(
  schema_version = "bdev.lineage_anchored_narrow_bdev_summary.v1", status = "complete_diagnostic_only",
  run_identity = cfg$run_identity, fixture_cells = nrow(decisions), sample_ids = sample_report,
  overall_status = overall_status, exact_pt_preserved = TRUE,
  pt_defined = sum(!is.na(decisions$pt)), pt_missing = sum(is.na(decisions$pt)),
  populations_competing = 64L, broad_gate_used = FALSE, Bdev_used_for_gene_selection = FALSE,
  bundle_fingerprint = manifest_after$current_state$fingerprint, bundle_registry_unchanged = TRUE,
  accepted_labels_emitted = FALSE, artifacts_promoted = FALSE)
write_json(summary, file.path(dirs[["reports"]], "02_bdev_summary.json"))
write_json(list(
  schema_version = "bdev.lineage_anchored_narrow_session_receipt.v1", status = "complete_without_promotion",
  session_id = session_id, score_run_id = score_result$run$run_id,
  calibration_run_id = calibration_result$run$run_id, execution_bundle = bundle_path,
  bundle_fingerprint = manifest_after$current_state$fingerprint, registry_unchanged = TRUE,
  scworkbench_artifacts_promoted = FALSE), file.path(dirs[["manifests"]], "03_session_receipt.json"))
