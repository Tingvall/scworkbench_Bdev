#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config")) stop("Usage: 01_reduce_and_compare.R --config <config.yaml>")
suppressPackageStartupMessages({ library(scworkbench); library(yaml); library(jsonlite); library(digest) })

cfg <- yaml::read_yaml(normalizePath(args[[2L]], mustWork = TRUE))
root <- cfg$storage$results_root
intermediate <- cfg$storage$intermediate_root
dirs <- setNames(file.path(root, c("01_tables", "03_reports", "04_logs", "06_manifests")), c("tables", "reports", "logs", "manifests"))
session_root <- file.path(intermediate, "05_objects", "01_sessions")
invisible(lapply(c(dirs, sessions = session_root), dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_raw_ucell_comparison.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
validate_signal <- function(x, cells, keys) {
  stopifnot(inherits(x, "wb_primitive"), identical(x$kind, "cell_vector_family"),
    length(unique(x$cell_id)) == cells, length(unique(x$key)) == keys,
    length(x$value) == cells * keys, all(is.finite(x$value)))
}
top_report <- function(decision, sample_field, unknown_label) {
  sample_ids <- sort(unique(as.character(decision[[sample_field]])))
  lapply(sample_ids, function(sample_id) {
    rows <- decision[as.character(decision[[sample_field]]) == sample_id, , drop = FALSE]
    assigned <- rows[rows$label != unknown_label, , drop = FALSE]
    counts <- sort(table(assigned$label), decreasing = TRUE)
    top <- head(counts, 3L)
    list(sampleId = sample_id, cells = nrow(rows), assigned = nrow(assigned), unassigned = nrow(rows) - nrow(assigned),
      top3 = lapply(seq_along(top), function(i) list(label = names(top)[[i]], cells = as.integer(top[[i]]),
        fraction_all = as.numeric(top[[i]]) / nrow(rows), fraction_assigned = as.numeric(top[[i]]) / max(1L, nrow(assigned)))))
  })
}

stopifnot(identical(sha256(cfg$inputs$user47_raw_signal$path), cfg$inputs$user47_raw_signal$sha256),
  identical(sha256(cfg$inputs$generated64_raw_signal$path), cfg$inputs$generated64_raw_signal$sha256),
  identical(sha256(cfg$inputs$generated64_module_map$path), cfg$inputs$generated64_module_map$sha256),
  identical(sha256(file.path(cfg$inputs$execution_bundle$path, "manifest.json")), cfg$inputs$execution_bundle$manifest_sha256))
fixture <- read.delim(cfg$inputs$fixture_cells$path, check.names = FALSE)
stopifnot(nrow(fixture) == as.integer(cfg$fixture$expected_cells), !anyDuplicated(fixture$cell_id))

user_signal <- readRDS(cfg$inputs$user47_raw_signal$path)
validate_signal(user_signal, nrow(fixture), 47L)
generated_all <- readRDS(cfg$inputs$generated64_raw_signal$path)
validate_signal(generated_all, nrow(fixture), 128L)
module_map <- read.delim(cfg$inputs$generated64_module_map$path, check.names = FALSE)
local_map <- module_map[module_map$module == "local", , drop = FALSE]
stopifnot(nrow(local_map) == 64L, !anyDuplicated(local_map$geneset_id), !anyDuplicated(local_map$target))
keep <- generated_all$key %in% local_map$geneset_id
generated_signal <- list(kind = "cell_vector_family",
  key = unname(local_map$target[match(generated_all$key[keep], local_map$geneset_id)]),
  cell_id = generated_all$cell_id[keep], value = generated_all$value[keep],
  value_type = generated_all$value_type, meta = list(source = "generated64_local_modules_only"))
class(generated_signal) <- "wb_primitive"
validate_signal(generated_signal, nrow(fixture), 64L)

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
metadata <- fixture_md[, c("cell_id", cfg$fixture$sample_field, cfg$fixture$pseudotime_field), drop = FALSE]
names(metadata)[names(metadata) == cfg$fixture$sample_field] <- cfg$fixture$reported_sample_field

session <- wb_session_create(bundle, workspace_root = session_root, session_id = "raw_ucell_comparison_v1")
reduce_signal <- function(session, signal, note) {
  wb_session_run_op(session, bundle, "scwbgeneprograms::cell_vector_family_reduce", inputs = list(
    signal = signal, cell_ids = fixture$cell_id, min_cells = 10L,
    rule = cfg$reduction$rule, min_score = as.numeric(cfg$reduction$minimum_raw_score),
    min_margin = as.numeric(cfg$reduction$minimum_raw_margin), unknown_label = cfg$reduction$unknown_label,
    ties = cfg$reduction$ties, tie_tol = as.numeric(cfg$reduction$tie_tolerance)), notes = note)
}
user_result <- reduce_signal(session, user_signal, "Raw-UCell reduction of the user-supplied 47 equal-size gene sets.")
generated_result <- reduce_signal(user_result$session, generated_signal, "Raw-UCell reduction of the preceding generated 64 local modules; anchor modules excluded.")

user_decision <- merge(user_result$run$outputs$decision_table$table, metadata, by = "cell_id", all.x = TRUE, sort = FALSE)
generated_decision <- merge(generated_result$run$outputs$decision_table$table, metadata, by = "cell_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(user_decision) == nrow(fixture), nrow(generated_decision) == nrow(fixture),
  !anyNA(user_decision[[cfg$fixture$reported_sample_field]]), !anyNA(generated_decision[[cfg$fixture$reported_sample_field]]))
write_tsv(user_decision, file.path(dirs[["tables"]], "01_user47_raw_ucell_decisions.tsv"))
write_tsv(generated_decision, file.path(dirs[["tables"]], "02_generated64_local_raw_ucell_decisions.tsv"))

user_report <- top_report(user_decision, cfg$fixture$reported_sample_field, cfg$reduction$unknown_label)
generated_report <- top_report(generated_decision, cfg$fixture$reported_sample_field, cfg$reduction$unknown_label)
comparison_rows <- do.call(rbind, lapply(seq_along(user_report), function(i) {
  u <- user_report[[i]]; g <- generated_report[[i]]
  data.frame(sampleId = u$sampleId, user47_assigned = u$assigned,
    user47_top = if (length(u$top3)) u$top3[[1L]]$label else NA_character_,
    user47_top_cells = if (length(u$top3)) u$top3[[1L]]$cells else 0L,
    generated64_assigned = g$assigned,
    generated64_top = if (length(g$top3)) g$top3[[1L]]$label else NA_character_,
    generated64_top_cells = if (length(g$top3)) g$top3[[1L]]$cells else 0L)
}))
write_tsv(comparison_rows, file.path(dirs[["tables"]], "03_top_winner_comparison.tsv"))
manifest_after <- wb_manifest(bundle)
stopifnot(identical(manifest_after$current_state$fingerprint, manifest_before$current_state$fingerprint))
summary <- list(schema_version = "bdev.raw_ucell_geneset_comparison_summary.v1", status = "complete_diagnostic_only",
  run_identity = cfg$run_identity, cells = nrow(fixture), minimum_raw_score = cfg$reduction$minimum_raw_score,
  minimum_raw_margin = cfg$reduction$minimum_raw_margin, percentile_calibration_used = FALSE,
  user47 = user_report, generated64_local = generated_report,
  broad_or_lineage_gate_used = FALSE, generated_anchor_modules_excluded = TRUE,
  bundle_fingerprint_unchanged = TRUE, artifacts_promoted = FALSE, accepted_labels_emitted = FALSE)
write_json(summary, file.path(dirs[["reports"]], "01_raw_ucell_comparison_summary.json"))
write_json(list(schema_version = "bdev.raw_ucell_geneset_comparison_session.v1", status = "complete_without_promotion",
  session_id = "raw_ucell_comparison_v1", user47_reduction_run_id = user_result$run$run_id,
  generated64_reduction_run_id = generated_result$run$run_id, artifacts_promoted = FALSE),
  file.path(dirs[["manifests"]], "01_session_receipt.json"))

