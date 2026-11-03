#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config")) stop("Usage: 02_score_with_scworkbench.R --config <config.yaml>")
suppressPackageStartupMessages({ library(scworkbench); library(yaml); library(jsonlite) })

cfg <- yaml::read_yaml(normalizePath(args[[2L]], mustWork = TRUE))
root <- cfg$storage$results_root
intermediate <- cfg$storage$intermediate_root
dirs <- setNames(file.path(root, c("01_tables", "03_reports", "04_logs", "06_manifests")), c("tables", "reports", "logs", "manifests"))
session_root <- file.path(intermediate, "05_objects", "02_sessions")
invisible(lapply(c(dirs, sessions = session_root), dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_scworkbench_labeling.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
run_ref <- function(run, output) wb_ref_run(run$run_id, output)
read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  out <- lapply(fields, function(x) unique(x[seq.int(3L, length(x))]))
  names(out) <- vapply(fields, `[[`, character(1), 1L)
  out
}
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
stopifnot(nrow(fixture) == as.integer(cfg$fixture$expected_cells), !anyDuplicated(fixture$cell_id))
gmt_path <- file.path(dirs[["tables"]], "01_user_narrow47_genesets.gmt")
gene_sets <- read_gmt(gmt_path)
stopifnot(length(gene_sets) == 47L, all(lengths(gene_sets) == 20L), !anyDuplicated(names(gene_sets)))

wb_operation_pack_load("scwbGenePrograms")
pack_record <- wb_operation_packs()
stopifnot(nrow(pack_record) == 1L, identical(as.character(pack_record$manifest_sha256[[1L]]), cfg$runtime$operation_pack_manifest_sha256))
for (op in c("feature_set_import", "scwbgeneprograms::feature_set_overlap_qc", "scwbgeneprograms::ucell", "scwbgeneprograms::cell_vector_family_qc", "scwbgeneprograms::cell_vector_family_calibrate", "scwbgeneprograms::cell_vector_family_reduce")) {
  write_json(wb_op_interface(op), file.path(dirs[["manifests"]], paste0("interface_", gsub("[^A-Za-z0-9]+", "_", op), ".json")))
}

bundle_path <- file.path(intermediate, "05_objects", "01_execution_bundle")
bundle <- wb_open_bundle(bundle_path, validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest_before <- wb_manifest(bundle)
fixture_md <- wb_meta(bundle, cells = fixture$cell_id)
stopifnot(identical(as.character(fixture_md$cell_id), fixture$cell_id))
grouping <- wb_prim_cell_vector(stats::setNames(as.character(fixture_md[[cfg$fixture$cluster_field]]), fixture_md$cell_id), "categorical")

registered_sets <- wb_feature_sets(bundle, matrix_id = cfg$scoring$matrix_id, namespace = cfg$scoring$namespace)
import_reused <- nrow(registered_sets) == 47L
if (import_reused) {
  import_manifest <- read.delim(file.path(dirs[["tables"]], "03_feature_set_import_manifest.tsv"), check.names = FALSE)
  stopifnot(nrow(import_manifest) == 47L, sum(import_manifest$status == "imported") == 47L)
  import_result <- NULL
} else {
  stopifnot(nrow(registered_sets) == 0L, identical(manifest_before$current_state$fingerprint, cfg$inputs$source_bundle$fingerprint))
  import_session <- wb_session_create(bundle, workspace_root = session_root, session_id = "import_user_narrow47_v1")
  import_result <- wb_session_run_op(import_session, bundle, "feature_set_import", inputs = list(
    matrix = cfg$scoring$matrix_id, cell_ids = fixture$cell_id, min_cells = 10L,
    source_path = gmt_path, source_format = "gmt", namespace = cfg$scoring$namespace,
    collection = "User-supplied 47 narrow population gene sets", species = "Mus musculus",
    gene_id_type = "gene_symbol", overwrite = TRUE, strict = FALSE,
    created_by = "scworkbench_Bdev exploratory development", source_id = cfg$run_identity,
    source_version = "20261102", origin = "user_supplied_exploratory", is_readonly = TRUE),
    notes = "User-supplied population sets; import is promoted only into the isolated execution clone for native overlap QC.")
  import_manifest <- import_result$run$outputs$import_manifest$table
  stopifnot(nrow(import_manifest) == 47L, sum(import_manifest$status == "imported") == 47L)
  write_tsv(import_manifest, file.path(dirs[["tables"]], "03_feature_set_import_manifest.tsv"))
  wb_session_promote_current(bundle, import_result$session,
    outputs = c("import_summary", "import_manifest", "rejected_sets"),
    name_prefix = "user_narrow47_import", description = "Exploratory user-supplied narrow gene sets for overlap QC",
    created_by = "scworkbench_Bdev exploratory development", fail_fast = TRUE)
}
import_run_id <- if (import_reused) "reused_required_import_from_interrupted_attempt" else import_result$run$run_id

session <- wb_session_create(bundle, workspace_root = session_root, session_id = "score_user_narrow47_v1")
overlap_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::feature_set_overlap_qc", inputs = list(
  matrix = cfg$scoring$matrix_id, cell_ids = fixture$cell_id, min_cells = 10L,
  namespace = cfg$scoring$namespace, min_set_size = as.integer(cfg$scoring$minimum_mapped_genes),
  min_detected_fraction = as.numeric(cfg$scoring$minimum_detected_fraction), deduplicate_sets = TRUE),
  notes = "Default feature-set overlap QC on the bounded Bdev fixture.")
session <- overlap_result$session
overlap <- overlap_result$run$outputs$overlap_qc$table
write_tsv(overlap, file.path(dirs[["tables"]], "04_feature_set_overlap_qc.tsv"))
write_tsv(overlap_result$run$outputs$namespace_summary$table, file.path(dirs[["tables"]], "05_feature_set_overlap_summary.tsv"))
included <- as.character(overlap$feature_set[overlap$screening_status == "include"])
included_names <- sub(paste0("^", cfg$scoring$namespace, ":"), "", included)
if (!all(included_names %in% names(gene_sets))) included_names <- sub(paste0("^", cfg$scoring$namespace, "_"), "", included)
stopifnot(length(included_names) > 1L, all(included_names %in% names(gene_sets)))
included_sets <- gene_sets[included_names]

score_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::ucell", inputs = list(
  matrix = cfg$scoring$matrix_id, cell_ids = fixture$cell_id, min_cells = 10L,
  features = included_sets, name = "USER_NARROW47", max_rank = as.integer(cfg$scoring$ucell_max_rank),
  ncores = as.integer(cfg$scoring$ncores), parallel_backend = cfg$scoring$parallel_backend,
  chunk_size = as.integer(cfg$scoring$chunk_size), missing_features_policy = cfg$scoring$missing_features_policy),
  notes = "Native UCell scoring of user-supplied narrow gene sets.")
session <- score_result$session
qc_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::cell_vector_family_qc", inputs = list(
  signal = run_ref(score_result$run, "signature_scores"), grouping = grouping,
  cell_ids = fixture$cell_id, min_cells = 10L, correlation_methods = c("pearson", "spearman"),
  tie_tol = as.numeric(cfg$scoring$tie_tolerance)), notes = "Default raw score-family QC by frozen ATAC cluster.")
session <- qc_result$session
calibration_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::cell_vector_family_calibrate", inputs = list(
  signal = run_ref(score_result$run, "signature_scores"), cell_ids = fixture$cell_id,
  min_cells = 10L, method = cfg$scoring$calibration_method,
  ties_method = cfg$scoring$calibration_ties_method, na_rm = TRUE),
  notes = "Default per-program percentile calibration.")
session <- calibration_result$session
reduction_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::cell_vector_family_reduce", inputs = list(
  signal = run_ref(calibration_result$run, "calibrated_signal"), cell_ids = fixture$cell_id,
  min_cells = 10L, rule = cfg$scoring$reduction_rule,
  min_score = as.numeric(cfg$scoring$minimum_calibrated_score),
  min_margin = as.numeric(cfg$scoring$minimum_calibrated_margin),
  unknown_label = cfg$scoring$unknown_label, ties = cfg$scoring$reduction_ties,
  tie_tol = as.numeric(cfg$scoring$tie_tolerance)),
  notes = "Default argmax reduction with calibrated score and margin abstention.")

write_tsv(score_result$run$outputs$signature_summary$table, file.path(dirs[["tables"]], "06_ucell_signature_summary.tsv"))
write_tsv(qc_result$run$outputs$key_summary$table, file.path(dirs[["tables"]], "07_raw_score_key_qc.tsv"))
write_tsv(qc_result$run$outputs$pairwise_summary$table, file.path(dirs[["tables"]], "08_raw_score_pairwise_qc.tsv"))
write_tsv(qc_result$run$outputs$group_summary$table, file.path(dirs[["tables"]], "09_raw_score_cluster_qc.tsv"))
write_tsv(calibration_result$run$outputs$calibration_summary$table, file.path(dirs[["tables"]], "10_calibration_summary.tsv"))
decision <- reduction_result$run$outputs$decision_table$table
decision$population <- sub("^USER_NARROW47_", "", as.character(decision$label))
decision$population[decision$label == cfg$scoring$unknown_label] <- cfg$scoring$unknown_label
metadata <- fixture_md[, c("cell_id", cfg$fixture$sample_field, cfg$fixture$pseudotime_field), drop = FALSE]
names(metadata)[names(metadata) == cfg$fixture$sample_field] <- cfg$fixture$reported_sample_field
decision <- merge(decision, metadata, by = "cell_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(decision) == nrow(fixture), !anyNA(decision[[cfg$fixture$reported_sample_field]]))
write_tsv(decision, file.path(dirs[["tables"]], "11_cell_label_decisions.tsv"))
summary_table <- count_summary(decision, c(cfg$fixture$reported_sample_field, "population"))
write_tsv(summary_table, file.path(dirs[["tables"]], "12_labels_by_sampleid.tsv"))

sample_ids <- sort(unique(as.character(decision[[cfg$fixture$reported_sample_field]])))
sample_report <- lapply(sample_ids, function(sample_id) {
  rows <- decision[as.character(decision[[cfg$fixture$reported_sample_field]]) == sample_id, , drop = FALSE]
  accepted <- rows[rows$population != cfg$scoring$unknown_label, , drop = FALSE]
  counts <- sort(table(accepted$population), decreasing = TRUE)
  top <- head(counts, 3L)
  list(sampleId = sample_id, cells = nrow(rows), assigned = nrow(accepted), unassigned = nrow(rows) - nrow(accepted),
    top3 = lapply(seq_along(top), function(i) list(label = names(top)[[i]], cells = as.integer(top[[i]]),
      fraction_all = as.numeric(top[[i]]) / nrow(rows), fraction_assigned = as.numeric(top[[i]]) / max(1L, nrow(accepted)))))
})
manifest_after <- wb_manifest(bundle)
report <- list(schema_version = "bdev.user_narrow47_default_labeling_summary.v1", status = "complete_diagnostic_only",
  run_identity = cfg$run_identity, cells = nrow(decision), populations_supplied = 47L,
  populations_included = length(included_names), sample_ids = sample_report,
  reduction = list(minimum_calibrated_score = cfg$scoring$minimum_calibrated_score,
    minimum_calibrated_margin = cfg$scoring$minimum_calibrated_margin),
  broad_or_lineage_gate_used = FALSE, source_bundle_mutated = FALSE,
  execution_clone_fingerprint_after_required_import = manifest_after$current_state$fingerprint,
  score_artifacts_promoted = FALSE, accepted_labels_emitted = FALSE)
write_json(report, file.path(dirs[["reports"]], "01_labeling_summary.json"))
write_json(list(schema_version = "bdev.user_narrow47_session_receipt.v1", status = "complete_without_score_promotion",
  required_import_promoted_to_execution_clone = TRUE, import_run_id = import_run_id,
  score_run_id = score_result$run$run_id, qc_run_id = qc_result$run$run_id,
  calibration_run_id = calibration_result$run$run_id, reduction_run_id = reduction_result$run$run_id,
  score_artifacts_promoted = FALSE), file.path(dirs[["manifests"]], "02_session_receipt.json"))
