#!/usr/bin/env Rscript

# Recovery-only finalizer for a state run whose atomic promotion completed but
# whose report serialization stopped on a non-JSON wb_session object. It never
# reruns analysis or promotion and refuses any bundle/report shape other than
# that exact recoverable state.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 07_finalize_signature_states_report.R BUNDLE_DIR REPORT_DIR SESSION_WORKSPACE", call. = FALSE)
}
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
report_dir <- normalizePath(args[[2L]], mustWork = TRUE)
session_workspace <- normalizePath(args[[3L]], mustWork = TRUE)

options(width = 180, warn = 1)
suppressPackageStartupMessages({ library(scworkbench); library(jsonlite) })

write_json_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  jsonlite::write_json(value, temporary_path, pretty = TRUE, auto_unbox = TRUE,
                       null = "null", na = "null", digits = NA)
  if (!file.rename(temporary_path, path)) stop("Could not publish JSON: ", path, call. = FALSE)
  invisible(path)
}

write_tsv_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  utils::write.table(value, temporary_path, sep = "\t", quote = TRUE,
                     row.names = FALSE, col.names = TRUE, na = "NA")
  if (!file.rename(temporary_path, path)) stop("Could not publish TSV: ", path, call. = FALSE)
  invisible(path)
}

main <- function() {
  stopifnot(
    file.exists(file.path(report_dir, "manifest_before.json")),
    file.exists(file.path(report_dir, "session_dag_before_promotion.json")),
    !file.exists(file.path(report_dir, "completion.json"))
  )
  bundle <- wb_open_bundle(bundle_dir, read_only = FALSE, validate = FALSE)
  on.exit(if (!is.null(bundle)) try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  imported_cells <- as.character(wb_cells(bundle, cell_set = "set:parse_cytokines"))
  analysis_cells <- as.character(wb_cells(bundle, cell_set = "set:parse_qc_control_single"))
  expected_sets <- sort(c("set:parse_cytokines", "set:parse_qc_accepted",
                          "set:parse_qc_control_single", "set:parse_qc_excluded"))
  stopifnot(
    length(imported_cells) == 251599L,
    length(analysis_cells) == 154209L,
    identical(sort(as.character(wb_cell_sets(bundle))), expected_sets),
    "RNA_log1p_cp10k" %in% wb_matrices(bundle)$matrix_id,
    "parse.signature_pca" %in% wb_spaces(bundle)$space_id
  )
  artifacts <- wb_art_list(bundle)
  promoted <- artifacts[grepl("parse_signature_states", artifacts$name, fixed = TRUE), , drop = FALSE]
  stopifnot(nrow(promoted) == 7L)
  cluster_row <- promoted[promoted$type == "cell_vector" & promoted$source_op == "leiden", , drop = FALSE]
  stopifnot(nrow(cluster_row) == 1L)
  cluster_artifact_id <- as.character(cluster_row$artifact_id[[1L]])

  metadata <- wb_meta(
    bundle, cells = analysis_cells,
    cols = list(wb_ref_metadata_field("core:cell_id"), "sample_id", "condition",
                "replicate", "treatment_category")
  )
  clusters <- wb_art_load(bundle, cluster_artifact_id, expected_kind = "cell_vector")
  index <- match(analysis_cells, as.character(clusters$cell_id))
  stopifnot(!anyNA(index))
  state <- as.character(clusters$value[index])
  stopifnot(!anyNA(state), length(unique(state)) >= 2L)
  state_by_sample <- as.data.frame(table(state = state, sample_id = metadata$sample_id),
                                   stringsAsFactors = FALSE)
  names(state_by_sample)[[3L]] <- "cells"
  state_by_sample <- state_by_sample[state_by_sample$cells > 0L, , drop = FALSE]
  sample_lookup <- unique(metadata[, c("sample_id", "condition", "replicate", "treatment_category")])
  state_by_sample <- merge(state_by_sample, sample_lookup, by = "sample_id", all.x = TRUE, sort = FALSE)
  state_by_sample <- state_by_sample[order(state_by_sample$state, state_by_sample$sample_id), , drop = FALSE]
  state_summary_raw <- aggregate(state_by_sample$cells, by = list(state = state_by_sample$state),
                                 FUN = function(x) c(cells = sum(x), samples = length(x),
                                                    min_per_sample = min(x), median_per_sample = median(x),
                                                    max_per_sample = max(x)))
  state_summary <- data.frame(
    state = state_summary_raw$state,
    cells = state_summary_raw$x[, "cells"],
    samples_present = state_summary_raw$x[, "samples"],
    min_cells_per_present_sample = state_summary_raw$x[, "min_per_sample"],
    median_cells_per_present_sample = state_summary_raw$x[, "median_per_sample"],
    max_cells_per_present_sample = state_summary_raw$x[, "max_per_sample"],
    stringsAsFactors = FALSE
  )
  state_condition <- as.data.frame(table(state = state, condition = metadata$condition),
                                   stringsAsFactors = FALSE)
  names(state_condition)[[3L]] <- "cells"
  stopifnot(sum(state_summary$cells) == length(analysis_cells))

  write_tsv_atomic(state_summary, file.path(report_dir, "state_summary.tsv"))
  write_tsv_atomic(state_by_sample, file.path(report_dir, "state_by_sample.tsv"))
  write_tsv_atomic(state_condition, file.path(report_dir, "state_by_condition.tsv"))
  write_json_atomic(list(
    recovery_status = "promotion_completed_report_serialization_recovered",
    promoted = promoted[, intersect(c("artifact_id", "type", "name", "source_op", "run_id", "snapshot_id"), names(promoted)), drop = FALSE],
    durable_session_workspace = session_workspace
  ), file.path(report_dir, "promotion.json"))
  write_json_atomic(list(
    state_method = "Leiden on 30-PC representation",
    hvg_learning_population = "QC-accepted control cells only",
    pca_and_clustering_population = "set:parse_qc_control_single",
    combined_cells_used = FALSE,
    bdev_used = FALSE,
    hvg_count = 2000L,
    pca_components = 30L,
    leiden_k = 30L,
    leiden_resolution = 0.4,
    seed = 1729L,
    normalized_matrix = "RNA_log1p_cp10k",
    pca_space = "parse.signature_pca",
    cluster_artifact_id = cluster_artifact_id,
    durable_session_workspace = session_workspace,
    compute_profile = list(profile = "cluster", scworkbench.max_dense_materialize_mb = 32768,
                           scworkbench.feature_plan_chunk_size = 2000L)
  ), file.path(report_dir, "state_analysis_parameters.json"))

  readiness <- wb_finalize_recipient_bundle(bundle)
  stopifnot(isTRUE(readiness$ready))
  write_json_atomic(readiness, file.path(report_dir, "recipient_readiness.json"))
  wb_close_bundle(bundle)
  bundle <- NULL
  validation <- wb_bundle_validation_report(bundle_dir, level = "recipient", progress = TRUE)
  stopifnot(isTRUE(validation$ok))
  write_json_atomic(validation, file.path(report_dir, "validation.json"))
  recipient <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
  manifest <- wb_manifest(recipient)
  stopifnot(
    identical(as.character(wb_cells(recipient, cell_set = "set:parse_cytokines")), imported_cells),
    identical(as.character(wb_cells(recipient, cell_set = "set:parse_qc_control_single")), analysis_cells),
    identical(sort(as.character(wb_cell_sets(recipient))), expected_sets)
  )
  wb_close_bundle(recipient)
  write_json_atomic(manifest, file.path(report_dir, "manifest.json"))
  write_json_atomic(list(
    status = "complete",
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    bundle_id = "parse_cytokines",
    imported_cells_retained = length(imported_cells),
    analysis_cells = length(analysis_cells),
    control_cells_for_hvg = sum(metadata$treatment_category == "control"),
    states = nrow(state_summary),
    cluster_artifact_id = cluster_artifact_id,
    normalized_matrix = "RNA_log1p_cp10k",
    pca_space = "parse.signature_pca",
    combined_cells_used = FALSE,
    bdev_used = FALSE,
    bundle_state_fingerprint = manifest$current_state$fingerprint,
    validation_status = validation$status,
    report_recovery = "promotion completed before wb_session JSON serialization error"
  ), file.path(report_dir, "completion.json"))
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
}

main()
