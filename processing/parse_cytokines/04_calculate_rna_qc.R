#!/usr/bin/env Rscript

# Calculate and promote factual per-cell RNA QC metrics for every cell already
# imported into the parse_cytokines bundle. This stage deliberately does not
# choose thresholds, define exclusions, or derive an accepted analysis set.
#
# Usage:
#   Rscript --vanilla 04_calculate_rna_qc.R BUNDLE_DIR REPORT_DIR

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 04_calculate_rna_qc.R BUNDLE_DIR REPORT_DIR", call. = FALSE)
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
report_dir <- normalizePath(args[[2L]], mustWork = FALSE)
if (file.exists(report_dir) || dir.exists(report_dir)) {
  stop("Report destination already exists; refusing to overwrite: ", report_dir, call. = FALSE)
}
if (!identical(basename(bundle_dir), "bundle") ||
    !identical(basename(dirname(bundle_dir)), "parse_cytokines")) {
  stop("Expected bundle destination .../parse_cytokines/bundle", call. = FALSE)
}

options(width = 180, warn = 1)
suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
})

matrix_id <- "RNA_counts"
imported_cell_set <- "set:parse_cytokines"
metric_name_prefix <- "parse_rna_qc"
mitochondrial_pattern <- "^mt-"
ribosomal_pattern <- "^Rp[sl]"
# This bundle's BPCells getter currently uses a guarded dense fallback for a
# feature-by-cell slice. At 78,298 features, 400 cells remain below the
# engine's 256 MiB materialization ceiling; this affects memory batching only.
chunk_size <- suppressWarnings(as.integer(Sys.getenv("PARSE_QC_CHUNK_SIZE", unset = "400")))
created_by <- Sys.getenv("SCWORKBENCH_CREATED_BY", unset = "jtg")

if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
  stop("PARSE_QC_CHUNK_SIZE must be a positive integer", call. = FALSE)
}

log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}

write_json_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  jsonlite::write_json(
    value, temporary_path, pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null", digits = NA
  )
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish JSON report: ", path, call. = FALSE)
  }
  invisible(path)
}

write_tsv_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  utils::write.table(
    value, file = temporary_path, sep = "\t", quote = TRUE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish TSV report: ", path, call. = FALSE)
  }
  invisible(path)
}

artifact_id_for <- function(promotion, output_name) {
  rows <- promotion$promoted[promotion$promoted$output == output_name, , drop = FALSE]
  if (nrow(rows) != 1L || is.na(rows$artifact_id[[1L]]) || !nzchar(rows$artifact_id[[1L]])) {
    stop("Expected exactly one promoted artifact for output: ", output_name, call. = FALSE)
  }
  as.character(rows$artifact_id[[1L]])
}

main <- function() {
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(report_dir)) stop("Could not create report directory: ", report_dir, call. = FALSE)

  bundle <- wb_open_bundle(bundle_dir, read_only = FALSE, validate = TRUE)
  on.exit(if (!is.null(bundle)) try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)

  manifest_before <- wb_manifest(bundle)
  if (!identical(manifest_before$bundle_id, "parse_cytokines")) {
    stop("Bundle identity is not parse_cytokines", call. = FALSE)
  }
  write_json_atomic(manifest_before, file.path(report_dir, "manifest_before_qc.json"))

  matrix_info <- wb_matrix_info(bundle, matrix_id)
  stopifnot(
    nrow(matrix_info) == 1L,
    identical(as.character(matrix_info$layer[[1L]]), "counts"),
    as.integer(matrix_info$n_features[[1L]]) == 78298L,
    as.integer(matrix_info$n_cells[[1L]]) == 251599L
  )

  imported_cells_before <- as.character(wb_cells(bundle, cell_set = imported_cell_set))
  all_cells_before <- as.character(wb_cells(bundle, cell_set = "all"))
  stopifnot(
    identical(imported_cells_before, all_cells_before),
    length(all_cells_before) == 251599L,
    !anyNA(all_cells_before),
    !anyDuplicated(all_cells_before)
  )
  cell_sets_before <- sort(as.character(wb_cell_sets(bundle)))
  if (any(grepl("qc_(accepted|excluded)", cell_sets_before, ignore.case = TRUE))) {
    stop("A QC accepted/excluded cell set already exists; refusing ambiguous rerun", call. = FALSE)
  }

  artifacts_before <- wb_art_list(bundle)
  existing_qc <- artifacts_before[
    artifacts_before$source_op == "rna_qc_metrics" &
      grepl(metric_name_prefix, artifacts_before$name, fixed = TRUE),
    , drop = FALSE
  ]
  if (nrow(existing_qc)) {
    expected_types <- c("cell_vector_family", "feature_table", "obs_table")
    if (nrow(existing_qc) != 3L || !setequal(existing_qc$type, expected_types) ||
        anyDuplicated(existing_qc$type)) {
      stop("Existing RNA QC artifacts are incomplete or ambiguous; refusing to recompute", call. = FALSE)
    }
  }

  features <- as.character(wb_feature_ids(bundle, matrix_id))
  stopifnot(length(features) == 78298L, !anyNA(features), !anyDuplicated(features))
  mitochondrial_features <- features[grepl(mitochondrial_pattern, features, perl = TRUE)]
  ribosomal_features <- features[grepl(ribosomal_pattern, features, perl = TRUE)]
  if (!length(mitochondrial_features) || !length(ribosomal_features)) {
    stop("One or more declared QC feature patterns matched no features", call. = FALSE)
  }
  feature_match_preview <- rbind(
    data.frame(
      feature_id = mitochondrial_features,
      metric_group = "mitochondrial",
      pattern = mitochondrial_pattern,
      stringsAsFactors = FALSE
    ),
    data.frame(
      feature_id = ribosomal_features,
      metric_group = "ribosomal",
      pattern = ribosomal_pattern,
      stringsAsFactors = FALSE
    )
  )
  write_tsv_atomic(feature_match_preview, file.path(report_dir, "qc_feature_match_preview.tsv"))

  interface <- wb_op_interface("rna_qc_metrics", bundle = bundle)
  write_json_atomic(interface, file.path(report_dir, "rna_qc_metrics_interface.json"))
  write_json_atomic(list(
    matrix_id = matrix_id,
    input_cell_set = imported_cell_set,
    cells = length(imported_cells_before),
    mitochondrial_pattern = mitochondrial_pattern,
    mitochondrial_features = length(mitochondrial_features),
    ribosomal_pattern = ribosomal_pattern,
    ribosomal_features = length(ribosomal_features),
    ribosomal_definition_note = paste(
      "Case-sensitive mouse-symbol Rp[sl] prefix; exact membership is retained",
      "for audit and is not asserted to be a curated structural-ribosome set."
    ),
    threshold_selection = "not_performed",
    exclusion_set_created = FALSE,
    accepted_analysis_set_created = FALSE
  ), file.path(report_dir, "qc_plan.json"))

  if (!nrow(existing_qc)) {
    log_event("Calculating factual RNA QC metrics for", length(imported_cells_before), "imported cells")
    result <- wb_session_run_op(
      wb_session_new(), bundle, "rna_qc_metrics",
      inputs = list(
        matrix = matrix_id,
        cell_set = imported_cell_set,
        min_cells = 1L,
        mitochondrial_pattern = mitochondrial_pattern,
        ribosomal_pattern = ribosomal_pattern,
        pattern_ignore_case = FALSE,
        missing_feature_policy = "error",
        require_integer_counts = TRUE,
        integer_tolerance = 1e-8,
        chunk_size = chunk_size
      ),
      notes = paste(
        "Calculate factual RNA QC metrics for every cell in the imported Parse cell set.",
        "No filtering threshold, exclusion set, or accepted analysis population is selected."
      )
    )
    write_json_atomic(
      wb_session_dag(result$session, bundle = bundle, include_artifacts = TRUE, include_active = TRUE),
      file.path(report_dir, "session_dag_before_promotion.json")
    )

    log_event("Promoting the metric family and its feature/summary audit artifacts")
    promotion <- wb_session_promote_current(
      bundle, result$session,
      outputs = c("metrics", "qc_features", "summary"),
      name_prefix = metric_name_prefix,
      description = paste(
        "Factual RNA QC metrics for all imported parse_cytokines cells;",
        "threshold selection and cell exclusion were not performed"
      ),
      created_by = created_by,
      fail_fast = TRUE
    )
    write_json_atomic(promotion, file.path(report_dir, "promotion.json"))
    artifact_ids <- c(
      metrics = artifact_id_for(promotion, "metrics"),
      qc_features = artifact_id_for(promotion, "qc_features"),
      summary = artifact_id_for(promotion, "summary")
    )
  } else {
    log_event("Resuming validation from the complete promoted RNA QC artifact family")
    id_for_type <- function(type) {
      as.character(existing_qc$artifact_id[existing_qc$type == type][[1L]])
    }
    artifact_ids <- c(
      metrics = id_for_type("cell_vector_family"),
      qc_features = id_for_type("feature_table"),
      summary = id_for_type("obs_table")
    )
    write_json_atomic(list(
      status = "reused_complete_promoted_family",
      reason = paste(
        "A prior run promoted all three artifacts and stopped only because",
        "post-promotion validation treated wrapped table primitives as data frames."
      ),
      artifact_ids = as.list(artifact_ids),
      recomputed = FALSE
    ), file.path(report_dir, "promotion_resume.json"))
  }

  metrics <- wb_art_load(bundle, artifact_ids[["metrics"]], expected_kind = "cell_vector_family")
  qc_features_primitive <- wb_art_load(bundle, artifact_ids[["qc_features"]], expected_kind = "feature_table")
  qc_summary_primitive <- wb_art_load(bundle, artifact_ids[["summary"]], expected_kind = "obs_table")
  qc_features <- qc_features_primitive$table
  qc_summary <- qc_summary_primitive$table

  expected_metric_keys <- c(
    "total_counts", "detected_features", "mitochondrial_counts",
    "mitochondrial_percent", "ribosomal_counts", "ribosomal_percent"
  )
  stopifnot(setequal(unique(as.character(metrics$key)), expected_metric_keys))
  for (key in expected_metric_keys) {
    rows <- metrics$key == key
    stopifnot(
      sum(rows) == length(imported_cells_before),
      identical(as.character(metrics$cell_id[rows]), imported_cells_before),
      all(is.finite(as.numeric(metrics$value[rows]))),
      all(as.numeric(metrics$value[rows]) >= 0)
    )
  }
  total_counts <- as.numeric(metrics$value[metrics$key == "total_counts"])
  stopifnot(all(total_counts >= 500))

  stopifnot(
    nrow(qc_features) == nrow(feature_match_preview),
    setequal(as.character(qc_features$feature_id), feature_match_preview$feature_id),
    setequal(as.character(qc_summary$metric), expected_metric_keys),
    all(as.integer(qc_summary$n_cells) == length(imported_cells_before)),
    identical(qc_summary_primitive$meta$threshold_selection, "not_performed")
  )
  write_tsv_atomic(qc_features, file.path(report_dir, "qc_feature_membership.tsv"))
  write_tsv_atomic(qc_summary, file.path(report_dir, "qc_metric_summary.tsv"))

  imported_cells_after <- as.character(wb_cells(bundle, cell_set = imported_cell_set))
  all_cells_after <- as.character(wb_cells(bundle, cell_set = "all"))
  cell_sets_after <- sort(as.character(wb_cell_sets(bundle)))
  stopifnot(
    identical(imported_cells_after, imported_cells_before),
    identical(all_cells_after, all_cells_before),
    identical(cell_sets_after, cell_sets_before)
  )

  certificates <- lapply(artifact_ids, function(artifact_id) {
    wb_entity_certificate(bundle, "artifact", artifact_id)
  })
  write_json_atomic(certificates, file.path(report_dir, "entity_certificates.json"))
  write_json_atomic(list(
    imported_cell_set = imported_cell_set,
    imported_cells_before = length(imported_cells_before),
    imported_cells_after = length(imported_cells_after),
    imported_membership_unchanged = TRUE,
    cell_sets_before = cell_sets_before,
    cell_sets_after = cell_sets_after,
    cell_set_registry_unchanged = TRUE,
    exclusion_set_created = FALSE,
    accepted_analysis_set_created = FALSE,
    threshold_selection = "not_performed"
  ), file.path(report_dir, "cell_universe_check.json"))
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))

  log_event("Finalizing the augmented bundle for recipient validation")
  readiness <- wb_finalize_recipient_bundle(bundle)
  stopifnot(isTRUE(readiness$ready))
  write_json_atomic(readiness, file.path(report_dir, "recipient_readiness.json"))
  wb_close_bundle(bundle)
  bundle <- NULL

  validation <- wb_bundle_validation_report(bundle_dir, level = "recipient", progress = TRUE)
  write_json_atomic(validation, file.path(report_dir, "validation.json"))
  stopifnot(isTRUE(validation$ok))

  # The immediately preceding recipient report already performed the complete
  # validation gate. Reopen without repeating that expensive full audit solely
  # to capture the finalized manifest and membership.
  recipient <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
  manifest_after <- wb_manifest(recipient)
  final_cells <- as.character(wb_cells(recipient, cell_set = imported_cell_set))
  final_sets <- sort(as.character(wb_cell_sets(recipient)))
  wb_close_bundle(recipient)
  stopifnot(
    identical(final_cells, imported_cells_before),
    identical(final_sets, cell_sets_before)
  )
  write_json_atomic(manifest_after, file.path(report_dir, "manifest.json"))

  completion <- list(
    status = "complete",
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    bundle_id = "parse_cytokines",
    bundle_dir = bundle_dir,
    matrix_id = matrix_id,
    imported_cell_set = imported_cell_set,
    imported_cells = length(imported_cells_before),
    artifact_ids = as.list(artifact_ids),
    metric_keys = expected_metric_keys,
    mitochondrial_pattern = mitochondrial_pattern,
    mitochondrial_features = length(mitochondrial_features),
    ribosomal_pattern = ribosomal_pattern,
    ribosomal_features = length(ribosomal_features),
    threshold_selection = "not_performed",
    exclusion_set_created = FALSE,
    accepted_analysis_set_created = FALSE,
    imported_membership_unchanged = TRUE,
    bundle_state_fingerprint = manifest_after$current_state$fingerprint,
    validation_status = validation$status
  )
  write_json_atomic(completion, file.path(report_dir, "completion.json"))
  log_event("Complete. No cells were removed and no analysis cell set was derived.")
}

main()
