#!/usr/bin/env Rscript

# Apply the user-selected Parse RNA QC rule to the promoted scWorkbench QC
# metrics. All imported cells remain in set:parse_cytokines and in the bundle;
# filtering is represented by separate durable exclusion and accepted-analysis
# cell sets.
#
# The exclusion preview, criterion overlaps, and per-sample effects are written
# before either durable set is promoted.
#
# Usage:
#   Rscript --vanilla 05_apply_rna_qc_filter.R BUNDLE_DIR REPORT_DIR

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 05_apply_rna_qc_filter.R BUNDLE_DIR REPORT_DIR", call. = FALSE)
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

imported_cell_set <- "set:parse_cytokines"
excluded_cell_set <- "set:parse_qc_excluded"
accepted_cell_set <- "set:parse_qc_accepted"
minimum_total_counts <- 1000
minimum_detected_features <- 500
maximum_mitochondrial_percent <- 10
created_by <- Sys.getenv("SCWORKBENCH_CREATED_BY", unset = "jtg")

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

write_tsv_gz_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  connection <- gzfile(temporary_path, open = "wt", compression = 6L)
  connection_open <- TRUE
  on.exit({
    if (isTRUE(connection_open)) try(close(connection), silent = TRUE)
    if (file.exists(temporary_path)) unlink(temporary_path)
  }, add = TRUE)
  utils::write.table(
    value, file = connection, sep = "\t", quote = TRUE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
  close(connection)
  connection_open <- FALSE
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish compressed TSV report: ", path, call. = FALSE)
  }
  invisible(path)
}

write_lines_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  writeLines(value, temporary_path, useBytes = TRUE)
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish text report: ", path, call. = FALSE)
  }
  invisible(path)
}

markdown_table <- function(value, digits = 6L) {
  if (!nrow(value)) return("_No rows._")
  format_value <- function(x) {
    if (is.numeric(x)) {
      out <- ifelse(
        is.na(x), "NA",
        format(signif(x, digits), trim = TRUE, scientific = FALSE, big.mark = ",")
      )
    } else {
      out <- ifelse(is.na(x), "NA", as.character(x))
    }
    gsub("|", "\\\\|", out, fixed = TRUE)
  }
  formatted <- lapply(value, format_value)
  header <- paste0("| ", paste(names(value), collapse = " | "), " |")
  divider <- paste0("| ", paste(rep("---", ncol(value)), collapse = " | "), " |")
  rows <- vapply(seq_len(nrow(value)), function(i) {
    paste0("| ", paste(vapply(formatted, `[[`, character(1L), i), collapse = " | "), " |")
  }, character(1L))
  paste(c(header, divider, rows), collapse = "\n")
}

artifact_row <- function(artifacts, type) {
  rows <- artifacts[
    artifacts$source_op == "rna_qc_metrics" &
      artifacts$type == type &
      grepl("parse_rna_qc", artifacts$name, fixed = TRUE),
    , drop = FALSE
  ]
  if (nrow(rows) != 1L) {
    stop("Expected exactly one current parse_rna_qc artifact of type ", type, "; found ", nrow(rows), call. = FALSE)
  }
  rows
}

artifact_id_for <- function(promotion, output_name) {
  rows <- promotion$promoted[promotion$promoted$output == output_name, , drop = FALSE]
  if (nrow(rows) != 1L || is.na(rows$artifact_id[[1L]]) || !nzchar(rows$artifact_id[[1L]])) {
    stop("Expected exactly one promoted artifact for output: ", output_name, call. = FALSE)
  }
  as.character(rows$artifact_id[[1L]])
}

metric_vector <- function(metrics, key, expected_cells) {
  rows <- metrics$key == key
  values <- as.numeric(metrics$value[rows])
  names(values) <- as.character(metrics$cell_id[rows])
  stopifnot(
    identical(names(values), expected_cells),
    length(values) == length(expected_cells),
    all(is.finite(values)),
    all(values >= 0)
  )
  values
}

summarize_population_metrics <- function(metric_values, populations) {
  probabilities <- c(0, 0.05, 0.25, 0.50, 0.75, 0.95, 1)
  labels <- c("min", "p05", "p25", "p50", "p75", "p95", "max")
  do.call(rbind, lapply(names(populations), function(population) {
    indices <- populations[[population]]
    do.call(rbind, lapply(names(metric_values), function(metric) {
      values <- metric_values[[metric]][indices]
      quantiles <- as.numeric(stats::quantile(values, probabilities, names = FALSE, type = 7))
      data.frame(
        population = population,
        cells = length(indices),
        metric = metric,
        stats::setNames(as.list(quantiles), labels),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }))
  }))
}

main <- function() {
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(report_dir)) stop("Could not create report directory: ", report_dir, call. = FALSE)

  bundle <- wb_open_bundle(bundle_dir, read_only = FALSE, validate = TRUE)
  on.exit(if (!is.null(bundle)) try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)

  manifest_before <- wb_manifest(bundle)
  stopifnot(identical(manifest_before$bundle_id, "parse_cytokines"))
  write_json_atomic(manifest_before, file.path(report_dir, "manifest_before_filtering.json"))
  write_json_atomic(
    wb_op_interface("cell_set_derive", bundle = bundle),
    file.path(report_dir, "cell_set_derive_interface.json")
  )

  imported_cells_before <- as.character(wb_cells(bundle, cell_set = imported_cell_set))
  all_cells_before <- as.character(wb_cells(bundle, cell_set = "all"))
  stopifnot(
    identical(imported_cells_before, all_cells_before),
    length(imported_cells_before) == 251599L,
    !anyNA(imported_cells_before),
    !anyDuplicated(imported_cells_before)
  )

  cell_sets_before <- sort(as.character(wb_cell_sets(bundle)))
  if (any(c(excluded_cell_set, accepted_cell_set) %in% cell_sets_before)) {
    stop("One or more target QC cell sets already exist; refusing to overwrite", call. = FALSE)
  }

  artifacts <- wb_art_list(bundle)
  metrics_row <- artifact_row(artifacts, "cell_vector_family")
  metrics_artifact_id <- as.character(metrics_row$artifact_id[[1L]])
  metrics_ref <- paste0("artifact:", metrics_artifact_id)
  metrics <- wb_art_load(bundle, metrics_artifact_id, expected_kind = "cell_vector_family")

  metric_values <- list(
    total_counts = metric_vector(metrics, "total_counts", imported_cells_before),
    detected_features = metric_vector(metrics, "detected_features", imported_cells_before),
    mitochondrial_percent = metric_vector(metrics, "mitochondrial_percent", imported_cells_before)
  )
  stopifnot(all(metric_values$total_counts >= 500))

  metadata <- wb_meta(
    bundle,
    cells = imported_cells_before,
    cols = list(wb_ref_metadata_field("core:cell_id"), "sample_id", "condition", "replicate")
  )
  stopifnot(
    identical(as.character(metadata$cell_id), imported_cells_before),
    !anyNA(metadata$sample_id),
    length(unique(metadata$sample_id)) == 96L
  )

  fail_total_counts <- metric_values$total_counts < minimum_total_counts
  fail_detected_features <- metric_values$detected_features < minimum_detected_features
  fail_mitochondrial_percent <- metric_values$mitochondrial_percent > maximum_mitochondrial_percent
  proposed_exclusion <- fail_total_counts | fail_detected_features | fail_mitochondrial_percent
  accepted_flag <- !proposed_exclusion
  if (!any(proposed_exclusion) || !any(accepted_flag)) {
    stop("Selected QC rule produced an empty excluded or accepted population", call. = FALSE)
  }

  accepted_selector <- wb_sel_all(bundle, cell_set = imported_cell_set)
  accepted_selector <- wb_sel_vector_family(
    accepted_selector, metrics_ref, "total_counts",
    min = minimum_total_counts, bundle = bundle
  )
  accepted_selector <- wb_sel_vector_family(
    accepted_selector, metrics_ref, "detected_features",
    min = minimum_detected_features, bundle = bundle
  )
  accepted_selector <- wb_sel_vector_family(
    accepted_selector, metrics_ref, "mitochondrial_percent",
    max = maximum_mitochondrial_percent, bundle = bundle
  )
  excluded_selector <- wb_sel_not(accepted_selector, bundle = bundle)
  accepted_cells <- as.character(accepted_selector)
  excluded_cells <- as.character(excluded_selector)
  stopifnot(
    identical(accepted_cells, imported_cells_before[accepted_flag]),
    identical(excluded_cells, imported_cells_before[proposed_exclusion]),
    !length(intersect(accepted_cells, excluded_cells)),
    identical(sort(c(accepted_cells, excluded_cells)), sort(imported_cells_before))
  )

  reason_for <- function(i) {
    paste(c(
      if (fail_total_counts[[i]]) "total_counts_below_1000",
      if (fail_detected_features[[i]]) "detected_features_below_500",
      if (fail_mitochondrial_percent[[i]]) "mitochondrial_percent_above_10"
    ), collapse = ";")
  }
  excluded_indices <- which(proposed_exclusion)
  proposed_exclusions <- data.frame(
    cell_id = imported_cells_before[excluded_indices],
    sample_id = as.character(metadata$sample_id[excluded_indices]),
    condition = as.character(metadata$condition[excluded_indices]),
    replicate = as.character(metadata$replicate[excluded_indices]),
    total_counts = unname(metric_values$total_counts[excluded_indices]),
    detected_features = unname(metric_values$detected_features[excluded_indices]),
    mitochondrial_percent = unname(metric_values$mitochondrial_percent[excluded_indices]),
    fail_total_counts = fail_total_counts[excluded_indices],
    fail_detected_features = fail_detected_features[excluded_indices],
    fail_mitochondrial_percent = fail_mitochondrial_percent[excluded_indices],
    exclusion_reason = vapply(excluded_indices, reason_for, character(1L)),
    stringsAsFactors = FALSE
  )

  criterion_summary <- data.frame(
    population = c(
      "imported", "fail_total_counts", "fail_detected_features",
      "fail_mitochondrial_percent", "excluded_union", "accepted_all_criteria"
    ),
    rule = c(
      "set:parse_cytokines",
      "total_counts < 1000",
      "detected_features < 500",
      "mitochondrial_percent > 10",
      "failure of one or more criteria",
      "total_counts >= 1000 AND detected_features >= 500 AND mitochondrial_percent <= 10"
    ),
    cells = c(
      length(imported_cells_before), sum(fail_total_counts), sum(fail_detected_features),
      sum(fail_mitochondrial_percent), length(excluded_cells), length(accepted_cells)
    ),
    stringsAsFactors = FALSE
  )
  criterion_summary$percent_of_imported <- 100 * criterion_summary$cells / length(imported_cells_before)

  overlap_observed <- aggregate(
    rep(1L, length(imported_cells_before)),
    by = list(
      fail_total_counts = fail_total_counts,
      fail_detected_features = fail_detected_features,
      fail_mitochondrial_percent = fail_mitochondrial_percent
    ),
    FUN = sum
  )
  names(overlap_observed)[[4L]] <- "cells"
  overlap_grid <- expand.grid(
    fail_total_counts = c(FALSE, TRUE),
    fail_detected_features = c(FALSE, TRUE),
    fail_mitochondrial_percent = c(FALSE, TRUE),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  overlap <- merge(
    overlap_grid, overlap_observed,
    by = c("fail_total_counts", "fail_detected_features", "fail_mitochondrial_percent"),
    all.x = TRUE, sort = FALSE
  )
  overlap$cells[is.na(overlap$cells)] <- 0L
  overlap$percent_of_imported <- 100 * overlap$cells / length(imported_cells_before)
  overlap$status <- ifelse(
    overlap$fail_total_counts | overlap$fail_detected_features | overlap$fail_mitochondrial_percent,
    "excluded", "accepted"
  )
  overlap <- overlap[order(overlap$status, -overlap$cells), , drop = FALSE]
  rownames(overlap) <- NULL

  audit_rows <- data.frame(
    sample_id = as.character(metadata$sample_id),
    condition = as.character(metadata$condition),
    replicate = as.character(metadata$replicate),
    imported = 1L,
    fail_total_counts = as.integer(fail_total_counts),
    fail_detected_features = as.integer(fail_detected_features),
    fail_mitochondrial_percent = as.integer(fail_mitochondrial_percent),
    excluded_union = as.integer(proposed_exclusion),
    accepted = as.integer(accepted_flag),
    stringsAsFactors = FALSE
  )
  per_sample <- aggregate(
    audit_rows[, c(
      "imported", "fail_total_counts", "fail_detected_features",
      "fail_mitochondrial_percent", "excluded_union", "accepted"
    )],
    by = audit_rows[, c("sample_id", "condition", "replicate")],
    FUN = sum
  )
  per_sample$percent_excluded <- 100 * per_sample$excluded_union / per_sample$imported
  per_sample$percent_accepted <- 100 * per_sample$accepted / per_sample$imported
  per_sample <- per_sample[order(per_sample$sample_id), , drop = FALSE]
  rownames(per_sample) <- NULL

  per_condition <- aggregate(
    audit_rows[, c(
      "imported", "fail_total_counts", "fail_detected_features",
      "fail_mitochondrial_percent", "excluded_union", "accepted"
    )],
    by = audit_rows[, "condition", drop = FALSE],
    FUN = sum
  )
  per_condition$percent_excluded <- 100 * per_condition$excluded_union / per_condition$imported
  per_condition$percent_accepted <- 100 * per_condition$accepted / per_condition$imported
  per_condition <- per_condition[order(per_condition$condition), , drop = FALSE]
  rownames(per_condition) <- NULL

  population_metric_summary <- summarize_population_metrics(
    metric_values,
    populations = list(
      imported = seq_along(imported_cells_before),
      excluded = excluded_indices,
      accepted = which(accepted_flag)
    )
  )
  rownames(population_metric_summary) <- NULL

  preview_generated_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  write_json_atomic(list(
    status = "selected_by_user_pending_cell_set_promotion",
    generated_utc = preview_generated_utc,
    input_cell_set = imported_cell_set,
    input_cells = length(imported_cells_before),
    metrics_artifact_id = metrics_artifact_id,
    criteria = list(
      total_counts = list(operator = ">=", threshold = minimum_total_counts),
      detected_features = list(operator = ">=", threshold = minimum_detected_features),
      mitochondrial_percent = list(operator = "<=", threshold = maximum_mitochondrial_percent)
    ),
    combination = "AND for acceptance; failure of any criterion is exclusion",
    excluded_cells = length(excluded_cells),
    accepted_cells = length(accepted_cells),
    durable_sets_created_at_preview_time = FALSE
  ), file.path(report_dir, "threshold_selection_and_preview.json"))
  write_tsv_atomic(criterion_summary, file.path(report_dir, "criterion_summary.tsv"))
  write_tsv_atomic(overlap, file.path(report_dir, "criterion_overlap.tsv"))
  write_tsv_atomic(per_sample, file.path(report_dir, "exclusions_by_sample.tsv"))
  write_tsv_atomic(per_condition, file.path(report_dir, "exclusions_by_condition.tsv"))
  write_tsv_atomic(population_metric_summary, file.path(report_dir, "population_metric_summary.tsv"))
  write_tsv_gz_atomic(proposed_exclusions, file.path(report_dir, "proposed_exclusions.tsv.gz"))

  most_affected_samples <- per_sample[order(-per_sample$percent_excluded, per_sample$sample_id), , drop = FALSE]
  most_affected_samples <- utils::head(most_affected_samples, 15L)
  preview_report <- paste0(
    "# Parse cytokine RNA QC exclusion preview\n\n",
    "Generated before durable QC cell-set promotion: `", preview_generated_utc, "`\n\n",
    "## Selected rule\n\n",
    "A cell is accepted only when all three user-selected inclusive criteria hold: ",
    "`total_counts >= 1000`, `detected_features >= 500`, and ",
    "`mitochondrial_percent <= 10`. The mitochondrial criterion uses the percentage metric, not raw mitochondrial counts.\n\n",
    markdown_table(criterion_summary), "\n\n",
    "## Criterion overlap\n\n",
    markdown_table(overlap), "\n\n",
    "## Samples with the largest excluded percentages\n\n",
    markdown_table(most_affected_samples), "\n\n",
    "The complete proposed-exclusion list is `proposed_exclusions.tsv.gz`; complete sample and condition summaries are retained beside this report. ",
    "At the time this preview was written, no QC-derived cell set had been promoted.\n"
  )
  write_lines_atomic(preview_report, file.path(report_dir, "exclusion_preview.md"))
  log_event(
    "Exclusion preview written before durable set derivation:",
    length(excluded_cells), "excluded and", length(accepted_cells), "accepted"
  )

  excluded_description <- paste(
    "Imported Parse cells failing one or more user-selected RNA QC criteria:",
    "total_counts >= 1000, detected_features >= 500, mitochondrial_percent <= 10."
  )
  excluded_run <- wb_session_run_op(
    wb_session_new(), bundle, "cell_set_derive",
    inputs = list(
      cell_ids = excluded_cells,
      min_cells = 1L,
      target_cell_set = excluded_cell_set,
      role = "exclusion_set",
      status = "accepted",
      source = "user_selected_parse_rna_qc_thresholds",
      created_by = created_by,
      description = excluded_description,
      recipe = wb_recipe(excluded_selector),
      overwrite = FALSE
    ),
    notes = paste(
      "Materialize the exclusion complement before the accepted analysis set.",
      "All imported cells remain retained in set:parse_cytokines."
    )
  )
  write_json_atomic(
    wb_session_dag(excluded_run$session, bundle = bundle, include_artifacts = TRUE, include_active = TRUE),
    file.path(report_dir, "excluded_session_dag_before_promotion.json")
  )
  excluded_promotion <- wb_session_promote_current(
    bundle, excluded_run$session,
    outputs = c("manifest", "rule_table"),
    name_prefix = "parse_qc_excluded",
    description = excluded_description,
    created_by = created_by,
    fail_fast = TRUE
  )
  write_json_atomic(excluded_promotion, file.path(report_dir, "excluded_promotion.json"))
  excluded_artifacts <- c(
    manifest = artifact_id_for(excluded_promotion, "manifest"),
    rule_table = artifact_id_for(excluded_promotion, "rule_table")
  )
  stopifnot(identical(as.character(wb_cells(bundle, cell_set = excluded_cell_set)), excluded_cells))

  accepted_description <- paste(
    "Accepted Parse RNA analysis population meeting all user-selected RNA QC criteria:",
    "total_counts >= 1000, detected_features >= 500, mitochondrial_percent <= 10."
  )
  accepted_run <- wb_session_run_op(
    wb_session_new(), bundle, "cell_set_derive",
    inputs = list(
      cell_ids = accepted_cells,
      min_cells = 1L,
      target_cell_set = accepted_cell_set,
      role = "analysis_universe",
      status = "accepted",
      source = "user_selected_parse_rna_qc_thresholds",
      created_by = created_by,
      description = accepted_description,
      recipe = wb_recipe(accepted_selector),
      overwrite = FALSE
    ),
    notes = paste(
      "Derive the accepted analysis universe after writing and promoting the exclusion audit.",
      "Do not replace or narrow the imported cell set."
    )
  )
  write_json_atomic(
    wb_session_dag(accepted_run$session, bundle = bundle, include_artifacts = TRUE, include_active = TRUE),
    file.path(report_dir, "accepted_session_dag_before_promotion.json")
  )
  accepted_promotion <- wb_session_promote_current(
    bundle, accepted_run$session,
    outputs = c("manifest", "rule_table"),
    name_prefix = "parse_qc_accepted",
    description = accepted_description,
    created_by = created_by,
    fail_fast = TRUE
  )
  write_json_atomic(accepted_promotion, file.path(report_dir, "accepted_promotion.json"))
  accepted_artifacts <- c(
    manifest = artifact_id_for(accepted_promotion, "manifest"),
    rule_table = artifact_id_for(accepted_promotion, "rule_table")
  )

  imported_cells_after <- as.character(wb_cells(bundle, cell_set = imported_cell_set))
  all_cells_after <- as.character(wb_cells(bundle, cell_set = "all"))
  excluded_cells_after <- as.character(wb_cells(bundle, cell_set = excluded_cell_set))
  accepted_cells_after <- as.character(wb_cells(bundle, cell_set = accepted_cell_set))
  cell_sets_after <- sort(as.character(wb_cell_sets(bundle)))
  stopifnot(
    identical(imported_cells_after, imported_cells_before),
    identical(all_cells_after, all_cells_before),
    identical(excluded_cells_after, excluded_cells),
    identical(accepted_cells_after, accepted_cells),
    !length(intersect(excluded_cells_after, accepted_cells_after)),
    identical(sort(c(excluded_cells_after, accepted_cells_after)), sort(imported_cells_before)),
    all(c(excluded_cell_set, accepted_cell_set) %in% cell_sets_after),
    all(metric_values$total_counts[accepted_cells_after] >= minimum_total_counts),
    all(metric_values$detected_features[accepted_cells_after] >= minimum_detected_features),
    all(metric_values$mitochondrial_percent[accepted_cells_after] <= maximum_mitochondrial_percent),
    all(
      metric_values$total_counts[excluded_cells_after] < minimum_total_counts |
        metric_values$detected_features[excluded_cells_after] < minimum_detected_features |
        metric_values$mitochondrial_percent[excluded_cells_after] > maximum_mitochondrial_percent
    )
  )

  excluded_recipe <- wb_get_cell_set_recipe(bundle, excluded_cell_set, parse = FALSE)
  accepted_recipe <- wb_get_cell_set_recipe(bundle, accepted_cell_set, parse = FALSE)
  stopifnot(
    nrow(excluded_recipe) == 1L,
    identical(as.character(excluded_recipe$role[[1L]]), "exclusion_set"),
    identical(as.character(excluded_recipe$status[[1L]]), "accepted"),
    nrow(accepted_recipe) == 1L,
    identical(as.character(accepted_recipe$role[[1L]]), "analysis_universe"),
    identical(as.character(accepted_recipe$status[[1L]]), "accepted")
  )
  write_tsv_atomic(excluded_recipe, file.path(report_dir, "excluded_cell_set_recipe.tsv"))
  write_tsv_atomic(accepted_recipe, file.path(report_dir, "accepted_cell_set_recipe.tsv"))

  certificates <- list(
    metrics_artifact = wb_entity_certificate(bundle, "artifact", metrics_artifact_id),
    excluded_cell_set = wb_entity_certificate(bundle, "cell_set", excluded_cell_set),
    accepted_cell_set = wb_entity_certificate(bundle, "cell_set", accepted_cell_set),
    excluded_manifest_artifact = wb_entity_certificate(bundle, "artifact", excluded_artifacts[["manifest"]]),
    accepted_manifest_artifact = wb_entity_certificate(bundle, "artifact", accepted_artifacts[["manifest"]])
  )
  write_json_atomic(certificates, file.path(report_dir, "entity_certificates.json"))
  write_json_atomic(list(
    imported_cell_set = imported_cell_set,
    imported_cells_before = length(imported_cells_before),
    imported_cells_after = length(imported_cells_after),
    imported_membership_unchanged = identical(imported_cells_after, imported_cells_before),
    excluded_cell_set = excluded_cell_set,
    excluded_cells = length(excluded_cells_after),
    accepted_cell_set = accepted_cell_set,
    accepted_cells = length(accepted_cells_after),
    excluded_and_accepted_disjoint = !length(intersect(excluded_cells_after, accepted_cells_after)),
    excluded_and_accepted_partition_imported = identical(
      sort(c(excluded_cells_after, accepted_cells_after)), sort(imported_cells_before)
    ),
    cell_sets_before = cell_sets_before,
    cell_sets_after = cell_sets_after
  ), file.path(report_dir, "cell_universe_check.json"))
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))

  log_event("Finalizing the bundle and running recipient validation")
  readiness <- wb_finalize_recipient_bundle(bundle)
  stopifnot(isTRUE(readiness$ready))
  write_json_atomic(readiness, file.path(report_dir, "recipient_readiness.json"))
  wb_close_bundle(bundle)
  bundle <- NULL

  validation <- wb_bundle_validation_report(bundle_dir, level = "recipient", progress = TRUE)
  write_json_atomic(validation, file.path(report_dir, "validation.json"))
  stopifnot(isTRUE(validation$ok))

  recipient <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
  manifest_after <- wb_manifest(recipient)
  final_imported <- as.character(wb_cells(recipient, cell_set = imported_cell_set))
  final_excluded <- as.character(wb_cells(recipient, cell_set = excluded_cell_set))
  final_accepted <- as.character(wb_cells(recipient, cell_set = accepted_cell_set))
  final_sets <- sort(as.character(wb_cell_sets(recipient)))
  wb_close_bundle(recipient)
  stopifnot(
    identical(final_imported, imported_cells_before),
    identical(final_excluded, excluded_cells),
    identical(final_accepted, accepted_cells),
    identical(final_sets, cell_sets_after)
  )
  write_json_atomic(manifest_after, file.path(report_dir, "manifest.json"))

  completion <- list(
    status = "complete",
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    bundle_id = "parse_cytokines",
    bundle_dir = bundle_dir,
    input_cell_set = imported_cell_set,
    imported_cells = length(imported_cells_before),
    imported_membership_unchanged = TRUE,
    metrics_artifact_id = metrics_artifact_id,
    criteria = list(
      total_counts = list(operator = ">=", threshold = minimum_total_counts),
      detected_features = list(operator = ">=", threshold = minimum_detected_features),
      mitochondrial_percent = list(operator = "<=", threshold = maximum_mitochondrial_percent)
    ),
    excluded = list(
      cell_set = excluded_cell_set,
      cells = length(excluded_cells),
      percent_of_imported = 100 * length(excluded_cells) / length(imported_cells_before),
      artifacts = as.list(excluded_artifacts)
    ),
    accepted = list(
      cell_set = accepted_cell_set,
      cells = length(accepted_cells),
      percent_of_imported = 100 * length(accepted_cells) / length(imported_cells_before),
      artifacts = as.list(accepted_artifacts)
    ),
    partition_verified = TRUE,
    preview_written_before_set_promotion = TRUE,
    bundle_state_fingerprint = manifest_after$current_state$fingerprint,
    validation_status = validation$status
  )
  write_json_atomic(completion, file.path(report_dir, "completion.json"))
  log_event(
    "Complete. Imported cells retained:", length(final_imported),
    "; accepted analysis cells:", length(final_accepted),
    "; excluded QC cells:", length(final_excluded)
  )
}

main()
