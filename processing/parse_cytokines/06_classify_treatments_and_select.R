#!/usr/bin/env Rscript

# Classify every imported Parse sample/cell as control, single, or combined,
# materialize the classification as provenance-backed scWorkbench metadata and
# a categorical cell-vector artifact, and derive a downstream analysis set from
# QC-accepted control and single-cytokine cells. Imported and QC-derived cell
# sets are preserved unchanged.
#
# Usage:
#   Rscript --vanilla 06_classify_treatments_and_select.R \
#     BUNDLE_DIR REPORT_DIR TREATMENT_MAP_TSV

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: 06_classify_treatments_and_select.R BUNDLE_DIR REPORT_DIR TREATMENT_MAP_TSV",
    call. = FALSE
  )
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
report_dir <- normalizePath(args[[2L]], mustWork = FALSE)
treatment_map_path <- normalizePath(args[[3L]], mustWork = TRUE)
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
qc_accepted_cell_set <- "set:parse_qc_accepted"
qc_excluded_cell_set <- "set:parse_qc_excluded"
downstream_cell_set <- "set:parse_qc_control_single"
category_vector_name <- "parse_treatment_category"
included_categories <- c("control", "single")
category_order <- c("control", "single", "combined")
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

artifact_id_for <- function(promotion, run_id, output_name) {
  rows <- promotion$promoted[
    promotion$promoted$run_id == run_id &
      promotion$promoted$output == output_name &
      promotion$promoted$status == "ok",
    , drop = FALSE
  ]
  if (nrow(rows) != 1L || is.na(rows$artifact_id[[1L]]) || !nzchar(rows$artifact_id[[1L]])) {
    stop("Expected one promoted artifact for run/output: ", run_id, "/", output_name, call. = FALSE)
  }
  as.character(rows$artifact_id[[1L]])
}

count_components <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  vapply(x, function(value) {
    if (!nzchar(value)) return(0L)
    length(strsplit(value, ";", fixed = TRUE)[[1L]])
  }, integer(1L))
}

main <- function() {
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(report_dir)) stop("Could not create report directory: ", report_dir, call. = FALSE)

  treatment_map <- utils::read.delim(
    treatment_map_path, header = TRUE, sep = "\t", quote = "",
    comment.char = "", stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = c("character", "character", "integer", "character", "character")
  )
  expected_columns <- c(
    "condition", "treatment_category", "cytokine_count",
    "cytokine_components", "classification_basis"
  )
  stopifnot(
    identical(names(treatment_map), expected_columns),
    nrow(treatment_map) == 32L,
    !anyNA(treatment_map$condition),
    !anyDuplicated(treatment_map$condition),
    all(treatment_map$treatment_category %in% category_order),
    identical(unname(count_components(treatment_map$cytokine_components)), treatment_map$cytokine_count),
    all(treatment_map$cytokine_count[treatment_map$treatment_category == "control"] == 0L),
    all(treatment_map$cytokine_count[treatment_map$treatment_category == "single"] == 1L),
    all(treatment_map$cytokine_count[treatment_map$treatment_category == "combined"] > 1L)
  )
  category_condition_counts <- table(factor(treatment_map$treatment_category, levels = category_order))
  stopifnot(identical(as.integer(category_condition_counts), c(1L, 26L, 5L)))

  bundle <- wb_open_bundle(bundle_dir, read_only = FALSE, validate = TRUE)
  on.exit(if (!is.null(bundle)) try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)

  manifest_before <- wb_manifest(bundle)
  stopifnot(identical(manifest_before$bundle_id, "parse_cytokines"))
  write_json_atomic(manifest_before, file.path(report_dir, "manifest_before_treatment_selection.json"))
  write_json_atomic(
    list(
      metadata_cell_vector = wb_op_interface("metadata_cell_vector", bundle = bundle),
      cell_set_derive = wb_op_interface("cell_set_derive", bundle = bundle)
    ),
    file.path(report_dir, "operation_interfaces.json")
  )

  imported_cells_before <- as.character(wb_cells(bundle, cell_set = imported_cell_set))
  all_cells_before <- as.character(wb_cells(bundle, cell_set = "all"))
  qc_accepted_before <- as.character(wb_cells(bundle, cell_set = qc_accepted_cell_set))
  qc_excluded_before <- as.character(wb_cells(bundle, cell_set = qc_excluded_cell_set))
  cell_sets_before <- sort(as.character(wb_cell_sets(bundle)))
  stopifnot(
    identical(imported_cells_before, all_cells_before),
    length(imported_cells_before) == 251599L,
    length(qc_accepted_before) == 194215L,
    length(qc_excluded_before) == 57384L,
    !anyDuplicated(imported_cells_before),
    !length(intersect(qc_accepted_before, qc_excluded_before)),
    identical(sort(c(qc_accepted_before, qc_excluded_before)), sort(imported_cells_before)),
    !(downstream_cell_set %in% cell_sets_before)
  )

  fields_before <- wb_metadata_fields(bundle, include_core = FALSE)
  if (any(fields_before$source_name %in% c("treatment_category", "cytokine_count"))) {
    stop("Treatment classification metadata already exists; refusing to overwrite", call. = FALSE)
  }
  artifacts_before <- wb_art_list(bundle)
  if (nrow(artifacts_before) && any(
    artifacts_before$type == "cell_vector" & artifacts_before$name == category_vector_name
  )) {
    stop("Treatment category vector already exists; refusing to overwrite", call. = FALSE)
  }

  metadata <- wb_meta(
    bundle,
    cells = imported_cells_before,
    cols = list(wb_ref_metadata_field("core:cell_id"), "sample_id", "condition", "replicate")
  )
  stopifnot(
    identical(as.character(metadata$cell_id), imported_cells_before),
    !anyNA(metadata$sample_id),
    !anyNA(metadata$condition),
    length(unique(metadata$sample_id)) == 96L,
    setequal(unique(as.character(metadata$condition)), treatment_map$condition)
  )

  mapping_index <- match(as.character(metadata$condition), treatment_map$condition)
  stopifnot(!anyNA(mapping_index))
  treatment_category <- treatment_map$treatment_category[mapping_index]
  cytokine_count <- treatment_map$cytokine_count[mapping_index]
  stopifnot(
    all(treatment_category %in% category_order),
    !anyNA(treatment_category),
    !anyNA(cytokine_count)
  )

  sample_identity <- unique(data.frame(
    sample_id = as.character(metadata$sample_id),
    condition = as.character(metadata$condition),
    replicate = as.character(metadata$replicate),
    stringsAsFactors = FALSE
  ))
  stopifnot(
    nrow(sample_identity) == 96L,
    !anyDuplicated(sample_identity$sample_id),
    all(table(sample_identity$condition) == 3L)
  )
  sample_map_index <- match(sample_identity$condition, treatment_map$condition)
  sample_identity$treatment_category <- treatment_map$treatment_category[sample_map_index]
  sample_identity$cytokine_count <- treatment_map$cytokine_count[sample_map_index]
  sample_identity$cytokine_components <- treatment_map$cytokine_components[sample_map_index]
  sample_identity$classification_basis <- treatment_map$classification_basis[sample_map_index]
  sample_category_counts <- table(factor(sample_identity$treatment_category, levels = category_order))
  stopifnot(identical(as.integer(sample_category_counts), c(3L, 78L, 15L)))

  qc_accepted_flag <- imported_cells_before %in% qc_accepted_before
  qc_excluded_flag <- imported_cells_before %in% qc_excluded_before
  downstream_flag <- qc_accepted_flag & treatment_category %in% included_categories
  treatment_scope_excluded_flag <- qc_accepted_flag & treatment_category == "combined"
  downstream_cells <- imported_cells_before[downstream_flag]
  treatment_scope_excluded_cells <- imported_cells_before[treatment_scope_excluded_flag]
  stopifnot(
    identical(qc_accepted_before, imported_cells_before[qc_accepted_flag]),
    identical(qc_excluded_before, imported_cells_before[qc_excluded_flag]),
    !length(intersect(downstream_cells, treatment_scope_excluded_cells)),
    identical(
      sort(c(downstream_cells, treatment_scope_excluded_cells)),
      sort(qc_accepted_before)
    ),
    all(treatment_category[downstream_flag] %in% included_categories),
    all(treatment_category[treatment_scope_excluded_flag] == "combined")
  )

  audit_rows <- data.frame(
    sample_id = as.character(metadata$sample_id),
    condition = as.character(metadata$condition),
    replicate = as.character(metadata$replicate),
    treatment_category = treatment_category,
    cytokine_count = cytokine_count,
    imported = 1L,
    qc_accepted = as.integer(qc_accepted_flag),
    qc_excluded = as.integer(qc_excluded_flag),
    proposed_downstream = as.integer(downstream_flag),
    excluded_by_treatment_scope = as.integer(treatment_scope_excluded_flag),
    stringsAsFactors = FALSE
  )
  count_columns <- c(
    "imported", "qc_accepted", "qc_excluded",
    "proposed_downstream", "excluded_by_treatment_scope"
  )
  per_sample <- aggregate(
    audit_rows[, count_columns, drop = FALSE],
    by = audit_rows[, c(
      "sample_id", "condition", "replicate",
      "treatment_category", "cytokine_count"
    ), drop = FALSE],
    FUN = sum
  )
  per_sample$included_in_downstream_scope <- per_sample$treatment_category %in% included_categories
  per_sample$percent_qc_accepted <- 100 * per_sample$qc_accepted / per_sample$imported
  per_sample <- per_sample[order(per_sample$sample_id), , drop = FALSE]
  rownames(per_sample) <- NULL
  stopifnot(
    nrow(per_sample) == 96L,
    all(per_sample$qc_accepted > 0L),
    sum(per_sample$included_in_downstream_scope) == 81L
  )

  condition_counts <- aggregate(
    audit_rows[, count_columns, drop = FALSE],
    by = audit_rows[, "condition", drop = FALSE],
    FUN = sum
  )
  condition_sample_counts <- aggregate(
    rep(1L, nrow(per_sample)),
    by = list(condition = per_sample$condition),
    FUN = sum
  )
  names(condition_sample_counts)[[2L]] <- "samples"
  condition_classification <- merge(
    treatment_map, condition_sample_counts, by = "condition", all.x = TRUE, sort = FALSE
  )
  condition_classification <- merge(
    condition_classification, condition_counts, by = "condition", all.x = TRUE, sort = FALSE
  )
  condition_classification <- condition_classification[
    match(treatment_map$condition, condition_classification$condition), , drop = FALSE
  ]
  condition_classification$included_in_downstream_scope <-
    condition_classification$treatment_category %in% included_categories
  rownames(condition_classification) <- NULL

  category_summary <- do.call(rbind, lapply(category_order, function(category) {
    cell_rows <- treatment_category == category
    sample_rows <- per_sample$treatment_category == category
    condition_rows <- treatment_map$treatment_category == category
    data.frame(
      treatment_category = category,
      cytokine_count_rule = switch(category, control = "0", single = "1", combined = ">1"),
      conditions = sum(condition_rows),
      samples = sum(sample_rows),
      imported_cells = sum(cell_rows),
      qc_accepted_cells = sum(qc_accepted_flag & cell_rows),
      qc_excluded_cells = sum(qc_excluded_flag & cell_rows),
      proposed_downstream_cells = sum(downstream_flag & cell_rows),
      excluded_by_treatment_scope = sum(treatment_scope_excluded_flag & cell_rows),
      downstream_status = if (category %in% included_categories) "included" else "excluded_combined",
      stringsAsFactors = FALSE
    )
  }))
  rownames(category_summary) <- NULL
  stopifnot(
    sum(category_summary$imported_cells) == length(imported_cells_before),
    sum(category_summary$qc_accepted_cells) == length(qc_accepted_before),
    sum(category_summary$qc_excluded_cells) == length(qc_excluded_before),
    sum(category_summary$proposed_downstream_cells) == length(downstream_cells),
    sum(category_summary$excluded_by_treatment_scope) == length(treatment_scope_excluded_cells)
  )

  cell_classification <- data.frame(
    cell_id = imported_cells_before,
    treatment_category = treatment_category,
    cytokine_count = cytokine_count,
    stringsAsFactors = FALSE
  )
  proposed_scope_exclusions <- data.frame(
    cell_id = treatment_scope_excluded_cells,
    sample_id = as.character(metadata$sample_id[treatment_scope_excluded_flag]),
    condition = as.character(metadata$condition[treatment_scope_excluded_flag]),
    replicate = as.character(metadata$replicate[treatment_scope_excluded_flag]),
    treatment_category = treatment_category[treatment_scope_excluded_flag],
    cytokine_count = cytokine_count[treatment_scope_excluded_flag],
    exclusion_reason = "combined_more_than_one_cytokine",
    stringsAsFactors = FALSE
  )

  write_tsv_atomic(treatment_map, file.path(report_dir, "treatment_category_map.tsv"))
  write_tsv_atomic(condition_classification, file.path(report_dir, "condition_treatment_classification.tsv"))
  write_tsv_atomic(per_sample, file.path(report_dir, "sample_treatment_classification.tsv"))
  write_tsv_atomic(category_summary, file.path(report_dir, "treatment_category_summary.tsv"))
  classification_path <- file.path(report_dir, "cell_treatment_classification.tsv")
  write_tsv_atomic(cell_classification, classification_path)
  write_tsv_gz_atomic(
    proposed_scope_exclusions,
    file.path(report_dir, "proposed_treatment_scope_exclusions.tsv.gz")
  )

  preview_generated_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  preview <- paste0(
    "# Parse cytokine treatment classification and selection preview\n\n",
    "Generated before treatment metadata import, categorical-vector promotion, or downstream cell-set derivation: `",
    preview_generated_utc, "`\n\n",
    "Every one of the 96 samples and every one of the 251,599 imported cells maps to exactly one of three categories. ",
    "`control` has zero cytokine stimuli, `single` has exactly one, and `combined` has more than one. ",
    "The source spellings are preserved; in particular, `Csf_1` and `Csf_2` are reviewed single-stimulus labels.\n\n",
    "## Category summary\n\n",
    markdown_table(category_summary), "\n\n",
    "## Complete condition map\n\n",
    markdown_table(condition_classification[, c(
      "condition", "treatment_category", "cytokine_count", "cytokine_components",
      "samples", "imported", "qc_accepted", "proposed_downstream",
      "excluded_by_treatment_scope", "included_in_downstream_scope"
    )]), "\n\n",
    "The proposed downstream population is the exact intersection of `", qc_accepted_cell_set,
    "` with categories `control` or `single`. QC-accepted `combined` cells are listed in ",
    "`proposed_treatment_scope_exclusions.tsv.gz`. All imported cells, including combined and QC-excluded cells, ",
    "remain in `", imported_cell_set, "`; the existing QC sets also remain unchanged.\n"
  )
  write_lines_atomic(preview, file.path(report_dir, "treatment_selection_preview.md"))
  write_json_atomic(list(
    status = "user_selected_pending_metadata_and_cell_set_promotion",
    generated_utc = preview_generated_utc,
    imported_cell_set = imported_cell_set,
    qc_accepted_cell_set = qc_accepted_cell_set,
    downstream_cell_set = downstream_cell_set,
    included_categories = included_categories,
    proposed_downstream_cells = length(downstream_cells),
    proposed_treatment_scope_exclusions = length(treatment_scope_excluded_cells),
    imported_membership_will_be_preserved = TRUE,
    durable_treatment_vector_created_at_preview_time = FALSE,
    durable_downstream_set_created_at_preview_time = FALSE
  ), file.path(report_dir, "treatment_selection_preview.json"))
  log_event(
    "Treatment preview written before bundle mutation:",
    length(downstream_cells), "proposed downstream cells and",
    length(treatment_scope_excluded_cells), "QC-accepted combined exclusions"
  )

  classification_source <- wb_source_input(
    path = classification_path,
    source_id = "parse_treatment_classification_v1",
    accession = "parse_cytokines_treatment_classification_v1",
    role = "cell_metadata",
    format = "tsv",
    scope = "external_file"
  )
  import_plan <- wb_metadata_import(
    bundle = bundle,
    input = classification_path,
    cell_id_col = "cell_id",
    source_identity = classification_source,
    fields = c("treatment_category", "cytokine_count"),
    unknown_cell_policy = "error",
    validate_only = TRUE,
    created_by = created_by,
    description = "Reviewed Parse cytokine treatment category and cytokine-count classification for every imported cell."
  )
  stopifnot(
    isTRUE(import_plan$valid),
    import_plan$counts$matched == length(imported_cells_before),
    import_plan$counts$unknown == 0L,
    import_plan$counts$unsupplied_known == 0L
  )
  write_json_atomic(unclass(import_plan), file.path(report_dir, "metadata_import_plan.json"))

  metadata_import <- wb_metadata_import(
    bundle = bundle,
    input = classification_path,
    cell_id_col = "cell_id",
    source_identity = classification_source,
    fields = c("treatment_category", "cytokine_count"),
    unknown_cell_policy = "error",
    validate_only = FALSE,
    created_by = created_by,
    description = "Reviewed Parse cytokine treatment category and cytokine-count classification for every imported cell."
  )
  write_json_atomic(unclass(metadata_import), file.path(report_dir, "metadata_import_result.json"))

  category_field <- metadata_import$fields[
    metadata_import$fields$source_name == "treatment_category", , drop = FALSE
  ]
  count_field <- metadata_import$fields[
    metadata_import$fields$source_name == "cytokine_count", , drop = FALSE
  ]
  stopifnot(nrow(category_field) == 1L, nrow(count_field) == 1L)
  category_field_id <- as.character(category_field$field_id[[1L]])
  count_field_id <- as.character(count_field$field_id[[1L]])
  imported_classification <- wb_meta(
    bundle,
    cells = imported_cells_before,
    cols = list(
      wb_ref_metadata_field("core:cell_id"),
      wb_ref_metadata_field(category_field_id),
      wb_ref_metadata_field(count_field_id)
    )
  )
  stopifnot(
    identical(as.character(imported_classification[[1L]]), imported_cells_before),
    identical(as.character(imported_classification[[2L]]), treatment_category),
    identical(as.integer(imported_classification[[3L]]), cytokine_count)
  )

  downstream_selector <- wb_sel_all(bundle, cell_set = qc_accepted_cell_set)
  downstream_selector <- wb_sel_meta(
    downstream_selector,
    col = wb_ref_metadata_field(category_field_id),
    values = included_categories,
    bundle = bundle
  )
  stopifnot(identical(as.character(downstream_selector), downstream_cells))

  vector_result <- wb_session_run_op(
    wb_session_new(), bundle, "metadata_cell_vector",
    inputs = list(
      source_column = wb_ref_metadata_field(category_field_id),
      vector_name = category_vector_name,
      value_type = "categorical",
      missing_policy = "error",
      cell_set = imported_cell_set,
      min_cells = 1L,
      overwrite = FALSE
    ),
    notes = paste(
      "Materialize the reviewed control/single/combined classification for all imported Parse cells.",
      "No cell filtering occurs in this operation."
    )
  )
  downstream_description <- paste(
    "Accepted Parse downstream analysis population: cells in set:parse_qc_accepted",
    "whose reviewed treatment category is control or single. Combined-stimulation",
    "cells remain retained in the imported and QC-derived cell sets."
  )
  downstream_result <- wb_session_run_op(
    vector_result$session, bundle, "cell_set_derive",
    inputs = list(
      cell_ids = downstream_cells,
      min_cells = 1L,
      target_cell_set = downstream_cell_set,
      role = "analysis_universe",
      status = "accepted",
      source = "user_selected_parse_treatment_scope",
      created_by = created_by,
      description = downstream_description,
      recipe = wb_recipe(downstream_selector),
      overwrite = FALSE
    ),
    notes = paste(
      "Derive the downstream universe only after writing the complete treatment-classification preview.",
      "Preserve the imported and QC-derived cell-set memberships."
    )
  )
  write_json_atomic(
    wb_session_dag(
      downstream_result$session, bundle = bundle,
      include_artifacts = TRUE, include_active = TRUE
    ),
    file.path(report_dir, "session_dag_before_promotion.json")
  )

  promotion <- wb_session_promote_current(
    bundle,
    downstream_result$session,
    outputs = list(
      wb_ref_run(vector_result$run$run_id, "vector"),
      wb_ref_run(downstream_result$run$run_id, "manifest"),
      wb_ref_run(downstream_result$run$run_id, "rule_table")
    ),
    description = downstream_description,
    created_by = created_by,
    fail_fast = TRUE
  )
  write_json_atomic(promotion, file.path(report_dir, "promotion.json"))
  category_vector_artifact_id <- artifact_id_for(
    promotion, vector_result$run$run_id, "vector"
  )
  downstream_manifest_artifact_id <- artifact_id_for(
    promotion, downstream_result$run$run_id, "manifest"
  )
  downstream_rule_table_artifact_id <- artifact_id_for(
    promotion, downstream_result$run$run_id, "rule_table"
  )

  imported_cells_after <- as.character(wb_cells(bundle, cell_set = imported_cell_set))
  all_cells_after <- as.character(wb_cells(bundle, cell_set = "all"))
  qc_accepted_after <- as.character(wb_cells(bundle, cell_set = qc_accepted_cell_set))
  qc_excluded_after <- as.character(wb_cells(bundle, cell_set = qc_excluded_cell_set))
  downstream_cells_after <- as.character(wb_cells(bundle, cell_set = downstream_cell_set))
  cell_sets_after <- sort(as.character(wb_cell_sets(bundle)))
  stopifnot(
    identical(imported_cells_after, imported_cells_before),
    identical(all_cells_after, all_cells_before),
    identical(qc_accepted_after, qc_accepted_before),
    identical(qc_excluded_after, qc_excluded_before),
    identical(downstream_cells_after, downstream_cells),
    setequal(cell_sets_after, c(cell_sets_before, downstream_cell_set)),
    all(downstream_cells_after %in% qc_accepted_after),
    !length(intersect(downstream_cells_after, treatment_scope_excluded_cells))
  )

  treatment_vector <- wb_art_load(
    bundle, category_vector_artifact_id, expected_kind = "cell_vector"
  )
  stopifnot(
    identical(as.character(treatment_vector$cell_id), imported_cells_before),
    identical(as.character(treatment_vector$value), treatment_category)
  )
  downstream_recipe <- wb_get_cell_set_recipe(bundle, downstream_cell_set, parse = FALSE)
  stopifnot(
    nrow(downstream_recipe) == 1L,
    identical(as.character(downstream_recipe$role[[1L]]), "analysis_universe"),
    identical(as.character(downstream_recipe$status[[1L]]), "accepted")
  )
  write_tsv_atomic(downstream_recipe, file.path(report_dir, "downstream_cell_set_recipe.tsv"))

  certificates <- list(
    treatment_category_metadata = wb_entity_certificate(bundle, "metadata_field", category_field_id),
    cytokine_count_metadata = wb_entity_certificate(bundle, "metadata_field", count_field_id),
    treatment_category_vector = wb_entity_certificate(bundle, "artifact", category_vector_artifact_id),
    downstream_cell_set = wb_entity_certificate(bundle, "cell_set", downstream_cell_set),
    downstream_manifest_artifact = wb_entity_certificate(
      bundle, "artifact", downstream_manifest_artifact_id
    )
  )
  write_json_atomic(certificates, file.path(report_dir, "entity_certificates.json"))
  write_json_atomic(list(
    imported_cell_set = imported_cell_set,
    imported_cells_before = length(imported_cells_before),
    imported_cells_after = length(imported_cells_after),
    imported_membership_unchanged = identical(imported_cells_after, imported_cells_before),
    qc_accepted_cell_set = qc_accepted_cell_set,
    qc_accepted_cells_before = length(qc_accepted_before),
    qc_accepted_cells_after = length(qc_accepted_after),
    qc_accepted_membership_unchanged = identical(qc_accepted_after, qc_accepted_before),
    qc_excluded_cell_set = qc_excluded_cell_set,
    qc_excluded_cells_before = length(qc_excluded_before),
    qc_excluded_cells_after = length(qc_excluded_after),
    qc_excluded_membership_unchanged = identical(qc_excluded_after, qc_excluded_before),
    downstream_cell_set = downstream_cell_set,
    downstream_cells = length(downstream_cells_after),
    treatment_scope_excluded_qc_accepted_combined_cells = length(treatment_scope_excluded_cells),
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
  final_qc_accepted <- as.character(wb_cells(recipient, cell_set = qc_accepted_cell_set))
  final_qc_excluded <- as.character(wb_cells(recipient, cell_set = qc_excluded_cell_set))
  final_downstream <- as.character(wb_cells(recipient, cell_set = downstream_cell_set))
  final_vector <- wb_art_load(
    recipient, category_vector_artifact_id, expected_kind = "cell_vector"
  )
  final_metadata <- wb_meta(
    recipient,
    cells = imported_cells_before,
    cols = list(
      wb_ref_metadata_field("core:cell_id"),
      wb_ref_metadata_field(category_field_id),
      wb_ref_metadata_field(count_field_id)
    )
  )
  final_sets <- sort(as.character(wb_cell_sets(recipient)))
  wb_close_bundle(recipient)
  stopifnot(
    identical(final_imported, imported_cells_before),
    identical(final_qc_accepted, qc_accepted_before),
    identical(final_qc_excluded, qc_excluded_before),
    identical(final_downstream, downstream_cells),
    identical(final_sets, cell_sets_after),
    identical(as.character(final_vector$cell_id), imported_cells_before),
    identical(as.character(final_vector$value), treatment_category),
    identical(as.character(final_metadata[[1L]]), imported_cells_before),
    identical(as.character(final_metadata[[2L]]), treatment_category),
    identical(as.integer(final_metadata[[3L]]), cytokine_count)
  )
  write_json_atomic(manifest_after, file.path(report_dir, "manifest.json"))

  completion <- list(
    status = "complete",
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    bundle_id = "parse_cytokines",
    bundle_dir = bundle_dir,
    treatment_map = treatment_map_path,
    treatment_categories = category_order,
    included_downstream_categories = included_categories,
    conditions = nrow(treatment_map),
    samples = nrow(per_sample),
    imported = list(
      cell_set = imported_cell_set,
      cells = length(imported_cells_before),
      membership_unchanged = TRUE
    ),
    qc_accepted = list(
      cell_set = qc_accepted_cell_set,
      cells = length(qc_accepted_before),
      membership_unchanged = TRUE
    ),
    qc_excluded = list(
      cell_set = qc_excluded_cell_set,
      cells = length(qc_excluded_before),
      membership_unchanged = TRUE
    ),
    downstream = list(
      cell_set = downstream_cell_set,
      cells = length(downstream_cells),
      samples = sum(per_sample$proposed_downstream > 0L),
      conditions = sum(condition_classification$proposed_downstream > 0L),
      role = "analysis_universe",
      lifecycle_status = "accepted",
      manifest_artifact_id = downstream_manifest_artifact_id,
      rule_table_artifact_id = downstream_rule_table_artifact_id
    ),
    excluded_by_treatment_scope = list(
      definition = "QC-accepted cells with treatment_category == combined",
      cells = length(treatment_scope_excluded_cells)
    ),
    classification = list(
      metadata_import_artifact_id = metadata_import$report_artifact_id,
      metadata_source_sha256 = metadata_import$source_sha256,
      treatment_category_field_id = category_field_id,
      cytokine_count_field_id = count_field_id,
      treatment_category_vector_artifact_id = category_vector_artifact_id
    ),
    category_summary = category_summary,
    checks = list(
      preview_written_before_bundle_mutation = TRUE,
      every_imported_cell_classified_once = TRUE,
      every_sample_classified_once = TRUE,
      combined_has_more_than_one_cytokine = TRUE,
      imported_membership_unchanged = TRUE,
      qc_accepted_membership_unchanged = TRUE,
      qc_excluded_membership_unchanged = TRUE,
      downstream_is_qc_accepted_control_or_single = TRUE,
      combined_absent_from_downstream = TRUE
    ),
    validation = list(
      ok = validation$ok,
      status = validation$status,
      bundle_state_fingerprint = manifest_after$current_state$fingerprint,
      recipient_validation_fingerprint = validation$post_validation_fingerprint
    )
  )
  write_json_atomic(completion, file.path(report_dir, "completion.json"))

  final_report <- paste0(
    "# Parse cytokine treatment classification and downstream selection\n\n",
    "All 96 samples and all 251,599 imported cells were classified as `control`, `single`, or `combined`. ",
    "Combined means more than one cytokine stimulus.\n\n",
    markdown_table(category_summary), "\n\n",
    "The accepted downstream set is `", downstream_cell_set, "` with ",
    format(length(downstream_cells), big.mark = ","), " QC-accepted control/single cells across ",
    sum(per_sample$proposed_downstream > 0L), " samples and ",
    sum(condition_classification$proposed_downstream > 0L), " conditions. ",
    format(length(treatment_scope_excluded_cells), big.mark = ","),
    " QC-accepted combined-stimulation cells are outside that downstream set.\n\n",
    "The imported and pre-existing QC cell sets are unchanged. Combined and QC-excluded cells remain in the bundle. ",
    "The complete condition and sample classifications, per-category counts, proposed treatment-scope exclusions, ",
    "cell-set recipe, provenance records, and recipient validation are retained beside this report.\n"
  )
  write_lines_atomic(final_report, file.path(report_dir, "treatment_selection_report.md"))
  log_event(
    "Treatment classification and downstream selection complete:",
    length(downstream_cells), "downstream cells; validation", validation$status
  )
}

main()
