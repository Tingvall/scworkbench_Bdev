#!/usr/bin/env Rscript

# Recover the successfully committed stage-09 UCell payload after the durable
# session receipt exceeded R's single-string serialization ceiling. Scores are
# read one signature at a time from scWorkbench's documented payload table;
# no score is recomputed and no source cell membership is changed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(paste(
    "Usage: 09_finalize_committed_ucell_report.R BUNDLE_DIR DISCOVERY_DIR REPORT_DIR",
    "SESSION_ROOT SCORE_ARTIFACT_ID SCORE_SUMMARY_ARTIFACT_ID",
    "QC_KEY_ARTIFACT_ID QC_PAIRWISE_ARTIFACT_ID COMMIT_ID"
  ), call. = FALSE)
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
discovery_dir <- normalizePath(args[[2L]], mustWork = TRUE)
report_dir <- normalizePath(args[[3L]], mustWork = FALSE)
session_root <- normalizePath(args[[4L]], mustWork = FALSE)
score_artifact_id <- args[[5L]]
score_summary_artifact_id <- args[[6L]]
qc_key_artifact_id <- args[[7L]]
qc_pairwise_artifact_id <- args[[8L]]
score_commit_id <- args[[9L]]

if (file.exists(report_dir) || dir.exists(report_dir)) {
  stop("Report destination already exists; refusing to overwrite: ", report_dir, call. = FALSE)
}

options(
  width = 180,
  warn = 1,
  error = function() {
    traceback(40L)
    quit(save = "no", status = 1L, runLast = FALSE)
  }
)
suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
  library(DBI)
})

analysis_set <- "set:parse_qc_control_single"
imported_set <- "set:parse_cytokines"
counts_matrix <- "RNA_counts"
created_by <- Sys.getenv("SCWORKBENCH_CREATED_BY", unset = "jtg")
seed <- 1729L

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

write_lines_atomic <- function(value, path) {
  temporary_path <- file.path(dirname(path), paste0(".", basename(path), ".tmp-", Sys.getpid()))
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  writeLines(value, temporary_path, useBytes = TRUE)
  if (!file.rename(temporary_path, path)) stop("Could not publish text: ", path, call. = FALSE)
  invisible(path)
}

artifact_id_for <- function(promotion, run_id, output_name) {
  rows <- promotion$promoted[
    promotion$promoted$run_id == run_id & promotion$promoted$output == output_name &
      promotion$promoted$status == "ok", , drop = FALSE
  ]
  if (nrow(rows) != 1L || is.na(rows$artifact_id[[1L]])) {
    stop("Expected one promoted artifact for ", run_id, "/", output_name, call. = FALSE)
  }
  as.character(rows$artifact_id[[1L]])
}

auc_binary <- function(score, positive) {
  keep <- is.finite(score) & !is.na(positive)
  score <- score[keep]
  positive <- as.logical(positive[keep])
  n_positive <- sum(positive)
  n_negative <- sum(!positive)
  if (!n_positive || !n_negative) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[positive]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

markdown_table <- function(value, digits = 4L) {
  if (!nrow(value)) return("_No rows._")
  formatted <- lapply(value, function(x) {
    if (is.numeric(x)) x <- ifelse(is.na(x), "NA", format(signif(x, digits), trim = TRUE))
    x <- ifelse(is.na(x), "NA", as.character(x))
    gsub("|", "\\\\|", x, fixed = TRUE)
  })
  rows <- vapply(seq_len(nrow(value)), function(i) {
    paste0("| ", paste(vapply(formatted, `[[`, character(1L), i), collapse = " | "), " |")
  }, character(1L))
  paste(c(
    paste0("| ", paste(names(value), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(value)), collapse = " | "), " |"),
    rows
  ), collapse = "\n")
}

group_quantiles <- function(values, group, signature_id) {
  do.call(rbind, lapply(c("control", "target", "other_single"), function(label) {
    local <- values[group == label]
    quantiles <- quantile(local, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                          names = FALSE, na.rm = TRUE)
    data.frame(
      signature_id = signature_id, group = label, cells = length(local),
      min = quantiles[[1L]], p05 = quantiles[[2L]], p25 = quantiles[[3L]],
      median = quantiles[[4L]], p75 = quantiles[[5L]], p95 = quantiles[[6L]],
      max = quantiles[[7L]], stringsAsFactors = FALSE
    )
  }))
}

draw_heatmap <- function(matrix_value, path, type = c("png", "pdf")) {
  type <- match.arg(type)
  if (type == "png") {
    grDevices::png(path, width = 3600, height = max(2200, 90 * nrow(matrix_value)), res = 220)
  } else {
    grDevices::pdf(path, width = 15, height = max(9, 0.32 * nrow(matrix_value)), onefile = TRUE)
  }
  limit <- stats::quantile(abs(matrix_value[is.finite(matrix_value)]), 0.98, names = FALSE)
  if (!is.finite(limit) || limit == 0) limit <- 1
  clipped <- matrix(
    pmax(-limit, pmin(limit, as.numeric(matrix_value))),
    nrow = nrow(matrix_value), ncol = ncol(matrix_value),
    dimnames = dimnames(matrix_value)
  )
  palette <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
  old <- par(mar = c(12, 20, 3, 2))
  on.exit({
    par(old)
    grDevices::dev.off()
  }, add = FALSE)
  image(seq_len(ncol(clipped)), seq_len(nrow(clipped)), t(clipped),
        col = palette, zlim = c(-limit, limit), axes = FALSE,
        xlab = "", ylab = "",
        main = "Oriented UCell sample-median shift from control")
  axis(1, at = seq_len(ncol(clipped)), labels = colnames(clipped), las = 2, cex.axis = 0.65)
  axis(2, at = seq_len(nrow(clipped)), labels = rownames(clipped), las = 2, cex.axis = 0.60)
  mtext("Single-cytokine condition", side = 1, line = 9.5)
  mtext("Signature", side = 2, line = 18)
  box()
}

draw_distributions <- function(sample_scores, signature_summary, path) {
  grDevices::pdf(path, width = 11, height = 8.5, onefile = TRUE)
  old <- par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))
  on.exit({
    par(old)
    grDevices::dev.off()
  }, add = FALSE)
  set.seed(seed)
  for (i in seq_len(nrow(signature_summary))) {
    info <- signature_summary[i, , drop = FALSE]
    local <- sample_scores[sample_scores$signature_id == info$signature_id, , drop = FALSE]
    oriented <- if (info$direction == "up") local$sample_median_ucell else -local$sample_median_ucell
    group <- ifelse(local$condition == info$cytokine, "target",
                    ifelse(local$treatment_category == "control", "control", "other single"))
    boxplot(split(oriented, factor(group, levels = c("control", "target", "other single"))),
            col = c("grey80", "#D73027", "#91BFDB"), outline = FALSE,
            ylab = "Oriented sample-median UCell score", xlab = "",
            main = paste(info$cytokine, info$direction, sep = " - "))
    stripchart(split(oriented, factor(group, levels = c("control", "target", "other single"))),
               vertical = TRUE, method = "jitter", add = TRUE, pch = 16,
               cex = 0.55, col = grDevices::adjustcolor("black", 0.55))
    legend("topleft", legend = paste0("final: ", info$final_status), bty = "n", cex = 0.8)
  }
}

read_committed_score <- function(bundle, artifact_id, signature_id, cells) {
  # wb_art_load() first materializes provenance_json. The current artifact's
  # valid 1.4 GB row-order descriptor exceeds R's 2^31-1-byte string ceiling
  # once scWorkbench/jsonlite expands it into a containing JSON document, so
  # use the same durable payload table that wb_art_load() reads after metadata.
  score <- DBI::dbGetQuery(
    bundle$con,
    paste(
      "SELECT cell_id, pseudotime FROM artifact_pseudotime",
      "WHERE artifact_id = ? AND lineage_id = ?;"
    ),
    params = list(artifact_id, signature_id)
  )
  if (nrow(score) != length(cells) || anyDuplicated(score$cell_id)) {
    stop("Committed UCell payload has invalid coverage for ", signature_id, call. = FALSE)
  }
  index <- match(cells, as.character(score$cell_id))
  if (anyNA(index)) stop("Committed UCell payload is missing analysis cells for ", signature_id, call. = FALSE)
  values <- as.numeric(score$pseudotime[index])
  if (anyNA(values) || any(!is.finite(values))) {
    stop("Committed UCell payload contains non-finite scores for ", signature_id, call. = FALSE)
  }
  values
}

main <- function() {
  discovery_completion <- jsonlite::read_json(file.path(discovery_dir, "completion.json"), simplifyVector = TRUE)
  discovery_parameters <- jsonlite::read_json(file.path(discovery_dir, "analysis_parameters.json"), simplifyVector = TRUE)
  stopifnot(
    identical(discovery_completion$status, "complete"),
    isFALSE(discovery_completion$bdev_used),
    isFALSE(discovery_completion$combined_cells_used)
  )
  candidate_sets <- utils::read.delim(
    file.path(discovery_dir, "validated_feature_sets.tsv"), sep = "\t", quote = "\"",
    comment.char = "", stringsAsFactors = FALSE, check.names = FALSE
  )
  discovery_summary <- utils::read.delim(
    file.path(discovery_dir, "signature_summary.tsv"), sep = "\t", quote = "\"",
    comment.char = "", stringsAsFactors = FALSE, check.names = FALSE
  )
  discovery_members <- utils::read.delim(
    file.path(discovery_dir, "selected_signature_members.tsv"), sep = "\t", quote = "\"",
    comment.char = "", stringsAsFactors = FALSE, check.names = FALSE
  )
  state_report_dir <- discovery_parameters$analysis_scope$state_report_dir
  state_summary <- utils::read.delim(
    file.path(state_report_dir, "state_summary.tsv"), sep = "\t", quote = "\"",
    comment.char = "", stringsAsFactors = FALSE
  )
  state_by_sample <- utils::read.delim(
    file.path(state_report_dir, "state_by_sample.tsv"), sep = "\t", quote = "\"",
    comment.char = "", stringsAsFactors = FALSE
  )
  state_by_sample$sample_analysis_cells <- ave(state_by_sample$cells, state_by_sample$sample_id, FUN = sum)
  state_by_sample$state_proportion <- state_by_sample$cells / state_by_sample$sample_analysis_cells
  state_composition <- aggregate(
    state_by_sample$state_proportion,
    by = state_by_sample[, c("state", "condition", "treatment_category"), drop = FALSE],
    FUN = function(x) c(samples = length(x), mean = mean(x), sd = stats::sd(x),
                        min = min(x), median = median(x), max = max(x))
  )
  state_composition <- data.frame(
    state = state_composition$state,
    condition = state_composition$condition,
    treatment_category = state_composition$treatment_category,
    samples_present = state_composition$x[, "samples"],
    mean_proportion = state_composition$x[, "mean"],
    sd_proportion = state_composition$x[, "sd"],
    min_proportion = state_composition$x[, "min"],
    median_proportion = state_composition$x[, "median"],
    max_proportion = state_composition$x[, "max"],
    stringsAsFactors = FALSE
  )

  candidate_ids <- unique(candidate_sets$feature_set)
  candidate_info <- discovery_summary[match(candidate_ids, discovery_summary$signature_id), , drop = FALSE]
  stopifnot(!anyNA(candidate_info$signature_id), all(candidate_info$status == "validated"))
  signature_list <- split(candidate_sets$feature_id, factor(candidate_sets$feature_set, levels = candidate_ids))
  signature_list <- lapply(signature_list, unique)

  dir.create(session_root, recursive = TRUE, showWarnings = FALSE)
  bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
  on.exit(if (!is.null(bundle)) try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  manifest_before <- wb_manifest(bundle)
  imported_cells_before <- as.character(wb_cells(bundle, cell_set = imported_set))
  analysis_cells <- as.character(wb_cells(bundle, cell_set = analysis_set))
  cell_sets_before <- sort(as.character(wb_cell_sets(bundle)))
  stopifnot(length(imported_cells_before) == 251599L, length(analysis_cells) == 154209L)

  metadata <- wb_meta(
    bundle, cells = analysis_cells,
    cols = list(wb_ref_metadata_field("core:cell_id"), "sample_id", "condition",
                "replicate", "treatment_category")
  )
  stopifnot(identical(as.character(metadata$cell_id), analysis_cells), !anyNA(metadata))
  sample_metadata <- unique(metadata[, c("sample_id", "condition", "replicate", "treatment_category")])
  control_condition <- unique(sample_metadata$condition[sample_metadata$treatment_category == "control"])
  single_conditions <- sort(unique(sample_metadata$condition[sample_metadata$treatment_category == "single"]))
  stopifnot(length(control_condition) == 1L, length(single_conditions) == 26L)

  artifact_audit <- DBI::dbGetQuery(
    bundle$con,
    paste(
      "SELECT artifact_id, type, name, source_op, cell_set,",
      "length(provenance_json) AS provenance_bytes",
      "FROM artifacts WHERE artifact_id IN (?, ?, ?, ?) ORDER BY artifact_id;"
    ),
    params = list(score_artifact_id, score_summary_artifact_id,
                  qc_key_artifact_id, qc_pairwise_artifact_id)
  )
  if (nrow(artifact_audit) != 4L ||
      !identical(artifact_audit$type[artifact_audit$artifact_id == score_artifact_id], "cell_vector_family")) {
    stop("The expected committed UCell artifacts are not present", call. = FALSE)
  }
  commit_audit <- DBI::dbGetQuery(
    bundle$con,
    "SELECT commit_id, created_at, created_by, state, selected_outputs_json, error_message FROM prov_commits WHERE commit_id = ?;",
    params = list(score_commit_id)
  )
  if (nrow(commit_audit) != 1L || !identical(commit_audit$state[[1L]], "complete")) {
    stop("The expected UCell promotion commit is not complete", call. = FALSE)
  }
  payload_audit <- DBI::dbGetQuery(
    bundle$con,
    paste(
      "SELECT lineage_id AS signature_id, COUNT(*) AS score_rows,",
      "COUNT(DISTINCT cell_id) AS distinct_cells, MIN(pseudotime) AS score_min,",
      "MAX(pseudotime) AS score_max",
      "FROM artifact_pseudotime WHERE artifact_id = ? GROUP BY lineage_id ORDER BY lineage_id;"
    ),
    params = list(score_artifact_id)
  )
  stopifnot(
    nrow(payload_audit) == length(candidate_ids),
    setequal(payload_audit$signature_id, candidate_ids),
    all(payload_audit$score_rows == length(analysis_cells)),
    all(payload_audit$distinct_cells == length(analysis_cells)),
    all(is.finite(payload_audit$score_min)),
    all(is.finite(payload_audit$score_max))
  )

  sample_score_rows <- cell_quantile_rows <- confirmation_rows <- condition_profile_rows <- list()
  for (i in seq_len(nrow(candidate_info))) {
    info <- candidate_info[i, , drop = FALSE]
    signature_id <- info$signature_id
    values <- read_committed_score(bundle, score_artifact_id, signature_id, analysis_cells)
    direction <- if (info$direction == "up") 1 else -1
    cell_group <- ifelse(metadata$condition == info$cytokine, "target",
                         ifelse(metadata$treatment_category == "control", "control", "other_single"))
    cell_quantile_rows[[i]] <- group_quantiles(values, cell_group, signature_id)

    sample_medians <- tapply(values, metadata$sample_id, median, na.rm = TRUE)
    local_metadata <- sample_metadata[match(names(sample_medians), sample_metadata$sample_id), , drop = FALSE]
    local <- data.frame(
      signature_id = signature_id,
      cytokine = info$cytokine,
      direction = info$direction,
      sample_id = names(sample_medians),
      condition = local_metadata$condition,
      replicate = local_metadata$replicate,
      treatment_category = local_metadata$treatment_category,
      sample_median_ucell = as.numeric(sample_medians),
      oriented_sample_median = direction * as.numeric(sample_medians),
      stringsAsFactors = FALSE
    )
    sample_score_rows[[i]] <- local

    response_rows <- local$condition %in% c(control_condition, info$cytokine)
    specificity_rows <- local$treatment_category == "single"
    response_auc <- auc_binary(local$oriented_sample_median[response_rows],
                               local$condition[response_rows] == info$cytokine)
    specificity_auc <- auc_binary(local$oriented_sample_median[specificity_rows],
                                  local$condition[specificity_rows] == info$cytokine)
    balanced_auc <- mean(c(response_auc, specificity_auc))
    confirmed <- response_auc >= 0.60 && specificity_auc >= 0.60 && balanced_auc >= 0.65
    confirmation_rows[[i]] <- data.frame(
      signature_id = signature_id,
      cytokine = info$cytokine,
      direction = info$direction,
      genes = length(signature_list[[signature_id]]),
      discovery_cv_response_auc = info$auc_response_vs_control,
      discovery_cv_specificity_auc = info$auc_specificity_vs_other_singles,
      discovery_cv_balanced_auc = info$balanced_auc,
      discovery_fold_gene_jaccard = info$fold_gene_jaccard,
      ucell_response_auc = response_auc,
      ucell_specificity_auc = specificity_auc,
      ucell_balanced_auc = balanced_auc,
      final_status = if (confirmed) "accepted" else "ucell_transfer_not_confirmed",
      stringsAsFactors = FALSE
    )

    control_median <- median(local$sample_median_ucell[local$treatment_category == "control"])
    condition_profile_rows[[i]] <- do.call(rbind, lapply(single_conditions, function(condition) {
      condition_values <- local$sample_median_ucell[local$condition == condition]
      data.frame(
        signature_id = signature_id,
        cytokine = info$cytokine,
        direction = info$direction,
        evaluated_condition = condition,
        is_target = condition == info$cytokine,
        condition_sample_median = median(condition_values),
        control_sample_median = control_median,
        oriented_shift_from_control = direction * (median(condition_values) - control_median),
        stringsAsFactors = FALSE
      )
    }))
    message("Recovered and evaluated UCell signature ", i, "/", nrow(candidate_info), ": ", signature_id)
  }

  sample_scores <- do.call(rbind, sample_score_rows)
  cell_quantiles <- do.call(rbind, cell_quantile_rows)
  confirmation <- do.call(rbind, confirmation_rows)
  condition_profile <- do.call(rbind, condition_profile_rows)
  accepted_ids <- confirmation$signature_id[confirmation$final_status == "accepted"]
  if (!length(accepted_ids)) stop("No candidate signature transferred successfully to UCell", call. = FALSE)

  accepted_sets <- candidate_sets[candidate_sets$feature_set %in% accepted_ids, , drop = FALSE]
  accepted_sets$description <- paste0(
    accepted_sets$description,
    " Final acceptance additionally required sample-level UCell response/specificity transfer."
  )

  existing_namespace <- wb_feature_sets(bundle, matrix_id = counts_matrix,
                                        namespace = "parse_cytokine_signatures")
  if (nrow(existing_namespace)) {
    stop("The parse_cytokine_signatures namespace already exists; refusing an ambiguous recovery import", call. = FALSE)
  }
  feature_set_interface <- list(
    preferred_op_name = "feature_set_import",
    interface_retrieved_in_original_stage_09 = TRUE,
    recovery_retrieval_skipped = TRUE,
    recovery_route = "read-only report export; bundle registration deferred",
    reason = paste(
      "Every supported writable open first refreshes the manifest and enumerates current artifacts.",
      "That path attempts to materialize the committed 1.398-GB provenance value before a",
      "feature-set import can begin."
    )
  )
  import_receipt <- data.frame(
    status = "deferred",
    namespace = "parse_cytokine_signatures",
    accepted_sets = length(accepted_ids),
    reason = paste(
      "Writable scWorkbench open is blocked by the existing 1.398-GB artifact provenance value;",
      "no direct SQLite registry mutation was attempted."
    ),
    stringsAsFactors = FALSE
  )
  registered_sets <- existing_namespace
  stopifnot(
    nrow(registered_sets) == 0L,
    identical(as.character(wb_cells(bundle, cell_set = imported_set)), imported_cells_before),
    identical(sort(as.character(wb_cell_sets(bundle))), cell_sets_before),
    identical(as.character(wb_cells(bundle, cell_set = analysis_set)), analysis_cells)
  )

  readiness <- wb_recipient_readiness(bundle_dir)
  stopifnot(isTRUE(readiness$ready))
  manifest_after <- wb_manifest(bundle)
  quick_check <- DBI::dbGetQuery(bundle$con, "PRAGMA quick_check;")[[1L]]
  stopifnot(identical(as.character(quick_check), "ok"))

  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  accepted_feature_path <- file.path(report_dir, "final_accepted_feature_sets.tsv")
  write_tsv_atomic(accepted_sets, accepted_feature_path)
  accepted_gmt <- vapply(unique(accepted_sets$feature_set), function(set_id) {
    rows <- accepted_sets$feature_set == set_id
    paste(c(set_id, accepted_sets$description[which(rows)[[1L]]], accepted_sets$feature_id[rows]), collapse = "\t")
  }, character(1L))
  write_lines_atomic(accepted_gmt, file.path(report_dir, "final_accepted_feature_sets.gmt"))

  complete_status <- merge(
    discovery_summary, confirmation,
    by = c("signature_id", "cytokine", "direction"),
    all.x = TRUE, sort = FALSE, suffixes = c("_discovery", "_ucell")
  )
  complete_status$final_status[is.na(complete_status$final_status)] <- "not_advanced_from_discovery"
  complete_status <- complete_status[order(match(complete_status$cytokine, single_conditions),
                                             match(complete_status$direction, c("up", "down"))), , drop = FALSE]
  top_gene_map <- tapply(discovery_members$feature_id, discovery_members$signature_id,
                         function(x) paste(utils::head(x, 8L), collapse = ", "))
  complete_status$top_genes <- unname(top_gene_map[complete_status$signature_id])

  write_tsv_atomic(sample_scores, file.path(report_dir, "ucell_sample_medians.tsv"))
  write_tsv_atomic(cell_quantiles, file.path(report_dir, "ucell_cell_distribution_quantiles.tsv"))
  write_tsv_atomic(condition_profile, file.path(report_dir, "ucell_condition_specificity_profile.tsv"))
  write_tsv_atomic(confirmation, file.path(report_dir, "ucell_transfer_validation.tsv"))
  write_tsv_atomic(complete_status, file.path(report_dir, "all_signature_status.tsv"))
  write_tsv_atomic(state_summary, file.path(report_dir, "parse_state_summary.tsv"))
  write_tsv_atomic(state_by_sample, file.path(report_dir, "parse_state_composition_by_sample.tsv"))
  write_tsv_atomic(state_composition, file.path(report_dir, "parse_state_composition_by_condition.tsv"))
  write_tsv_atomic(artifact_audit, file.path(report_dir, "committed_ucell_artifacts.tsv"))
  write_tsv_atomic(payload_audit, file.path(report_dir, "committed_ucell_payload_audit.tsv"))
  write_tsv_atomic(commit_audit, file.path(report_dir, "committed_ucell_commit.tsv"))
  write_tsv_atomic(registered_sets, file.path(report_dir, "registered_feature_sets.tsv"))
  write_json_atomic(feature_set_interface, file.path(report_dir, "feature_set_import_interface.json"))
  write_tsv_atomic(import_receipt, file.path(report_dir, "feature_set_import_receipt.tsv"))
  write_json_atomic(readiness, file.path(report_dir, "recipient_readiness.json"))
  write_json_atomic(manifest_before, file.path(report_dir, "manifest_before.json"))
  write_json_atomic(manifest_after, file.path(report_dir, "manifest.json"))
  write_json_atomic(list(
    status = "bounded_validation_passed",
    sqlite_quick_check = quick_check,
    imported_cell_set_preserved = TRUE,
    analysis_cell_set_preserved = TRUE,
    cell_set_registry_preserved = TRUE,
    score_artifact_registry_rows = nrow(artifact_audit),
    score_signatures = nrow(payload_audit),
    scores_per_signature = unique(payload_audit$score_rows),
    score_cells_match_analysis_set = TRUE,
    final_feature_sets_exported = length(accepted_ids),
    final_feature_sets_registered_readonly = 0L,
    feature_set_registration_status = "deferred_pending_oversized_artifact_remediation",
    recipient_ready = readiness$ready,
    full_recipient_validation_run = FALSE,
    full_recipient_validation_reason = paste(
      "The committed score artifact contains a valid 1.398-GB provenance row-order descriptor.",
      "Serializing the enclosing scWorkbench manifest/receipt exceeds R's 2^31-1-byte",
      "character-string ceiling, so the standard validator cannot inspect this artifact."
    )
  ), file.path(report_dir, "validation.json"))

  ordered_candidates <- candidate_info[order(match(candidate_info$cytokine, single_conditions),
                                                match(candidate_info$direction, c("up", "down"))), , drop = FALSE]
  profile_matrix <- matrix(NA_real_, nrow = nrow(ordered_candidates), ncol = length(single_conditions),
                           dimnames = list(ordered_candidates$signature_id, single_conditions))
  for (i in seq_len(nrow(condition_profile))) {
    profile_matrix[condition_profile$signature_id[[i]], condition_profile$evaluated_condition[[i]]] <-
      condition_profile$oriented_shift_from_control[[i]]
  }
  draw_heatmap(profile_matrix, file.path(report_dir, "ucell_specificity_heatmap.png"), "png")
  draw_heatmap(profile_matrix, file.path(report_dir, "ucell_specificity_heatmap.pdf"), "pdf")
  plot_info <- merge(candidate_info, confirmation[, c("signature_id", "final_status")],
                     by = "signature_id", all.x = TRUE, sort = FALSE)
  draw_distributions(sample_scores, plot_info, file.path(report_dir, "ucell_signature_distributions.pdf"))

  accepted_summary <- confirmation[confirmation$final_status == "accepted", , drop = FALSE]
  accepted_summary <- accepted_summary[order(match(accepted_summary$cytokine, single_conditions),
                                               match(accepted_summary$direction, c("up", "down"))), , drop = FALSE]
  status_counts <- as.data.frame(table(final_status = complete_status$final_status), stringsAsFactors = FALSE)
  report <- paste0(
    "# Individual-cytokine transcriptional signatures from Parse\n\n",
    "Generated: `", format(Sys.time(), tz = "UTC", usetz = TRUE), "`\n\n",
    "## Result\n\n",
    "The workflow produced **", nrow(accepted_summary), " final accepted signatures** across **",
    length(unique(accepted_summary$cytokine)), " of 26 individual cytokines**. Of these, ",
    sum(accepted_summary$direction == "up"), " are induced and ",
    sum(accepted_summary$direction == "down"), " are repressed. Unsupported cytokine/direction pairs were ",
    "reported rather than forced into the library.\n\n",
    markdown_table(status_counts), "\n\n",
    "## Population boundaries\n\n",
    "All **251,599 imported cells** remain in `", imported_set, "`. Discovery and scoring used only the ",
    "**154,209 QC-accepted control/single cells** in `", analysis_set, "` (81 samples: 3 control and ",
    "78 single-cytokine samples). Combined-stimulation cells remain retained but were not used. ",
    "Bdev was not read or used in state construction, discovery, validation, or scoring.\n\n",
    "## Analysis strategy\n\n",
    "Parse-only states were learned from control-derived HVGs with unlabeled PCA/Leiden clustering. Raw counts ",
    "were summed by biological sample, then limma-voom tested target-versus-control response and target-versus-",
    "other-singles specificity contrasts. Induced and repressed programs were separate. Genes had to pass both ",
    "contrasts, effect-size cutoffs, FDR 0.05, and direction stability in all three leave-one-replicate-out fits. ",
    "Compact set size was selected from held-out sample-level rank performance and fold membership stability. ",
    "State composition is reported separately, and within-state concordance contributed to gene ranking.\n\n",
    "The 31 discovery candidates were scored per cell by scWorkbench UCell on raw `RNA_counts` ",
    "(`max_rank=1500`) without graph smoothing or cross-signature calibration. Biological-sample medians—not ",
    "cells as replicates—were the confirmation unit. Final acceptance required response and specificity UCell ",
    "AUC >=0.60 and balanced AUC >=0.65. Scores are compared within signatures, not across signatures.\n\n",
    "## Final accepted signatures\n\n",
    markdown_table(accepted_summary[, c(
      "cytokine", "direction", "genes", "discovery_cv_response_auc",
      "discovery_cv_specificity_auc", "ucell_response_auc", "ucell_specificity_auc"
    )]), "\n\n",
    "## Specificity diagnostics\n\n",
    "![Oriented UCell specificity heatmap](ucell_specificity_heatmap.png)\n\n",
    "Each heatmap value is the direction-oriented shift of a condition's median sample-level UCell score from ",
    "control. The multi-page PDF contains the sample-level distributions. Exact sample medians, cell-level ",
    "descriptive quantiles, off-target profiles, state-composition summaries, and every accepted or rejected ",
    "status are supplied as TSV files.\n\n",
    "## Durable outputs and validation boundary\n\n",
    "Final gene sets are in `final_accepted_feature_sets.tsv` and `.gmt`. Bundle registration in namespace ",
    "`parse_cytokine_signatures` is deferred: every supported writable scWorkbench open refreshes the manifest and ",
    "therefore reaches the oversized provenance value before import can begin. No direct SQLite registry mutation ",
    "was attempted. The committed UCell family is `", score_artifact_id, "` (commit `", score_commit_id, "`). ",
    "Recipient readiness and SQLite `quick_check` passed, all 31 score keys ",
    "have exactly 154,209 distinct analysis cells, and imported/analysis memberships were unchanged.\n\n",
    "The standard full recipient validator could not be run because this scWorkbench build stores a row-order ",
    "token for every family row in `provenance_json`; the 31-signature artifact's valid descriptor is 1.398 GB, ",
    "and serializing its enclosing manifest/receipt exceeds R's 2^31-1-byte string ceiling. The scores themselves ",
    "are in the normal ",
    "`artifact_pseudotime` payload and were read one signature at a time. `validation.json` records the bounded ",
    "checks and this limitation; no artifact was deleted or rewritten during recovery.\n"
  )
  write_lines_atomic(report, file.path(report_dir, "REPORT.md"))
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
  write_json_atomic(list(
    status = "report_complete_bundle_registration_blocked_by_scworkbench_json_string_limit",
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    final_accepted_signatures = nrow(accepted_summary),
    cytokines_with_final_signature = length(unique(accepted_summary$cytokine)),
    final_up_signatures = sum(accepted_summary$direction == "up"),
    final_down_signatures = sum(accepted_summary$direction == "down"),
    score_artifact_id = score_artifact_id,
    score_summary_artifact_id = score_summary_artifact_id,
    qc_key_artifact_id = qc_key_artifact_id,
    qc_pairwise_artifact_id = qc_pairwise_artifact_id,
    score_commit_id = score_commit_id,
    feature_set_import_route = "deferred_no_direct_database_mutation",
    imported_cells_retained = length(imported_cells_before),
    analysis_cells = length(analysis_cells),
    combined_cells_used = FALSE,
    bdev_used = FALSE,
    bundle_state_fingerprint = manifest_after$current_state$fingerprint,
    validation_status = "bounded_validation_passed_full_validator_blocked_by_r_json_string_limit"
  ), file.path(report_dir, "completion.json"))

  wb_close_bundle(bundle)
  bundle <- NULL
}

main()
