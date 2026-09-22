#!/usr/bin/env Rscript

# Score held-out-validated Parse cytokine signatures with scWorkbench UCell,
# confirm score transfer at the biological-sample level, register only the
# final accepted feature sets, and generate the final audit/report package.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4L, 5L)) {
  stop("Usage: 09_score_and_report_signatures.R BUNDLE_DIR DISCOVERY_DIR REPORT_DIR SESSION_ROOT [RESUME_WORKSPACE]", call. = FALSE)
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
discovery_dir <- normalizePath(args[[2L]], mustWork = TRUE)
report_dir <- normalizePath(args[[3L]], mustWork = FALSE)
session_root <- normalizePath(args[[4L]], mustWork = FALSE)
resume_workspace <- if (length(args) == 5L) normalizePath(args[[5L]], mustWork = TRUE) else NULL
if (file.exists(report_dir) || dir.exists(report_dir)) {
  stop("Report destination already exists; refusing to overwrite: ", report_dir, call. = FALSE)
}

options(width = 180, warn = 1)
suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
})

analysis_set <- "set:parse_qc_control_single"
imported_set <- "set:parse_cytokines"
counts_matrix <- "RNA_counts"
threads <- as.integer(Sys.getenv("PARSE_THREADS", unset = "8"))
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
  (sum(ranks[positive]) - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative)
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

score_vector <- function(score_family, signature, cells) {
  rows <- score_family$key == signature
  values <- as.numeric(score_family$value[rows])
  score_cells <- as.character(score_family$cell_id[rows])
  index <- match(cells, score_cells)
  if (anyNA(index)) stop("UCell score missing cells for signature ", signature, call. = FALSE)
  values[index]
}

group_quantiles <- function(values, group, signature_id) {
  groups <- c("control", "target", "other_single")
  do.call(rbind, lapply(groups, function(label) {
    local <- values[group == label]
    quantiles <- quantile(local, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                          names = FALSE, na.rm = TRUE)
    data.frame(
      signature_id = signature_id,
      group = label,
      cells = length(local),
      min = quantiles[[1L]], p05 = quantiles[[2L]], p25 = quantiles[[3L]],
      median = quantiles[[4L]], p75 = quantiles[[5L]], p95 = quantiles[[6L]],
      max = quantiles[[7L]],
      stringsAsFactors = FALSE
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
    rows <- sample_scores$signature_id == info$signature_id
    local <- sample_scores[rows, , drop = FALSE]
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
  state_summary <- utils::read.delim(
    file.path(discovery_parameters$analysis_scope$state_report_dir, "state_summary.tsv"),
    sep = "\t", quote = "\"", comment.char = "", stringsAsFactors = FALSE
  )
  state_by_sample <- utils::read.delim(
    file.path(discovery_parameters$analysis_scope$state_report_dir, "state_by_sample.tsv"),
    sep = "\t", quote = "\"", comment.char = "", stringsAsFactors = FALSE
  )
  state_by_sample$sample_analysis_cells <- ave(state_by_sample$cells, state_by_sample$sample_id, FUN = sum)
  state_by_sample$state_proportion <- state_by_sample$cells / state_by_sample$sample_analysis_cells
  state_composition_by_condition <- aggregate(
    state_by_sample$state_proportion,
    by = state_by_sample[, c("state", "condition", "treatment_category"), drop = FALSE],
    FUN = function(x) c(samples = length(x), mean = mean(x), sd = stats::sd(x),
                        min = min(x), median = median(x), max = max(x))
  )
  state_composition_by_condition <- data.frame(
    state = state_composition_by_condition$state,
    condition = state_composition_by_condition$condition,
    treatment_category = state_composition_by_condition$treatment_category,
    samples_present = state_composition_by_condition$x[, "samples"],
    mean_proportion = state_composition_by_condition$x[, "mean"],
    sd_proportion = state_composition_by_condition$x[, "sd"],
    min_proportion = state_composition_by_condition$x[, "min"],
    median_proportion = state_composition_by_condition$x[, "median"],
    max_proportion = state_composition_by_condition$x[, "max"],
    stringsAsFactors = FALSE
  )
  candidate_ids <- unique(candidate_sets$feature_set)
  candidate_info <- discovery_summary[match(candidate_ids, discovery_summary$signature_id), , drop = FALSE]
  stopifnot(!anyNA(candidate_info$signature_id), all(candidate_info$status == "validated"))
  signature_list <- split(candidate_sets$feature_id, factor(candidate_sets$feature_set, levels = candidate_ids))
  signature_list <- lapply(signature_list, unique)

  dir.create(session_root, recursive = TRUE, showWarnings = FALSE)
  # Stage 07 has already completed recipient validation and this stage performs
  # another full recipient gate after both mutations. Avoid paying for a third
  # full content audit before UCell while retaining exact manifest/set checks.
  bundle <- wb_open_bundle(bundle_dir, read_only = FALSE, validate = FALSE)
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

  interfaces <- list(
    ucell = wb_op_interface("ucell", bundle = bundle),
    cell_vector_family_qc = wb_op_interface("cell_vector_family_qc", bundle = bundle),
    feature_set_import = wb_op_interface("feature_set_import", bundle = bundle)
  )

  scoring_session <- if (is.null(resume_workspace)) {
    wb_session_create(
      bundle, workspace_root = session_root,
      session_id = paste0("parse_signature_ucell_", format(Sys.time(), "%Y%m%dT%H%M%S"))
    )
  } else {
    wb_session_resume(bundle, resume_workspace)
  }
  signature_list_up <- signature_list[candidate_info$signature_id[candidate_info$direction == "up"]]
  signature_list_down <- signature_list[candidate_info$signature_id[candidate_info$direction == "down"]]
  stopifnot(length(signature_list_up) > 0L, length(signature_list_down) > 0L)
  ucell_up <- wb_session_run_op(
    scoring_session, bundle, "ucell",
    inputs = list(
      matrix = counts_matrix,
      cell_set = analysis_set,
      min_cells = 10L,
      features = signature_list_up,
      name = "parse_cytokine_up",
      max_rank = 1500L,
      ncores = threads,
      parallel_backend = "auto",
      chunk_size = 400L,
      missing_features_policy = "error"
    ),
    notes = paste(
      "Rank-based per-cell scoring of replicate-validated induced candidate signatures on raw RNA counts.",
      "No graph smoothing and no cross-signature calibration are applied.",
      "A 400-cell chunk bounds the current BPCells dense fallback below 256 MiB without changing score values.",
      "Induced and repressed families are separate to bound in-session memory; per-cell families are not promoted."
    )
  )
  ucell_down <- wb_session_run_op(
    ucell_up$session, bundle, "ucell",
    inputs = list(
      matrix = counts_matrix,
      cell_set = analysis_set,
      min_cells = 10L,
      features = signature_list_down,
      name = "parse_cytokine_down",
      max_rank = 1500L,
      ncores = threads,
      parallel_backend = "auto",
      chunk_size = 400L,
      missing_features_policy = "error"
    ),
    notes = paste(
      "Rank-based per-cell scoring of replicate-validated repressed candidate signatures on raw RNA counts.",
      "No graph smoothing and no cross-signature calibration are applied.",
      "A 400-cell chunk bounds the current BPCells dense fallback below 256 MiB without changing score values.",
      "Induced and repressed families are separate to bound in-session memory; per-cell families are not promoted."
    )
  )
  score_qc_up <- wb_session_run_op(
    ucell_down$session, bundle, "cell_vector_family_qc",
    inputs = list(
      signal = wb_ref_run(ucell_up$run$run_id, "signature_scores"),
      cell_set = analysis_set,
      min_cells = 10L,
      correlation_methods = c("pearson", "spearman"),
      tie_tol = 1e-8
    ),
    notes = "QC the raw induced UCell family without smoothing or cross-signature calibration."
  )
  score_qc_down <- wb_session_run_op(
    score_qc_up$session, bundle, "cell_vector_family_qc",
    inputs = list(
      signal = wb_ref_run(ucell_down$run$run_id, "signature_scores"),
      cell_set = analysis_set,
      min_cells = 10L,
      correlation_methods = c("pearson", "spearman"),
      tie_tol = 1e-8
    ),
    notes = "QC the raw repressed UCell family without smoothing or cross-signature calibration."
  )
  scoring_dag <- wb_session_dag(score_qc_down$session, bundle = bundle,
                                include_artifacts = TRUE, include_active = TRUE)
  scoring_promotion <- wb_session_promote_current(
    bundle, score_qc_down$session,
    outputs = list(
      wb_ref_run(ucell_up$run$run_id, "signature_summary"),
      wb_ref_run(ucell_down$run$run_id, "signature_summary"),
      wb_ref_run(score_qc_up$run$run_id, "key_summary"),
      wb_ref_run(score_qc_up$run$run_id, "pairwise_summary"),
      wb_ref_run(score_qc_down$run$run_id, "key_summary"),
      wb_ref_run(score_qc_down$run$run_id, "pairwise_summary")
    ),
    name_prefix = "parse_cytokine_candidate_ucell",
    description = paste(
      "UCell scores for held-out-validated individual-cytokine signature candidates",
      "within the accepted Parse control-plus-single analysis population."
    ),
    created_by = created_by,
    fail_fast = TRUE
  )
  score_summary_artifact_ids <- c(
    up = artifact_id_for(scoring_promotion, ucell_up$run$run_id, "signature_summary"),
    down = artifact_id_for(scoring_promotion, ucell_down$run$run_id, "signature_summary")
  )
  score_family <- rbind(
    ucell_up$run$outputs$signature_scores,
    ucell_down$run$outputs$signature_scores
  )
  stopifnot(setequal(unique(as.character(score_family$key)), candidate_ids))

  sample_score_rows <- list()
  cell_quantile_rows <- list()
  confirmation_rows <- list()
  condition_profile_rows <- list()
  sample_counter <- quantile_counter <- confirmation_counter <- profile_counter <- 0L
  for (i in seq_len(nrow(candidate_info))) {
    info <- candidate_info[i, , drop = FALSE]
    signature_id <- info$signature_id
    values <- score_vector(score_family, signature_id, analysis_cells)
    direction <- if (info$direction == "up") 1 else -1
    cell_group <- ifelse(metadata$condition == info$cytokine, "target",
                         ifelse(metadata$treatment_category == "control", "control", "other_single"))
    quantile_counter <- quantile_counter + 1L
    cell_quantile_rows[[quantile_counter]] <- group_quantiles(values, cell_group, signature_id)

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
    sample_counter <- sample_counter + 1L
    sample_score_rows[[sample_counter]] <- local

    response_rows <- local$condition %in% c(control_condition, info$cytokine)
    specificity_rows <- local$treatment_category == "single"
    response_auc <- auc_binary(local$oriented_sample_median[response_rows],
                               local$condition[response_rows] == info$cytokine)
    specificity_auc <- auc_binary(local$oriented_sample_median[specificity_rows],
                                  local$condition[specificity_rows] == info$cytokine)
    balanced_auc <- mean(c(response_auc, specificity_auc))
    confirmed <- response_auc >= 0.60 && specificity_auc >= 0.60 && balanced_auc >= 0.65
    confirmation_counter <- confirmation_counter + 1L
    confirmation_rows[[confirmation_counter]] <- data.frame(
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
    for (condition in single_conditions) {
      profile_counter <- profile_counter + 1L
      condition_values <- local$sample_median_ucell[local$condition == condition]
      condition_profile_rows[[profile_counter]] <- data.frame(
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
    }
  }

  sample_scores <- do.call(rbind, sample_score_rows)
  cell_quantiles <- do.call(rbind, cell_quantile_rows)
  confirmation <- do.call(rbind, confirmation_rows)
  condition_profile <- do.call(rbind, condition_profile_rows)
  accepted_ids <- confirmation$signature_id[confirmation$final_status == "accepted"]
  if (!length(accepted_ids)) stop("No candidate signature transferred successfully to UCell", call. = FALSE)

  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(report_dir)) stop("Could not create report directory", call. = FALSE)
  accepted_sets <- candidate_sets[candidate_sets$feature_set %in% accepted_ids, , drop = FALSE]
  accepted_sets$description <- paste0(
    accepted_sets$description,
    " Final acceptance additionally required sample-level UCell response/specificity transfer."
  )
  accepted_feature_path <- file.path(report_dir, "final_accepted_feature_sets.tsv")
  write_tsv_atomic(accepted_sets, accepted_feature_path)
  accepted_gmt <- vapply(unique(accepted_sets$feature_set), function(set_id) {
    rows <- accepted_sets$feature_set == set_id
    paste(c(set_id, accepted_sets$description[which(rows)[[1L]]], accepted_sets$feature_id[rows]), collapse = "\t")
  }, character(1L))
  write_lines_atomic(accepted_gmt, file.path(report_dir, "final_accepted_feature_sets.gmt"))

  import_session <- wb_session_create(
    bundle, workspace_root = session_root,
    session_id = paste0("parse_signature_import_", format(Sys.time(), "%Y%m%dT%H%M%S"))
  )
  imported <- wb_session_run_op(
    import_session, bundle, "feature_set_import",
    inputs = list(
      matrix = counts_matrix,
      cell_set = analysis_set,
      min_cells = 10L,
      source_path = accepted_feature_path,
      source_format = "table",
      namespace = "parse_cytokine_signatures",
      collection = "individual_cytokine",
      species = "Mus musculus",
      organism = "Mus musculus",
      gene_id_type = "gene_symbol",
      set_id_col = "feature_set",
      set_name_col = "set_name",
      gene_col = "feature_id",
      description_col = "description",
      split_by_collection = FALSE,
      overwrite = FALSE,
      strict = TRUE,
      created_by = created_by,
      source_id = "parse_cytokine_signature_workflow_v1",
      source_version = "2026-09-21",
      origin = "derived",
      is_readonly = TRUE
    ),
    notes = paste(
      "Register only signatures that passed replicate-aware held-out discovery",
      "and sample-level transfer to scWorkbench UCell scoring."
    )
  )
  import_dag <- wb_session_dag(imported$session, bundle = bundle,
                               include_artifacts = TRUE, include_active = TRUE)
  import_promotion <- wb_session_promote_current(
    bundle, imported$session,
    outputs = list(
      wb_ref_run(imported$run$run_id, "import_summary"),
      wb_ref_run(imported$run$run_id, "import_manifest"),
      wb_ref_run(imported$run$run_id, "rejected_sets")
    ),
    name_prefix = "parse_cytokine_final_signatures",
    description = "Final accepted individual-cytokine transcriptional signature feature sets.",
    created_by = created_by,
    fail_fast = TRUE
  )
  import_manifest_artifact_id <- artifact_id_for(import_promotion, imported$run$run_id, "import_manifest")

  stopifnot(
    identical(as.character(wb_cells(bundle, cell_set = imported_set)), imported_cells_before),
    identical(sort(as.character(wb_cell_sets(bundle))), cell_sets_before),
    identical(as.character(wb_cells(bundle, cell_set = analysis_set)), analysis_cells)
  )

  complete_status <- merge(
    discovery_summary,
    confirmation,
    by = c("signature_id", "cytokine", "direction"),
    all.x = TRUE, sort = FALSE, suffixes = c("_discovery", "_ucell")
  )
  complete_status$final_status[is.na(complete_status$final_status)] <- "not_advanced_from_discovery"
  complete_status <- complete_status[order(match(complete_status$cytokine, single_conditions),
                                             match(complete_status$direction, c("up", "down"))), , drop = FALSE]
  top_gene_map <- tapply(
    discovery_members$feature_id,
    discovery_members$signature_id,
    function(x) paste(utils::head(x, 8L), collapse = ", ")
  )
  complete_status$top_genes <- unname(top_gene_map[complete_status$signature_id])

  write_tsv_atomic(sample_scores, file.path(report_dir, "ucell_sample_medians.tsv"))
  write_tsv_atomic(cell_quantiles, file.path(report_dir, "ucell_cell_distribution_quantiles.tsv"))
  write_tsv_atomic(condition_profile, file.path(report_dir, "ucell_condition_specificity_profile.tsv"))
  write_tsv_atomic(confirmation, file.path(report_dir, "ucell_transfer_validation.tsv"))
  write_tsv_atomic(complete_status, file.path(report_dir, "all_signature_status.tsv"))
  write_tsv_atomic(state_summary, file.path(report_dir, "parse_state_summary.tsv"))
  write_tsv_atomic(state_by_sample, file.path(report_dir, "parse_state_composition_by_sample.tsv"))
  write_tsv_atomic(state_composition_by_condition,
                   file.path(report_dir, "parse_state_composition_by_condition.tsv"))
  write_json_atomic(interfaces, file.path(report_dir, "operation_interfaces.json"))
  write_json_atomic(scoring_dag, file.path(report_dir, "ucell_session_dag.json"))
  write_json_atomic(list(
    promoted = scoring_promotion$promoted,
    snapshot_id = if (is.null(scoring_promotion$snapshot_id)) NA_character_ else scoring_promotion$snapshot_id,
    commit_id = if (is.null(scoring_promotion$commit_id)) NA_character_ else scoring_promotion$commit_id,
    durable_session_state = if (is.null(scoring_promotion$session$durable$state)) "promoted" else scoring_promotion$session$durable$state
  ), file.path(report_dir, "ucell_promotion.json"))
  write_json_atomic(import_dag, file.path(report_dir, "feature_set_import_session_dag.json"))
  write_json_atomic(list(
    promoted = import_promotion$promoted,
    snapshot_id = if (is.null(import_promotion$snapshot_id)) NA_character_ else import_promotion$snapshot_id,
    commit_id = if (is.null(import_promotion$commit_id)) NA_character_ else import_promotion$commit_id,
    durable_session_state = if (is.null(import_promotion$session$durable$state)) "promoted" else import_promotion$session$durable$state
  ), file.path(report_dir, "feature_set_import_promotion.json"))

  ordered_candidates <- candidate_info[order(match(candidate_info$cytokine, single_conditions),
                                                match(candidate_info$direction, c("up", "down"))), , drop = FALSE]
  row_order <- ordered_candidates$signature_id
  profile_matrix <- matrix(NA_real_, nrow = length(row_order), ncol = length(single_conditions),
                           dimnames = list(row_order, single_conditions))
  for (i in seq_len(nrow(condition_profile))) {
    profile_matrix[condition_profile$signature_id[[i]], condition_profile$evaluated_condition[[i]]] <-
      condition_profile$oriented_shift_from_control[[i]]
  }
  draw_heatmap(profile_matrix, file.path(report_dir, "ucell_specificity_heatmap.png"), "png")
  draw_heatmap(profile_matrix, file.path(report_dir, "ucell_specificity_heatmap.pdf"), "pdf")
  plot_info <- merge(candidate_info, confirmation[, c("signature_id", "final_status")],
                     by = "signature_id", all.x = TRUE, sort = FALSE)
  draw_distributions(sample_scores, plot_info, file.path(report_dir, "ucell_signature_distributions.pdf"))

  readiness <- wb_finalize_recipient_bundle(bundle)
  stopifnot(isTRUE(readiness$ready))
  write_json_atomic(readiness, file.path(report_dir, "recipient_readiness.json"))
  wb_close_bundle(bundle)
  bundle <- NULL

  validation <- wb_bundle_validation_report(bundle_dir, level = "recipient", progress = TRUE)
  stopifnot(isTRUE(validation$ok))
  write_json_atomic(validation, file.path(report_dir, "validation.json"))
  recipient <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
  manifest_after <- wb_manifest(recipient)
  stopifnot(
    identical(as.character(wb_cells(recipient, cell_set = imported_set)), imported_cells_before),
    identical(as.character(wb_cells(recipient, cell_set = analysis_set)), analysis_cells),
    identical(sort(as.character(wb_cell_sets(recipient))), cell_sets_before)
  )
  wb_close_bundle(recipient)
  write_json_atomic(manifest_before, file.path(report_dir, "manifest_before.json"))
  write_json_atomic(manifest_after, file.path(report_dir, "manifest.json"))

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
    sum(accepted_summary$direction == "up"), " are induced signatures and ",
    sum(accepted_summary$direction == "down"), " are repressed signatures. No unsupported cytokine/direction ",
    "was forced into the accepted library.\n\n",
    markdown_table(status_counts), "\n\n",
    "## Population boundaries\n\n",
    "All **251,599 imported cells** remain in `", imported_set, "`. Discovery and scoring used only the ",
    "**154,209 QC-accepted control/single cells** in `", analysis_set, "` (81 samples: 3 control and ",
    "78 single-cytokine samples). Combined-stimulation cells were retained in the bundle but were not used. ",
    "Bdev was not read or used anywhere in discovery, validation, or scoring.\n\n",
    "## Analysis and validation\n\n",
    "Parse-only states were constructed from control-derived HVGs and unlabeled PCA/Leiden clustering. Raw counts ",
    "were summed by biological sample, then limma-voom tested two contrasts per cytokine: target versus control ",
    "(response) and target versus the equal-weight mean of the other 25 single cytokines (specificity). Up and down ",
    "programs were treated separately. Genes had to pass both contrasts, effect-size cutoffs, FDR 0.05, and ",
    "direction stability in all three leave-one-replicate-out fits. Compact size was selected by held-out, ",
    "sample-level rank performance and fold membership stability. State composition is reported separately in ",
    "`parse_state_composition_by_sample.tsv` and `parse_state_composition_by_condition.tsv`; within-state response/",
    "specificity concordance is retained for every candidate and used in ranking, so composition shifts are not ",
    "silently presented as expression evidence.\n\n",
    "Candidates were then scored per cell by scWorkbench UCell on raw RNA counts (`max_rank=1500`) without graph ",
    "smoothing. Biological-sample medians—not cells as replicates—were the primary confirmation unit. Raw scores ",
    "were compared only within each signature; they were not treated as directly comparable across signatures. ",
    "Final acceptance required both response and specificity UCell AUC >=0.60 and balanced AUC >=0.65.\n\n",
    "## Final accepted signatures\n\n",
    markdown_table(accepted_summary[, c(
      "cytokine", "direction", "genes", "discovery_cv_response_auc",
      "discovery_cv_specificity_auc", "ucell_response_auc", "ucell_specificity_auc"
    )]), "\n\n",
    "## Specificity diagnostics\n\n",
    "![Oriented UCell specificity heatmap](ucell_specificity_heatmap.png)\n\n",
    "Each heatmap value is the direction-oriented shift of the condition's median sample-level UCell score from ",
    "control. A specific signature should peak in its own cytokine column. Multi-page sample-median distributions ",
    "are in `ucell_signature_distributions.pdf`; exact sample medians, cell-level descriptive quantiles, off-target ",
    "profiles, and all rejected/unsupported statuses are provided as TSV files.\n\n",
    "## Durable outputs\n\n",
    "The final accepted gene sets are `final_accepted_feature_sets.tsv` and `.gmt` and are registered read-only in ",
    "the scWorkbench namespace `parse_cytokine_signatures`. ",
    "Per-cell score families remain in the durable scWorkbench session rather than being promoted: promoting ",
    "millions of family rows makes this scWorkbench build embed an unsafe row-order list in artifact provenance. ",
    "Promoted score summaries are `", paste(score_summary_artifact_ids, collapse = "`, `"),
    "`; the final feature-set import manifest is `", import_manifest_artifact_id, "`. ",
    "The source discovery directory is `", discovery_dir, "`. Bundle recipient validation finished with status `",
    validation$status, "`.\n"
  )
  write_lines_atomic(report, file.path(report_dir, "REPORT.md"))
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
  write_json_atomic(list(
    status = "complete",
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    final_accepted_signatures = nrow(accepted_summary),
    cytokines_with_final_signature = length(unique(accepted_summary$cytokine)),
    final_up_signatures = sum(accepted_summary$direction == "up"),
    final_down_signatures = sum(accepted_summary$direction == "down"),
    per_cell_score_families_promoted = FALSE,
    score_summary_artifact_ids = as.list(score_summary_artifact_ids),
    feature_set_import_manifest_artifact_id = import_manifest_artifact_id,
    imported_cells_retained = length(imported_cells_before),
    analysis_cells = length(analysis_cells),
    combined_cells_used = FALSE,
    bdev_used = FALSE,
    bundle_state_fingerprint = manifest_after$current_state$fingerprint,
    validation_status = validation$status
  ), file.path(report_dir, "completion.json"))
}

main()
