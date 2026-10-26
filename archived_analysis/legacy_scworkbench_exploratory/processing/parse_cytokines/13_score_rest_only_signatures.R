#!/usr/bin/env Rscript
# Native scWorkbench scoring and staged feature-set import. No promotion.
source("src/parse_rest_only_workflow.R")
main <- function() {
  stopifnot(file.exists(file.path(report, "discovery.json")), file.exists(file.path(report, "preparation.json")),
    !file.exists(file.path(report, "scoring.json")), isFALSE(cfg$publication$promote_results))
  parent_check()
  discovery_path <- file.path(discovery_dir, "discovery.rds")
  discovery <- readRDS(discovery_path)
  info <- discovery$summary
  members <- discovery$members
  ids <- info$signature_id[info$discovery_status == "discovery_pass"]
  signatures <- setNames(lapply(ids, function(id) members$feature_id[members$signature_id == id]), ids)
  stopifnot(length(ids) == 47L, all(lengths(signatures) >= cfg$analysis$minimum_genes),
    !anyDuplicated(members[c("signature_id", "feature_id")]))
  b <- wb_open_bundle(cfg$bundle_dir, validate = FALSE)
  on.exit(wb_close_bundle(b), add = TRUE)
  before <- wb_manifest(b)$current_state$fingerprint
  cells <- wb_cells(b, "all")
  md <- wb_meta(b, cells = cells, cols = list(wb_ref_metadata_field("core:cell_id"),
    "sample_id", "condition", "replicate", "treatment_category"))
  stopifnot(length(cells) == 154209L, identical(md$cell_id, cells), !anyNA(md))
  write_json(wb_op_interface("ucell", bundle = b), "native_ucell_interface.json")
  workspace <- file.path(session_root, "rest_only_ucell")
  s <- if (dir.exists(workspace)) wb_session_resume(b, workspace) else
    wb_session_create(b, workspace_root = session_root, session_id = "rest_only_ucell")
  on.exit(try(wb_close_bundle(wb_session_view(s)), silent = TRUE), add = TRUE)
  runs <- s$runs
  score_runs <- names(runs)[vapply(runs, function(r) identical(r$op_name, "ucell"), logical(1))]
  # A completed scoring checkpoint can also be recognized by its typed output;
  # this avoids relying on administrative run-name spelling across runtimes.
  if (!length(score_runs) && length(runs)) score_runs <- names(runs)[vapply(runs,
    function(r) "signature_scores" %in% names(r$outputs), logical(1))]
  if (!length(score_runs)) {
    log_event("Native UCell scoring", length(ids), "new gene sets over", length(cells), "cells")
    rr <- wb_session_run_op(s, b, "ucell", inputs = list(matrix = "RNA_counts", cell_set = "all",
      min_cells = 10L, features = signatures, name = cfg$namespace, max_rank = cfg$analysis$max_rank,
      ncores = 4L, parallel_backend = "auto", chunk_size = cfg$analysis$chunk_size,
      missing_features_policy = "error"),
      notes = "Rest-only Parse candidate scoring. No treatment-control comparison, graph smoothing or cross-signature calibration. Outputs remain in this durable session until user review.")
    s <- rr$session
    score_run <- rr$run
  } else {
    stopifnot(length(score_runs) == 1L)
    log_event("Reusing completed native UCell checkpoint")
    score_run <- s$runs[[score_runs]]
  }
  score <- score_run$outputs$signature_scores
  native_summary <- score_run$outputs$signature_summary$table
  stopifnot(length(score$value) == length(cells) * length(ids),
    setequal(unique(score$key), ids), all(is.finite(score$value)), all(score$value >= 0 & score$value <= 1),
    setequal(native_summary$signature, ids), all(native_summary$n_features_missing == 0L),
    all(native_summary$n_cells == length(cells)))
  for (id in ids) {
    q <- native_summary[native_summary$signature == id, ]
    stopifnot(q$n_features_matched == length(signatures[[id]]),
      setequal(strsplit(q$matched_features, ";", fixed = TRUE)[[1]], signatures[[id]]))
  }
  score_path <- file.path(work, "ucell_scores.rds")
  saveRDS(score, score_path, compress = "gzip")
  write_tsv(native_summary, "native_ucell_summary.tsv")
  sample_md <- unique(md[, c("sample_id", "condition", "replicate", "treatment_category")])
  stopifnot(nrow(sample_md) == 81L, all(table(sample_md$condition) == 3L),
    sum(md$treatment_category == "single") == 149393L)
  singles <- sort(unique(sample_md$condition[sample_md$treatment_category == "single"]))
  sample_rows <- profile_rows <- quantile_rows <- confirmation_rows <- list()
  for (id in ids) {
    entry <- info[info$signature_id == id, ]
    sign <- if (entry$direction == "up") 1 else -1
    rows <- score$key == id
    score_cells <- score$cell_id[rows]
    stopifnot(!anyDuplicated(score_cells), setequal(score_cells, cells))
    values <- as.numeric(score$value[rows])[match(cells, score_cells)]
    stopifnot(length(values) == length(cells), all(is.finite(values)))
    medians <- tapply(values, md$sample_id, median)
    local <- sample_md[match(names(medians), sample_md$sample_id), ]
    stopifnot(!anyNA(local), nrow(local) == 81L)
    local$signature_id <- id
    local$cytokine <- entry$cytokine
    local$direction <- entry$direction
    local$sample_median_ucell <- as.numeric(medians)
    local$oriented_sample_median <- sign * local$sample_median_ucell
    sample_rows[[length(sample_rows) + 1L]] <- local
    single <- local$treatment_category == "single"
    target <- single & local$condition == entry$cytokine
    rest <- single & !target
    stopifnot(sum(target) == 3L, sum(rest) == 75L)
    rest_auc <- auc(local$oriented_sample_median[single], local$condition[single] == entry$cytokine)
    confirmation_rows[[length(confirmation_rows) + 1L]] <- data.frame(signature_id = id,
      ucell_rest_auc = rest_auc, ucell_check = if (rest_auc >= cfg$analysis$ucell_rest_auc_minimum) "PASS" else "FAIL",
      target_sample_median = median(local$sample_median_ucell[target]),
      target_sample_iqr = IQR(local$sample_median_ucell[target]),
      rest_sample_median = median(local$sample_median_ucell[rest]),
      rest_sample_iqr = IQR(local$sample_median_ucell[rest]),
      zero_sample_medians_of_78_singles = sum(local$sample_median_ucell[single] == 0))
    for (condition in singles) {
      use <- local$condition == condition
      profile_rows[[length(profile_rows) + 1L]] <- data.frame(signature_id = id,
        cytokine = entry$cytokine, direction = entry$direction, condition = condition,
        is_target = condition == entry$cytokine, raw_median = median(local$sample_median_ucell[use]),
        sample_iqr = IQR(local$sample_median_ucell[use]),
        oriented_difference_from_other_single_median = sign * (median(local$sample_median_ucell[use]) - median(local$sample_median_ucell[rest])))
    }
    for (group in c("target", "other_single")) {
      use <- if (group == "target") md$condition == entry$cytokine else
        md$treatment_category == "single" & md$condition != entry$cytokine
      q <- quantile(values[use], c(.05, .25, .5, .75, .95), names = FALSE)
      quantile_rows[[length(quantile_rows) + 1L]] <- data.frame(signature_id = id, group = group,
        cells = sum(use), p05 = q[1], p25 = q[2], median = q[3], p75 = q[4], p95 = q[5], iqr = q[4] - q[2])
    }
  }
  confirmation <- do.call(rbind, confirmation_rows)
  final <- merge(info, confirmation, by = "signature_id", all.x = TRUE, sort = FALSE)
  final <- final[match(info$signature_id, final$signature_id), ]
  final$final_status <- ifelse(final$discovery_status != "discovery_pass", final$discovery_status,
    ifelse(final$ucell_check == "PASS", "passes_rest_only_checks", "ucell_not_confirmed"))
  stopifnot(!anyNA(final$final_status))
  selected_ids <- final$signature_id[final$final_status == "passes_rest_only_checks"]
  write_tsv(final, "all_signature_status.tsv")
  write_tsv(do.call(rbind, sample_rows), "ucell_sample_medians.tsv")
  write_tsv(do.call(rbind, profile_rows), "condition_profiles.tsv")
  write_tsv(do.call(rbind, quantile_rows), "cell_distribution_quantiles.tsv")
  write_sets(members, selected_ids, "rest_only_gene_sets")

  # Register the review-ready collection only in the private session overlay.
  # User instructions require a report before promoting any results.
  csv_path <- file.path(report, "rest_only_gene_sets.csv")
  export <- read_tsv(file.path(report, "rest_only_gene_sets.tsv"))
  write.csv(export, csv_path, row.names = FALSE)
  write_json(wb_op_interface("feature_set_import", bundle = b), "native_feature_set_import_interface.json")
  stage <- wb_session_view(s)
  registered <- wb_feature_sets(stage, namespace = cfg$namespace, matrix_id = "RNA_counts")
  if (!nrow(registered) && length(selected_ids)) {
    rr <- wb_session_run_op(s, b, "feature_set_import", inputs = list(matrix = "RNA_counts", cell_set = "all",
      source_path = csv_path, source_format = "table", namespace = cfg$namespace,
      collection = "individual_cytokine_rest_only", species = "Mus musculus", organism = "Mus musculus",
      gene_id_type = "gene_symbol", set_id_col = "feature_set", set_name_col = "set_name",
      gene_col = "feature_id", description_col = "description", split_by_collection = FALSE,
      overwrite = FALSE, strict = TRUE, created_by = "jtg", source_id = "parse_rest_only_signatures",
      source_version = cfg$run_id, origin = "derived", is_readonly = TRUE),
      notes = "Stage new rest-only gene sets passing discovery and sample-median UCell checks. Await report review; do not promote.")
    s <- rr$session
  }
  stage <- wb_session_view(s)
  registered <- wb_feature_sets(stage, namespace = cfg$namespace, matrix_id = "RNA_counts")
  stopifnot(setequal(registered$feature_set, selected_ids))
  for (id in selected_ids) stopifnot(identical(wb_feature_set_get(stage, id, "RNA_counts"),
    export$feature_id[export$feature_set == id]))
  stopifnot(identical(before, wb_manifest(b)$current_state$fingerprint),
    nrow(wb_feature_sets(b, namespace = cfg$namespace)) == 0L)
  write_json(wb_session_dag(s, bundle = b, include_artifacts = FALSE, include_active = TRUE), "native_session_dag.json")
  save_provenance("scoring", c(discovery_path, file.path(report, "discovery_gene_sets.tsv"),
    csv_path, "processing/parse_cytokines/13_score_rest_only_signatures.R"),
    list(scored_sets = length(ids), scored_cells_per_set = length(cells), scores = length(score$value),
      validation_single_samples = 78L, control_auc_computed = FALSE, combined_auc_computed = FALSE,
      passed_sets = length(selected_ids), up = sum(final$direction == "up" & final$final_status == "passes_rest_only_checks"),
      down = sum(final$direction == "down" & final$final_status == "passes_rest_only_checks"),
      cytokines = length(unique(final$cytokine[final$final_status == "passes_rest_only_checks"])),
      score_rds = score_path, score_sha256 = wb_source_sha256(score_path),
      session_workspace = workspace, staged_bundle = stage$dir, namespace_staged = TRUE,
      result_promotions = 0L, input_bundle_state_unchanged = TRUE,
      raw_down_score_interpretation = "Higher expression of genes lower in target than other stimuli; score is not inverted except for AUC and oriented differences."))
  log_event("Scoring complete;", length(selected_ids), "sets pass; results staged for report review")
}
main()
