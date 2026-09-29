#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 04_promote_flattened_rank_scores.R",
  "<bdev_bundle> <breadth_intensity_report> <promotion_report_dir>",
  "<expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
score_report <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
expected_fingerprint <- args[[4L]]
if (dir.exists(output_dir)) stop("Promotion report already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete promotion report already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
})

write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n"); flush.console()
}
artifact_id_for <- function(promotion, run_id, output_name) {
  rows <- promotion$promoted[
    promotion$promoted$run_id == run_id &
      promotion$promoted$output == output_name &
      promotion$promoted$status == "ok", , drop = FALSE]
  if (nrow(rows) != 1L || is.na(rows$artifact_id[[1L]]) ||
      !nzchar(rows$artifact_id[[1L]])) {
    stop("Expected one promoted artifact for ", run_id, "/", output_name)
  }
  as.character(rows$artifact_id[[1L]])
}
as_family <- function(x, component, cells, source_sha) {
  stopifnot(is.matrix(x), identical(colnames(x), cells),
    !is.null(rownames(x)), !anyDuplicated(rownames(x)),
    all(is.finite(x)), all(x >= 0), all(x <= 1))
  wb_prim_cell_vector_family(
    key = rep(rownames(x), each = ncol(x)),
    cell_id = rep(colnames(x), times = nrow(x)),
    value = as.numeric(t(x)), value_type = "continuous",
    meta = list(
      score_method = "flattened_rank_breadth_intensity_v1",
      component = component,
      source_report = score_report,
      source_score_rds_sha256 = source_sha,
      rank_cutoff = 1500L, alpha = 0.5,
      detection_matrix = "RNA_counts", rank_matrix = "RNA_counts",
      rank_universe = "all genes in RNA_counts",
      descending_rank = TRUE, tie_method = "average"))
}
assert_family_equal <- function(observed, expected, label) {
  stopifnot(inherits(observed, "wb_primitive"),
    identical(observed$kind, "cell_vector_family"))
  observed_df <- data.frame(key = as.character(observed$key),
    cell_id = as.character(observed$cell_id), value = as.numeric(observed$value))
  expected_df <- data.frame(key = as.character(expected$key),
    cell_id = as.character(expected$cell_id), value = as.numeric(expected$value))
  observed_df <- observed_df[order(observed_df$key, observed_df$cell_id), ]
  expected_df <- expected_df[order(expected_df$key, expected_df$cell_id), ]
  rownames(observed_df) <- NULL; rownames(expected_df) <- NULL
  if (!identical(observed_df, expected_df)) stop(label, " values changed")
  TRUE
}

required <- c("completion.json", "checks.json",
  "breadth_intensity_cell_scores.rds", "score_summary.tsv")
stopifnot(all(file.exists(file.path(score_report, required))))
source_completion <- read_json(file.path(score_report, "completion.json"),
  simplifyVector = TRUE)
source_checks <- read_json(file.path(score_report, "checks.json"),
  simplifyVector = TRUE)
stopifnot(identical(source_completion$status, "complete"),
  identical(source_completion$promoted, FALSE),
  identical(as.integer(source_completion$signatures), 38L),
  identical(as.integer(source_completion$cells), 26723L),
  identical(as.integer(source_completion$parameters$rank_cutoff), 1500L),
  identical(as.numeric(source_completion$parameters$alpha), 0.5),
  all(unlist(source_checks)))
score_path <- file.path(score_report, "breadth_intensity_cell_scores.rds")
source_sha <- wb_source_sha256(score_path)
source_files <- data.frame(
  role = c("score_values", "completion", "checks", "score_summary"),
  path = normalizePath(file.path(score_report, c(
    "breadth_intensity_cell_scores.rds", "completion.json", "checks.json",
    "score_summary.tsv")), mustWork = TRUE), stringsAsFactors = FALSE)
source_files$sha256 <- vapply(source_files$path, wb_source_sha256, character(1))
write_tsv(source_files, file.path(stage_dir, "input_manifest.tsv"))

score_object <- readRDS(score_path)
stopifnot(identical(score_object$parameters$rank_cutoff, 1500L),
  identical(score_object$parameters$alpha, 0.5),
  setequal(names(score_object$scores), c("combined", "breadth", "intensity")))

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest_before <- wb_manifest(bundle)
stopifnot(identical(manifest_before$bundle_id, "bdev_reference"),
  identical(manifest_before$current_state$fingerprint, expected_fingerprint),
  identical(source_completion$bundle_fingerprint, expected_fingerprint))
cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "pt"))
stopifnot(identical(as.character(metadata$cell_id), cells),
  length(cells) == 26723L, sum(is.finite(metadata$pt)) == 12642L,
  sum(is.na(metadata$pt)) == 14081L)

families <- list(
  combined = as_family(score_object$scores$combined, "combined_default", cells,
    source_sha),
  breadth = as_family(score_object$scores$breadth, "breadth_component", cells,
    source_sha),
  intensity = as_family(score_object$scores$intensity, "intensity_component",
    cells, source_sha))

# Identity calibration is a supported native terminal operation. With z-score
# centering and scaling both disabled it returns every input value unchanged.
interface <- wb_op_interface("cell_vector_family_calibrate", bundle = bundle)
write_json(interface, file.path(stage_dir, "identity_operation_interface.json"))
qc_interface <- wb_op_interface("cell_vector_family_qc", bundle = bundle)
write_json(qc_interface, file.path(stage_dir, "qc_operation_interface.json"))

session_root <- file.path(dirname(bundle_dir), "sessions",
  "bdev_default_flattened_rank_scores_v1")
session_id <- "bdev_default_flattened_rank_scores_v1"
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) wb_session_resume(bundle, workspace) else
  wb_session_create(bundle, workspace_root = session_root, session_id = session_id)
if (!is.null(session$durable$promotion_receipt)) {
  stop("Durable session is already promoted but no finalized project report exists")
}
if (length(session$runs)) stop("Promotion session contains unexpected prior runs")

identity_runs <- list(); qc_runs <- list()
for (component in c("combined", "breadth", "intensity")) {
  default_note <- if (component == "combined")
    "Accepted default gene-set score" else "Retained component score"
  log_event("Staging exact", component, "cell-vector family")
  result <- wb_session_run_op(session, bundle, "cell_vector_family_calibrate",
    inputs = list(signal = families[[component]], cell_set = "all",
      min_cells = 1L, method = "zscore", center = FALSE, scale = FALSE,
      na_rm = TRUE),
    notes = paste0("[identity_", component, "] ", default_note,
      "; identity materialization of user-approved flattened-rank ", component,
      " values from ", source_sha,
      "; zscore method with center=FALSE and scale=FALSE."))
  session <- result$session; identity_runs[[component]] <- result$run
  stopifnot(assert_family_equal(result$run$outputs$calibrated_signal,
    families[[component]], component))
  log_event("Running native QC for", component)
  result <- wb_session_run_op(session, bundle, "cell_vector_family_qc",
    inputs = list(
      signal = wb_ref_run(identity_runs[[component]]$run_id,
        "calibrated_signal"),
      cell_set = "all", min_cells = 1L,
      correlation_methods = c("pearson", "spearman"), tie_tol = 1e-8),
    notes = paste0("[qc_", component,
      "] Native QC for the accepted flattened-rank ", component,
      " score family before promotion."))
  session <- result$session; qc_runs[[component]] <- result$run
}

write_json(wb_session_dag(session, bundle = bundle, include_artifacts = TRUE,
  include_active = TRUE), file.path(stage_dir, "session_dag_before_promotion.json"))
for (component in names(identity_runs)) {
  write_tsv(identity_runs[[component]]$outputs$calibration_summary$table,
    file.path(stage_dir, paste0(component, "_identity_summary.tsv")))
  write_tsv(qc_runs[[component]]$outputs$key_summary$table,
    file.path(stage_dir, paste0(component, "_key_qc.tsv")))
  write_tsv(qc_runs[[component]]$outputs$pairwise_summary$table,
    file.path(stage_dir, paste0(component, "_pairwise_qc.tsv")))
}

promote_outputs <- unlist(lapply(c("combined", "breadth", "intensity"),
  function(component) list(
    wb_ref_run(identity_runs[[component]]$run_id, "calibrated_signal"),
    wb_ref_run(identity_runs[[component]]$run_id, "calibration_summary"),
    wb_ref_run(qc_runs[[component]]$run_id, "key_summary"),
    wb_ref_run(qc_runs[[component]]$run_id, "pairwise_summary"))),
  recursive = FALSE)
log_event("Promoting default Combined plus Breadth and Intensity components")
promotion <- wb_session_promote_current(bundle, session,
  outputs = promote_outputs,
  name_prefix = "bdev_flattened_rank",
  description = paste(
    "Accepted Bdev gene-set scores for 38 reviewed cytokine signatures.",
    "Combined flattened rank (R=1500, alpha=0.5) is the project default;",
    "Breadth and Intensity are retained as component cell-vector families.",
    "Raw RNA counts, descending within-cell average ranks, exact values from",
    source_sha),
  created_by = "scworkbench_Bdev flattened-rank promotion",
  fail_fast = TRUE)

artifact_ids <- list()
for (component in names(identity_runs)) {
  artifact_ids[[component]] <- artifact_id_for(promotion,
    identity_runs[[component]]$run_id, "calibrated_signal")
}
promotion_record <- list(
  promoted = promotion$promoted,
  snapshot_id = promotion$snapshot_id,
  commit_id = promotion$commit_id,
  durable_session_state = promotion$session$durable$state,
  session_workspace = workspace,
  default_score = "combined_flattened_rank",
  default_artifact_id = artifact_ids$combined,
  component_artifact_ids = list(breadth = artifact_ids$breadth,
    intensity = artifact_ids$intensity))
saveRDS(promotion_record, file.path(stage_dir, "promotion_receipt.rds"))
write_json(promotion_record, file.path(stage_dir, "promotion_receipt.json"))

for (component in names(artifact_ids)) {
  promoted_family <- wb_art_load(bundle, artifact_ids[[component]],
    expected_kind = "cell_vector_family")
  stopifnot(assert_family_equal(promoted_family, families[[component]],
    paste("promoted", component)))
}
manifest_after <- wb_manifest(bundle)
fingerprint_after <- manifest_after$current_state$fingerprint
stopifnot(!identical(fingerprint_after, expected_fingerprint))

write_json(list(
  status = "promoted", promoted = TRUE,
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  bundle = bundle_dir,
  bundle_fingerprint_before = expected_fingerprint,
  bundle_fingerprint_after = fingerprint_after,
  source_report = score_report, source_score_sha256 = source_sha,
  default_score = list(
    name = "combined_flattened_rank", artifact_id = artifact_ids$combined,
    formula = "0.5 * breadth + 0.5 * intensity", rank_cutoff = 1500L),
  retained_components = list(
    breadth = artifact_ids$breadth, intensity = artifact_ids$intensity),
  signatures = 38L, cells = 26723L,
  exact_pt_defined = 12642L, exact_pt_missing = 14081L,
  identity_materialization = list(operation = "cell_vector_family_calibrate",
    method = "zscore", center = FALSE, scale = FALSE),
  all_promoted_values_exact = TRUE,
  promotion_snapshot_id = promotion$snapshot_id,
  promotion_commit_id = promotion$commit_id,
  output_dir = output_dir), file.path(stage_dir, "completion.json"))

wb_close_bundle(bundle); bundle <- NULL
log_event("Running recipient-level bundle validation")
validation <- wb_bundle_validation_report(bundle_dir, level = "recipient",
  progress = TRUE)
write_json(validation, file.path(stage_dir, "recipient_validation.json"))
stopifnot(validation$status %in% c("ok", "ok_with_warnings"),
  identical(validation$bundle_fingerprint_before,
    validation$bundle_fingerprint_after))

if (!file.rename(stage_dir, output_dir)) stop("Could not finalize promotion report")
log_event("Flattened-rank score promotion complete at", output_dir)
