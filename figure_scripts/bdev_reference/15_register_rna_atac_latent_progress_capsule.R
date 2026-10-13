#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 15_register_rna_atac_latent_progress_capsule.R",
  "<bdev_bundle> <analysis_script> <output_dir>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
analysis_script <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
})

required <- c("window_distance_summary.tsv", "window_membership.tsv.gz",
  "bb1_bb2_highlights.tsv", "trajectory_cells_ordered.tsv",
  "rna_vs_atac_relative_latent_distance.png",
  "rna_vs_atac_relative_latent_distance.pdf", "diagnostics.json", "REPORT.md")
stopifnot(all(file.exists(file.path(output_dir, required))))
diagnostics <- read_json(file.path(output_dir, "diagnostics.json"),
  simplifyVector = TRUE)
stopifnot(identical(diagnostics$status, "complete"),
  identical(diagnostics$bundle_id, "bdev_reference"),
  identical(as.integer(diagnostics$trajectory_cells), 12642L),
  identical(as.integer(diagnostics$windows$size), 200L),
  identical(as.integer(diagnostics$windows$stride), 20L))

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
stopifnot(identical(wb_manifest(bundle)$current_state$fingerprint,
  diagnostics$bundle_fingerprint))

catalog <- wb_operation_catalog_identity()
capsule <- wb_capsule_register(
  bundle = bundle,
  name = "Bdev RNA versus ATAC relative latent progression v1",
  description = paste(
    "Two-hundred-cell sliding-window comparison of cumulative RNA-PCA and",
    "ATAC-LSI centroid path distances along authoritative Bdev pt."),
  inputs = list(
    rna_pca = wb_ref_entity("space", "pca"),
    atac_lsi = wb_ref_entity("space", "lsi"),
    pseudotime = wb_ref_entity("artifact",
      "cell_vector_20260914_162517_7q4y8xbb"),
    bb_up = wb_ref_metadata_field("meta:000053")),
  script = analysis_script,
  outputs = list(
    window_metrics = list(path = file.path(output_dir,
      "window_distance_summary.tsv"), role = "metrics_table",
      media_type = "text/tab-separated-values"),
    window_membership = list(path = file.path(output_dir,
      "window_membership.tsv.gz"), role = "membership_table",
      media_type = "application/gzip"),
    bb_highlights = list(path = file.path(output_dir,
      "bb1_bb2_highlights.tsv"), role = "metrics_table",
      media_type = "text/tab-separated-values"),
    trajectory_cells = list(path = file.path(output_dir,
      "trajectory_cells_ordered.tsv"), role = "input_audit_table",
      media_type = "text/tab-separated-values"),
    figure_png = list(path = file.path(output_dir,
      "rna_vs_atac_relative_latent_distance.png"), role = "figure",
      media_type = "image/png"),
    figure_pdf = list(path = file.path(output_dir,
      "rna_vs_atac_relative_latent_distance.pdf"), role = "figure",
      media_type = "application/pdf"),
    diagnostics = list(path = file.path(output_dir, "diagnostics.json"),
      role = "diagnostics", media_type = "application/json"),
    report = list(path = file.path(output_dir, "REPORT.md"), role = "report",
      media_type = "text/markdown")),
  parameters = list(window_size = 200L, window_stride = 20L,
    trajectory_rule = "non-missing authoritative pt",
    cell_order = "pt ascending; cell_id tie-breaker",
    rna_space = "pca", rna_dimensions = 1:50,
    atac_space = "lsi", atac_dimensions = 1:50,
    centroid_summary = "arithmetic mean per latent dimension",
    distance_metric = "Euclidean",
    distance_path = "cumulative consecutive-window centroid distance",
    standardization = "divide by modality endpoint cumulative distance",
    bb_highlight = "nearest window median pt to bb_up up1/up2 median pt"),
  command = c("Rscript", "--vanilla", basename(analysis_script),
    "<bdev_bundle>", "<output_dir>", "<expected_bundle_fingerprint>",
    "200", "20"),
  runtime = list(R = R.version.string,
    scworkbench = as.character(packageVersion("scworkbench")),
    container = "scworkbench-bdev-signac.sif"),
  executed_at = format(file.info(file.path(output_dir, "diagnostics.json"))$mtime,
    "%Y-%m-%dT%H:%M:%S%z", tz = "UTC"),
  created_by = "Codex for jtg",
  purpose = "custom_visualization",
  rationale = paste(
    "No registered operation or plot contract computes a sliding-window",
    "comparison of cumulative centroid path length across two latent spaces;",
    "the result is a dataset-specific terminal figure and summary."),
  catalog_review = list(
    catalog_identity = catalog$scheme,
    catalog_version = paste0("sha256:", catalog$sha256),
    reviewed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", tz = "UTC"),
    candidate_op_names = c("space_path_define", "space_path_score",
      "trajectory_space_path_segment"),
    decision = "retain_as_capsule"),
  family_hint = "bdev-cross-modality-latent-window-progress",
  limitations = list(
    "All 50 imported dimensions are used without coordinate rescaling.",
    paste("Relative path length depends on the imported representations,",
      "pseudotime ordering, window size, and stride."),
    paste("The comparison is descriptive and does not establish causal lead-lag",
      "between modalities.")))

receipt <- list(status = "registered", capsule_id = capsule$capsule_id,
  capsule_name = capsule$name,
  registered_at = capsule$registered_at,
  bundle_fingerprint_after = wb_manifest(bundle)$current_state$fingerprint)
write_json(receipt, file.path(output_dir, "capsule_registration.json"),
  pretty = TRUE, auto_unbox = TRUE)
cat("Registered capsule", capsule$capsule_id, "\n")
