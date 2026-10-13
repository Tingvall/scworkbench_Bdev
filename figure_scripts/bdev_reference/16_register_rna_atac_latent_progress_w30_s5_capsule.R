#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 16_register_rna_atac_latent_progress_w30_s5_capsule.R",
  "<bdev_bundle> <analysis_script> <output_dir>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
analysis_script <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
})

required <- c("window_distance_and_slope.tsv", "window_membership.tsv.gz",
  "bb1_bb2_regions.tsv", "trajectory_cells_ordered.tsv",
  "rna_vs_atac_relative_latent_distance_w30_s5.png",
  "rna_vs_atac_relative_latent_distance_w30_s5.pdf",
  "rna_vs_atac_local_slope_w30_s5.png",
  "rna_vs_atac_local_slope_w30_s5.pdf", "diagnostics.json", "REPORT.md")
stopifnot(all(file.exists(file.path(output_dir, required))))
diagnostics <- read_json(file.path(output_dir, "diagnostics.json"),
  simplifyVector = TRUE)
stopifnot(identical(diagnostics$status, "complete"),
  identical(diagnostics$bundle_id, "bdev_reference"),
  identical(as.integer(diagnostics$trajectory_cells), 12642L),
  identical(as.integer(diagnostics$windows$size), 30L),
  identical(as.integer(diagnostics$windows$stride), 5L))

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
stopifnot(identical(wb_manifest(bundle)$current_state$fingerprint,
  diagnostics$bundle_fingerprint))
catalog <- wb_operation_catalog_identity()

capsule <- wb_capsule_register(
  bundle = bundle,
  name = "Bdev RNA versus ATAC relative latent progression w30 s5 v1",
  description = paste(
    "Thirty-cell, five-cell-stride comparison of cumulative RNA-PCA and",
    "ATAC-LSI centroid path distances and their local slope along Bdev pt."),
  inputs = list(
    rna_pca = wb_ref_entity("space", "pca"),
    atac_lsi = wb_ref_entity("space", "lsi"),
    pseudotime = wb_ref_entity("artifact",
      "cell_vector_20260914_162517_7q4y8xbb"),
    bb_up = wb_ref_metadata_field("meta:000053")),
  script = analysis_script,
  outputs = list(
    window_metrics = list(path = file.path(output_dir,
      "window_distance_and_slope.tsv"), role = "metrics_table",
      media_type = "text/tab-separated-values"),
    window_membership = list(path = file.path(output_dir,
      "window_membership.tsv.gz"), role = "membership_table",
      media_type = "application/gzip"),
    bb_regions = list(path = file.path(output_dir, "bb1_bb2_regions.tsv"),
      role = "metrics_table", media_type = "text/tab-separated-values"),
    trajectory_cells = list(path = file.path(output_dir,
      "trajectory_cells_ordered.tsv"), role = "input_audit_table",
      media_type = "text/tab-separated-values"),
    distance_png = list(path = file.path(output_dir,
      "rna_vs_atac_relative_latent_distance_w30_s5.png"), role = "figure",
      media_type = "image/png"),
    distance_pdf = list(path = file.path(output_dir,
      "rna_vs_atac_relative_latent_distance_w30_s5.pdf"), role = "figure",
      media_type = "application/pdf"),
    slope_png = list(path = file.path(output_dir,
      "rna_vs_atac_local_slope_w30_s5.png"), role = "figure",
      media_type = "image/png"),
    slope_pdf = list(path = file.path(output_dir,
      "rna_vs_atac_local_slope_w30_s5.pdf"), role = "figure",
      media_type = "application/pdf"),
    diagnostics = list(path = file.path(output_dir, "diagnostics.json"),
      role = "diagnostics", media_type = "application/json"),
    report = list(path = file.path(output_dir, "REPORT.md"), role = "report",
      media_type = "text/markdown")),
  parameters = list(window_size = 30L, window_stride = 5L,
    trajectory_rule = "non-missing authoritative pt",
    cell_order = "pt ascending; cell_id tie-breaker",
    rna_space = "pca", rna_dimensions = 1:50,
    atac_space = "lsi", atac_dimensions = 1:50,
    centroid_summary = "arithmetic mean per latent dimension",
    distance_metric = "Euclidean",
    distance_path = "cumulative consecutive-window centroid distance",
    standardization = "divide by modality endpoint cumulative distance",
    slope = "centered delta RNA relative distance / centered delta ATAC relative distance; one-sided endpoints",
    bb_regions = "source pt min/max and median for bb_up up1/up2"),
  command = c("Rscript", "--vanilla", basename(analysis_script),
    "<bdev_bundle>", "<output_dir>", "<expected_bundle_fingerprint>"),
  runtime = list(R = R.version.string,
    scworkbench = as.character(packageVersion("scworkbench")),
    container = "scworkbench-bdev-signac.sif"),
  executed_at = format(file.info(file.path(output_dir, "diagnostics.json"))$mtime,
    "%Y-%m-%dT%H:%M:%S%z", tz = "UTC"),
  created_by = "Codex for jtg", purpose = "custom_visualization",
  rationale = paste(
    "No registered operation or plot contract computes cross-modality",
    "sliding-window cumulative centroid distance and its local derivative;",
    "the outputs are dataset-specific terminal figures and summaries."),
  catalog_review = list(catalog_identity = catalog$scheme,
    catalog_version = paste0("sha256:", catalog$sha256),
    reviewed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", tz = "UTC"),
    candidate_op_names = c("space_path_define", "space_path_score",
      "trajectory_space_path_segment"), decision = "retain_as_capsule"),
  family_hint = "bdev-cross-modality-latent-window-progress-w30-s5",
  limitations = list(
    "All 50 imported dimensions are used without coordinate rescaling.",
    "Local finite-difference slopes depend on representation geometry, window size, and stride.",
    "The comparison is descriptive and does not establish causal lead-lag."))

receipt <- list(status = "registered", capsule_id = capsule$capsule_id,
  capsule_name = capsule$name, registered_at = capsule$registered_at,
  bundle_fingerprint_after = wb_manifest(bundle)$current_state$fingerprint)
write_json(receipt, file.path(output_dir, "capsule_registration.json"),
  pretty = TRUE, auto_unbox = TRUE)
cat("Registered capsule", capsule$capsule_id, "\n")
