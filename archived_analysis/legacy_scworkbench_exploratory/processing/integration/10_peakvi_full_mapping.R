#!/usr/bin/env Rscript

# Fit PeakVI on the complete Bdev reference, project all accepted KO query
# cells through the frozen model, create a model-bound union, and fit a display
# UMAP. All outputs remain in a durable, unpromoted review session.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 10_peakvi_full_mapping.R BUNDLE_DIR REPORT_DIR SESSION_ROOT",
       call. = FALSE)
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
report_dir <- normalizePath(args[[2L]], mustWork = FALSE)
session_root <- normalizePath(args[[3L]], mustWork = FALSE)
if (file.exists(report_dir) || dir.exists(report_dir)) {
  stop("Full mapping report directory already exists; refusing to overwrite: ",
       report_dir, call. = FALSE)
}

options(width = 160, warn = 1)
suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
})

write_json <- function(value, name) {
  jsonlite::write_json(
    value, file.path(report_dir, name), pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null", digits = NA
  )
}

completed_run <- function(session, op_name) {
  ids <- names(session$runs)[vapply(
    session$runs, function(x) identical(x$op_name, op_name), logical(1)
  )]
  if (length(ids) > 1L) stop("Multiple runs found for ", op_name, call. = FALSE)
  if (!length(ids)) return(NULL)
  list(session = session, run = session$runs[[ids]])
}

main <- function() {
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(session_root, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    if (!file.exists(file.path(report_dir, "mapping.json"))) {
      unlink(report_dir, recursive = TRUE)
    }
  }, add = TRUE)

  bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
  on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  manifest_before <- wb_manifest(bundle)
  stopifnot(identical(manifest_before$bundle_id, "bdev_integrated_ko_mapping_v1"))
  reference_cells <- wb_cells(bundle, "set:bdev_reference")
  query_cells <- wb_cells(bundle, "set:ko_query")
  all_mapping_cells <- c(reference_cells, query_cells)
  stopifnot(length(reference_cells) == 26723L, length(query_cells) == 23252L,
            !anyDuplicated(all_mapping_cells))

  for (op in c("peakvi_fit", "peakvi_project", "model_bound_latent_union", "umap_fit")) {
    write_json(wb_op_interface(op, bundle = bundle), paste0("interface_", op, ".json"))
  }
  workspace <- file.path(session_root, "peakvi_full_mapping_v1")
  session <- if (dir.exists(workspace)) {
    wb_session_resume(bundle, workspace)
  } else {
    wb_session_create(bundle, workspace_root = session_root,
                      session_id = "peakvi_full_mapping_v1")
  }

  fit <- completed_run(session, "peakvi_fit")
  if (is.null(fit)) {
    fit <- wb_session_run_op(
      session, bundle, "peakvi_fit",
      inputs = list(
        matrix = "ATAC_bdev_mm10_peaks", cell_ids = reference_cells,
        min_cells = 10L, n_latent = 20L, n_layers = 2L, n_hidden = 256L,
        dropout = 0.2, max_epochs = 500L, lr = 0.001,
        early_stopping = TRUE, seed = 1729L, use_gpu = TRUE,
        space_id_out = "peakvi.full.reference"
      ),
      notes = paste(
        "Full Bdev reference PeakVI fit for frozen-model KO projection;",
        "GPU, early stopping, retained in unpromoted review session."
      )
    )
  }
  session <- fit$session

  projected <- completed_run(session, "peakvi_project")
  if (is.null(projected)) {
    projected <- wb_session_run_op(
      session, bundle, "peakvi_project",
      inputs = list(
        matrix = "ATAC_bdev_mm10_peaks",
        model = wb_ref_run(fit$run$run_id, "model"),
        cell_ids = query_cells, min_cells = 10L, use_gpu = TRUE,
        space_id_out = "peakvi.full.ko_projected"
      ),
      notes = "All accepted KO cells projected through the frozen Bdev PeakVI model; no refit."
    )
  }
  session <- projected$session

  joined <- completed_run(session, "model_bound_latent_union")
  if (is.null(joined)) {
    joined <- wb_session_run_op(
      session, bundle, "model_bound_latent_union",
      inputs = list(
        reference_space = "peakvi.full.reference",
        query_space = "peakvi.full.ko_projected",
        model = wb_ref_run(fit$run$run_id, "model"),
        projection_summary = wb_ref_run(projected$run$run_id, "fit_summary"),
        cell_ids = all_mapping_cells, min_cells = 10L,
        space_id_out = "peakvi.full.joint"
      ),
      notes = "Model-bound reference/query row union without coordinate transformation."
    )
  }
  session <- joined$session

  display <- completed_run(session, "umap_fit")
  if (is.null(display)) {
    display <- wb_session_run_op(
      session, bundle, "umap_fit",
      inputs = list(
        space = "peakvi.full.joint", cell_ids = all_mapping_cells,
        min_cells = 10L, dims = 1:20, n_neighbors = 30L,
        min_dist = 0.3, metric = "euclidean", nn_method = "annoy",
        seed = 1729L, space_id_out = "umap.peakvi.full.joint"
      ),
      notes = "Full mapping display UMAP retained for review; no promotion."
    )
  }
  session <- display$session

  staged <- wb_open_bundle(file.path(workspace, "stage"), validate = FALSE)
  on.exit(try(wb_close_bundle(staged), silent = TRUE), add = TRUE)
  reference_space <- wb_space(staged, "peakvi.full.reference", cells = reference_cells)
  query_space <- wb_space(staged, "peakvi.full.ko_projected", cells = query_cells)
  joint_space <- wb_space(staged, "peakvi.full.joint", cells = all_mapping_cells)
  umap_space <- wb_space(staged, "umap.peakvi.full.joint", cells = all_mapping_cells)
  project_summary <- projected$run$outputs$fit_summary$table
  union_summary <- joined$run$outputs$union_summary$table
  checks <- list(
    reference_dimensions = identical(dim(reference_space), c(26723L, 20L)),
    query_dimensions = identical(dim(query_space), c(23252L, 20L)),
    joint_dimensions = identical(dim(joint_space), c(49975L, 20L)),
    umap_dimensions = identical(dim(umap_space[, c("x", "y"), drop = FALSE]),
                                c(49975L, 2L)),
    reference_cell_order = identical(rownames(reference_space), reference_cells),
    query_cell_order = identical(rownames(query_space), query_cells),
    joint_cell_order = identical(rownames(joint_space), all_mapping_cells),
    umap_cell_order = identical(as.character(umap_space$cell_id), all_mapping_cells),
    reference_finite = all(is.finite(as.matrix(reference_space))),
    query_finite = all(is.finite(as.matrix(query_space))),
    joint_finite = all(is.finite(as.matrix(joint_space))),
    umap_finite = all(is.finite(as.matrix(umap_space[, c("x", "y")]))),
    projection_without_refit = identical(project_summary$projection_mode,
                                         "transform_without_refit"),
    projection_cells = identical(project_summary$n_cells, 23252L),
    union_without_transformation = identical(union_summary$coordinate_operation,
                                             "row_bind_without_transformation"),
    session_unpromoted = is.null(session$durable$promotion_receipt)
  )
  stopifnot(all(unlist(checks)))
  manifest_after <- wb_manifest(bundle)
  stopifnot(identical(manifest_before$current_state$fingerprint,
                      manifest_after$current_state$fingerprint))

  write_json(wb_session_dag(session, bundle = bundle, include_artifacts = TRUE,
                            include_active = TRUE), "session_dag.json")
  write_json(list(
    status = "complete_unpromoted",
    scope = "Full Bdev-trained PeakVI fit and frozen-model KO projection",
    promoted = FALSE,
    runs = list(fit = fit$run$run_id, project = projected$run$run_id,
                union = joined$run$run_id, umap = display$run$run_id),
    parameters = fit$run$params,
    checks = checks,
    bundle_fingerprint_unchanged = TRUE,
    session_workspace = workspace
  ), "mapping.json")
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
}

main()
