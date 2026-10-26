#!/usr/bin/env Rscript

# Runtime and contract smoke test only. Train for one epoch on deterministic
# small reference/query subsets, project the query, create a model-bound union,
# and fit a display UMAP. Nothing is promoted.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 08_peakvi_mapping_smoke.R BUNDLE_DIR REPORT_DIR SESSION_ROOT", call. = FALSE)
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
report_dir <- normalizePath(args[[2L]], mustWork = FALSE)
session_root <- normalizePath(args[[3L]], mustWork = FALSE)
if (file.exists(report_dir) || dir.exists(report_dir)) {
  stop("Smoke report directory already exists; refusing to overwrite: ", report_dir, call. = FALSE)
}

options(width = 160, warn = 1)
suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
})

log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}

write_json <- function(value, name) {
  jsonlite::write_json(
    value, file.path(report_dir, name), pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null", digits = NA
  )
}

even_sample <- function(x, n) {
  x[unique(as.integer(round(seq.int(1L, length(x), length.out = min(n, length(x))))))]
}

completed_run <- function(session, op_name) {
  ids <- names(session$runs)[vapply(
    session$runs, function(x) identical(x$op_name, op_name), logical(1)
  )]
  if (length(ids) > 1L) {
    stop("Multiple completed runs found for operation ", op_name, call. = FALSE)
  }
  if (!length(ids)) return(NULL)
  list(session = session, run = session$runs[[ids]])
}

main <- function() {
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(session_root, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    if (!file.exists(file.path(report_dir, "smoke.json"))) {
      unlink(report_dir, recursive = TRUE)
    }
  }, add = TRUE)

  bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
  on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  manifest_before <- wb_manifest(bundle)
  stopifnot(identical(manifest_before$bundle_id, "bdev_integrated_ko_mapping_v1"))
  reference_cells <- wb_cells(bundle, "set:bdev_reference")
  query_cells <- wb_cells(bundle, "set:ko_query")
  metadata <- wb_meta(bundle, cells = c(reference_cells, query_cells))
  ref_smoke <- even_sample(reference_cells, 256L)

  query_md <- metadata[match(query_cells, metadata$cell_id), , drop = FALSE]
  query_groups <- interaction(query_md$sample_id, query_md$genotype,
                              query_md$knockout_target, drop = TRUE)
  query_smoke <- unlist(lapply(split(query_cells, query_groups), even_sample, n = 64L),
                        use.names = FALSE)
  query_smoke <- unique(query_smoke)
  stopifnot(length(ref_smoke) == 256L, length(query_smoke) > 0L,
            !anyDuplicated(c(ref_smoke, query_smoke)),
            all(ref_smoke %in% reference_cells), all(query_smoke %in% query_cells))

  for (op in c("peakvi_fit", "peakvi_project", "model_bound_latent_union", "umap_fit")) {
    write_json(wb_op_interface(op, bundle = bundle), paste0("interface_", op, ".json"))
  }
  writeLines(ref_smoke, file.path(report_dir, "reference_smoke_cells.txt"))
  writeLines(query_smoke, file.path(report_dir, "query_smoke_cells.txt"))

  workspace <- file.path(session_root, "peakvi_smoke_v1")
  session <- if (dir.exists(workspace)) {
    log_event("Resuming durable smoke session", workspace)
    wb_session_resume(bundle, workspace)
  } else {
    wb_session_create(bundle, workspace_root = session_root,
                      session_id = "peakvi_smoke_v1")
  }
  fit <- completed_run(session, "peakvi_fit")
  if (is.null(fit)) {
    log_event("Running one-epoch CPU PeakVI fit on", length(ref_smoke), "Bdev cells")
    fit <- wb_session_run_op(
      session, bundle, "peakvi_fit",
      inputs = list(
        matrix = "ATAC_bdev_mm10_peaks", cell_ids = ref_smoke, min_cells = 10L,
        n_latent = 5L, n_layers = 1L, n_hidden = 64L, dropout = 0.1,
        max_epochs = 1L, lr = 0.001, early_stopping = FALSE, seed = 1729L,
        use_gpu = FALSE, space_id_out = "peakvi.smoke.reference"
      ),
      notes = "Runtime smoke only: deterministic 256-cell Bdev subset, one epoch, CPU, no promotion."
    )
  } else {
    log_event("Reusing completed PeakVI fit", fit$run$run_id)
  }
  session <- fit$session

  projected <- completed_run(session, "peakvi_project")
  if (is.null(projected)) {
    log_event("Projecting", length(query_smoke), "KO smoke cells through the frozen smoke model")
    projected <- wb_session_run_op(
      session, bundle, "peakvi_project",
      inputs = list(
        matrix = "ATAC_bdev_mm10_peaks", model = wb_ref_run(fit$run$run_id, "model"),
        cell_ids = query_smoke, min_cells = 10L, use_gpu = FALSE,
        space_id_out = "peakvi.smoke.ko_projected"
      ),
      notes = "Runtime smoke only: query projection without reference refit; no promotion."
    )
  } else {
    log_event("Reusing completed KO projection", projected$run$run_id)
  }
  session <- projected$session

  smoke_cells <- c(ref_smoke, query_smoke)
  joined <- completed_run(session, "model_bound_latent_union")
  if (is.null(joined)) {
    log_event("Combining unchanged smoke reference and projected coordinates")
    joined <- wb_session_run_op(
      session, bundle, "model_bound_latent_union",
      inputs = list(
        reference_space = "peakvi.smoke.reference",
        query_space = "peakvi.smoke.ko_projected",
        model = wb_ref_run(fit$run$run_id, "model"),
        projection_summary = wb_ref_run(projected$run$run_id, "fit_summary"),
        cell_ids = smoke_cells, min_cells = 10L,
        space_id_out = "peakvi.smoke.joint"
      ),
      notes = "Runtime smoke only: model-bound row union; no coordinate transformation or promotion."
    )
  } else {
    log_event("Reusing completed model-bound union", joined$run$run_id)
  }
  session <- joined$session

  display <- completed_run(session, "umap_fit")
  if (is.null(display)) {
    log_event("Fitting smoke display UMAP")
    display <- wb_session_run_op(
      session, bundle, "umap_fit",
      inputs = list(
        space = "peakvi.smoke.joint", cell_ids = smoke_cells, min_cells = 10L,
        dims = 1:5, n_neighbors = 15L, min_dist = 0.3, metric = "euclidean",
        nn_method = "annoy", seed = 1729L,
        space_id_out = "umap.peakvi.smoke.joint"
      ),
      notes = "Runtime smoke display only; not biological mapping evidence and not promoted."
    )
  } else {
    log_event("Reusing completed smoke UMAP", display$run$run_id)
  }
  session <- display$session

  # Unpromoted outputs live in the durable session's staged bundle. Validate
  # that overlay directly; the committed parent bundle must not expose them.
  staged <- wb_open_bundle(file.path(workspace, "stage"), validate = FALSE)
  on.exit(try(wb_close_bundle(staged), silent = TRUE), add = TRUE)
  reference_space <- wb_space(staged, "peakvi.smoke.reference", cells = ref_smoke)
  query_space <- wb_space(staged, "peakvi.smoke.ko_projected", cells = query_smoke)
  joint_space <- wb_space(staged, "peakvi.smoke.joint", cells = smoke_cells)
  umap_space <- wb_space(staged, "umap.peakvi.smoke.joint", cells = smoke_cells)
  project_summary <- projected$run$outputs$fit_summary$table
  union_summary <- joined$run$outputs$union_summary$table
  checks <- list(
    reference_latent_dimensions = identical(dim(reference_space), c(length(ref_smoke), 5L)),
    query_latent_dimensions = identical(dim(query_space), c(length(query_smoke), 5L)),
    joint_latent_dimensions = identical(dim(joint_space), c(length(smoke_cells), 5L)),
    umap_dimensions = identical(dim(umap_space[, c("x", "y"), drop = FALSE]),
                                c(length(smoke_cells), 2L)),
    reference_cell_order = identical(rownames(reference_space), ref_smoke),
    query_cell_order = identical(rownames(query_space), query_smoke),
    joint_cell_order = identical(rownames(joint_space), smoke_cells),
    umap_cell_order = identical(as.character(umap_space$cell_id), smoke_cells),
    reference_values_finite = all(is.finite(as.matrix(reference_space))),
    query_values_finite = all(is.finite(as.matrix(query_space))),
    joint_values_finite = all(is.finite(as.matrix(joint_space))),
    umap_values_finite = all(is.finite(as.matrix(umap_space[, c("x", "y")]))),
    projection_without_refit = identical(project_summary$projection_mode,
                                         "transform_without_refit"),
    projection_summary_dimensions = identical(
      c(project_summary$n_cells, project_summary$n_features, project_summary$n_latent),
      c(length(query_smoke), 97788L, 5L)
    ),
    union_without_transformation = identical(union_summary$coordinate_operation,
                                             "row_bind_without_transformation"),
    union_summary_dimensions = identical(
      c(union_summary$n_reference_cells, union_summary$n_query_cells,
        union_summary$n_cells, union_summary$n_dimensions),
      c(length(ref_smoke), length(query_smoke), length(smoke_cells), 5L)
    ),
    session_unpromoted = is.null(session$durable$promotion_receipt)
  )
  stopifnot(all(unlist(checks)))
  manifest_after <- wb_manifest(bundle)
  stopifnot(identical(manifest_before$current_state$fingerprint,
                      manifest_after$current_state$fingerprint))

  write_json(wb_session_dag(session, bundle = bundle, include_artifacts = TRUE,
                            include_active = TRUE), "session_dag.json")
  write_json(list(
    status = "passed",
    scope = "PeakVI runtime/contract smoke test; not a scientific mapping result",
    promoted = FALSE,
    matrix = "ATAC_bdev_mm10_peaks",
    parameters = list(
      reference_cells = length(ref_smoke), query_cells = length(query_smoke),
      n_latent = 5L, n_layers = 1L, n_hidden = 64L, max_epochs = 1L,
      use_gpu = FALSE, seed = 1729L
    ),
    runs = list(fit = fit$run$run_id, project = projected$run$run_id,
                union = joined$run$run_id, umap = display$run$run_id),
    checks = checks,
    bundle_fingerprint_unchanged = TRUE,
    session_workspace = workspace
  ), "smoke.json")
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
  log_event("PeakVI smoke test passed without promotion")
}

main()
