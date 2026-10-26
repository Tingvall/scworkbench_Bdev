#!/usr/bin/env Rscript

# Five-epoch, full-reference PeakVI benchmark. This measures the real staging
# and training path on every Bdev reference cell but does not project query
# cells or promote any result.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 09_peakvi_full_reference_benchmark.R BUNDLE_DIR REPORT_DIR SESSION_ROOT",
       call. = FALSE)
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
report_dir <- normalizePath(args[[2L]], mustWork = FALSE)
session_root <- normalizePath(args[[3L]], mustWork = FALSE)
if (file.exists(report_dir) || dir.exists(report_dir)) {
  stop("Benchmark report directory already exists; refusing to overwrite: ",
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

main <- function() {
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(session_root, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    if (!file.exists(file.path(report_dir, "benchmark.json"))) {
      unlink(report_dir, recursive = TRUE)
    }
  }, add = TRUE)

  bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
  on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  manifest_before <- wb_manifest(bundle)
  stopifnot(identical(manifest_before$bundle_id, "bdev_integrated_ko_mapping_v1"))
  reference_cells <- wb_cells(bundle, "set:bdev_reference")
  stopifnot(length(reference_cells) == 26723L)

  write_json(wb_op_interface("peakvi_fit", bundle = bundle),
             "interface_peakvi_fit.json")
  workspace <- file.path(session_root, "peakvi_full_reference_benchmark_v1")
  session <- if (dir.exists(workspace)) {
    wb_session_resume(bundle, workspace)
  } else {
    wb_session_create(bundle, workspace_root = session_root,
                      session_id = "peakvi_full_reference_benchmark_v1")
  }
  matching <- names(session$runs)[vapply(
    session$runs, function(x) identical(x$op_name, "peakvi_fit"), logical(1)
  )]
  if (length(matching) > 1L) stop("Multiple benchmark PeakVI fits found", call. = FALSE)

  started <- Sys.time()
  if (length(matching)) {
    fit <- list(session = session, run = session$runs[[matching]])
    stopifnot(
      identical(fit$run$params$n_cells, 26723L),
      identical(fit$run$params$n_features, 97788L),
      identical(fit$run$params$n_latent, 20L),
      identical(fit$run$params$n_layers, 2L),
      identical(fit$run$params$n_hidden, 256L),
      identical(fit$run$params$max_epochs, 5L),
      isTRUE(fit$run$params$use_gpu)
    )
  } else {
    fit <- wb_session_run_op(
      session, bundle, "peakvi_fit",
      inputs = list(
        matrix = "ATAC_bdev_mm10_peaks",
        cell_ids = reference_cells,
        min_cells = 10L,
        n_latent = 20L,
        n_layers = 2L,
        n_hidden = 256L,
        dropout = 0.2,
        max_epochs = 5L,
        lr = 0.001,
        early_stopping = FALSE,
        seed = 1729L,
        use_gpu = TRUE,
        space_id_out = "peakvi.benchmark.reference"
      ),
      notes = paste(
        "Resource benchmark only: all Bdev reference cells, five epochs, GPU;",
        "no query projection and no promotion."
      )
    )
  }
  finished <- Sys.time()
  session <- fit$session

  staged <- wb_open_bundle(file.path(workspace, "stage"), validate = FALSE)
  on.exit(try(wb_close_bundle(staged), silent = TRUE), add = TRUE)
  latent <- wb_space(staged, "peakvi.benchmark.reference", cells = reference_cells)
  checks <- list(
    latent_dimensions = identical(dim(latent), c(26723L, 20L)),
    cell_order = identical(rownames(latent), reference_cells),
    values_finite = all(is.finite(as.matrix(latent))),
    gpu_requested = isTRUE(fit$run$params$use_gpu),
    five_epochs_requested = identical(fit$run$params$max_epochs, 5L),
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
    scope = "Full-reference five-epoch GPU resource benchmark; not a retained scientific model",
    promoted = FALSE,
    run_id = fit$run$run_id,
    parameters = fit$run$params,
    invocation_elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
    checks = checks,
    bundle_fingerprint_unchanged = TRUE,
    session_workspace = workspace
  ), "benchmark.json")
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
}

main()
