#!/usr/bin/env Rscript

# Build a Parse-only expression representation and coarse cell states for
# cytokine-signature discovery. Highly variable genes are learned only from
# QC-accepted control cells. PCA and Leiden then operate on the already
# accepted control-plus-single population. Bdev and combined-stimulation cells
# are not inputs to feature learning, PCA, or clustering.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L)) {
  stop("Usage: 07_build_signature_states.R BUNDLE_DIR REPORT_DIR SESSION_ROOT [RESUME_WORKSPACE]", call. = FALSE)
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
report_dir <- normalizePath(args[[2L]], mustWork = FALSE)
session_root <- normalizePath(args[[3L]], mustWork = FALSE)
resume_workspace <- if (length(args) == 4L) normalizePath(args[[4L]], mustWork = TRUE) else NULL
if (file.exists(report_dir) || dir.exists(report_dir)) {
  stop("Report destination already exists; refusing to overwrite: ", report_dir, call. = FALSE)
}

options(width = 180, warn = 1)
suppressPackageStartupMessages({
  library(scworkbench)
  library(jsonlite)
})

imported_set <- "set:parse_cytokines"
analysis_set <- "set:parse_qc_control_single"
counts_matrix <- "RNA_counts"
normalized_matrix <- "RNA_log1p_cp10k"
pca_space <- "parse.signature_pca"
created_by <- Sys.getenv("SCWORKBENCH_CREATED_BY", unset = "jtg")

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

artifact_id_for <- function(promotion, run_id, output_name) {
  rows <- promotion$promoted[
    promotion$promoted$run_id == run_id &
      promotion$promoted$output == output_name &
      promotion$promoted$status == "ok", , drop = FALSE
  ]
  if (nrow(rows) != 1L || is.na(rows$artifact_id[[1L]])) {
    stop("Expected one promoted artifact for ", run_id, "/", output_name, call. = FALSE)
  }
  as.character(rows$artifact_id[[1L]])
}

main <- function() {
  dir.create(session_root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(session_root)) stop("Could not create session root", call. = FALSE)

  compute_profile <- wb_compute_profile("cluster")
  bundle <- wb_open_bundle(bundle_dir, read_only = FALSE, validate = is.null(resume_workspace))
  on.exit(if (!is.null(bundle)) try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)

  manifest_before <- wb_manifest(bundle)
  stopifnot(identical(manifest_before$bundle_id, "parse_cytokines"))
  cells_before <- as.character(wb_cells(bundle, cell_set = imported_set))
  analysis_cells <- as.character(wb_cells(bundle, cell_set = analysis_set))
  sets_before <- sort(as.character(wb_cell_sets(bundle)))
  matrices_before <- wb_matrices(bundle)$matrix_id
  spaces_before <- wb_spaces(bundle)$space_id
  stopifnot(
    length(cells_before) == 251599L,
    length(analysis_cells) == 154209L,
    !(normalized_matrix %in% matrices_before),
    !(pca_space %in% spaces_before)
  )

  metadata <- wb_meta(
    bundle, cells = analysis_cells,
    cols = list(wb_ref_metadata_field("core:cell_id"), "sample_id", "condition",
                "replicate", "treatment_category")
  )
  stopifnot(
    identical(as.character(metadata$cell_id), analysis_cells),
    !anyNA(metadata),
    all(metadata$treatment_category %in% c("control", "single")),
    length(unique(metadata$sample_id)) == 81L,
    length(unique(metadata$condition)) == 27L
  )
  control_cells <- analysis_cells[metadata$treatment_category == "control"]
  stopifnot(length(control_cells) > 0L, length(unique(metadata$sample_id[metadata$treatment_category == "control"])) == 3L)

  interfaces <- lapply(c("matrix_log_normalize", "hvg_select", "pca_fit", "leiden"), function(op) {
    wb_op_interface(op, bundle = bundle)
  })
  names(interfaces) <- c("matrix_log_normalize", "hvg_select", "pca_fit", "leiden")

  session_id <- paste0("parse_signature_states_", format(Sys.time(), "%Y%m%dT%H%M%S"))
  session <- if (is.null(resume_workspace)) {
    wb_session_create(bundle, workspace_root = session_root, session_id = session_id)
  } else {
    wb_session_resume(bundle, resume_workspace)
  }
  session_workspace <- session$durable$workspace

  completed_result <- function(current_session, op_name) {
    matches <- Filter(function(run) identical(run$op_name, op_name), current_session$runs)
    if (!length(matches)) return(NULL)
    if (length(matches) != 1L) stop("Ambiguous completed run for operation: ", op_name, call. = FALSE)
    list(session = current_session, run = matches[[1L]])
  }

  normalized <- completed_result(session, "matrix_log_normalize")
  if (is.null(normalized)) normalized <- wb_session_run_op(
    session, bundle, "matrix_log_normalize",
    inputs = list(
      matrix = counts_matrix,
      cell_set = "all",
      min_cells = 1L,
      output_matrix_id = normalized_matrix,
      scale_factor = 10000,
      overwrite = FALSE
    ),
    notes = paste(
      "Per-cell CP10K log normalization is materialized for the full imported matrix",
      "because scWorkbench registers full matrices. Downstream feature learning and",
      "state construction are explicitly restricted to control or control-plus-single cells."
    )
  )

  hvg <- completed_result(normalized$session, "hvg_select")
  if (is.null(hvg)) hvg <- wb_session_run_op(
    normalized$session, bundle, "hvg_select",
    inputs = list(
      matrix = wb_ref_run(normalized$run$run_id, "matrix"),
      cell_ids = control_cells,
      min_cells = 10L,
      n_features = 2000L,
      method = "vst",
      seed = 1729L,
      chunk_size = 256L
    ),
    notes = paste(
      "Select 2,000 HVGs exclusively from the three QC-accepted control samples",
      "to reduce direct cytokine-label leakage into the state feature space."
    )
  )

  pca <- completed_result(hvg$session, "pca_fit")
  if (is.null(pca)) pca <- wb_session_run_op(
    hvg$session, bundle, "pca_fit",
    inputs = list(
      matrix = wb_ref_run(normalized$run$run_id, "matrix"),
      features = wb_ref_run(hvg$run$run_id, "selected_features"),
      cell_set = analysis_set,
      min_cells = 10L,
      n_components = 30L,
      center = TRUE,
      scale = TRUE,
      backend = "irlba",
      seed = 1729L,
      space_id_out = pca_space,
      overwrite = FALSE,
      max_dense_gb = 8
    ),
    notes = "Fit PCA on the accepted control-plus-single Parse population using control-derived HVGs."
  )

  leiden <- completed_result(pca$session, "leiden")
  if (is.null(leiden)) leiden <- wb_session_run_op(
    pca$session, bundle, "leiden",
    inputs = list(
      space = wb_ref_run(pca$run$run_id, "space"),
      cell_set = analysis_set,
      min_cells = 10L,
      dims = 1:30,
      k = 30L,
      resolution = 0.4,
      n_iterations = 2L,
      seed = 1729L,
      weight = TRUE
    ),
    notes = paste(
      "Construct coarse Parse-only states without supplying treatment labels.",
      "Resolution 0.4, k=30, and seed 1729 were fixed before inspecting clusters."
    )
  )

  session_dag <- wb_session_dag(leiden$session, bundle = bundle,
                                include_artifacts = TRUE, include_active = TRUE)
  promotion <- wb_session_promote_current(
    bundle, leiden$session,
    outputs = list(
      wb_ref_run(normalized$run$run_id, "summary"),
      wb_ref_run(hvg$run$run_id, "selected_features"),
      wb_ref_run(hvg$run$run_id, "fit_summary"),
      wb_ref_run(pca$run$run_id, "loadings"),
      wb_ref_run(pca$run$run_id, "variance"),
      wb_ref_run(pca$run$run_id, "fit_summary"),
      wb_ref_run(leiden$run$run_id, "clusters")
    ),
    name_prefix = "parse_signature_states",
    description = paste(
      "Parse-only state representation for individual-cytokine signature discovery;",
      "control-derived HVGs and the accepted control-plus-single analysis population."
    ),
    created_by = created_by,
    fail_fast = TRUE
  )

  cluster_artifact_id <- artifact_id_for(promotion, leiden$run$run_id, "clusters")
  clusters <- wb_art_load(bundle, cluster_artifact_id, expected_kind = "cell_vector")
  cluster_index <- match(analysis_cells, as.character(clusters$cell_id))
  stopifnot(!anyNA(cluster_index))
  state <- as.character(clusters$value[cluster_index])
  stopifnot(!anyNA(state), length(unique(state)) >= 2L)

  state_by_sample <- as.data.frame(table(
    state = state,
    sample_id = as.character(metadata$sample_id)
  ), stringsAsFactors = FALSE)
  names(state_by_sample)[[3L]] <- "cells"
  state_by_sample <- state_by_sample[state_by_sample$cells > 0L, , drop = FALSE]
  sample_lookup <- unique(metadata[, c("sample_id", "condition", "replicate", "treatment_category")])
  state_by_sample <- merge(state_by_sample, sample_lookup, by = "sample_id", all.x = TRUE, sort = FALSE)
  state_by_sample <- state_by_sample[order(state_by_sample$state, state_by_sample$sample_id), , drop = FALSE]
  state_summary <- aggregate(state_by_sample$cells, by = list(state = state_by_sample$state),
                             FUN = function(x) c(cells = sum(x), samples = length(x),
                                                min_per_sample = min(x), median_per_sample = median(x),
                                                max_per_sample = max(x)))
  state_summary <- data.frame(
    state = state_summary$state,
    cells = state_summary$x[, "cells"],
    samples_present = state_summary$x[, "samples"],
    min_cells_per_present_sample = state_summary$x[, "min_per_sample"],
    median_cells_per_present_sample = state_summary$x[, "median_per_sample"],
    max_cells_per_present_sample = state_summary$x[, "max_per_sample"],
    stringsAsFactors = FALSE
  )
  state_condition <- as.data.frame(table(state = state, condition = as.character(metadata$condition)),
                                   stringsAsFactors = FALSE)
  names(state_condition)[[3L]] <- "cells"

  stopifnot(
    identical(as.character(wb_cells(bundle, cell_set = imported_set)), cells_before),
    identical(sort(as.character(wb_cell_sets(bundle))), sets_before),
    normalized_matrix %in% wb_matrices(bundle)$matrix_id,
    pca_space %in% wb_spaces(bundle)$space_id,
    sum(state_summary$cells) == length(analysis_cells)
  )

  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(report_dir)) stop("Could not create report directory", call. = FALSE)
  write_json_atomic(manifest_before, file.path(report_dir, "manifest_before.json"))
  write_json_atomic(interfaces, file.path(report_dir, "operation_interfaces.json"))
  write_json_atomic(session_dag, file.path(report_dir, "session_dag_before_promotion.json"))
  write_json_atomic(list(
    promoted = promotion$promoted,
    snapshot_id = if (is.null(promotion$snapshot_id)) NA_character_ else promotion$snapshot_id,
    commit_id = if (is.null(promotion$commit_id)) NA_character_ else promotion$commit_id,
    durable_session_state = if (is.null(promotion$session$durable$state)) "promoted" else promotion$session$durable$state
  ), file.path(report_dir, "promotion.json"))
  write_tsv_atomic(state_summary, file.path(report_dir, "state_summary.tsv"))
  write_tsv_atomic(state_by_sample, file.path(report_dir, "state_by_sample.tsv"))
  write_tsv_atomic(state_condition, file.path(report_dir, "state_by_condition.tsv"))
  write_json_atomic(list(
    state_method = "Leiden on 30-PC representation",
    hvg_learning_population = "QC-accepted control cells only",
    pca_and_clustering_population = analysis_set,
    combined_cells_used = FALSE,
    bdev_used = FALSE,
    hvg_count = 2000L,
    pca_components = 30L,
    leiden_k = 30L,
    leiden_resolution = 0.4,
    seed = 1729L,
    normalized_matrix = normalized_matrix,
    pca_space = pca_space,
    cluster_artifact_id = cluster_artifact_id,
    durable_session_workspace = session_workspace,
    resumed_workspace = resume_workspace,
    compute_profile = compute_profile
  ), file.path(report_dir, "state_analysis_parameters.json"))

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
    identical(as.character(wb_cells(recipient, cell_set = imported_set)), cells_before),
    identical(sort(as.character(wb_cell_sets(recipient))), sets_before)
  )
  wb_close_bundle(recipient)
  write_json_atomic(manifest_after, file.path(report_dir, "manifest.json"))
  write_json_atomic(list(
    status = "complete",
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    bundle_id = "parse_cytokines",
    imported_cells_retained = length(cells_before),
    analysis_cells = length(analysis_cells),
    control_cells_for_hvg = length(control_cells),
    states = nrow(state_summary),
    cluster_artifact_id = cluster_artifact_id,
    normalized_matrix = normalized_matrix,
    pca_space = pca_space,
    combined_cells_used = FALSE,
    bdev_used = FALSE,
    bundle_state_fingerprint = manifest_after$current_state$fingerprint,
    validation_status = validation$status
  ), file.path(report_dir, "completion.json"))
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
}

main()
