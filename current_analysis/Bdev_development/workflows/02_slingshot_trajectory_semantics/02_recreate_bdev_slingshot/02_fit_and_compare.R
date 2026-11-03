#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config"))
  stop("Usage: 02_fit_and_compare.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbTrajectory)
  library(yaml)
  library(jsonlite)
  library(digest)
})

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
results <- cfg$storage$results_root
dirs <- setNames(file.path(results, c("01_tables", "02_figures", "03_reports", "04_logs", "06_manifests")),
  c("tables", "figures", "reports", "logs", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
log_path <- file.path(dirs[["logs"]], "01_fit_and_compare.log")
sink(log_path, split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
write_tsv <- function(x, path, compress = FALSE) {
  connection <- if (compress) gzfile(path, "wt") else path
  on.exit(if (compress) close(connection), add = TRUE)
  write.table(x, connection, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
pt_binary_sha256 <- function(cell_id, pt) {
  ordering <- order(cell_id, method = "radix")
  cell_id <- as.character(cell_id[ordering]); pt <- as.numeric(pt[ordering])
  connection <- rawConnection(raw(), open = "wb")
  on.exit(close(connection))
  for (i in seq_along(cell_id)) {
    encoded <- charToRaw(enc2utf8(cell_id[[i]]))
    writeBin(as.integer(length(encoded)), connection, size = 4L, endian = "big")
    writeBin(encoded, connection)
    if (is.na(pt[[i]])) writeBin(as.raw(0L), connection) else {
      writeBin(as.raw(1L), connection)
      writeBin(pt[[i]], connection, size = 8L, endian = "big")
    }
  }
  digest::digest(rawConnectionValue(connection), algo = "sha256", serialize = FALSE)
}
numeric_summary <- function(x) c(n = sum(is.finite(x)), min = min(x, na.rm = TRUE),
  median = median(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
safe_cor <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = method))
}
draw_values <- function(xy, values, main, curves = NULL) {
  finite <- is.finite(values)
  plot(xy[, 1L], xy[, 2L], pch = 16, cex = 0.18, col = "grey85",
    xlab = "ATAC UMAP 1", ylab = "ATAC UMAP 2", main = main, asp = 1)
  if (any(finite)) {
    palette <- hcl.colors(256L, "Viridis")
    limits <- range(values[finite])
    scaled <- if (diff(limits) == 0) rep(128L, sum(finite)) else
      pmax(1L, pmin(256L, 1L + floor(255 * (values[finite] - limits[[1L]]) / diff(limits))))
    points(xy[finite, 1L], xy[finite, 2L], pch = 16, cex = 0.22,
      col = palette[scaled])
  }
  if (!is.null(curves)) for (curve in curves) {
    coordinates <- as.matrix(curve$s)
    lines(coordinates[, 1L], coordinates[, 2L], lwd = 1.4, col = "black")
  }
}

stopifnot(identical(as.character(packageVersion("scworkbench")), cfg$runtime$scworkbench$version),
  identical(as.character(packageVersion("scwbTrajectory")), cfg$runtime$scwbTrajectory$version),
  identical(as.character(packageVersion("slingshot")), cfg$runtime$slingshot_version))
pack_manifest <- system.file("scworkbench", "operation-pack.json", package = "scwbTrajectory")
stopifnot(identical(digest::digest(pack_manifest, algo = "sha256", file = TRUE),
  cfg$runtime$scwbTrajectory$operation_pack_manifest_sha256))
wb_operation_pack_load("scwbTrajectory")
interface <- wb_op_interface("scwbtrajectory::slingshot_fit")
parameters <- setNames(interface$params, vapply(interface$params, `[[`, "", "name"))
stopifnot(identical(parameters$omega$type, "bool_or_dbl"),
  identical(parameters$omega_scale$type, "dbl"))
write_json(interface, file.path(dirs[["manifests"]], "02_slingshot_fit_interface.json"))

source_manifest <- file.path(cfg$inputs$source_bundle$path, "manifest.json")
stopifnot(identical(digest::digest(source_manifest, algo = "sha256", file = TRUE),
  cfg$inputs$source_bundle$manifest_sha256))
log_event("Opening isolated execution clone and applying clone-only schema migration")
bundle <- wb_open_bundle(cfg$storage$execution_clone, validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
all_cells <- as.character(wb_cells(bundle, "all"))
stopifnot(length(all_cells) == cfg$inputs$source_bundle$expected_cells)
metadata <- wb_meta(bundle, cells = all_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$cluster_field,
  cfg$inputs$preserved_pseudotime$field))
stopifnot(identical(as.character(metadata$cell_id), all_cells))
preserved_pt <- as.numeric(metadata[[cfg$inputs$preserved_pseudotime$field]])
stopifnot(sum(is.finite(preserved_pt)) == cfg$inputs$preserved_pseudotime$defined_cells,
  sum(is.na(preserved_pt)) == cfg$inputs$preserved_pseudotime$missing_cells,
  identical(pt_binary_sha256(metadata$cell_id, preserved_pt),
    cfg$inputs$preserved_pseudotime$binary_sha256))
pt_before <- pt_binary_sha256(metadata$cell_id, preserved_pt)
clusters <- wb_prim_cell_vector(stats::setNames(
  as.character(metadata[[cfg$inputs$cluster_field]]), metadata$cell_id),
  value_type = "categorical")

session_workspace <- file.path(cfg$storage$sessions,
  "fit_published_bdev_slingshot_recipe")
if (dir.exists(session_workspace)) {
  log_event("Resuming the durable completed Slingshot fit")
  session <- wb_session_resume(bundle, session_workspace)
  stopifnot(length(session$runs) == 1L)
  fit <- list(session = session, run = session$runs[[1L]])
} else {
  session <- wb_session_create(bundle, workspace_root = cfg$storage$sessions,
    session_id = "fit_published_bdev_slingshot_recipe")
  log_event("Fitting Slingshot on", cfg$inputs$fit_space, "with omega=TRUE and omega_scale=5")
  fit <- wb_session_run_op(session, bundle, "scwbtrajectory::slingshot_fit",
    inputs = list(
      fit_space = cfg$inputs$fit_space,
      clusters = clusters,
      start.clus = cfg$fit$start_cluster,
      exclude_clusters = unlist(cfg$fit$exclude_clusters, use.names = FALSE),
      stretch = as.numeric(cfg$fit$stretch),
      omega = isTRUE(cfg$fit$omega),
      omega_scale = as.numeric(cfg$fit$omega_scale)),
    notes = "Development-only recreation of the published Bdev ATAC-UMAP Slingshot recipe; outputs remain staged and are not promoted.")
}

family <- fit$run$outputs$pseudotime
weights <- fit$run$outputs$lineage_weights
keys <- unique(as.character(family$key))
stopifnot(cfg$fit$selected_lineage %in% keys)
write_tsv(data.frame(lineage = family$key, cell_id = family$cell_id,
  pseudotime = family$value), file.path(dirs[["tables"]],
  "01_fitted_lineage_pseudotime.tsv.gz"), compress = TRUE)
write_tsv(data.frame(lineage = weights$key, cell_id = weights$cell_id,
  weight = weights$value), file.path(dirs[["tables"]],
  "02_fitted_lineage_weights.tsv.gz"), compress = TRUE)

model_ref <- fit$run$outputs$model
model_path <- file.path(fit$session$staging_bundle$dir,
  model_ref$model_tree_root, model_ref$model_primary_path)
model <- readRDS(model_path)
lineages <- slingshot::slingLineages(model)
model_params <- slingshot::slingParams(model)
lineage_paths <- do.call(rbind, lapply(seq_along(lineages), function(i)
  data.frame(lineage = paste0("L", i), order = seq_along(lineages[[i]]),
    cluster = as.character(lineages[[i]]))))
write_tsv(lineage_paths, file.path(dirs[["tables"]], "03_lineage_cluster_paths.tsv"))

selected_index <- family$key == cfg$fit$selected_lineage
selected_values <- stats::setNames(as.numeric(family$value[selected_index]),
  as.character(family$cell_id[selected_index]))
fitted_full <- unname(selected_values[metadata$cell_id])
negative_fitted <- -fitted_full
comparison <- data.frame(cell_id = metadata$cell_id,
  cluster = as.character(metadata[[cfg$inputs$cluster_field]]),
  fitted_selected_lineage = fitted_full,
  negative_fitted_selected_lineage = negative_fitted,
  preserved_pt = preserved_pt)
comparison$membership <- ifelse(is.finite(fitted_full) & is.finite(preserved_pt), "both_defined",
  ifelse(is.finite(fitted_full), "fitted_only",
    ifelse(is.finite(preserved_pt), "preserved_only", "neither_defined")))
comparison$difference_negative_fit_minus_preserved <- negative_fitted - preserved_pt
write_tsv(comparison, file.path(dirs[["tables"]],
  "04_selected_lineage_vs_preserved_pt.tsv.gz"), compress = TRUE)

lineage_summary <- do.call(rbind, lapply(keys, function(key) {
  current <- family$key == key
  values <- stats::setNames(as.numeric(family$value[current]),
    as.character(family$cell_id[current]))
  aligned <- unname(values[metadata$cell_id])
  shared <- is.finite(aligned) & is.finite(preserved_pt)
  data.frame(lineage = key, cells = sum(is.finite(aligned)), shared_with_preserved_pt = sum(shared),
    spearman_raw = safe_cor(aligned, preserved_pt, "spearman"),
    spearman_negative = safe_cor(-aligned, preserved_pt, "spearman"),
    pearson_negative = safe_cor(-aligned, preserved_pt, "pearson"))
}))
write_tsv(lineage_summary, file.path(dirs[["tables"]], "05_lineage_comparison_summary.tsv"))

cluster_levels <- sort(unique(as.character(metadata[[cfg$inputs$cluster_field]])))
cluster_summary <- do.call(rbind, lapply(cluster_levels, function(cluster) {
  index <- as.character(metadata[[cfg$inputs$cluster_field]]) == cluster
  data.frame(cluster = cluster, all_cells = sum(index),
    selected_lineage_cells = sum(index & is.finite(fitted_full)),
    preserved_pt_defined = sum(index & is.finite(preserved_pt)),
    both_defined = sum(index & is.finite(fitted_full) & is.finite(preserved_pt)))
}))
write_tsv(cluster_summary, file.path(dirs[["tables"]], "06_cluster_coverage_summary.tsv"))

space_values <- wb_space(bundle, cfg$inputs$fit_space, cells = metadata$cell_id)
if (is.data.frame(space_values)) {
  stopifnot(all(c("cell_id", "x", "y") %in% names(space_values)))
  xy <- cbind(x = as.numeric(space_values$x), y = as.numeric(space_values$y))
  rownames(xy) <- as.character(space_values$cell_id)
} else {
  xy <- as.matrix(space_values)
}
stopifnot(identical(rownames(xy), metadata$cell_id), ncol(xy) == 2L)
png(file.path(dirs[["figures"]], "01_atac_umap_recreated_L4_vs_preserved_pt.png"),
  width = 2400, height = 1200, res = 180)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
draw_values(xy, negative_fitted, "Recreated -L4", slingshot::slingCurves(model))
draw_values(xy, preserved_pt, "Preserved Bdev pt")
dev.off()

png(file.path(dirs[["figures"]], "02_atac_umap_fitted_lineage_family.png"),
  width = 1800, height = 1800, res = 180)
grid_n <- ceiling(sqrt(length(keys)))
par(mfrow = c(grid_n, grid_n), mar = c(2.5, 2.5, 2.5, 0.5))
for (key in keys) {
  current <- family$key == key
  values <- stats::setNames(as.numeric(family$value[current]),
    as.character(family$cell_id[current]))
  draw_values(xy, unname(values[metadata$cell_id]), key)
}
dev.off()

shared <- is.finite(negative_fitted) & is.finite(preserved_pt)
membership_counts <- as.list(table(factor(comparison$membership,
  levels = c("both_defined", "fitted_only", "preserved_only", "neither_defined"))))
names(membership_counts) <- c("both_defined", "fitted_only", "preserved_only", "neither_defined")
fit_summary <- list(
  schema_version = "bdev.slingshot_trajectory_recreation_result.v1",
  status = "complete_unpromoted_development_validation",
  run_identity = cfg$run_identity,
  runtime = list(scworkbench = as.character(packageVersion("scworkbench")),
    scwbTrajectory = as.character(packageVersion("scwbTrajectory")),
    slingshot = as.character(packageVersion("slingshot"))),
  recipe = cfg$fit,
  cells = list(total = length(all_cells), fitted = length(model_ref$cells),
    excluded_from_fit = length(all_cells) - length(model_ref$cells)),
  lineages = list(count = length(keys), keys = keys,
    cluster_paths = lapply(lineages, as.character)),
  selected_lineage = list(key = cfg$fit$selected_lineage,
    fitted_cells = sum(is.finite(fitted_full)), membership_counts = membership_counts,
    shared_cells = sum(shared),
    spearman_negative_fit_vs_preserved = safe_cor(negative_fitted, preserved_pt, "spearman"),
    pearson_negative_fit_vs_preserved = safe_cor(negative_fitted, preserved_pt, "pearson"),
    mean_absolute_difference = mean(abs(negative_fitted[shared] - preserved_pt[shared])),
    maximum_absolute_difference = max(abs(negative_fitted[shared] - preserved_pt[shared])),
    exact_numeric_matches = sum(negative_fitted[shared] == preserved_pt[shared])),
  preserved_pt = list(defined = sum(is.finite(preserved_pt)), missing = sum(is.na(preserved_pt)),
    binary_sha256_before = pt_before,
    binary_sha256_after = pt_binary_sha256(metadata$cell_id, preserved_pt),
    overwritten = FALSE),
  slingshot_parameters = list(omega = model_params$omega,
    omega_scale = model_params$omega_scale, stretch = model_params$stretch,
    start_cluster = cfg$fit$start_cluster),
  scworkbench = list(operation = "scwbtrajectory::slingshot_fit",
    run_id = fit$run$run_id, session_workspace = fit$session$durable$workspace,
    outputs_promoted = FALSE, staged_session_retained = TRUE),
  source_bundle_mutated = FALSE)
stopifnot(identical(fit_summary$preserved_pt$binary_sha256_before,
  fit_summary$preserved_pt$binary_sha256_after),
  identical(digest::digest(source_manifest, algo = "sha256", file = TRUE),
    cfg$inputs$source_bundle$manifest_sha256))
write_json(fit_summary, file.path(dirs[["reports"]], "01_trajectory_recreation_summary.json"))
write_json(list(schema_version = "bdev.slingshot_trajectory_session_receipt.v1",
  run_identity = cfg$run_identity, operation = "scwbtrajectory::slingshot_fit",
  run_id = fit$run$run_id, session_workspace = fit$session$durable$workspace,
  staged_outputs = names(fit$run$outputs), promoted = FALSE,
  source_bundle_mutated = FALSE, pt_overwritten = FALSE),
  file.path(dirs[["manifests"]], "03_scworkbench_session_receipt.json"))
log_event("Complete: retained staged lineage family; no outputs promoted; preserved pt unchanged")
