#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config"))
  stop("Usage: 01_plot_up_geneset_gsea_umaps.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(yaml)
  library(jsonlite)
  library(digest)
})

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "06_manifests")), c("tables", "figures", "reports", "logs",
  "manifests"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_up_geneset_gsea_umaps.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
for (entry in list(cfg$inputs$parent_inventory, cfg$inputs$fixture,
    cfg$inputs$window_enrichment, cfg$inputs$enrichment_intervals))
  stopifnot(identical(sha256(entry$path), entry$sha256))
stopifnot(identical(sha256(file.path(cfg$inputs$execution_bundle$path,
  "manifest.json")), cfg$inputs$execution_bundle$manifest_sha256),
  identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256),
  identical(as.character(packageVersion("scworkbench")),
    as.character(cfg$runtime$scworkbench_version)))

read_table <- function(path) read.delim(path, check.names = FALSE,
  quote = "\"", na.strings = "NA")
fixture <- read_table(cfg$inputs$fixture$path)
local <- read_table(cfg$inputs$window_enrichment$path)
intervals <- read_table(cfg$inputs$enrichment_intervals$path)
stopifnot(!anyDuplicated(fixture$cell_id), all(is.finite(fixture$pseudotime)),
  all(local$lineage_id == cfg$inputs$lineage),
  all(intervals$lineage_id == cfg$inputs$lineage))

is_supported <- if (is.logical(local[[cfg$selection$supported_field]]))
  local[[cfg$selection$supported_field]] else
    tolower(as.character(local[[cfg$selection$supported_field]])) == "true"
local$supported_flag <- is_supported
local$window_center_fraction <- (local$window_start_fraction +
  local$window_end_fraction) / 2
valid_scale <- is.finite(local$window_start_fraction) &
  is.finite(local$window_end_fraction) &
  local$window_end_fraction > local$window_start_fraction
pt_span <- stats::median((local$window_end[valid_scale] -
  local$window_start[valid_scale]) /
  (local$window_end_fraction[valid_scale] -
    local$window_start_fraction[valid_scale]))
pt_origin <- stats::median(local$window_start[valid_scale] -
  local$window_start_fraction[valid_scale] * pt_span)
stopifnot(is.finite(pt_span), pt_span > 0, is.finite(pt_origin))
fixture$analysis_fraction <- (fixture$pseudotime - pt_origin) / pt_span
fixture$display_fraction <- 1 - fixture$analysis_fraction

suffix <- cfg$selection$feature_set_suffix
feature_sets <- sort(unique(local$feature_set[endsWith(local$feature_set,
  suffix)]))
stopifnot(length(feature_sets) == cfg$acceptance$expected_up_signatures)
supported_up <- local$supported_flag & local$feature_set %in% feature_sets
nes_limit <- max(abs(local$normalized_enrichment_score[supported_up]),
  na.rm = TRUE)
if (!is.finite(nes_limit) || nes_limit == 0) nes_limit <- 1
nes_limit <- max(2, ceiling(nes_limit * 2) / 2)
palette <- grDevices::colorRampPalette(
  c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(256L)
map_color <- function(value) {
  scaled <- (pmax(-nes_limit, pmin(nes_limit, value)) + nes_limit) /
    (2 * nes_limit)
  palette[round(scaled * (length(palette) - 1L)) + 1L]
}
human_label <- function(feature_set) {
  value <- sub("^parse_cytokine_rest_only_refinement_all_top100_v1[.]", "",
    feature_set)
  value <- sub("[.]up$", "", value)
  toupper(gsub("_", "-", value, fixed = TRUE))
}
slug <- function(feature_set) {
  value <- tolower(human_label(feature_set))
  value <- gsub("[^a-z0-9]+", "_", value)
  gsub("(^_+|_+$)", "", value)
}

bundle <- wb_open_bundle(cfg$inputs$execution_bundle$path, validate = FALSE,
  read_only = TRUE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
space <- wb_space(bundle, cfg$inputs$fit_space)
stopifnot(!anyDuplicated(space$cell_id), all(fixture$cell_id %in% space$cell_id))

empty_assignments <- function() data.frame(
  feature_set = character(), interval_id = character(),
  enrichment_direction = character(), cell_id = character(),
  pseudotime = numeric(), analysis_fraction = numeric(),
  display_fraction = numeric(), ranking_id = character(),
  window_index = integer(), window_center_fraction = numeric(),
  normalized_enrichment_score = numeric(), assignment_distance = numeric(),
  stringsAsFactors = FALSE)

assign_feature <- function(feature_set) {
  feature_intervals <- intervals[intervals$feature_set == feature_set,
    , drop = FALSE]
  if (!nrow(feature_intervals)) return(empty_assignments())
  rows <- list()
  out_index <- 0L
  for (ii in seq_len(nrow(feature_intervals))) {
    interval <- feature_intervals[ii, , drop = FALSE]
    ranking_ids <- strsplit(interval$constituent_ranking_ids[[1L]], ";",
      fixed = TRUE)[[1L]]
    windows <- local[local$feature_set == feature_set &
      local$ranking_id %in% ranking_ids & local$supported_flag &
      local$enrichment_direction == interval$enrichment_direction[[1L]],
      , drop = FALSE]
    stopifnot(nrow(windows) == interval$n_constituent_windows[[1L]])
    selected <- fixture$analysis_fraction >=
      interval$interval_start_fraction[[1L]] - 1e-12 &
      fixture$analysis_fraction <=
      interval$interval_end_fraction[[1L]] + 1e-12
    cells <- fixture[selected, , drop = FALSE]
    if (!nrow(cells)) next
    nearest <- vapply(cells$analysis_fraction, function(value)
      which.min(abs(windows$window_center_fraction - value)), integer(1L))
    chosen <- windows[nearest, , drop = FALSE]
    out_index <- out_index + 1L
    rows[[out_index]] <- data.frame(
      feature_set = feature_set,
      interval_id = interval$interval_id[[1L]],
      enrichment_direction = interval$enrichment_direction[[1L]],
      cell_id = cells$cell_id,
      pseudotime = cells$pseudotime,
      analysis_fraction = cells$analysis_fraction,
      display_fraction = cells$display_fraction,
      ranking_id = chosen$ranking_id,
      window_index = chosen$window_index,
      window_center_fraction = chosen$window_center_fraction,
      normalized_enrichment_score = chosen$normalized_enrichment_score,
      assignment_distance = abs(chosen$window_center_fraction -
        cells$analysis_fraction), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(empty_assignments())
  output <- do.call(rbind, rows)
  output <- output[order(output$cell_id, output$assignment_distance,
    output$interval_id), , drop = FALSE]
  output <- output[!duplicated(output$cell_id), , drop = FALSE]
  rownames(output) <- NULL
  output
}

draw_colorbar <- function() {
  par(mar = c(5.2, 0.5, 4.5, 4.2))
  plot(NA, xlim = c(0, 1), ylim = c(-nes_limit, nes_limit), axes = FALSE,
    xlab = "", ylab = "", type = "n")
  edges <- seq(-nes_limit, nes_limit, length.out = length(palette) + 1L)
  rect(0.18, edges[-length(edges)], 0.72, edges[-1L], col = palette,
    border = NA)
  axis(4, at = seq(-nes_limit, nes_limit, length.out = 5L), las = 1,
    cex.axis = 0.8)
  mtext("Local GSEA NES", side = 4, line = 2.8, cex = 0.9)
  box()
}

draw_map <- function(feature_set, assignments, label) {
  layout(matrix(c(1, 2), nrow = 1L), widths = c(7.2, 1.15))
  par(mar = c(5.2, 5.0, 4.5, 0.7))
  plot(space$x, space$y, pch = 16, cex = 0.15, col = "#D2D2D270",
    asp = 1, xlab = "ATAC UMAP 1", ylab = "ATAC UMAP 2",
    main = paste0(label, " up: supported local-GSEA regions"))
  if (nrow(assignments)) {
    matched <- match(assignments$cell_id, space$cell_id)
    stopifnot(!anyNA(matched))
    foreground <- data.frame(x = space$x[matched], y = space$y[matched],
      nes = assignments$normalized_enrichment_score)
    foreground <- foreground[order(abs(foreground$nes)), , drop = FALSE]
    points(foreground$x, foreground$y, pch = 16, cex = 0.62,
      col = map_color(foreground$nes))
    intervals_n <- length(unique(assignments$interval_id))
    subtitle <- sprintf("%d colored L4 fixture cells; %d supported interval%s",
      nrow(assignments), intervals_n, if (intervals_n == 1L) "" else "s")
  } else {
    text(mean(range(space$x)), mean(range(space$y)),
      "No supported local-GSEA region", cex = 1.15, font = 2,
      col = "#555555")
    subtitle <- "No cells colored"
  }
  mtext(subtitle, side = 3, line = 0.5, cex = 0.85)
  mtext("Gray: all Bdev cells; color: cells in supported merged L4 intervals",
    side = 1, line = 3.8, cex = 0.72)
  draw_colorbar()
}

draw_devices <- function(stem, plot_fun) {
  grDevices::png(file.path(dirs[["figures"]], paste0(stem, ".png")),
    width = 8.6 * 180, height = 7.0 * 180, res = 180, type = "cairo-png")
  plot_fun(); grDevices::dev.off()
  grDevices::cairo_pdf(file.path(dirs[["figures"]], paste0(stem, ".pdf")),
    width = 8.6, height = 7.0, onefile = TRUE)
  plot_fun(); grDevices::dev.off()
}

all_assignments <- list()
index_rows <- list()
for (ii in seq_along(feature_sets)) {
  feature_set <- feature_sets[[ii]]
  label <- human_label(feature_set)
  assignments <- assign_feature(feature_set)
  if (nrow(assignments)) all_assignments[[length(all_assignments) + 1L]] <-
    assignments
  figure_number <- sprintf("%02d", ii)
  stem <- paste0(figure_number, "_", slug(feature_set),
    "_up_gsea_significant_regions_umap")
  draw_devices(stem, function() draw_map(feature_set, assignments, label))
  feature_local <- local[local$feature_set == feature_set &
    local$supported_flag, , drop = FALSE]
  feature_intervals <- intervals[intervals$feature_set == feature_set,
    , drop = FALSE]
  index_rows[[ii]] <- data.frame(
    figure_number = figure_number, feature_set = feature_set,
    display_label = paste(label, "up"),
    supported_windows = nrow(feature_local),
    supported_intervals = nrow(feature_intervals),
    colored_cells = nrow(assignments),
    minimum_nes = if (nrow(assignments))
      min(assignments$normalized_enrichment_score) else NA_real_,
    maximum_nes = if (nrow(assignments))
      max(assignments$normalized_enrichment_score) else NA_real_,
    png = paste0(stem, ".png"), pdf = paste0(stem, ".pdf"),
    stringsAsFactors = FALSE)
}
index <- do.call(rbind, index_rows)
assignments <- if (length(all_assignments)) do.call(rbind, all_assignments) else
  empty_assignments()
write.table(index, file.path(dirs[["tables"]],
  "01_up_geneset_gsea_umap_index.tsv"), sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write.table(assignments, file.path(dirs[["tables"]],
  "02_up_geneset_gsea_cell_assignments.tsv"), sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")

with_regions <- sum(index$supported_intervals > 0L)
stopifnot(with_regions ==
    cfg$acceptance$expected_signatures_with_supported_intervals,
  nrow(index) == cfg$acceptance$expected_up_signatures,
  length(list.files(dirs[["figures"]], pattern = "[.]png$")) == nrow(index),
  length(list.files(dirs[["figures"]], pattern = "[.]pdf$")) == nrow(index))
report <- list(
  schema_version = "bdev.parse38_up_gsea_umaps_summary.v1",
  status = "complete_development_visualization",
  run_identity = cfg$run_identity,
  parent_run_identity = cfg$scope$parent_run_identity,
  up_signatures = nrow(index),
  signatures_with_supported_intervals = with_regions,
  signatures_without_supported_intervals = nrow(index) - with_regions,
  shared_nes_limit = nes_limit,
  colored_cell_assignments = nrow(assignments),
  assignment_rule = cfg$selection$nes_assignment,
  background = cfg$selection$background_cells,
  foreground = cfg$selection$foreground_cells,
  parent_run_mutated = FALSE,
  outputs_promoted = FALSE)
jsonlite::write_json(report, file.path(dirs[["reports"]],
  "01_up_geneset_gsea_umap_summary.json"), pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
receipt <- c(report, list(config_path = cfg_path, config_sha256 = sha256(cfg_path),
  parent_inventory_sha256 = cfg$inputs$parent_inventory$sha256,
  execution_bundle_manifest_sha256 =
    cfg$inputs$execution_bundle$manifest_sha256))
jsonlite::write_json(receipt, file.path(dirs[["manifests"]],
  "01_plotting_receipt.json"), pretty = TRUE, auto_unbox = TRUE,
  digits = NA, null = "null", na = "null")
cat("Created", nrow(index), "up-signature GSEA-region UMAPs in",
  dirs[["figures"]], "\n")
