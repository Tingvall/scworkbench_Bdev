#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 35_preliminary_residual_ratio_umap.R",
  "<bdev_bundle> <window_report> <ratio_depth_report> <output_dir>",
  "<expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
ratio_depth_report <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_fingerprint <- args[[5L]]
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(jsonlite)
})
read_tsv <- function(path) read.delim(path, check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = "NA")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
windows <- read_tsv(file.path(window_report, "window_distance_and_slope.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
ratio <- read_tsv(file.path(ratio_depth_report,
  "ratio_total_signal_fitted_and_residual_windows.tsv"))
stopifnot(nrow(trajectory) == 12642L, !anyDuplicated(trajectory$cell_id),
  nrow(windows) == 2524L, nrow(ratio) == nrow(windows),
  identical(windows$window_id, ratio$window_id),
  !anyDuplicated(ratio$window_id),
  all(table(membership$window_id) == 30L),
  all(membership$window_id %in% ratio$window_id),
  all(is.finite(ratio$full_residual)))

# Retain the original membership-mean projection for comparison, but use a
# continuous linear interpolation from window-center trajectory ranks for the
# UMAP color. Membership sets change only every five cells and otherwise
# quantize the map into visible stripes.
window_index <- match(membership$window_id, ratio$window_id)
stopifnot(!anyNA(window_index))
membership$residual_log2_rna_atac <- ratio$full_residual[window_index]
cell_ratio <- aggregate(residual_log2_rna_atac ~ cell_id, membership, mean)
coverage <- as.data.frame(table(membership$cell_id), stringsAsFactors = FALSE)
names(coverage) <- c("cell_id", "contributing_windows")
cell_ratio <- merge(cell_ratio, coverage, by = "cell_id", sort = FALSE)
cell_ratio <- cell_ratio[match(trajectory$cell_id, cell_ratio$cell_id), ]
stopifnot(identical(cell_ratio$cell_id, trajectory$cell_id),
  all(is.finite(cell_ratio$residual_log2_rna_atac)),
  all(cell_ratio$contributing_windows >= 1L))
window_center_rank <- (windows$start_rank + windows$end_rank) / 2
stopifnot(all(diff(window_center_rank) > 0))
cell_ratio$interpolated_residual_log2_rna_atac <- approx(
  x = window_center_rank, y = ratio$full_residual,
  xout = seq_len(nrow(trajectory)), method = "linear", rule = 2,
  ties = "ordered")$y
stopifnot(all(is.finite(cell_ratio$interpolated_residual_log2_rna_atac)))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))
cells <- as.character(wb_cells(bundle, "all"))
space <- wb_space(bundle, "umap.atac", cells = cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells),
    all(c("x", "y") %in% names(space)))
  umap <- data.frame(cell_id = cells, x = as.numeric(space$x),
    y = as.numeric(space$y), stringsAsFactors = FALSE)
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  umap <- data.frame(cell_id = cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]), stringsAsFactors = FALSE)
}
umap$residual_ratio <- cell_ratio$interpolated_residual_log2_rna_atac[
  match(umap$cell_id, cell_ratio$cell_id)]
umap$contributing_windows <- cell_ratio$contributing_windows[
  match(umap$cell_id, cell_ratio$cell_id)]
background <- umap[is.na(umap$residual_ratio), ]
colored <- umap[!is.na(umap$residual_ratio), ]
color_limit <- as.numeric(quantile(abs(colored$residual_ratio), 0.99,
  names = FALSE, type = 7))
full_color_limit <- max(abs(colored$residual_ratio))
stopifnot(nrow(colored) == nrow(trajectory), all(is.finite(colored$x)),
  all(is.finite(colored$y)), is.finite(color_limit), color_limit > 0,
  is.finite(full_color_limit), full_color_limit >= color_limit)
# Draw values nearest zero first and extremes last so strong residual regions
# remain visible rather than being hidden by overlapping points.
colored <- colored[order(abs(colored$residual_ratio)), ]

make_umap <- function(limit, clipped) {
  endpoint_labels <- if (clipped) c(
    paste0("≤ ", sprintf("%.2f", -limit)), "0",
    paste0("≥ ", sprintf("%.2f", limit))) else
    c(sprintf("%.2f", -limit), "0", sprintf("%.2f", limit))
  caption <- if (clipped) paste0(
    "Window residuals were linearly interpolated to exact cell trajectory rank; ",
    "gray cells are outside the trajectory. Color clipped symmetrically at ",
    "the absolute p99 (", sprintf("%.3f", limit), ").") else paste0(
    "Window residuals were linearly interpolated to exact cell trajectory rank; ",
    "gray cells are outside the trajectory. No clipping; symmetric limits use ",
    "the observed maximum absolute residual (", sprintf("%.3f", limit), ").")
  ggplot() +
    geom_point(data = background, aes(x, y), color = "#D8D8D8",
      size = 0.24, alpha = 0.65, stroke = 0) +
    geom_point(data = colored, aes(x, y, color = residual_ratio),
      size = 0.36, alpha = 0.95, stroke = 0) +
    scale_color_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
      midpoint = 0, limits = c(-limit, limit), oob = scales::squish,
      breaks = c(-limit, 0, limit), labels = endpoint_labels,
      name = "Residual log2\nRNA / ATAC") +
    coord_equal() +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 8.5, hjust = 0),
      legend.position = "right", legend.key.height = grid::unit(25, "mm"),
      plot.margin = margin(12, 16, 12, 16)) +
    labs(title = "Total-signal-adjusted RNA/ATAC balance on Bdev UMAP",
      subtitle = paste("Red: RNA-faster than total-signal expectation;",
        "blue: ATAC-faster; white: expected balance; rank-interpolated display"),
      caption = caption)
}
residual_umap <- make_umap(color_limit, clipped = TRUE)
residual_umap_unclipped <- make_umap(full_color_limit, clipped = FALSE)
ggsave(file.path(stage_dir, "residual_rna_atac_ratio_umap.png"),
  residual_umap, width = 8.5, height = 7.5, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "residual_rna_atac_ratio_umap.pdf"),
  residual_umap, width = 8.5, height = 7.5, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "residual_rna_atac_ratio_umap_unclipped.png"),
  residual_umap_unclipped, width = 8.5, height = 7.5, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "residual_rna_atac_ratio_umap_unclipped.pdf"),
  residual_umap_unclipped, width = 8.5, height = 7.5,
  device = cairo_pdf, bg = "white")

write_tsv(cell_ratio, file.path(stage_dir,
  "trajectory_cell_residual_rna_atac_ratio.tsv"))
diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  residual_source = ratio_depth_report,
  projection = paste("linear interpolation of full-path residual ratio from",
    "window center rank to exact cell trajectory rank"),
  comparison_projection = "mean residual ratio across containing windows",
  trajectory_cells = nrow(colored), background_cells = nrow(background),
  color_scale = "blue-white-red centered at zero",
  symmetric_absolute_p99_limit = color_limit,
  symmetric_unclipped_maximum_absolute_limit = full_color_limit,
  bundle_mutated = FALSE, capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Residual RNA/ATAC ratio on Bdev UMAP", "",
  paste("The full-path total-signal-adjusted log2 RNA/ATAC ratio was linearly",
    "interpolated from window-center ranks to every cell's exact trajectory",
    "rank, removing the five-cell quantization of the membership-mean map."),
  paste("The exported table retains both the interpolated display value and",
    "the original mean across all windows containing each cell."),
  paste("Red denotes RNA progression faster than expected from total RNA and",
    "ATAC signal; blue denotes ATAC-faster residual progression; white denotes",
    "the expected balance. Cells outside the pt trajectory are gray."),
  paste0("The symmetric color scale is clipped at absolute p99 = ",
    sprintf("%.4f", color_limit), ". A second figure is unclipped and uses",
    " the observed maximum absolute residual as both symmetric endpoints: ",
    sprintf("%.4f", full_color_limit), ".")),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed residual-ratio UMAP at", output_dir, "\n")
