#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 22_preliminary_lsi1_sensitivity.R",
  "<bdev_bundle> <w30_s5_report> <smoothed_report> <turning_report>",
  "<output_dir> <expected_bundle_fingerprint> <loess_span>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
smoothed_report <- normalizePath(args[[3L]], mustWork = TRUE)
turning_report <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
expected_fingerprint <- args[[6L]]
loess_span <- as.numeric(args[[7L]])
stopifnot(is.finite(loess_span), loess_span > 0, loess_span <= 1)
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

read_tsv <- function(path) read.delim(path, check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = "NA")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

original <- read_tsv(file.path(window_report, "window_distance_and_slope.tsv"))
smoothed_original <- read_tsv(file.path(smoothed_report,
  "smoothed_slope_windows.tsv"))
trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
turning <- read_tsv(file.path(turning_report, "selected_turning_regions.tsv"))
stopifnot(nrow(original) == 2524L,
  identical(original$window_id, smoothed_original$window_id),
  nrow(trajectory) == 12642L, nrow(turning) == 4L,
  all(original$local_slope_rna_per_atac > 0))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))
lsi <- wb_space(bundle, "lsi", cells = trajectory$cell_id)
stopifnot(is.matrix(lsi), identical(rownames(lsi), trajectory$cell_id),
  ncol(lsi) == 50L, all(is.finite(lsi)))

window_starts <- original$start_rank
window_ends <- original$end_rank
centroid_matrix <- function(space) {
  out <- t(vapply(seq_len(nrow(original)), function(i)
    colMeans(space[window_starts[[i]]:window_ends[[i]], , drop = FALSE]),
    numeric(ncol(space))))
  rownames(out) <- original$window_id
  out
}
path_metrics <- function(centroids) {
  delta <- centroids[-1L, , drop = FALSE] -
    centroids[-nrow(centroids), , drop = FALSE]
  step <- c(0, sqrt(rowSums(delta^2)))
  cumulative <- cumsum(step)
  list(step = step, cumulative = cumulative,
    relative = cumulative / tail(cumulative, 1L),
    total = tail(cumulative, 1L))
}
centered_change <- function(x) {
  n <- length(x)
  c(x[[2L]] - x[[1L]], x[3L:n] - x[1L:(n - 2L)],
    x[[n]] - x[[n - 1L]])
}
smooth_log_slope <- function(values) {
  fit_data <- data.frame(x = original$window_index, response = log10(values))
  fit <- loess(response ~ x, data = fit_data, span = loess_span,
    degree = 2L, family = "symmetric", normalize = TRUE,
    control = loess.control(surface = "direct", statistics = "approximate"))
  predicted <- as.numeric(predict(fit,
    newdata = data.frame(x = original$window_index)))
  stopifnot(all(is.finite(predicted)))
  10^predicted
}

lsi_without_first <- lsi[, -1L, drop = FALSE]
excluded_path <- path_metrics(centroid_matrix(lsi_without_first))
atac_change_excluded <- centered_change(excluded_path$relative)
rna_change <- original$centered_rna_relative_change
excluded_slope <- rna_change / atac_change_excluded
stopifnot(all(is.finite(excluded_slope)), all(excluded_slope > 0),
  all(diff(excluded_path$relative) > 0),
  isTRUE(all.equal(original$rna_relative_distance,
    smoothed_original$rna_relative_distance, tolerance = 0)))

comparison <- data.frame(
  window_id = original$window_id,
  window_index = original$window_index,
  pt_median = original$pt_median,
  rna_relative_distance = original$rna_relative_distance,
  atac_relative_full_lsi = original$atac_relative_distance,
  atac_relative_lsi2_50 = excluded_path$relative,
  delta_atac_relative_lsi2_50_minus_full =
    excluded_path$relative - original$atac_relative_distance,
  slope_full_lsi = original$local_slope_rna_per_atac,
  slope_lsi2_50 = excluded_slope,
  smoothed_slope_full_lsi = smoothed_original$smoothed_local_slope,
  smoothed_slope_lsi2_50 = smooth_log_slope(excluded_slope),
  stringsAsFactors = FALSE)
comparison$log2_smoothed_slope_ratio_lsi2_50_vs_full <-
  log2(comparison$smoothed_slope_lsi2_50 /
    comparison$smoothed_slope_full_lsi)

method_levels <- c("Full LSI (1–50)", "LSI2–50 (LSI1 excluded)")
method_colors <- c(`Full LSI (1–50)` = "#4B1D82",
  `LSI2–50 (LSI1 excluded)` = "#D55E00")
xy_long <- rbind(
  data.frame(window_id = comparison$window_id,
    window_index = comparison$window_index, pt_median = comparison$pt_median,
    method = method_levels[[1L]],
    atac_relative_distance = comparison$atac_relative_full_lsi,
    rna_relative_distance = comparison$rna_relative_distance),
  data.frame(window_id = comparison$window_id,
    window_index = comparison$window_index, pt_median = comparison$pt_median,
    method = method_levels[[2L]],
    atac_relative_distance = comparison$atac_relative_lsi2_50,
    rna_relative_distance = comparison$rna_relative_distance))
xy_long$method <- factor(xy_long$method, levels = method_levels)
slope_long <- rbind(
  data.frame(window_id = comparison$window_id,
    pt_median = comparison$pt_median, method = method_levels[[1L]],
    raw_slope = comparison$slope_full_lsi,
    smoothed_slope = comparison$smoothed_slope_full_lsi),
  data.frame(window_id = comparison$window_id,
    pt_median = comparison$pt_median, method = method_levels[[2L]],
    raw_slope = comparison$slope_lsi2_50,
    smoothed_slope = comparison$smoothed_slope_lsi2_50))
slope_long$method <- factor(slope_long$method, levels = method_levels)

turning_xy <- merge(xy_long,
  turning[, c("window_id", "region_id", "pt_min", "pt_max")],
  by = "window_id", all = FALSE, sort = FALSE)
turning_xy$method <- factor(turning_xy$method, levels = method_levels)
turning_bands <- turning[, c("region_id", "pt_min", "pt_max")]

summary_table <- data.frame(
  metric = c("ATAC total path distance: full LSI",
    "ATAC total path distance: LSI2-50",
    "LSI2-50 / full total ATAC path",
    "median absolute shift in ATAC relative distance",
    "maximum absolute shift in ATAC relative distance",
    "Pearson correlation of smoothed log10 slopes",
    "median absolute log2 smoothed-slope ratio",
    "maximum absolute log2 smoothed-slope ratio"),
  value = c(tail(original$atac_cumulative_distance, 1L),
    excluded_path$total,
    excluded_path$total / tail(original$atac_cumulative_distance, 1L),
    median(abs(comparison$delta_atac_relative_lsi2_50_minus_full)),
    max(abs(comparison$delta_atac_relative_lsi2_50_minus_full)),
    cor(log10(comparison$smoothed_slope_full_lsi),
      log10(comparison$smoothed_slope_lsi2_50), method = "pearson"),
    median(abs(comparison$log2_smoothed_slope_ratio_lsi2_50_vs_full)),
    max(abs(comparison$log2_smoothed_slope_ratio_lsi2_50_vs_full))),
  stringsAsFactors = FALSE)
write_tsv(comparison, file.path(stage_dir, "window_lsi1_comparison.tsv"))
write_tsv(summary_table, file.path(stage_dir, "comparison_summary.tsv"))
write_tsv(turning_xy,
  file.path(stage_dir, "turning_points_lsi1_comparison.tsv"))

xy_plot <- ggplot(xy_long,
    aes(atac_relative_distance, rna_relative_distance, color = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
    linewidth = 0.75, color = "#555555") +
  geom_path(linewidth = 0.75, alpha = 0.82) +
  geom_point(data = xy_long[xy_long$window_index %% 10L == 1L, ],
    size = 0.75, alpha = 0.65) +
  geom_point(data = turning_xy, size = 3, shape = 21, fill = "white",
    stroke = 1.05) +
  geom_text(data = turning_xy,
    aes(label = region_id), fontface = "bold", size = 3,
    vjust = -0.9, show.legend = FALSE) +
  scale_color_manual(values = method_colors, name = "ATAC dimensions") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom") +
  labs(title = "A  RNA-versus-ATAC relative progression",
    subtitle = "RNA coordinates are identical; colored paths differ only in ATAC dimensions",
    x = "Relative ATAC latent path distance",
    y = "Relative RNA latent path distance")

region_colors <- c(T1 = "#0072B2", T2 = "#D55E00", T3 = "#009E73",
  T4 = "#CC79A7")
slope_plot <- ggplot(slope_long,
    aes(pt_median, smoothed_slope, color = method)) +
  geom_rect(data = turning_bands,
    aes(xmin = pt_min, xmax = pt_max, fill = region_id),
    ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.12) +
  geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.7,
    color = "#555555") +
  geom_line(aes(y = raw_slope), linewidth = 0.25, alpha = 0.08) +
  geom_line(linewidth = 1.05) +
  scale_y_log10() +
  scale_color_manual(values = method_colors, name = "ATAC dimensions") +
  scale_fill_manual(values = region_colors, guide = "none") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom") +
  labs(title = "B  Smoothed local RNA-versus-ATAC slope",
    subtitle = paste0("Robust LOESS on log10 slope, span ", loess_span,
      "; faint lines show raw slopes"),
    x = "Bdev pseudotime (window median pt)",
    y = "ΔRNA / ΔATAC (log10 scale)")

combined <- xy_plot / slope_plot +
  plot_layout(heights = c(1.1, 1), guides = "collect") +
  plot_annotation(
    title = "Sensitivity of Bdev modality progression to excluding LSI1",
    subtitle = "Same 12,642 trajectory cells and 30-cell windows with 5-cell stride",
    caption = paste0(
      "Median |ATAC relative-distance shift| = ",
      sprintf("%.4f", summary_table$value[[4L]]),
      "; max = ", sprintf("%.4f", summary_table$value[[5L]]),
      ". Smoothed log-slope correlation = ",
      sprintf("%.3f", summary_table$value[[6L]]), "."),
    theme = theme(plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10.5),
      plot.caption = element_text(size = 8.5))) &
  theme(legend.position = "bottom")

ggsave(file.path(stage_dir, "lsi1_inclusion_xy_comparison.png"), xy_plot,
  width = 8.5, height = 7.2, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "lsi1_inclusion_xy_comparison.pdf"), xy_plot,
  width = 8.5, height = 7.2, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "lsi1_inclusion_slope_comparison.png"), slope_plot,
  width = 10, height = 5.8, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "lsi1_inclusion_slope_comparison.pdf"), slope_plot,
  width = 10, height = 5.8, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "lsi1_inclusion_combined_comparison.png"), combined,
  width = 10, height = 12.5, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "lsi1_inclusion_combined_comparison.pdf"), combined,
  width = 10, height = 12.5, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id,
  bundle_fingerprint = expected_fingerprint,
  full_atac_dimensions = colnames(lsi),
  sensitivity_atac_dimensions = colnames(lsi_without_first),
  excluded_dimension = colnames(lsi)[[1L]],
  rna_coordinates = "reused unchanged from accepted 30-cell/5-cell report",
  windows = nrow(original), loess_span = loess_span,
  metrics = setNames(as.list(summary_table$value), summary_table$metric),
  turning_windows = "existing T1-T4 windows marked for paired comparison; not re-called",
  bundle_mutated = FALSE, capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Preliminary LSI1 sensitivity comparison", "",
  "The accepted full-LSI trajectory is compared with an otherwise identical",
  "calculation using LSI dimensions 2–50. RNA distances, cells, pseudotime,",
  "window membership, stride, and smoothing are unchanged.", "",
  paste0("Median absolute ATAC relative-distance shift: ",
    sprintf("%.4f", summary_table$value[[4L]]), "."),
  paste0("Maximum absolute ATAC relative-distance shift: ",
    sprintf("%.4f", summary_table$value[[5L]]), "."),
  paste0("Smoothed log10-slope Pearson correlation: ",
    sprintf("%.3f", summary_table$value[[6L]]), "."),
  "T1–T4 mark the existing windows for direct comparison and were not re-called.",
  "The bundle was opened read-only and no capsule was registered."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed preliminary LSI1 sensitivity report at", output_dir, "\n")
