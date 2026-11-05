#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config"))
  stop("Usage: 04_plot_cytokine_trajectory_methods.R --config <config.yaml>")

suppressPackageStartupMessages(library(yaml))

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
root <- cfg$storage$results_root
tables <- file.path(root, "01_tables")
figures <- file.path(root, "02_figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

read_table <- function(name) {
  read.delim(file.path(tables, name), check.names = FALSE, quote = "\"",
    na.strings = "NA")
}

signal <- read_table("09_dynamics_signal_table.tsv")
dynamics <- read_table("10_dynamics_dynamics_table.tsv")
interval_global <- read_table("14_intervals_global_test.tsv")
intervals <- read_table("16_intervals_signal_intervals.tsv")
gsea_global <- read_table("19_gsea_global_enrichment.tsv")
gsea_local <- read_table("20_gsea_window_enrichment.tsv")

stopifnot(length(unique(dynamics$signal_key)) == 38L,
  all(signal$lineage_id == "L4"),
  all(interval_global$permutation_mode == "unblocked"),
  all(intervals$permutation_mode == "unblocked"),
  all(gsea_local$lineage_id == "L4"))

short_label <- function(x) {
  x <- sub("^parse_cytokine_rest_only_refinement_all_top100_v1\\.", "", x)
  direction <- sub("^.*\\.", "", x)
  cytokine <- sub("\\.(up|down)$", "", x)
  cytokine <- toupper(gsub("_", "-", cytokine, fixed = TRUE))
  paste0(cytokine, " [", direction, "]")
}

keys <- dynamics$signal_key[order(dynamics$rank)]
plot_keys <- rev(keys)
labels <- stats::setNames(short_label(keys), keys)
y <- seq_along(plot_keys)

palette_diverging <- grDevices::colorRampPalette(
  c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(201L)
map_color <- function(value, limit) {
  clipped <- pmax(-limit, pmin(limit, value))
  index <- round((clipped + limit) / (2 * limit) * (length(palette_diverging) - 1L)) + 1L
  palette_diverging[index]
}

draw_devices <- function(stem, width, height, plot_fun) {
  png_path <- file.path(figures, paste0(stem, ".png"))
  pdf_path <- file.path(figures, paste0(stem, ".pdf"))
  grDevices::png(png_path, width = width * 180, height = height * 180,
    res = 180, type = "cairo-png")
  plot_fun()
  grDevices::dev.off()
  grDevices::cairo_pdf(pdf_path, width = width, height = height,
    onefile = TRUE)
  plot_fun()
  grDevices::dev.off()
}

draw_color_key <- function(limit, label) {
  par(mar = c(3.0, 11.5, 0.2, 2.0))
  values <- seq(-limit, limit, length.out = length(palette_diverging))
  image(values, 1, matrix(values, ncol = 1L), col = palette_diverging,
    axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  axis(1, at = c(-limit, 0, limit), labels = c(sprintf("-%.1f", limit), "0", sprintf("+%.1f", limit)))
  mtext(label, side = 1, line = 1.8, cex = 0.8)
  box()
}

# Figure 1: row-standardized UCell dynamics across the common 20-bin scale.
bins <- seq_len(cfg$methods$dynamics$n_bins)
dynamics_matrix <- matrix(NA_real_, nrow = length(keys), ncol = length(bins),
  dimnames = list(keys, bins))
for (ii in seq_len(nrow(signal))) {
  row <- match(signal$signal_key[[ii]], keys)
  col <- match(signal$bin[[ii]], bins)
  if (!is.na(row) && !is.na(col)) dynamics_matrix[row, col] <- signal$value[[ii]]
}
z_matrix <- t(apply(dynamics_matrix, 1L, function(value) {
  observed <- is.finite(value)
  if (sum(observed) < 2L || stats::sd(value[observed]) == 0) {
    value[observed] <- 0
  } else {
    value[observed] <- (value[observed] - mean(value[observed])) / stats::sd(value[observed])
  }
  value
}))
rownames(z_matrix) <- keys
constant_keys <- rownames(z_matrix)[apply(dynamics_matrix, 1L, function(value) {
  observed <- value[is.finite(value)]
  length(observed) > 0L && diff(range(observed)) == 0
})]
z_limit <- 2.5

draw_dynamics <- function() {
  layout(matrix(c(1, 2), ncol = 1L), heights = c(12, 0.8))
  par(mar = c(4.2, 11.5, 3.0, 2.0))
  image((bins - 0.5) / length(bins), y, t(z_matrix[plot_keys, , drop = FALSE]),
    col = palette_diverging, zlim = c(-z_limit, z_limit), axes = FALSE,
    xlab = "", ylab = "", useRaster = TRUE)
  abline(h = seq(0.5, length(plot_keys) + 0.5, by = 1), col = "#D8D8D8", lwd = 0.35)
  constant_y <- match(constant_keys, plot_keys)
  text(rep(0.5, length(constant_y)), constant_y, "constant zero signal",
    col = "#777777", cex = 0.48)
  axis(1, at = seq(0, 1, by = 0.2), labels = sprintf("%.1f", seq(0, 1, by = 0.2)))
  axis(2, at = y, labels = labels[plot_keys], las = 1, tick = FALSE, cex.axis = 0.58)
  mtext("Native L4 trajectory position", side = 1, line = 2.5)
  title("Parse38 cytokine dynamics along fitted Slingshot L4",
    sub = "Median UCell signal in 20 bins; signatures ordered by dynamics score",
    cex.main = 1.15, cex.sub = 0.82)
  box()
  draw_color_key(z_limit, "Within-signature standardized median UCell signal")
}
draw_devices("01_dynamics_ucell_heatmap", 11, 12.5, draw_dynamics)

# Infer the native-pseudotime range from GSEA's recorded fractional windows.
valid_window <- is.finite(gsea_local$window_start_fraction) &
  is.finite(gsea_local$window_end_fraction) &
  gsea_local$window_end_fraction > gsea_local$window_start_fraction
pt_span <- stats::median((gsea_local$window_end[valid_window] - gsea_local$window_start[valid_window]) /
  (gsea_local$window_end_fraction[valid_window] - gsea_local$window_start_fraction[valid_window]))
pt_origin <- stats::median(gsea_local$window_start[valid_window] -
  gsea_local$window_start_fraction[valid_window] * pt_span)
to_fraction <- function(value) pmax(0, pmin(1, (value - pt_origin) / pt_span))

# Figure 2: supported signal intervals, using the same signature order.
draw_intervals <- function() {
  par(mar = c(4.2, 11.5, 3.0, 2.0))
  plot(NA, xlim = c(-0.065, 1), ylim = c(0.5, length(plot_keys) + 0.5),
    axes = FALSE, xlab = "", ylab = "")
  abline(h = y, col = "#ECECEC", lwd = 0.7)
  axis(1, at = seq(0, 1, by = 0.2), labels = sprintf("%.1f", seq(0, 1, by = 0.2)))
  axis(2, at = y, labels = labels[plot_keys], las = 1, tick = FALSE, cex.axis = 0.58)
  mtext("Native L4 trajectory position", side = 1, line = 2.5)
  effect_limit <- max(abs(intervals$effect_size), na.rm = TRUE)
  for (ii in seq_len(nrow(intervals))) {
    yy <- match(intervals$signal_key[[ii]], plot_keys)
    color <- if (intervals$direction[[ii]] == "elevated") "#B2182B" else "#2166AC"
    width <- 1.5 + 3.5 * abs(intervals$effect_size[[ii]]) / effect_limit
    segments(to_fraction(intervals$interval_start[[ii]]), yy,
      to_fraction(intervals$interval_end[[ii]]), yy, col = color, lwd = width,
      lend = "butt")
    points(to_fraction(intervals$peak_pseudotime[[ii]]), yy, pch = 21,
      bg = "white", col = color, cex = 0.55, lwd = 0.9)
  }
  globally_supported <- interval_global$signal_key[interval_global$global_significant]
  points(rep(-0.035, length(globally_supported)), match(globally_supported, plot_keys),
    pch = 19, cex = 0.6, col = "#333333")
  no_local_interval <- setdiff(plot_keys, intervals$signal_key)
  text(rep(0.5, length(no_local_interval)), match(no_local_interval, plot_keys),
    "no supported local interval", col = "#888888", cex = 0.47)
  legend("bottomright", inset = 0.01,
    legend = c("Elevated interval", "Reduced interval", "Global q <= 0.05"),
    col = c("#B2182B", "#2166AC", "#333333"),
    lwd = c(4, 4, NA), pch = c(NA, NA, 19), bty = "n", cex = 0.78)
  title("Permutation-supported Parse38 signal intervals",
    sub = "Unblocked null; 5% windows at 2.5% steps; point marks interval peak",
    cex.main = 1.15, cex.sub = 0.82)
  box()
}
draw_devices("02_signal_interval_map", 11, 12, draw_intervals)

# Figure 3: local NES heatmap with aligned global NES summary.
window_centers <- sort(unique((gsea_local$window_start_fraction + gsea_local$window_end_fraction) / 2))
nes_matrix <- matrix(NA_real_, nrow = length(keys), ncol = length(window_centers),
  dimnames = list(keys, sprintf("%.6f", window_centers)))
supported_matrix <- matrix(FALSE, nrow = length(keys), ncol = length(window_centers),
  dimnames = dimnames(nes_matrix))
for (ii in seq_len(nrow(gsea_local))) {
  row <- match(gsea_local$feature_set[[ii]], keys)
  center <- (gsea_local$window_start_fraction[[ii]] + gsea_local$window_end_fraction[[ii]]) / 2
  col <- match(center, window_centers)
  if (!is.na(row) && !is.na(col)) {
    nes_matrix[row, col] <- gsea_local$normalized_enrichment_score[[ii]]
    supported_matrix[row, col] <- isTRUE(gsea_local$supported[[ii]])
  }
}
global_nes <- stats::setNames(gsea_global$normalized_enrichment_score, gsea_global$feature_set)
global_supported <- stats::setNames(gsea_global$supported, gsea_global$feature_set)
nes_limit <- max(2, ceiling(max(abs(nes_matrix), na.rm = TRUE) * 2) / 2)

draw_gsea <- function() {
  layout(matrix(c(1, 2, 3, 3), nrow = 2L, byrow = TRUE),
    widths = c(5.2, 1.4), heights = c(12, 0.8))
  par(mar = c(4.2, 11.5, 3.0, 1.0))
  image(window_centers, y, t(nes_matrix[plot_keys, , drop = FALSE]),
    col = palette_diverging, zlim = c(-nes_limit, nes_limit), axes = FALSE,
    xlab = "", ylab = "", useRaster = TRUE)
  missing_gsea <- plot_keys[apply(nes_matrix[plot_keys, , drop = FALSE], 1L,
    function(value) all(!is.finite(value)))]
  if (length(missing_gsea)) {
    missing_y <- match(missing_gsea, plot_keys)
    half_step <- stats::median(diff(window_centers)) / 2
    rect(min(window_centers) - half_step, missing_y - 0.5,
      max(window_centers) + half_step, missing_y + 0.5,
      col = "#E2E2E2", border = NA)
    text(rep(0.5, length(missing_y)), missing_y, "below minimum set size",
      col = "#666666", cex = 0.47)
  }
  abline(h = seq(0.5, length(plot_keys) + 0.5, by = 1), col = "#FFFFFF66", lwd = 0.3)
  axis(1, at = seq(0, 1, by = 0.2), labels = sprintf("%.1f", seq(0, 1, by = 0.2)))
  axis(2, at = y, labels = labels[plot_keys], las = 1, tick = FALSE, cex.axis = 0.58)
  supported <- which(supported_matrix[plot_keys, , drop = FALSE], arr.ind = TRUE)
  if (nrow(supported))
    points(window_centers[supported[, "col"]], supported[, "row"],
      pch = 21, bg = NA, col = "#111111", cex = 0.42, lwd = 0.7)
  mtext("Native L4 trajectory position", side = 1, line = 2.5)
  title("Local Parse38 preranked GSEA along fitted Slingshot L4",
    sub = "5% windows at 2.5% steps; circles mark global-family q <= 0.05",
    cex.main = 1.12, cex.sub = 0.80)
  box()

  par(mar = c(4.2, 0.7, 3.0, 1.2))
  xlimit <- max(2, ceiling(max(abs(global_nes), na.rm = TRUE) * 2) / 2)
  plot(NA, xlim = c(-xlimit, xlimit), ylim = c(0.5, length(plot_keys) + 0.5),
    axes = FALSE, xlab = "", ylab = "")
  abline(v = 0, col = "#888888", lty = 2)
  segments(0, y, global_nes[plot_keys], y,
    col = map_color(global_nes[plot_keys], xlimit), lwd = 1.4)
  points(global_nes[plot_keys], y,
    pch = ifelse(global_supported[plot_keys], 19, 1),
    col = map_color(global_nes[plot_keys], xlimit), cex = 0.75)
  axis(1, at = pretty(c(-xlimit, xlimit), n = 3))
  mtext("Global NES", side = 1, line = 2.5, cex = 0.85)
  title("Global", cex.main = 0.95)
  box()

  draw_color_key(nes_limit, "Local normalized enrichment score (NES)")
}
draw_devices("03_preranked_gsea_heatmap", 12.5, 12.5, draw_gsea)

cat("Created method figures in", figures, "\n")
