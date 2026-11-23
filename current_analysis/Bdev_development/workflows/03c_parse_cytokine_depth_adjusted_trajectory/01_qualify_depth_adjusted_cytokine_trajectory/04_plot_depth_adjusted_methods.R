#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config"))
  stop("Usage: 04_plot_depth_adjusted_methods.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(slingshot)
  library(yaml)
})

cfg <- yaml::read_yaml(normalizePath(args[[2L]], mustWork = TRUE))
root <- cfg$storage$results_root
tables <- file.path(root, "01_tables")
figures <- file.path(root, "02_figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
read_table <- function(name) read.delim(file.path(tables, name),
  check.names = FALSE, quote = "\"", na.strings = "NA")
draw_devices <- function(stem, width, height, plot_fun) {
  grDevices::png(file.path(figures, paste0(stem, ".png")),
    width = width * 180, height = height * 180, res = 180, type = "cairo-png")
  plot_fun(); grDevices::dev.off()
  grDevices::cairo_pdf(file.path(figures, paste0(stem, ".pdf")),
    width = width, height = height, onefile = TRUE)
  plot_fun(); grDevices::dev.off()
}

scores <- read_table("10_depth_adjusted_scores.tsv")
dynamics_curve <- read_table("13_dynamics_signal_table.tsv")
dynamics_summary <- read_table("14_dynamics_dynamics_table.tsv")
interval_curve <- read_table("17_intervals_signal_curve.tsv")
interval_global <- read_table("18_intervals_global_test.tsv")
signal_intervals <- read_table("20_intervals_signal_intervals.tsv")
gsea_global <- read_table("23_gsea_global_enrichment.tsv")
gsea_local <- read_table("24_gsea_window_enrichment.tsv")
gsea_intervals <- read_table("25_gsea_enrichment_intervals.tsv")

stopifnot(length(unique(scores$signal_key)) == 38L,
  length(unique(dynamics_summary$signal_key)) == 38L,
  all(dynamics_curve$lineage_id == cfg$fit$lineage),
  all(interval_global$permutation_mode == "unblocked"),
  all(gsea_local$lineage_id == cfg$fit$lineage))

feature_key <- function(signal_key) sub(paste0("^", cfg$fit$lineage, "::"),
  "", signal_key)
short_label <- function(x) {
  x <- feature_key(x)
  x <- sub("^parse_cytokine_rest_only_refinement_all_top100_v1[.]", "", x)
  direction <- sub("^.*[.]", "", x)
  cytokine <- sub("[.](up|down)$", "", x)
  cytokine <- toupper(gsub("_", "-", cytokine, fixed = TRUE))
  paste0(cytokine, " [", direction, "]")
}
development_axis <- function() axis(1, at = c(0, 0.25, 0.5, 0.75, 1),
  labels = c("LSK\n0", "0.25", "0.50", "0.75", "1\nCD19"))

diverging_palette <- grDevices::colorRampPalette(
  c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(256L)
map_diverging <- function(value, limit) {
  scaled <- (pmax(-limit, pmin(limit, value)) + limit) / (2 * limit)
  diverging_palette[round(scaled * (length(diverging_palette) - 1L)) + 1L]
}

keys <- dynamics_summary$signal_key[order(dynamics_summary$rank)]
plot_keys <- rev(keys)
labels <- stats::setNames(short_label(keys), keys)
y_axis <- seq_along(plot_keys)
bins <- seq_len(cfg$methods$dynamics$n_bins)
dynamics_matrix <- matrix(NA_real_, nrow = length(keys), ncol = length(bins),
  dimnames = list(keys, bins))
for (ii in seq_len(nrow(dynamics_curve))) {
  row <- match(dynamics_curve$signal_key[[ii]], keys)
  col <- match(dynamics_curve$bin[[ii]], bins)
  if (!is.na(row) && !is.na(col))
    dynamics_matrix[row, col] <- dynamics_curve$value[[ii]]
}
# Native L4 is reversed only for biological display.
display_matrix <- dynamics_matrix[, rev(bins), drop = FALSE]
heat_limit <- max(2, ceiling(max(abs(display_matrix), na.rm = TRUE) * 2) / 2)

draw_dynamics <- function() {
  par(mar = c(5.0, 11.5, 3.2, 2.0))
  image((bins - 0.5) / length(bins), y_axis,
    t(display_matrix[plot_keys, , drop = FALSE]), col = diverging_palette,
    zlim = c(-heat_limit, heat_limit), axes = FALSE, xlab = "", ylab = "",
    useRaster = TRUE)
  development_axis()
  axis(2, at = y_axis, labels = labels[plot_keys], las = 1, tick = FALSE,
    cex.axis = 0.58)
  abline(h = seq(0.5, length(plot_keys) + 0.5, by = 1),
    col = "#FFFFFF66", lwd = 0.35)
  mtext("B-cell development (LSK -> CD19)", side = 1, line = 3.2)
  title("Depth-adjusted Parse38 score dynamics",
    sub = "RNA log1p-CP10K; per-gene nCount_RNA adjustment with ns(L4, df=8); 20 median bins")
  box()
  legend("topright", inset = 0.01, bty = "n", cex = 0.72,
    legend = c(sprintf("blue <= -%.1f", heat_limit), "white = 0",
      sprintf("red >= +%.1f", heat_limit)),
    text.col = c("#2166AC", "#666666", "#B2182B"))
}
draw_devices("01_depth_adjusted_dynamics_heatmap", 11, 12, draw_dynamics)

valid_window <- is.finite(gsea_local$window_start_fraction) &
  is.finite(gsea_local$window_end_fraction) &
  gsea_local$window_end_fraction > gsea_local$window_start_fraction
pt_span <- stats::median((gsea_local$window_end[valid_window] -
  gsea_local$window_start[valid_window]) /
  (gsea_local$window_end_fraction[valid_window] -
    gsea_local$window_start_fraction[valid_window]))
pt_origin <- stats::median(gsea_local$window_start[valid_window] -
  gsea_local$window_start_fraction[valid_window] * pt_span)
native_fraction <- function(value) pmax(0, pmin(1, (value - pt_origin) / pt_span))
display_fraction <- function(value) 1 - native_fraction(value)
reverse_interval <- function(start, end)
  c(1 - native_fraction(end), 1 - native_fraction(start))

draw_intervals <- function() {
  par(mar = c(5.0, 11.5, 3.2, 2.0))
  plot(NA, xlim = c(-0.065, 1), ylim = c(0.5, length(plot_keys) + 0.5),
    axes = FALSE, xlab = "", ylab = "")
  abline(h = y_axis, col = "#ECECEC", lwd = 0.7)
  development_axis()
  axis(2, at = y_axis, labels = labels[plot_keys], las = 1, tick = FALSE,
    cex.axis = 0.58)
  mtext("B-cell development (LSK -> CD19)", side = 1, line = 3.2)
  effect_limit <- max(abs(signal_intervals$effect_size), na.rm = TRUE)
  if (!is.finite(effect_limit) || effect_limit == 0) effect_limit <- 1
  for (ii in seq_len(nrow(signal_intervals))) {
    yy <- match(signal_intervals$signal_key[[ii]], plot_keys)
    color <- if (signal_intervals$direction[[ii]] == "elevated")
      "#B2182B" else "#2166AC"
    limits <- reverse_interval(signal_intervals$interval_start[[ii]],
      signal_intervals$interval_end[[ii]])
    width <- 1.5 + 3.5 * abs(signal_intervals$effect_size[[ii]]) / effect_limit
    segments(limits[[1L]], yy, limits[[2L]], yy, col = color, lwd = width,
      lend = "butt")
    points(display_fraction(signal_intervals$peak_pseudotime[[ii]]), yy,
      pch = 21, bg = "white", col = color, cex = 0.55, lwd = 0.9)
  }
  globally_supported <- interval_global$signal_key[
    interval_global$global_significant]
  points(rep(-0.035, length(globally_supported)),
    match(globally_supported, plot_keys), pch = 19, cex = 0.6,
    col = "#333333")
  no_local <- setdiff(plot_keys, signal_intervals$signal_key)
  if (length(no_local)) text(rep(0.5, length(no_local)),
    match(no_local, plot_keys), "no supported local interval",
    col = "#888888", cex = 0.47)
  legend("bottomright", inset = 0.01, bty = "n", cex = 0.76,
    legend = c("Elevated adjusted score", "Reduced adjusted score",
      "Global q <= 0.05"), col = c("#B2182B", "#2166AC", "#333333"),
    lwd = c(4, 4, NA), pch = c(NA, NA, 19))
  title("Permutation-supported intervals from the same adjusted scores",
    sub = "5% windows at 2.5% steps; no second depth adjustment")
  box()
}
draw_devices("02_depth_adjusted_signal_interval_map", 11, 12, draw_intervals)

feature_sets <- feature_key(keys)
window_centers <- sort(unique((gsea_local$window_start_fraction +
  gsea_local$window_end_fraction) / 2))
nes_matrix <- matrix(NA_real_, nrow = length(keys), ncol = length(window_centers),
  dimnames = list(keys, sprintf("%.6f", window_centers)))
supported_matrix <- matrix(FALSE, nrow = nrow(nes_matrix),
  ncol = ncol(nes_matrix), dimnames = dimnames(nes_matrix))
for (ii in seq_len(nrow(gsea_local))) {
  signal_key <- paste(cfg$fit$lineage, gsea_local$feature_set[[ii]], sep = "::")
  row <- match(signal_key, keys)
  center <- (gsea_local$window_start_fraction[[ii]] +
    gsea_local$window_end_fraction[[ii]]) / 2
  col <- match(center, window_centers)
  if (!is.na(row) && !is.na(col)) {
    nes_matrix[row, col] <- gsea_local$normalized_enrichment_score[[ii]]
    supported_matrix[row, col] <- isTRUE(gsea_local$supported[[ii]])
  }
}
# Reverse native window order for LSK-to-CD19 display.
nes_matrix <- nes_matrix[, rev(seq_len(ncol(nes_matrix))), drop = FALSE]
supported_matrix <- supported_matrix[, rev(seq_len(ncol(supported_matrix))), drop = FALSE]
display_centers <- sort(1 - window_centers)
nes_limit <- max(2, ceiling(max(abs(nes_matrix), na.rm = TRUE) * 2) / 2)

draw_gsea <- function() {
  par(mar = c(5.0, 11.5, 3.2, 2.0))
  image(display_centers, y_axis, t(nes_matrix[plot_keys, , drop = FALSE]),
    col = diverging_palette, zlim = c(-nes_limit, nes_limit), axes = FALSE,
    xlab = "", ylab = "", useRaster = TRUE)
  development_axis()
  axis(2, at = y_axis, labels = labels[plot_keys], las = 1, tick = FALSE,
    cex.axis = 0.58)
  supported <- which(supported_matrix[plot_keys, , drop = FALSE], arr.ind = TRUE)
  if (nrow(supported)) points(display_centers[supported[, "col"]],
    supported[, "row"], pch = 21, bg = NA, col = "#111111", cex = 0.42)
  mtext("B-cell development (LSK -> CD19)", side = 1, line = 3.2)
  title("SCT-based local preranked GSEA",
    sub = "5% windows at 2.5% steps; depth adjustment off; circles mark supported windows")
  box()
  legend("topright", inset = 0.01, bty = "n", cex = 0.72,
    legend = c(sprintf("NES <= -%.1f", nes_limit), "NES = 0",
      sprintf("NES >= +%.1f", nes_limit)),
    text.col = c("#2166AC", "#666666", "#B2182B"))
}
draw_devices("03_sct_preranked_gsea_heatmap", 11, 12, draw_gsea)

targets <- data.frame(
  feature_set = c(
    "parse_cytokine_rest_only_refinement_all_top100_v1.flt3l.up",
    "parse_cytokine_rest_only_refinement_all_top100_v1.il_7.up"),
  label = c("FLT3L up", "IL-7 up"),
  stem = c("04_flt3l_up_depth_adjusted_four_view_summary",
    "05_il7_up_depth_adjusted_four_view_summary"))
targets$signal_key <- paste(cfg$fit$lineage, targets$feature_set, sep = "::")

bundle <- wb_open_bundle(cfg$storage$execution_clone, validate = FALSE,
  read_only = TRUE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
atac_umap <- wb_space(bundle, cfg$inputs$fit_space)
fit <- readRDS(file.path(root, "05_objects", "01_authenticated_l4_model.rds"))
curves <- slingshot::slingCurves(fit)
stopifnot("Lineage4" %in% names(curves))
native_curve <- curves[["Lineage4"]]
curve_xy <- native_curve$s[native_curve$ord, , drop = FALSE]
display_curve <- curve_xy[nrow(curve_xy):1L, , drop = FALSE]

panel_label <- function(label) mtext(label, side = 3, adj = 0, line = 0.45,
  font = 2, cex = 0.92)
display_curve_class <- function(value) {
  replacements <- c(turns_on = "turns_off", turns_off = "turns_on",
    increasing = "decreasing", decreasing = "increasing",
    early_peak = "late_peak", late_peak = "early_peak")
  if (value %in% names(replacements)) unname(replacements[[value]]) else value
}

plot_umap <- function(signal_key) {
  target <- scores[scores$signal_key == signal_key, , drop = FALSE]
  matched <- match(target$cell_id, atac_umap$cell_id)
  stopifnot(!anyNA(matched))
  foreground <- data.frame(x = atac_umap$x[matched], y = atac_umap$y[matched],
    value = target$value)
  limit <- as.numeric(stats::quantile(abs(foreground$value), 0.99,
    na.rm = TRUE, names = FALSE))
  if (!is.finite(limit) || limit == 0) limit <- max(abs(foreground$value))
  foreground <- foreground[order(abs(foreground$value)), , drop = FALSE]
  par(mar = c(4.2, 3.8, 2.2, 1.0))
  plot(atac_umap$x, atac_umap$y, pch = 16, cex = 0.12, col = "#D6D6D680",
    asp = 1, xlab = "ATAC UMAP 1", ylab = "ATAC UMAP 2")
  points(foreground$x, foreground$y, pch = 16, cex = 0.42,
    col = map_diverging(foreground$value, limit))
  lines(display_curve[, 1L], display_curve[, 2L], col = "#FFFFFFD9", lwd = 4.8)
  lines(display_curve[, 1L], display_curve[, 2L], col = "#101010", lwd = 2.0)
  n_curve <- nrow(display_curve)
  arrows(display_curve[max(1L, n_curve - 12L), 1L],
    display_curve[max(1L, n_curve - 12L), 2L], display_curve[n_curve, 1L],
    display_curve[n_curve, 2L], length = 0.13, angle = 24, lwd = 2,
    col = "#101010")
  text(display_curve[1L, 1L], display_curve[1L, 2L], "LSK (0)", pos = 2,
    font = 2, cex = 0.78)
  text(display_curve[n_curve, 1L], display_curve[n_curve, 2L], "CD19 (1)",
    pos = 4, font = 2, cex = 0.78)
  legend("bottomleft", bty = "n", cex = 0.68,
    legend = c("all Bdev cells", sprintf("adjusted score <= -%.2f", limit),
      "adjusted score = 0", sprintf("adjusted score >= +%.2f", limit),
      "fitted L4: LSK -> CD19"),
    pch = c(16, 16, 16, 16, NA), lty = c(NA, NA, NA, NA, 1),
    lwd = c(NA, NA, NA, NA, 2),
    col = c("#BDBDBD", "#2166AC", "#F7F7F7", "#B2182B", "#101010"))
  panel_label("A  ATAC-UMAP: adjusted score and fitted L4")
}

plot_dynamics <- function(signal_key) {
  curve <- dynamics_curve[dynamics_curve$signal_key == signal_key, , drop = FALSE]
  summary <- dynamics_summary[dynamics_summary$signal_key == signal_key, , drop = FALSE]
  x <- 1 - (curve$bin - 0.5) / cfg$methods$dynamics$n_bins
  ord <- order(x); x <- x[ord]; curve <- curve[ord, , drop = FALSE]
  par(mar = c(4.9, 4.2, 2.2, 1.0))
  plot(x, curve$value, type = "n", xlim = c(0, 1),
    ylim = range(c(0, curve$value)), xaxt = "n",
    xlab = "B-cell development (LSK -> CD19)",
    ylab = "Median adjusted score")
  development_axis(); abline(h = 0, col = "#999999", lty = 2)
  lines(x, curve$value, col = "#6A3D9A", lwd = 2.5)
  points(x, curve$value, pch = 21, bg = "white", col = "#6A3D9A", cex = 0.75)
  abline(v = 1 - (summary$peak_bin - 0.5) / cfg$methods$dynamics$n_bins,
    col = "#6A3D9A", lty = 3)
  legend("topright", bty = "n", cex = 0.72,
    legend = c(display_curve_class(summary$curve_class),
      sprintf("LSK->CD19 rho = %.2f", -summary$spearman_rho),
      sprintf("amplitude = %.3f", summary$amplitude)))
  panel_label("B  Adjusted score dynamics: 20 bins")
}

plot_intervals <- function(signal_key) {
  curve <- interval_curve[interval_curve$signal_key == signal_key, , drop = FALSE]
  global <- interval_global[interval_global$signal_key == signal_key, , drop = FALSE]
  regions <- signal_intervals[signal_intervals$signal_key == signal_key, , drop = FALSE]
  x <- display_fraction(curve$pseudotime); ord <- order(x)
  x <- x[ord]; curve <- curve[ord, , drop = FALSE]
  par(mar = c(4.9, 4.2, 2.2, 1.0))
  plot(x, curve$value, type = "n", xlim = c(0, 1),
    ylim = range(c(0, curve$value), finite = TRUE), xaxt = "n",
    xlab = "B-cell development (LSK -> CD19)",
    ylab = "Window median adjusted score")
  development_axis()
  if (nrow(regions)) for (ii in seq_len(nrow(regions))) {
    color <- if (regions$direction[[ii]] == "elevated") "#B2182B" else "#2166AC"
    limits <- reverse_interval(regions$interval_start[[ii]], regions$interval_end[[ii]])
    rect(limits[[1L]], par("usr")[[3L]], limits[[2L]], par("usr")[[4L]],
      col = grDevices::adjustcolor(color, alpha.f = 0.16), border = NA)
  }
  abline(h = global$baseline, col = "#777777", lty = 2)
  lines(x, curve$value, col = "#333333", lwd = 2.2)
  points(x, curve$value, pch = 16,
    col = ifelse(curve$local_supported,
      ifelse(curve$local_direction == "elevated", "#B2182B", "#2166AC"),
      "#777777"), cex = ifelse(curve$local_supported, 0.85, 0.45))
  legend("topright", bty = "n", cex = 0.72,
    legend = c(sprintf("global q = %.3f", global$q_global),
      sprintf("%d supported intervals", nrow(regions)),
      "already-adjusted input"))
  panel_label("C  Adjusted score intervals: 5% / 2.5%")
}

plot_gsea <- function(feature_set) {
  local <- gsea_local[gsea_local$feature_set == feature_set, , drop = FALSE]
  global <- gsea_global[gsea_global$feature_set == feature_set, , drop = FALSE]
  regions <- gsea_intervals[gsea_intervals$feature_set == feature_set, , drop = FALSE]
  x <- 1 - (local$window_start_fraction + local$window_end_fraction) / 2
  ord <- order(x); x <- x[ord]; local <- local[ord, , drop = FALSE]
  limit <- max(2.5,
    ceiling(max(abs(local$normalized_enrichment_score), na.rm = TRUE) * 2) / 2)
  par(mar = c(4.9, 4.2, 2.2, 1.0))
  plot(x, local$normalized_enrichment_score, type = "n", xlim = c(0, 1),
    ylim = c(-limit, limit), xaxt = "n",
    xlab = "B-cell development (LSK -> CD19)", ylab = "Local GSEA NES")
  development_axis()
  if (nrow(regions)) for (ii in seq_len(nrow(regions))) {
    color <- if (regions$enrichment_direction[[ii]] == "elevated")
      "#B2182B" else "#2166AC"
    rect(1 - regions$interval_end_fraction[[ii]], -limit,
      1 - regions$interval_start_fraction[[ii]], limit,
      col = grDevices::adjustcolor(color, alpha.f = 0.12), border = NA)
  }
  abline(h = 0, col = "#777777", lty = 2)
  lines(x, local$normalized_enrichment_score, col = "#777777", lwd = 1.2)
  colors <- map_diverging(local$normalized_enrichment_score, limit)
  points(x, local$normalized_enrichment_score, pch = 21, bg = colors,
    col = ifelse(local$supported, "#111111", colors),
    lwd = ifelse(local$supported, 1.2, 0.5),
    cex = ifelse(local$supported, 0.9, 0.65))
  legend("topright", bty = "n", cex = 0.72,
    legend = c(sprintf("LSK->CD19 global NES = %.2f",
      -global$normalized_enrichment_score),
      sprintf("global q = %.3f", global$q_global_family),
      sprintf("%d supported intervals", nrow(regions)), "SCT; no added depth fit"))
  panel_label("D  SCT-based preranked GSEA")
}

draw_focus <- function(signal_key, feature_set, label, stem) {
  draw <- function() {
    par(oma = c(0.2, 0.2, 2.3, 0.2), family = "sans")
    layout(matrix(seq_len(4L), nrow = 2L, byrow = TRUE))
    plot_umap(signal_key); plot_dynamics(signal_key)
    plot_intervals(signal_key); plot_gsea(feature_set)
    mtext(paste0(label, ": depth-adjusted score views plus SCT GSEA; LSK (0) -> CD19 (1)"),
      outer = TRUE, side = 3, line = 0.55, font = 2, cex = 1.15)
  }
  draw_devices(stem, 12, 10.5, draw)
}
for (ii in seq_len(nrow(targets)))
  draw_focus(targets$signal_key[[ii]], targets$feature_set[[ii]],
    targets$label[[ii]], targets$stem[[ii]])

cat("Created 03c method and focused figures in", figures, "\n")
