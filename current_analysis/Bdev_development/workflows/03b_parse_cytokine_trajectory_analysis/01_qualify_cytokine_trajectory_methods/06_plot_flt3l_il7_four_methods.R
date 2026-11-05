#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config"))
  stop("Usage: 06_plot_flt3l_il7_four_methods.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(yaml)
})

cfg <- yaml::read_yaml(normalizePath(args[[2L]], mustWork = TRUE))
root <- cfg$storage$results_root
tables <- file.path(root, "01_tables")
figures <- file.path(root, "02_figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

read_table <- function(name) {
  read.delim(file.path(tables, name), check.names = FALSE, quote = "\"",
    na.strings = "NA")
}

dynamics_curve <- read_table("09_dynamics_signal_table.tsv")
dynamics_summary <- read_table("10_dynamics_dynamics_table.tsv")
interval_curve <- read_table("13_intervals_signal_curve.tsv")
interval_global <- read_table("14_intervals_global_test.tsv")
signal_intervals <- read_table("16_intervals_signal_intervals.tsv")
gsea_global <- read_table("19_gsea_global_enrichment.tsv")
gsea_local <- read_table("20_gsea_window_enrichment.tsv")
gsea_intervals <- read_table("21_gsea_enrichment_intervals.tsv")

targets <- data.frame(
  key = c(
    "parse_cytokine_rest_only_refinement_all_top100_v1.flt3l.up",
    "parse_cytokine_rest_only_refinement_all_top100_v1.il_7.up"
  ),
  label = c("FLT3L up", "IL-7 up"),
  stem = c("04_flt3l_up_four_method_summary", "05_il7_up_four_method_summary")
)

score_payload <- file.path(cfg$storage$sessions,
  "parse38_cytokine_trajectory_qualification", "stage", "operation_pack_v2",
  "transition_gwdcodtyyl9dkv1x", "outputs", "signature_scores", "payload.rds")
scores <- readRDS(score_payload)
stopifnot(inherits(scores, "wb_primitive"), identical(scores$kind, "cell_vector_family"))
score_table <- data.frame(key = scores$key, cell_id = scores$cell_id,
  value = scores$value)
score_table <- score_table[score_table$key %in% targets$key, , drop = FALSE]
stopifnot(nrow(score_table) == 2L * 1810L,
  all(table(score_table$key) == 1810L),
  all(is.finite(score_table$value)))

bundle <- wb_open_bundle(cfg$storage$execution_clone, validate = FALSE,
  read_only = TRUE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
atac_umap <- wb_space(bundle, cfg$inputs$fit_space)
stopifnot(is.data.frame(atac_umap), nrow(atac_umap) == cfg$inputs$source_bundle$expected_cells,
  all(c("cell_id", "x", "y") %in% names(atac_umap)))

pt_rows <- is.finite(gsea_local$window_start_fraction) &
  is.finite(gsea_local$window_end_fraction) &
  gsea_local$window_end_fraction > gsea_local$window_start_fraction
pt_span <- stats::median((gsea_local$window_end[pt_rows] - gsea_local$window_start[pt_rows]) /
  (gsea_local$window_end_fraction[pt_rows] - gsea_local$window_start_fraction[pt_rows]))
pt_origin <- stats::median(gsea_local$window_start[pt_rows] -
  gsea_local$window_start_fraction[pt_rows] * pt_span)
to_fraction <- function(value) pmax(0, pmin(1, (value - pt_origin) / pt_span))

umap_palette <- grDevices::colorRampPalette(
  c("#D9D9D9", "#440154", "#31688E", "#35B779", "#FDE725"))(256L)
diverging_palette <- grDevices::colorRampPalette(
  c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(256L)
elevated_color <- "#B2182B"
reduced_color <- "#2166AC"
dynamics_color <- "#6A3D9A"

score_cap <- as.numeric(stats::quantile(score_table$value, 0.99, na.rm = TRUE,
  names = FALSE))
if (!is.finite(score_cap) || score_cap <= 0) score_cap <- max(score_table$value)

sequential_colors <- function(value) {
  scaled <- pmax(0, pmin(1, value / score_cap))
  umap_palette[round(scaled * (length(umap_palette) - 1L)) + 1L]
}

diverging_colors <- function(value, limit) {
  scaled <- (pmax(-limit, pmin(limit, value)) + limit) / (2 * limit)
  diverging_palette[round(scaled * (length(diverging_palette) - 1L)) + 1L]
}

panel_label <- function(label) {
  mtext(label, side = 3, adj = 0, line = 0.45, font = 2, cex = 0.95)
}

plot_umap <- function(key) {
  target_scores <- score_table[score_table$key == key, c("cell_id", "value")]
  matched <- match(target_scores$cell_id, atac_umap$cell_id)
  stopifnot(!anyNA(matched))
  foreground <- cbind(atac_umap[matched, c("cell_id", "x", "y")],
    value = target_scores$value)
  foreground <- foreground[order(foreground$value), , drop = FALSE]
  par(mar = c(4.2, 3.8, 2.2, 1.0))
  plot(atac_umap$x, atac_umap$y, pch = 16, cex = 0.12, col = "#D6D6D680",
    asp = 1, xlab = "ATAC UMAP 1", ylab = "ATAC UMAP 2", axes = TRUE)
  points(foreground$x, foreground$y, pch = 16, cex = 0.42,
    col = sequential_colors(foreground$value))
  legend("bottomleft", inset = 0.01, bty = "n", cex = 0.72,
    legend = c("all Bdev cells", "score = 0", sprintf("score >= %.3f", score_cap)),
    pch = 16, col = c("#BDBDBD", umap_palette[[1L]], tail(umap_palette, 1L)))
  panel_label("A  ATAC-UMAP: UCell program score")
}

plot_dynamics <- function(key) {
  curve <- dynamics_curve[dynamics_curve$signal_key == key, , drop = FALSE]
  summary <- dynamics_summary[dynamics_summary$signal_key == key, , drop = FALSE]
  stopifnot(nrow(curve) > 0L, nrow(summary) == 1L)
  x <- (curve$bin - 0.5) / cfg$methods$dynamics$n_bins
  par(mar = c(4.2, 4.2, 2.2, 1.0))
  plot(x, curve$value, type = "n", xlim = c(0, 1), ylim = range(c(0, curve$value)),
    xlab = "Native L4 trajectory position", ylab = "Median UCell score")
  polygon(c(x, rev(x)), c(rep(0, length(x)), rev(curve$value)),
    border = NA, col = grDevices::adjustcolor(dynamics_color, alpha.f = 0.15))
  lines(x, curve$value, col = dynamics_color, lwd = 2.5)
  points(x, curve$value, pch = 21, bg = "white", col = dynamics_color, cex = 0.75)
  abline(v = (summary$peak_bin - 0.5) / cfg$methods$dynamics$n_bins,
    col = dynamics_color, lty = 3)
  legend("topright", bty = "n", cex = 0.75,
    legend = c(summary$curve_class,
      sprintf("rho = %.2f", summary$spearman_rho),
      sprintf("amplitude = %.3f", summary$amplitude)))
  panel_label("B  Descriptive dynamics: 20 bins")
}

plot_intervals <- function(key) {
  curve <- interval_curve[interval_curve$signal_key == key, , drop = FALSE]
  global <- interval_global[interval_global$signal_key == key, , drop = FALSE]
  regions <- signal_intervals[signal_intervals$signal_key == key, , drop = FALSE]
  stopifnot(nrow(curve) > 0L, nrow(global) == 1L)
  x <- to_fraction(curve$pseudotime)
  par(mar = c(4.2, 4.2, 2.2, 1.0))
  yrange <- range(c(0, curve$value), finite = TRUE)
  plot(x, curve$value, type = "n", xlim = c(0, 1), ylim = yrange,
    xlab = "Native L4 trajectory position", ylab = "Window median UCell score")
  if (nrow(regions)) {
    for (ii in seq_len(nrow(regions))) {
      color <- if (regions$direction[[ii]] == "elevated") elevated_color else reduced_color
      rect(to_fraction(regions$interval_start[[ii]]), par("usr")[[3L]],
        to_fraction(regions$interval_end[[ii]]), par("usr")[[4L]],
        col = grDevices::adjustcolor(color, alpha.f = 0.16), border = NA)
    }
  }
  abline(h = global$baseline, col = "#777777", lty = 2)
  lines(x, curve$value, col = "#333333", lwd = 2.2)
  points(x, curve$value, pch = 16,
    col = ifelse(curve$local_supported,
      ifelse(curve$local_direction == "elevated", elevated_color, reduced_color),
      "#777777"), cex = ifelse(curve$local_supported, 0.85, 0.45))
  legend("topright", bty = "n", cex = 0.74,
    legend = c(sprintf("global q = %.3f", global$q_global),
      sprintf("%d supported intervals", nrow(regions)), "unblocked permutations"))
  panel_label("C  Signal intervals: 5% width / 2.5% step")
}

plot_gsea <- function(key) {
  local <- gsea_local[gsea_local$feature_set == key, , drop = FALSE]
  global <- gsea_global[gsea_global$feature_set == key, , drop = FALSE]
  regions <- gsea_intervals[gsea_intervals$feature_set == key, , drop = FALSE]
  stopifnot(nrow(local) > 0L, nrow(global) == 1L)
  x <- (local$window_start_fraction + local$window_end_fraction) / 2
  nes_limit <- max(2.5, ceiling(max(abs(local$normalized_enrichment_score), na.rm = TRUE) * 2) / 2)
  par(mar = c(4.2, 4.2, 2.2, 1.0))
  plot(x, local$normalized_enrichment_score, type = "n", xlim = c(0, 1),
    ylim = c(-nes_limit, nes_limit), xlab = "Native L4 trajectory position",
    ylab = "Local normalized enrichment score")
  if (nrow(regions)) {
    for (ii in seq_len(nrow(regions))) {
      color <- if (regions$enrichment_direction[[ii]] == "elevated") elevated_color else reduced_color
      rect(regions$interval_start_fraction[[ii]], -nes_limit,
        regions$interval_end_fraction[[ii]], nes_limit,
        col = grDevices::adjustcolor(color, alpha.f = 0.12), border = NA)
    }
  }
  abline(h = 0, col = "#777777", lty = 2)
  lines(x, local$normalized_enrichment_score, col = "#777777", lwd = 1.2)
  point_colors <- diverging_colors(local$normalized_enrichment_score, nes_limit)
  points(x, local$normalized_enrichment_score, pch = 21,
    bg = point_colors, col = ifelse(local$supported, "#111111", point_colors),
    lwd = ifelse(local$supported, 1.2, 0.5), cex = ifelse(local$supported, 0.9, 0.65))
  global_color <- if (global$normalized_enrichment_score >= 0) elevated_color else reduced_color
  legend("topright", bty = "n", cex = 0.72,
    legend = c(sprintf("global NES = %.2f", global$normalized_enrichment_score),
      sprintf("global q = %.3f", global$q_global_family),
      sprintf("%d supported intervals", nrow(regions))),
    text.col = c(global_color, "#333333", "#333333"))
  panel_label("D  Preranked GSEA: local and global")
}

draw_figure <- function(key, label, stem) {
  plot_all <- function() {
    par(oma = c(0.2, 0.2, 2.3, 0.2), family = "sans")
    layout(matrix(seq_len(4L), nrow = 2L, byrow = TRUE))
    plot_umap(key)
    plot_dynamics(key)
    plot_intervals(key)
    plot_gsea(key)
    mtext(paste0(label, " cytokine program across four qualified views"),
      outer = TRUE, side = 3, line = 0.55, font = 2, cex = 1.25)
  }
  grDevices::png(file.path(figures, paste0(stem, ".png")), width = 2400,
    height = 2100, res = 200, type = "cairo-png")
  plot_all()
  grDevices::dev.off()
  grDevices::cairo_pdf(file.path(figures, paste0(stem, ".pdf")),
    width = 12, height = 10.5, onefile = TRUE)
  plot_all()
  grDevices::dev.off()
}

for (ii in seq_len(nrow(targets)))
  draw_figure(targets$key[[ii]], targets$label[[ii]], targets$stem[[ii]])

cat("Created FLT3L-up and IL-7-up four-method figures in", figures, "\n")
