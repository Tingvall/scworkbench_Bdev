#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 34_preliminary_residual_ratio_bb_gene_stack.R",
  "<bdev_bundle> <window_report> <ratio_depth_report> <output_dir>",
  "<expected_bundle_fingerprint> <smoothing_fraction>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
ratio_depth_report <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_fingerprint <- args[[5L]]
smoothing_fraction <- as.numeric(args[[6L]])
stopifnot(is.finite(smoothing_fraction), smoothing_fraction > 0,
  smoothing_fraction <= 1)
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(Matrix)
  library(matrixStats)
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

trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
windows <- read_tsv(file.path(window_report, "window_distance_and_slope.tsv"))
bb <- read_tsv(file.path(window_report, "bb1_bb2_regions.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
ratio <- read_tsv(file.path(ratio_depth_report,
  "ratio_total_signal_fitted_and_residual_windows.tsv"))
derivative <- read_tsv(file.path(ratio_depth_report,
  "ratio_derivative_total_signal_fitted_and_residual_windows.tsv"))
stopifnot(nrow(trajectory) == 12642L, nrow(windows) == 2524L,
  nrow(ratio) == nrow(windows), nrow(derivative) == nrow(windows) - 2L,
  identical(windows$window_id, ratio$window_id),
  identical(bb$boundary, c("BB1", "BB2")),
  all(bb$pt_min < bb$pt_median), all(bb$pt_median < bb$pt_max),
  all(table(membership$window_id) == 30L))

cell_index <- match(membership$cell_id, trajectory$cell_id)
window_index <- match(membership$window_id, windows$window_id)
stopifnot(!anyNA(cell_index), !anyNA(window_index))
window_design <- sparseMatrix(i = cell_index, j = window_index,
  x = rep(1 / 30, nrow(membership)),
  dims = c(nrow(trajectory), nrow(windows)))
membership_split <- split(cell_index, factor(membership$window_id,
  levels = windows$window_id))
smoothing_half_width <- floor(nrow(windows) * smoothing_fraction / 2)
smoothing_rows <- unlist(lapply(seq_len(nrow(windows)), function(j)
  seq.int(max(1L, j - smoothing_half_width),
    min(nrow(windows), j + smoothing_half_width))), use.names = FALSE)
smoothing_columns <- rep(seq_len(nrow(windows)),
  vapply(seq_len(nrow(windows)), function(j)
    min(nrow(windows), j + smoothing_half_width) -
      max(1L, j - smoothing_half_width) + 1L, numeric(1)))
smoothing_counts <- tabulate(smoothing_columns, nbins = nrow(windows))
smoothing_design <- sparseMatrix(i = smoothing_rows, j = smoothing_columns,
  x = 1 / smoothing_counts[smoothing_columns],
  dims = c(nrow(windows), nrow(windows)))

genes <- c("Ebf1", "Ezh2", "Nipbl", "Wapl")
bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint),
  all(genes %in% wb_feature_ids(bundle, "SCT_data")))
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", genes,
  trajectory$cell_id, missing = "NA"))
stopifnot(identical(rownames(expression), genes),
  identical(colnames(expression), trajectory$cell_id),
  all(is.finite(expression)))
raw_mean <- as.matrix(expression %*% window_design)
raw_q25 <- raw_q75 <- matrix(NA_real_, nrow = length(genes),
  ncol = nrow(windows), dimnames = list(genes, NULL))
for (j in seq_len(nrow(windows))) {
  q <- rowQuantiles(expression[, membership_split[[j]], drop = FALSE],
    probs = c(0.25, 0.75), drop = FALSE)
  raw_q25[, j] <- q[, 1L]
  raw_q75[, j] <- q[, 2L]
}
gene_mean <- as.matrix(raw_mean %*% smoothing_design)
gene_q25 <- as.matrix(raw_q25 %*% smoothing_design)
gene_q75 <- as.matrix(raw_q75 %*% smoothing_design)
rownames(gene_mean) <- rownames(gene_q25) <- rownames(gene_q75) <- genes
gene_windows <- do.call(rbind, lapply(genes, function(gene)
  data.frame(window_id = windows$window_id, window_index = seq_len(nrow(windows)),
    pt = windows$pt_median, gene = gene, mean = gene_mean[gene, ],
    q25 = gene_q25[gene, ], q75 = gene_q75[gene, ])))
write_tsv(gene_windows, file.path(stage_dir,
  "ebf1_ezh2_nipbl_wapl_smoothed_expression_iqr.tsv"))

bb_fill <- c(BB1 = "#2C7FB8", BB2 = "#D95F0E")
gene_colors <- c(Ebf1 = "#0072B2", Ezh2 = "#CC79A7",
  Nipbl = "#009E73", Wapl = "#E69F00")
shared_x <- range(windows$pt_median)
bb_layers <- list(
  geom_rect(data = bb,
    aes(xmin = pt_min, xmax = pt_max, ymin = -Inf, ymax = Inf,
      fill = boundary), inherit.aes = FALSE, alpha = 0.15, color = NA),
  geom_vline(data = bb, aes(xintercept = pt_median, color = boundary),
    inherit.aes = FALSE, linetype = "dashed", linewidth = 0.45,
    show.legend = FALSE),
  scale_fill_manual(values = bb_fill, name = "Source bb_up region"),
  scale_color_manual(values = bb_fill, guide = "none"))
base_theme <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
    plot.margin = margin(3, 8, 3, 12))

make_gene_plot <- function(gene, show_x = FALSE) {
  d <- gene_windows[gene_windows$gene == gene, ]
  p <- ggplot(d, aes(pt, mean)) + bb_layers +
    geom_ribbon(aes(ymin = q25, ymax = q75),
      fill = gene_colors[[gene]], alpha = 0.17, color = NA) +
    geom_line(color = gene_colors[[gene]], linewidth = 0.85) +
    scale_x_continuous(limits = shared_x, expand = expansion(mult = 0)) +
    base_theme +
    labs(title = gene, x = if (show_x) "Bdev pseudotime" else NULL,
      y = "Mean SCT_data")
  if (!show_x) p <- p + theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank())
  p
}

make_residual_plot <- function(kind) {
  if (kind == "level") {
    d <- data.frame(pt = ratio$pt, value = ratio$full_residual)
    title <- "RNA/ATAC residual after total-signal adjustment"
    subtitle <- paste("Observed log2 RNA/ATAC minus the full-path fit from",
      "log-total RNA and ATAC")
    ylabel <- "Residual log2\nRNA / ATAC"
  } else {
    d <- data.frame(pt = derivative$pt, value = derivative$residual)
    title <- "RNA/ATAC ratio-change residual after total-signal adjustment"
    subtitle <- paste("Observed ratio derivative minus the total-signal",
      "derivative fit learned in T1-T3")
    ylabel <- "Residual central\nchange / window"
  }
  ggplot(d, aes(pt, value)) + bb_layers +
    geom_hline(yintercept = 0, color = "#777777", linetype = "dotted",
      linewidth = 0.4) +
    geom_line(color = "#008837", linewidth = 0.95) +
    geom_label(data = bb, aes(pt_median, Inf, label = boundary,
        fill = boundary), inherit.aes = FALSE, color = "white", linewidth = 0,
      fontface = "bold", size = 3, vjust = 1.15, show.legend = FALSE) +
    scale_x_continuous(limits = shared_x, expand = expansion(mult = 0)) +
    base_theme +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
      axis.ticks.x = element_blank(), plot.title = element_text(size = 13),
      plot.subtitle = element_text(size = 9)) +
    labs(title = title, subtitle = subtitle, y = ylabel)
}

gene_plots <- lapply(seq_along(genes), function(i)
  make_gene_plot(genes[[i]], show_x = i == length(genes)))
assemble <- function(kind) {
    make_residual_plot(kind) / gene_plots[[1L]] / gene_plots[[2L]] /
    gene_plots[[3L]] / gene_plots[[4L]] +
    plot_layout(heights = c(1.25, 1, 1, 1, 1), guides = "collect") +
    plot_annotation(
      caption = paste(
        "Full pseudotime; exact source BB1/BB2 ranges and medians; gene ribbons",
        "show smoothed within-window IQR. Bdev pt is preserved; bundle unchanged."),
      theme = theme(plot.caption = element_text(size = 8.5))) &
    theme(legend.position = "bottom")
}
level_stack <- assemble("level")
change_stack <- assemble("change")
ggsave(file.path(stage_dir, "residual_ratio_bb_ebf1_ezh2_nipbl_wapl.png"),
  level_stack, width = 11.5, height = 13, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "residual_ratio_bb_ebf1_ezh2_nipbl_wapl.pdf"),
  level_stack, width = 11.5, height = 13, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir,
  "residual_ratio_change_bb_ebf1_ezh2_nipbl_wapl.png"),
  change_stack, width = 11.5, height = 13, dpi = 300, bg = "white")
ggsave(file.path(stage_dir,
  "residual_ratio_change_bb_ebf1_ezh2_nipbl_wapl.pdf"),
  change_stack, width = 11.5, height = 13, device = cairo_pdf, bg = "white")

write_tsv(bb, file.path(stage_dir, "bb1_bb2_regions.tsv"))
diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  ratio_source = ratio_depth_report, genes = as.list(genes),
  expression_matrix = "SCT_data", trajectory_cells = nrow(trajectory),
  windows = nrow(windows), window_size = 30L, stride = 5L,
  smoothing_fraction = smoothing_fraction,
  displayed_range = "entire pseudotime",
  highlighted_regions = "source bb_up BB1 and BB2 pt min-to-max",
  bundle_mutated = FALSE, capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Residual RNA/ATAC ratio with BB regions and genes", "",
  paste("Two vertically aligned full-pseudotime figures show either the",
    "total-signal-adjusted ratio level or its adjusted local-change curve",
    "above Ebf1, Ezh2, Nipbl, and Wapl SCT_data trajectories."),
  paste("BB1 and BB2 shading uses the exact source bb_up pseudotime ranges;",
    "dashed lines mark their medians. Gene ribbons are smoothed exact",
    "within-window cell-level IQR.")), file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed residual-ratio BB gene stacks at", output_dir, "\n")
