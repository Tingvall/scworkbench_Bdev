#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 20_preliminary_modality_ratio_umap_five_genes.R",
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

distance <- read_tsv(file.path(window_report, "window_distance_and_slope.tsv"))
slope <- read_tsv(file.path(smoothed_report, "smoothed_slope_windows.tsv"))
trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
turning <- read_tsv(file.path(turning_report, "selected_turning_regions.tsv"))
turning_cells <- read_tsv(file.path(turning_report,
  "selected_turning_region_cells.tsv"))
stopifnot(identical(distance$window_id, slope$window_id),
  nrow(slope) == 2524L, nrow(trajectory) == 12642L,
  !anyDuplicated(trajectory$cell_id), nrow(turning) == 4L,
  all(table(membership$window_id) == 30L),
  all(table(turning_cells$region_id) == 30L),
  all(slope$smoothed_local_slope > 0))

region_colors <- c(T1 = "#0072B2", T2 = "#D55E00", T3 = "#009E73",
  T4 = "#CC79A7")
turning$region_id <- factor(turning$region_id, levels = names(region_colors))
stopifnot(!anyNA(turning$region_id))
turning_bands <- turning[, c("region_id", "pt_min", "pt_max")]
turning_xy <- merge(turning[, c("window_id", "region_id", "pt")],
  distance[, c("window_id", "atac_relative_distance",
    "rna_relative_distance")], by = "window_id", sort = FALSE)
turning_xy <- turning_xy[match(names(region_colors),
  as.character(turning_xy$region_id)), ]
turning_xy$region_id <- factor(turning_xy$region_id,
  levels = names(region_colors))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))

genes <- c("Ebf1", "Ezh2", "Nipbl", "Mau2", "Wapl", "Ctcf", "Eed",
  "Tcf7", "Lef1")
available <- wb_feature_ids(bundle, "SCT_data")
stopifnot(all(genes %in% available))
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", genes,
  trajectory$cell_id, missing = "NA"))
stopifnot(identical(rownames(expression), genes),
  identical(colnames(expression), trajectory$cell_id),
  all(is.finite(expression)))

cell_index <- match(membership$cell_id, trajectory$cell_id)
window_index <- match(membership$window_id, slope$window_id)
stopifnot(!anyNA(cell_index), !anyNA(window_index))
membership_split <- split(cell_index, factor(membership$window_id,
  levels = slope$window_id))
stopifnot(length(membership_split) == nrow(slope),
  all(lengths(membership_split) == 30L))

smoothing_fallback_genes <- character()
smooth_window_values <- function(values, gene) {
  fit_data <- data.frame(x = slope$window_index, response = values)
  fit <- tryCatch(loess(response ~ x, data = fit_data, span = loess_span,
      degree = 2L, family = "symmetric", normalize = TRUE,
      control = loess.control(surface = "direct", statistics = "approximate")),
    error = function(e) {
      smoothing_fallback_genes <<- c(smoothing_fallback_genes, gene)
      loess(response ~ x, data = fit_data, span = loess_span,
        degree = 2L, family = "gaussian", normalize = TRUE,
        control = loess.control(surface = "direct", statistics = "approximate"))
    })
  predicted <- as.numeric(predict(fit,
    newdata = data.frame(x = slope$window_index)))
  stopifnot(length(predicted) == length(values), all(is.finite(predicted)))
  predicted
}

expression_rows <- lapply(seq_along(genes), function(g) {
  raw <- vapply(membership_split, function(ix) mean(expression[g, ix]),
    numeric(1))
  data.frame(window_id = slope$window_id,
    window_index = slope$window_index, pt_median = slope$pt_median,
    gene = genes[[g]], mean_sct_data = raw,
    smoothed_mean_sct_data = smooth_window_values(raw, genes[[g]]),
    stringsAsFactors = FALSE)
})
expression_windows <- do.call(rbind, expression_rows)
expression_windows$gene <- factor(expression_windows$gene, levels = genes)
stopifnot(nrow(expression_windows) == length(genes) * nrow(slope),
  all(is.finite(expression_windows$smoothed_mean_sct_data)))

# A cell participates in up to six overlapping 30-cell windows. Assign its
# modality-balance value as the mean smoothed log2 slope over those memberships.
# Thus zero means matched RNA/ATAC progress, positive means RNA-faster, and
# negative means ATAC-faster.
membership$smoothed_log2_rna_atac_ratio <-
  log2(slope$smoothed_local_slope[window_index])
cell_ratio <- aggregate(smoothed_log2_rna_atac_ratio ~ cell_id, membership,
  mean)
names(cell_ratio)[[2L]] <- "mean_smoothed_log2_rna_atac_ratio"
coverage <- as.data.frame(table(membership$cell_id), stringsAsFactors = FALSE)
names(coverage) <- c("cell_id", "contributing_windows")
cell_ratio <- merge(cell_ratio, coverage, by = "cell_id", sort = FALSE)
cell_ratio <- cell_ratio[match(trajectory$cell_id, cell_ratio$cell_id), ]
stopifnot(identical(cell_ratio$cell_id, trajectory$cell_id),
  all(is.finite(cell_ratio$mean_smoothed_log2_rna_atac_ratio)))

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
umap$log2_ratio <- cell_ratio$mean_smoothed_log2_rna_atac_ratio[
  match(umap$cell_id, cell_ratio$cell_id)]
umap_background <- umap[is.na(umap$log2_ratio), ]
umap_trajectory <- umap[!is.na(umap$log2_ratio), ]
turning_highlight <- merge(umap, turning_cells[, c("cell_id", "region_id")],
  by = "cell_id", all = FALSE, sort = FALSE)
turning_highlight$region_id <- factor(turning_highlight$region_id,
  levels = names(region_colors))
centroids <- aggregate(cbind(x, y) ~ region_id, turning_highlight, median)
color_limit <- max(abs(umap_trajectory$log2_ratio))
stopifnot(nrow(umap_trajectory) == nrow(trajectory),
  nrow(turning_highlight) == 120L, is.finite(color_limit), color_limit > 0)

shared_x <- range(slope$pt_median)
base_pt_theme <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
    plot.margin = margin(2, 7, 2, 3))

slope_plot <- ggplot(slope, aes(pt_median, smoothed_local_slope)) +
  geom_rect(data = turning_bands,
    aes(xmin = pt_min, xmax = pt_max, fill = region_id),
    ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.20) +
  geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.6,
    color = "#4A4A4A") +
  geom_line(color = "#4B1D82", linewidth = 0.95) +
  geom_point(data = turning, aes(pt, smoothed_slope, color = region_id),
    inherit.aes = FALSE, size = 2) +
  geom_text(data = turning,
    aes(pt, smoothed_slope, label = region_id, color = region_id),
    inherit.aes = FALSE, vjust = -0.8, size = 2.8, fontface = "bold",
    show.legend = FALSE) +
  scale_x_continuous(limits = shared_x) + scale_y_log10() +
  scale_fill_manual(values = region_colors, guide = "none") +
  scale_color_manual(values = region_colors, guide = "none") +
  base_pt_theme +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  labs(title = "A  RNA/ATAC local slope", y = "ΔRNA / ΔATAC\n(log10)")

gene_colors <- c(Ebf1 = "#0072B2", Ezh2 = "#CC79A7", Nipbl = "#E69F00",
  Mau2 = "#009E73", Wapl = "#D55E00", Ctcf = "#6A3D9A",
  Eed = "#8C564B", Tcf7 = "#17BECF", Lef1 = "#7F7F7F")
panel_letters <- setNames(LETTERS[2:(length(genes) + 1L)], genes)
make_gene_plot <- function(gene, show_x = FALSE) {
  gene_data <- expression_windows[expression_windows$gene == gene, ]
  p <- ggplot(gene_data, aes(pt_median, smoothed_mean_sct_data)) +
    geom_rect(data = turning_bands,
      aes(xmin = pt_min, xmax = pt_max, fill = region_id),
      ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.20) +
    geom_line(color = gene_colors[[gene]], linewidth = 0.9) +
    scale_x_continuous(limits = shared_x) +
    scale_fill_manual(values = region_colors, guide = "none") +
    base_pt_theme +
    labs(title = paste0(panel_letters[[gene]], "  ", gene),
      x = if (show_x) "Bdev pseudotime (window median pt)" else NULL,
      y = "Mean SCT_data")
  if (!show_x) p <- p + theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank())
  p
}
gene_plots <- lapply(seq_along(genes), function(i)
  make_gene_plot(genes[[i]], show_x = i == length(genes)))

xy_plot <- ggplot(distance,
    aes(atac_relative_distance, rna_relative_distance)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
    linewidth = 0.7, color = "#555555") +
  geom_path(color = "#777777", linewidth = 0.65, alpha = 0.8) +
  geom_point(data = turning_xy, aes(color = region_id), size = 3) +
  geom_text(data = turning_xy, aes(label = region_id, color = region_id),
    vjust = -1, fontface = "bold", size = 3.2, show.legend = FALSE) +
  scale_color_manual(values = region_colors, guide = "none") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11)) +
  labs(title = paste0(LETTERS[[length(genes) + 2L]],
      "  Relative progression trajectory"),
    subtitle = "Dotted line: equal RNA and ATAC progress",
    x = "ATAC relative distance", y = "RNA relative distance")

umap_plot <- ggplot() +
  geom_point(data = umap_background, aes(x, y), color = "#D8D8D8",
    size = 0.20, alpha = 0.6, stroke = 0) +
  geom_point(data = umap_trajectory, aes(x, y, color = log2_ratio),
    size = 0.30, alpha = 0.9, stroke = 0) +
  geom_point(data = turning_highlight, aes(x, y), shape = 21, fill = NA,
    color = "#1A1A1A", size = 0.9, stroke = 0.25) +
  geom_label(data = centroids, aes(x, y, label = region_id), fill = "white",
    color = "#1A1A1A", linewidth = 0.18, fontface = "bold", size = 3) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-color_limit, color_limit),
    oob = scales::squish, name = expression(log[2](Delta*RNA/Delta*ATAC))) +
  coord_equal() + theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9), legend.position = "right") +
  labs(title = paste0(LETTERS[[length(genes) + 3L]],
      "  Local modality balance on Bdev umap.atac"),
    subtitle = paste(
      "Red: RNA-faster; blue: ATAC-faster; white: matched; gray: outside trajectory"))

left_column <- wrap_plots(c(list(slope_plot), gene_plots), ncol = 1,
  heights = c(1.08, rep(1, length(genes))))
right_column <- xy_plot / umap_plot + plot_layout(heights = c(0.8, 1.2))
combined <- (left_column | right_column) +
  plot_layout(widths = c(1.18, 1), guides = "collect") +
  plot_annotation(
    title = "Bdev RNA–ATAC progression, turning windows, and regulator expression",
    subtitle = paste(
      "Slope and nine expression curves share pseudotime;",
      "UMAP color is the membership-weighted smoothed local modality ratio"),
    caption = paste(
      "UMAP zero denotes matched RNA/ATAC progress (slope = 1).",
      "Black outlines and T1–T4 labels mark the exact turning-window cells.",
      "Preliminary and read-only."),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 8.2)))

ggsave(file.path(stage_dir,
  "modality_ratio_umap_nine_genes.png"), combined,
  width = 14.5, height = 22, dpi = 320, bg = "white")
ggsave(file.path(stage_dir,
  "modality_ratio_umap_nine_genes.pdf"), combined,
  width = 14.5, height = 22, device = cairo_pdf, bg = "white")

write_tsv(expression_windows,
  file.path(stage_dir, "nine_gene_window_expression.tsv"))
write_tsv(cell_ratio,
  file.path(stage_dir, "trajectory_cell_modality_log2_ratio.tsv"))
write_tsv(turning_xy,
  file.path(stage_dir, "selected_turning_points_rna_atac.tsv"))
diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id,
  bundle_fingerprint = expected_fingerprint,
  genes = genes, expression_matrix = "SCT_data",
  loess_span = loess_span,
  robust_loess_fallback_genes = unique(smoothing_fallback_genes),
  fallback_method = "Gaussian LOESS only if robust LOESS fails numerically",
  trajectory_cells = nrow(trajectory),
  cell_modality_value = paste(
    "mean across overlapping memberships of log2(smoothed local RNA/ATAC slope)"),
  zero_definition = "matched local RNA and ATAC progress; slope = 1",
  umap_color_limit = color_limit, nontrajectory_color = "gray",
  turning_regions = nrow(turning), bundle_mutated = FALSE,
  capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Preliminary modality-ratio UMAP and nine-gene summary", "",
  "The UMAP color is log2(delta RNA / delta ATAC), averaged over every",
  "smoothed 30-cell window containing each trajectory cell. Red is RNA-faster,",
  "blue is ATAC-faster, and white is matched progression. Non-trajectory cells",
  "are gray; T1–T4 cells are outlined and labeled.", "",
  paste(paste(genes, collapse = ", "), "are displayed as LOESS-smoothed"),
  paste0("Gaussian fallback after robust-LOESS numerical failure: ",
    if (length(smoothing_fallback_genes))
      paste(unique(smoothing_fallback_genes), collapse = ", ") else "none", "."),
  "30-cell window means from SCT_data. The bundle was opened read-only."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed preliminary modality-ratio report at", output_dir, "\n")
