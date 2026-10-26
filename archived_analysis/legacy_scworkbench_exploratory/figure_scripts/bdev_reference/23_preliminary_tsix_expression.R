#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 23_preliminary_tsix_expression.R",
  "<bdev_bundle> <w30_s5_report> <turning_report> <output_dir>",
  "<expected_bundle_fingerprint> <loess_span>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
turning_report <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_fingerprint <- args[[5L]]
loess_span <- as.numeric(args[[6L]])
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

trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
windows <- read_tsv(file.path(window_report, "window_distance_and_slope.tsv"))
bb <- read_tsv(file.path(window_report, "bb1_bb2_regions.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
turning <- read_tsv(file.path(turning_report, "selected_turning_regions.tsv"))
turning_cells <- read_tsv(file.path(turning_report,
  "selected_turning_region_cells.tsv"))
stopifnot(nrow(trajectory) == 12642L, !anyDuplicated(trajectory$cell_id),
  nrow(windows) == 2524L, !anyDuplicated(windows$window_id),
  identical(bb$boundary, c("BB1", "BB2")),
  all(bb$pt_min < bb$pt_median), all(bb$pt_median < bb$pt_max),
  all(table(membership$window_id) == 30L), nrow(turning) == 4L,
  all(table(turning_cells$region_id) == 30L))

region_colors <- c(T1 = "#0072B2", T2 = "#D55E00", T3 = "#009E73",
  T4 = "#CC79A7")
turning$region_id <- factor(turning$region_id, levels = names(region_colors))
turning_cells$region_id <- factor(turning_cells$region_id,
  levels = names(region_colors))
stopifnot(!anyNA(turning$region_id), !anyNA(turning_cells$region_id))
turning_bands <- turning[, c("region_id", "pt_min", "pt_max")]

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint),
  "Tsix" %in% wb_feature_ids(bundle, "SCT_data"))

cells <- as.character(wb_cells(bundle, "all"))
stopifnot(!anyDuplicated(cells), all(trajectory$cell_id %in% cells),
  all(turning_cells$cell_id %in% cells))
tsix_matrix <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", "Tsix",
  cells, missing = "NA"))
stopifnot(identical(rownames(tsix_matrix), "Tsix"),
  identical(colnames(tsix_matrix), cells), all(is.finite(tsix_matrix)))
tsix <- as.numeric(tsix_matrix[1L, ])

cell_index <- match(membership$cell_id, cells)
window_index <- match(membership$window_id, windows$window_id)
stopifnot(!anyNA(cell_index), !anyNA(window_index))
membership_split <- split(cell_index, factor(membership$window_id,
  levels = windows$window_id))
stopifnot(length(membership_split) == nrow(windows),
  all(lengths(membership_split) == 30L))
window_mean <- vapply(membership_split, function(ix) mean(tsix[ix]), numeric(1))
fit_data <- data.frame(x = seq_len(nrow(windows)), response = window_mean)
robust_fallback <- FALSE
fit <- tryCatch(loess(response ~ x, data = fit_data, span = loess_span,
    degree = 2L, family = "symmetric", normalize = TRUE,
    control = loess.control(surface = "direct", statistics = "approximate")),
  error = function(e) {
    robust_fallback <<- TRUE
    loess(response ~ x, data = fit_data, span = loess_span,
      degree = 2L, family = "gaussian", normalize = TRUE,
      control = loess.control(surface = "direct", statistics = "approximate"))
  })
smoothed <- as.numeric(predict(fit, newdata = data.frame(x = fit_data$x)))
stopifnot(length(smoothed) == nrow(windows), all(is.finite(smoothed)))
expression_windows <- data.frame(window_id = windows$window_id,
  window_index = seq_len(nrow(windows)), pt_median = windows$pt_median,
  mean_sct_data = window_mean, smoothed_mean_sct_data = smoothed)
turning_expression <- merge(turning[, c("window_id", "region_id", "pt")],
  expression_windows[, c("window_id", "smoothed_mean_sct_data")],
  by = "window_id", sort = FALSE)
turning_expression <- turning_expression[match(names(region_colors),
  as.character(turning_expression$region_id)), ]
turning_expression$region_id <- factor(turning_expression$region_id,
  levels = names(region_colors))

space <- wb_space(bundle, "umap.atac", cells = cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells),
    all(c("x", "y") %in% names(space)))
  umap <- data.frame(cell_id = cells, x = as.numeric(space$x),
    y = as.numeric(space$y), sct_data = tsix, stringsAsFactors = FALSE)
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  umap <- data.frame(cell_id = cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]), sct_data = tsix, stringsAsFactors = FALSE)
}
positive <- umap$sct_data[umap$sct_data > 0]
stopifnot(length(positive) > 0L, all(is.finite(umap$x)), all(is.finite(umap$y)))
positive_p99 <- as.numeric(quantile(positive, 0.99, names = FALSE, type = 7))
stopifnot(is.finite(positive_p99), positive_p99 > 0)
contrast_quantiles <- as.numeric(quantile(positive, c(0, 0.75, 0.90, 0.99),
  names = FALSE, type = 7))
stopifnot(all(is.finite(contrast_quantiles)),
  all(diff(contrast_quantiles) > 0))
umap$relative_expression <- ifelse(umap$sct_data > 0,
  pmin(umap$sct_data, positive_p99) / positive_p99, NA_real_)
umap <- umap[order(ifelse(is.na(umap$relative_expression), -Inf,
  umap$relative_expression)), ]
highlight <- merge(umap[, c("cell_id", "x", "y", "sct_data")],
  turning_cells[, c("cell_id", "region_id")], by = "cell_id", all = FALSE,
  sort = FALSE)
highlight$region_id <- factor(highlight$region_id, levels = names(region_colors))
centroids <- aggregate(cbind(x, y) ~ region_id, highlight, median)
stopifnot(nrow(highlight) == 120L, nrow(centroids) == 4L)

umap_plot <- ggplot(umap, aes(x, y, fill = relative_expression)) +
  geom_point(shape = 21, color = "#00000000", size = 0.24, alpha = 0.90,
    stroke = 0) +
  geom_point(data = highlight, aes(x, y, color = region_id), inherit.aes = FALSE,
    shape = 21, fill = NA, size = 0.95, stroke = 0.28) +
  geom_label(data = centroids, aes(x, y, label = region_id,
      color = region_id), inherit.aes = FALSE, fill = "white",
    linewidth = 0.18, fontface = "bold", size = 3.2, show.legend = FALSE) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1),
    oob = scales::squish, na.value = "#E2E2E2",
    breaks = c(0, 0.5, 1), labels = c("low", "mid", "p99+"),
    name = "Positive-cell\nSCT_data") +
  scale_color_manual(values = region_colors, guide = "none") +
  coord_equal() + theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9), legend.position = "bottom",
    legend.key.width = grid::unit(20, "mm")) +
  labs(title = "A  Tsix expression on Bdev umap.atac",
    subtitle = "Gray: SCT_data ≤ 0; colored outlines: exact T1–T4 cells")

pt_plot <- ggplot(expression_windows,
    aes(pt_median, smoothed_mean_sct_data)) +
  geom_rect(data = turning_bands,
    aes(xmin = pt_min, xmax = pt_max, fill = region_id),
    ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.20) +
  geom_line(color = "#8E0152", linewidth = 1.0) +
  geom_point(data = turning_expression,
    aes(pt, smoothed_mean_sct_data, color = region_id), inherit.aes = FALSE,
    size = 2.4) +
  geom_text(data = turning_expression,
    aes(pt, smoothed_mean_sct_data, label = region_id, color = region_id),
    inherit.aes = FALSE, vjust = -0.9, size = 3.2, fontface = "bold",
    show.legend = FALSE) +
  scale_fill_manual(values = region_colors, guide = "none") +
  scale_color_manual(values = region_colors, guide = "none") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9)) +
  labs(title = "B  Tsix expression along Bdev pseudotime",
    subtitle = "Robust-LOESS of 30-cell means at 5-cell stride",
    x = "Bdev pseudotime (window median pt)", y = "Mean SCT_data")

bb_fill <- c(BB1 = "#2C7FB8", BB2 = "#D95F0E")
bb$label_y <- max(expression_windows$smoothed_mean_sct_data) * 1.055
bb_y_max <- max(expression_windows$smoothed_mean_sct_data) * 1.12
bb_pt_plot <- ggplot(expression_windows,
    aes(pt_median, smoothed_mean_sct_data)) +
  geom_rect(data = bb, aes(xmin = pt_min, xmax = pt_max, fill = boundary),
    ymin = 0, ymax = Inf, inherit.aes = FALSE, alpha = 0.17) +
  geom_line(color = "#8E0152", linewidth = 1.05) +
  geom_vline(data = bb, aes(xintercept = pt_median, color = boundary),
    linetype = "dashed", linewidth = 0.65, show.legend = FALSE) +
  geom_label(data = bb, aes(pt_median, label_y, label = boundary,
      fill = boundary), inherit.aes = FALSE, color = "white", linewidth = 0,
    fontface = "bold", size = 3.3, show.legend = FALSE) +
  scale_fill_manual(values = bb_fill, name = "Source bb_up region") +
  scale_color_manual(values = bb_fill, guide = "none") +
  scale_y_continuous(limits = c(0, bb_y_max),
    expand = expansion(mult = c(0, 0))) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10), legend.position = "right") +
  labs(title = "Tsix expression along Bdev pseudotime",
    subtitle = paste(
      "BB1/BB2 bands span the source bb_up pseudotime ranges;",
      "dashed lines mark their medians"),
    x = "Bdev pseudotime (window median pt)", y = "Mean SCT_data")

combined <- umap_plot | pt_plot
combined <- combined + plot_layout(widths = c(1, 1.12)) +
  plot_annotation(title = "Tsix across B-cell development",
    subtitle = "The established RNA–ATAC turning windows are marked as T1–T4",
    caption = paste0("UMAP positive expression is clipped at p99 = ",
      signif(positive_p99, 4), ". Bdev bundle opened read-only."),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 8.5)))

# A separate, uncluttered UMAP uses a conventional linear white-to-crimson
# expression scale. Values above the positive-cell p99 are clipped only for
# display. No turning-window layer is included.
umap$display_sct_data <- pmax(umap$sct_data, 0)
umap <- umap[order(umap$display_sct_data), ]
high_contrast_umap <- ggplot(umap,
    aes(x, y, color = display_sct_data)) +
  geom_point(size = 0.34, alpha = 1, stroke = 0) +
  scale_color_gradient(low = "#FFFFFF", high = "#C51B5A",
    limits = c(0, contrast_quantiles[[4L]]), oob = scales::squish,
    breaks = c(0, contrast_quantiles[2:4]),
    labels = c("0", "p75", "p90", "p99+"), name = NULL) +
  coord_equal() + theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 8.5, hjust = 0),
    plot.margin = margin(10, 12, 14, 12), legend.position = "bottom",
    legend.key.width = grid::unit(24, "mm")) +
  labs(title = "Tsix expression on Bdev umap.atac",
    subtitle = "High-expression contrast; no T1–T4 overlay",
    caption = paste(
      "White: zero/non-positive. Crimson: positive-cell p99 or higher;",
      "values above p99 are clipped."))

ggsave(file.path(stage_dir, "tsix_umap_expression_and_pseudotime.png"),
  combined, width = 12.5, height = 5.8, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "tsix_umap_expression_and_pseudotime.pdf"),
  combined, width = 12.5, height = 5.8, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "tsix_umap_high_expression_contrast.png"),
  high_contrast_umap, width = 7.2, height = 7.2, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "tsix_umap_high_expression_contrast.pdf"),
  high_contrast_umap, width = 7.2, height = 7.2, device = cairo_pdf,
  bg = "white")
ggsave(file.path(stage_dir, "tsix_pseudotime_bb_regions_zero_baseline.png"),
  bb_pt_plot, width = 8.4, height = 5.3, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "tsix_pseudotime_bb_regions_zero_baseline.pdf"),
  bb_pt_plot, width = 8.4, height = 5.3, device = cairo_pdf, bg = "white")
write_tsv(expression_windows,
  file.path(stage_dir, "tsix_window_expression.tsv"))
write_tsv(turning_expression,
  file.path(stage_dir, "tsix_turning_point_expression.tsv"))
cell_expression_connection <- gzfile(file.path(stage_dir,
  "tsix_cell_expression.tsv.gz"), open = "wt")
write.table(data.frame(cell_id = cells, sct_data = tsix),
  cell_expression_connection, sep = "\t", quote = TRUE, row.names = FALSE,
  na = "NA")
close(cell_expression_connection)
diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  gene = "Tsix", expression_matrix = "SCT_data", embedding = "umap.atac",
  trajectory_cells = nrow(trajectory), all_bundle_cells = length(cells),
  positive_cells = length(positive), positive_cell_p99 = positive_p99,
  standalone_umap_color_knots = as.list(setNames(contrast_quantiles,
    c("positive_min", "positive_p75", "positive_p90", "positive_p99"))),
  standalone_umap_color_scale = "linear white-to-crimson from zero to positive-cell p99",
  standalone_umap_turning_overlay = FALSE,
  standalone_pseudotime_regions = c("BB1", "BB2"),
  standalone_pseudotime_region_rule = "source pt min-to-max for bb_up up1/up2; dashed line at median",
  standalone_pseudotime_y_min = 0,
  window_size = 30L, stride = 5L, loess_span = loess_span,
  robust_loess_fallback = robust_fallback, turning_regions = nrow(turning),
  highlighted_cells_per_region = 30L, bundle_mutated = FALSE,
  capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Preliminary Tsix expression figure", "",
  "Tsix SCT_data expression is shown on the existing Bdev `umap.atac` and",
  "along the accepted Bdev pseudotime. The pseudotime curve is a LOESS-smoothed",
  "series of the same 30-cell window means at 5-cell stride used in the RNA–ATAC",
  "distance analysis. T1–T4 are fixed from the existing turning-window report;",
  "they were not re-called.", "",
  "A separate UMAP omits T1–T4 and uses a linear white-to-crimson scale from",
  "zero to the positive-cell p99; values above p99 are clipped for display.", "",
  "A separate pseudotime figure replaces T1–T4 with the source BB1/BB2",
  "pseudotime ranges and medians, and fixes the expression y-axis minimum at zero.",
  "The bundle was opened read-only."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed preliminary Tsix report at", output_dir, "\n")
