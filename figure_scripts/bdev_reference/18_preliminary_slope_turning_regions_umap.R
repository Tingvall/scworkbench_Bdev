#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 18_preliminary_slope_turning_regions_umap.R",
  "<bdev_bundle> <w30_s5_report> <smoothed_report> <output_dir>",
  "<expected_bundle_fingerprint> <minimum_prominence_fold>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
smoothed_report <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_fingerprint <- args[[5L]]
minimum_prominence_fold <- as.numeric(args[[6L]])
stopifnot(is.finite(minimum_prominence_fold), minimum_prominence_fold > 1)
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
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}

required_smoothed <- c("smoothed_slope_windows.tsv",
  "ebf1_pax5_ezh2_window_expression.tsv", "diagnostics.json")
required_window <- c("window_membership.tsv.gz", "trajectory_cells_ordered.tsv",
  "diagnostics.json")
stopifnot(all(file.exists(file.path(smoothed_report, required_smoothed))),
  all(file.exists(file.path(window_report, required_window))))
smooth_diagnostics <- read_json(file.path(smoothed_report, "diagnostics.json"),
  simplifyVector = TRUE)
window_diagnostics <- read_json(file.path(window_report, "diagnostics.json"),
  simplifyVector = TRUE)
stopifnot(identical(smooth_diagnostics$status,
    "complete_preliminary_unregistered"),
  identical(window_diagnostics$status, "complete"),
  identical(as.integer(window_diagnostics$windows$size), 30L),
  identical(as.integer(window_diagnostics$windows$stride), 5L))

windows <- read_tsv(file.path(smoothed_report, "smoothed_slope_windows.tsv"))
expression_windows <- read_tsv(file.path(smoothed_report,
  "ebf1_pax5_ezh2_window_expression.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
stopifnot(nrow(windows) == 2524L, nrow(trajectory) == 12642L,
  !anyDuplicated(trajectory$cell_id),
  all(c("pt_min", "pt_median", "pt_max", "local_slope_rna_per_atac",
    "smoothed_local_slope") %in% names(windows)),
  all(windows$smoothed_local_slope > 0),
  all(expression_windows$gene %in% c("Ebf1", "Pax5", "Ezh2")))

# Find zero-crossings of the centered derivative of the already-smoothed
# log10 RNA/ATAC slope. Window index is used as the trajectory coordinate so
# tied pseudotime values cannot create zero denominators; crossings are reported
# in the original, unchanged Bdev pt coordinate.
y <- log10(windows$smoothed_local_slope)
n <- length(y)
derivative <- c(y[[2L]] - y[[1L]],
  (y[3L:n] - y[1L:(n - 2L)]) / 2,
  y[[n]] - y[[n - 1L]])
stopifnot(length(derivative) == n, all(is.finite(derivative)))
turn_index <- which(diff(sign(derivative)) != 0) + 1L
stopifnot(length(turn_index) > 0L)

# Prominence is the smaller log-slope excursion to the immediately adjacent
# extrema (or series endpoint). A 1.10-fold cutoff removes shallow wiggles.
extrema_with_ends <- c(1L, turn_index, n)
turning <- lapply(seq_along(turn_index), function(j) {
  i <- turn_index[[j]]
  is_minimum <- derivative[[i - 1L]] < 0 && derivative[[i + 1L]] > 0
  left_i <- extrema_with_ends[[j]]
  right_i <- extrema_with_ends[[j + 2L]]
  prominence_log10 <- if (is_minimum) {
    min(y[[left_i]] - y[[i]], y[[right_i]] - y[[i]])
  } else {
    min(y[[i]] - y[[left_i]], y[[i]] - y[[right_i]])
  }
  data.frame(window_index = i, window_id = windows$window_id[[i]],
    pt = windows$pt_median[[i]], pt_min = windows$pt_min[[i]],
    pt_max = windows$pt_max[[i]],
    smoothed_slope = windows$smoothed_local_slope[[i]],
    derivative_before = derivative[[i - 1L]],
    derivative_after = derivative[[i + 1L]],
    extremum = if (is_minimum) "minimum" else "maximum",
    interpretation = if (is_minimum) "RNA catch-up begins" else
      "ATAC catch-up begins",
    adjacent_left_pt = windows$pt_median[[left_i]],
    adjacent_right_pt = windows$pt_median[[right_i]],
    prominence_log10 = prominence_log10,
    prominence_fold = 10^prominence_log10,
    stringsAsFactors = FALSE)
})
turning <- do.call(rbind, turning)
turning$selected <- turning$prominence_fold >= minimum_prominence_fold
selected <- turning[turning$selected, , drop = FALSE]
selected$region_id <- paste0("T", seq_len(nrow(selected)))
stopifnot(nrow(selected) > 0L, !anyDuplicated(selected$window_id))

selected_membership <- merge(
  membership[membership$window_id %in% selected$window_id, , drop = FALSE],
  selected[, c("window_id", "region_id", "interpretation")],
  by = "window_id", all.x = TRUE, sort = FALSE)
selected_membership <- selected_membership[order(
  match(selected_membership$region_id, selected$region_id),
  selected_membership$position_in_window), ]
stopifnot(nrow(selected_membership) == 30L * nrow(selected),
  all(table(selected_membership$region_id) == 30L))

turning$region_id <- NA_character_
turning$region_id[turning$selected] <- selected$region_id
write_tsv(turning, file.path(stage_dir, "all_slope_turning_points.tsv"))
write_tsv(selected, file.path(stage_dir, "selected_turning_regions.tsv"))
write_tsv(selected_membership,
  file.path(stage_dir, "selected_turning_region_cells.tsv"))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))
cells <- as.character(wb_cells(bundle, "all"))
stopifnot(length(cells) == 26723L, !anyDuplicated(cells),
  all(selected_membership$cell_id %in% cells))
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
stopifnot(all(is.finite(umap$x)), all(is.finite(umap$y)))
highlight <- merge(umap, selected_membership[, c("cell_id", "region_id",
  "interpretation")], by = "cell_id", all = FALSE, sort = FALSE)
highlight$region_id <- factor(highlight$region_id, levels = selected$region_id)
centroids <- aggregate(cbind(x, y) ~ region_id, highlight, median)

region_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#F0E442", "#000000")
stopifnot(nrow(selected) <= length(region_colors))
region_colors <- setNames(region_colors[seq_len(nrow(selected))],
  selected$region_id)
selected$region_id <- factor(selected$region_id, levels = names(region_colors))

slope_plot <- ggplot(windows,
    aes(pt_median, smoothed_local_slope)) +
  geom_rect(data = selected,
    aes(xmin = pt_min, xmax = pt_max, fill = region_id),
    ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.20) +
  geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.65,
    color = "#4A4A4A") +
  geom_line(color = "#4B1D82", linewidth = 1.05) +
  geom_point(data = selected, aes(pt, smoothed_slope, color = region_id),
    inherit.aes = FALSE, size = 2.5) +
  geom_text(data = selected, aes(pt, smoothed_slope, label = region_id,
      color = region_id), inherit.aes = FALSE, vjust = -0.9,
    fontface = "bold", size = 3.4, show.legend = FALSE) +
  scale_y_log10() +
  scale_fill_manual(values = region_colors, name = "Turning window") +
  scale_color_manual(values = region_colors, name = "Turning window") +
  theme_classic(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"),
    legend.position = "right") +
  labs(title = "A  Persistent zero-crossings in slope-of-slope",
    subtitle = paste0("Selected extrema have ≥", minimum_prominence_fold,
      "-fold prominence; bands are their exact 30-cell windows"),
    x = NULL, y = "Smoothed local slope\nΔRNA / ΔATAC (log10)")

gene_order <- c("Ebf1", "Pax5", "Ezh2")
expression_windows$gene <- factor(expression_windows$gene, levels = gene_order)
gene_colors <- c(Ebf1 = "#0072B2", Pax5 = "#009E73", Ezh2 = "#CC79A7")
expression_plot <- ggplot(expression_windows,
    aes(pt_median, smoothed_mean_sct_data, group = gene)) +
  geom_rect(data = selected,
    aes(xmin = pt_min, xmax = pt_max, fill = region_id),
    ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.18) +
  geom_line(aes(color = gene), linewidth = 0.95) +
  facet_wrap(~ gene, nrow = 1L, scales = "free_y") +
  scale_color_manual(values = gene_colors, guide = "none") +
  scale_fill_manual(values = region_colors, guide = "none") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10.5),
    plot.title = element_text(face = "bold")) +
  labs(title = "B–D  Expression across the same turning windows",
    x = "Bdev pseudotime (window median pt)",
    y = "Smoothed mean SCT_data expression")

umap_plot <- ggplot(umap, aes(x, y)) +
  geom_point(color = "#D9D9D9", size = 0.20, alpha = 0.65, stroke = 0) +
  geom_point(data = highlight, aes(color = region_id), size = 0.75,
    alpha = 0.95, stroke = 0) +
  geom_label(data = centroids, aes(label = region_id, color = region_id),
    fill = "white", label.size = 0.18, fontface = "bold", size = 3,
    show.legend = FALSE) +
  scale_color_manual(values = region_colors, name = "Turning window") +
  coord_equal() +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9), legend.position = "right") +
  labs(title = "E  Turning-window cells on Bdev umap.atac",
    subtitle = "Gray: all Bdev cells; color: exact cells in each 30-cell window")

combined <- slope_plot / expression_plot / umap_plot +
  plot_layout(heights = c(1.08, 0.95, 1.12), guides = "collect") +
  plot_annotation(
    title = "Where RNA and ATAC progression begin to catch up",
    subtitle = paste0("Robust-LOESS slope turning points along Bdev pt; ",
      nrow(selected), " persistent regions retained"),
    caption = paste(
      "A minimum means the RNA/ATAC slope turns upward (RNA catch-up begins);",
      "a maximum means it turns downward (ATAC catch-up begins).",
      "Preliminary, read-only, and unregistered."),
    theme = theme(plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10.5),
      plot.caption = element_text(size = 8.5))) &
  theme(legend.position = "right")

ggsave(file.path(stage_dir,
  "slope_turning_regions_genes_umap_w30_s5.png"), combined,
  width = 13, height = 15, dpi = 320, bg = "white")
ggsave(file.path(stage_dir,
  "slope_turning_regions_genes_umap_w30_s5.pdf"), combined,
  width = 13, height = 15, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id,
  bundle_fingerprint = expected_fingerprint,
  source_window_report = window_report,
  source_smoothed_report = smoothed_report,
  trajectory_cells = nrow(trajectory), all_bdev_cells = length(cells),
  derivative = list(response = "log10(smoothed local RNA/ATAC slope)",
    coordinate = "ordered sliding-window index",
    method = "centered finite difference",
    candidate_zero_crossings = nrow(turning)),
  persistence_filter = list(method = paste(
    "minimum excursion to immediately adjacent extrema, on log10 scale"),
    minimum_prominence_fold = minimum_prominence_fold,
    selected_regions = nrow(selected)),
  region_definition = "exact 30-cell source window centered on selected extremum",
  embedding = "umap.atac", expression_matrix = "SCT_data",
  genes = gene_order, bundle_mutated = FALSE, capsule_registered = FALSE,
  intended_use = "preliminary sensitivity testing")
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))

region_lines <- apply(selected, 1L, function(z) paste0("- ", z[["region_id"]],
  ": pt ", sprintf("%.3f", as.numeric(z[["pt_min"]])), " to ",
  sprintf("%.3f", as.numeric(z[["pt_max"]])), "; center ",
  sprintf("%.3f", as.numeric(z[["pt"]])), "; ", z[["interpretation"]],
  "; prominence ", sprintf("%.3f", as.numeric(z[["prominence_fold"]])),
  "-fold."))
report <- c("# Preliminary slope-turning regions", "",
  paste0("Centered finite differences of the robust-LOESS log10 slope yielded ",
    nrow(turning), " candidate extrema. A ", minimum_prominence_fold,
    "-fold local-prominence filter retained ", nrow(selected),
    " turning windows:"), "", region_lines, "",
  "Each highlighted region is the exact 30-cell source window around the extremum.",
  "The same regions are shown on the Ebf1/Pax5/Ezh2 curves and Bdev umap.atac.",
  "The bundle was opened read-only and no capsule was registered.")
writeLines(report, file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
log_event("Completed preliminary slope-turning report at", output_dir)
