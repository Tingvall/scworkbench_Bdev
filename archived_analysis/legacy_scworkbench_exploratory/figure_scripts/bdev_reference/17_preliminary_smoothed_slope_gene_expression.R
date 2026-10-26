#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 17_preliminary_smoothed_slope_gene_expression.R",
  "<bdev_bundle> <w30_s5_report> <output_dir>",
  "<expected_bundle_fingerprint> <loess_span>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
source_report <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
expected_fingerprint <- args[[4L]]
loess_span <- as.numeric(args[[5L]])
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
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}

required <- c("window_distance_and_slope.tsv", "window_membership.tsv.gz",
  "bb1_bb2_regions.tsv", "trajectory_cells_ordered.tsv", "diagnostics.json")
stopifnot(all(file.exists(file.path(source_report, required))))
source_diagnostics <- read_json(file.path(source_report, "diagnostics.json"),
  simplifyVector = TRUE)
stopifnot(identical(source_diagnostics$status, "complete"),
  identical(as.integer(source_diagnostics$windows$size), 30L),
  identical(as.integer(source_diagnostics$windows$stride), 5L),
  identical(source_diagnostics$bundle_id, "bdev_reference"))

windows <- read_tsv(file.path(source_report, "window_distance_and_slope.tsv"))
membership_connection <- gzfile(file.path(source_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
bb <- read_tsv(file.path(source_report, "bb1_bb2_regions.tsv"))
trajectory <- read_tsv(file.path(source_report, "trajectory_cells_ordered.tsv"))
stopifnot(nrow(windows) == 2524L,
  identical(as.character(windows$window_id), unique(membership$window_id)),
  all(table(membership$window_id) == 30L),
  nrow(trajectory) == 12642L, !anyDuplicated(trajectory$cell_id),
  all(is.finite(windows$local_slope_rna_per_atac)),
  all(windows$local_slope_rna_per_atac > 0))

smooth_by_window <- function(values, log10_response = FALSE) {
  response <- if (log10_response) log10(values) else values
  fit_data <- data.frame(x = windows$window_index, response = response)
  fit <- loess(response ~ x, data = fit_data, span = loess_span,
    degree = 2L, family = "symmetric", normalize = TRUE,
    control = loess.control(surface = "direct", statistics = "approximate"))
  predicted <- as.numeric(predict(fit,
    newdata = data.frame(x = windows$window_index)))
  stopifnot(length(predicted) == length(values), all(is.finite(predicted)))
  if (log10_response) 10^predicted else predicted
}

windows$smoothed_local_slope <- smooth_by_window(
  windows$local_slope_rna_per_atac, log10_response = TRUE)
stopifnot(all(is.finite(windows$smoothed_local_slope),
  all(windows$smoothed_local_slope > 0)))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))

genes <- c("Ebf1", "Pax5", "Ezh2")
available <- wb_feature_ids(bundle, "SCT_data")
stopifnot(all(genes %in% available))
log_event("Reading SCT_data for", paste(genes, collapse = ", "))
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", genes,
  trajectory$cell_id, missing = "NA"))
stopifnot(identical(rownames(expression), genes),
  identical(colnames(expression), trajectory$cell_id),
  all(is.finite(expression)))

cell_index <- match(membership$cell_id, trajectory$cell_id)
stopifnot(!anyNA(cell_index))
membership_split <- split(cell_index, factor(membership$window_id,
  levels = windows$window_id))
stopifnot(length(membership_split) == nrow(windows),
  all(lengths(membership_split) == 30L))

expression_rows <- lapply(seq_along(genes), function(g) {
  raw <- vapply(membership_split, function(ix) mean(expression[g, ix]),
    numeric(1))
  data.frame(window_id = windows$window_id,
    window_index = windows$window_index, pt_median = windows$pt_median,
    gene = genes[[g]], mean_sct_data = raw,
    smoothed_mean_sct_data = smooth_by_window(raw),
    stringsAsFactors = FALSE)
})
expression_windows <- do.call(rbind, expression_rows)
expression_windows$gene <- factor(expression_windows$gene, levels = genes)
stopifnot(nrow(expression_windows) == length(genes) * nrow(windows),
  all(is.finite(expression_windows$mean_sct_data)),
  all(is.finite(expression_windows$smoothed_mean_sct_data)))

write_tsv(windows, file.path(stage_dir, "smoothed_slope_windows.tsv"))
write_tsv(expression_windows,
  file.path(stage_dir, "ebf1_pax5_ezh2_window_expression.tsv"))
write_tsv(bb, file.path(stage_dir, "bb1_bb2_regions.tsv"))

bb_fill <- c(BB1 = "#2C7FB8", BB2 = "#D95F0E")
slope_plot <- ggplot(windows, aes(pt_median, local_slope_rna_per_atac)) +
  geom_rect(data = bb, aes(xmin = pt_min, xmax = pt_max, fill = boundary),
    ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.13) +
  geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.7,
    color = "#4A4A4A") +
  geom_line(color = "#8F79AE", linewidth = 0.28, alpha = 0.28) +
  geom_line(aes(y = smoothed_local_slope), color = "#4B1D82",
    linewidth = 1.15) +
  geom_vline(data = bb, aes(xintercept = pt_median, color = boundary),
    linetype = "dashed", linewidth = 0.65, show.legend = FALSE) +
  scale_y_log10() +
  scale_fill_manual(values = bb_fill, name = "Source bb_up region") +
  scale_color_manual(values = bb_fill) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
    legend.position = "right") +
  labs(title = "A  Smoothed local RNA-versus-ATAC progression slope",
    subtitle = paste0("Robust LOESS on log10(slope), span = ", loess_span,
      "; faint line shows unsmoothed values"),
    x = NULL,
    y = "Local slope\nΔRNA / ΔATAC (log10 scale)")

gene_colors <- c(Ebf1 = "#0072B2", Pax5 = "#009E73", Ezh2 = "#CC79A7")
expression_plot <- ggplot(expression_windows,
  aes(pt_median, mean_sct_data, group = gene)) +
  geom_rect(data = bb, aes(xmin = pt_min, xmax = pt_max, fill = boundary),
    ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.10) +
  geom_line(aes(color = gene), linewidth = 0.25, alpha = 0.22) +
  geom_line(aes(y = smoothed_mean_sct_data, color = gene), linewidth = 1) +
  geom_vline(data = bb, aes(xintercept = pt_median), linetype = "dashed",
    color = "#555555", linewidth = 0.5) +
  facet_wrap(~ gene, nrow = 1L, scales = "free_y") +
  scale_color_manual(values = gene_colors, guide = "none") +
  scale_fill_manual(values = bb_fill, guide = "none") +
  theme_classic(base_size = 10.5) +
  theme(strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold")) +
  labs(title = "B–D  Sliding-window expression in Bdev SCT_data",
    subtitle = paste0("Window means with robust LOESS, span = ", loess_span,
      "; panels have independent expression ranges"),
    x = "Bdev pseudotime (window median pt)",
    y = "Mean SCT_data expression")

combined <- slope_plot / expression_plot +
  plot_layout(heights = c(1.12, 1)) +
  plot_annotation(
    title = "Smoothed modality progression and B-cell regulator expression",
    subtitle = "Preliminary 30-cell windows with 5-cell stride; BB1/BB2 source ranges are shaded",
    caption = paste(
      "Preliminary, unregistered sensitivity display. Source Bdev pt values and",
      "missingness are unchanged; smoothing affects display summaries only."),
    theme = theme(plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      plot.caption = element_text(size = 9)))

ggsave(file.path(stage_dir,
  "smoothed_slope_ebf1_pax5_ezh2_w30_s5.png"), combined,
  width = 13, height = 10, dpi = 320, bg = "white")
ggsave(file.path(stage_dir,
  "smoothed_slope_ebf1_pax5_ezh2_w30_s5.pdf"), combined,
  width = 13, height = 10, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id,
  bundle_fingerprint = expected_fingerprint,
  source_report_bundle_fingerprint = source_diagnostics$bundle_fingerprint,
  source_report = source_report,
  source_window_size = 30L, source_window_stride = 5L,
  source_window_count = nrow(windows),
  trajectory_cells = nrow(trajectory),
  slope_smoothing = list(method = "robust LOESS",
    response_scale = "log10", span = loess_span, degree = 2L,
    surface = "direct"),
  expression = list(matrix = "SCT_data", genes = genes,
    window_summary = "arithmetic mean", smoothing = "robust LOESS",
    span = loess_span, degree = 2L),
  bb_regions = "source pt min/max and median for bb_up up1/up2",
  bundle_mutated = FALSE, capsule_registered = FALSE,
  intended_use = "preliminary sensitivity testing")
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))

report <- c(
  "# Preliminary smoothed slope and regulator expression",
  "",
  paste0("This read-only sensitivity display reuses the accepted 30-cell/5-cell",
    " windows and applies robust LOESS with span ", loess_span, "."),
  "",
  "- Slope is smoothed on the log10 scale; raw slope remains in the output table.",
  "- Ebf1, Pax5, and Ezh2 are window means from Bdev `SCT_data` and are smoothed separately.",
  "- BB1/BB2 shading and dashed medians use the existing `bb_up` source annotations.",
  "- Source pseudotime, latent coordinates, expression, and missingness are unchanged.",
  "",
  "This is explicitly preliminary and is not registered as durable scWorkbench evidence.")
writeLines(report, file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
log_event("Completed preliminary smoothed slope/expression figure at", output_dir)
