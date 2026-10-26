#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 24_preliminary_log2_modality_ratio.R",
  "<bdev_bundle> <smoothed_report> <turning_report> <output_dir>",
  "<expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
smoothed_report <- normalizePath(args[[2L]], mustWork = TRUE)
turning_report <- normalizePath(args[[3L]], mustWork = TRUE)
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

slope <- read_tsv(file.path(smoothed_report, "smoothed_slope_windows.tsv"))
turning <- read_tsv(file.path(turning_report, "selected_turning_regions.tsv"))
stopifnot(nrow(slope) == 2524L, !anyDuplicated(slope$window_id),
  all(is.finite(slope$smoothed_local_slope)),
  all(slope$smoothed_local_slope > 0), nrow(turning) == 4L,
  all(turning$window_id %in% slope$window_id))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))

region_colors <- c(T1 = "#0072B2", T2 = "#D55E00", T3 = "#009E73",
  T4 = "#CC79A7")
turning$region_id <- factor(turning$region_id, levels = names(region_colors))
stopifnot(!anyNA(turning$region_id))
turning_bands <- turning[, c("region_id", "pt_min", "pt_max")]

ratio <- data.frame(window_id = slope$window_id,
  window_index = slope$window_index, pt_median = slope$pt_median,
  smoothed_rna_per_atac = slope$smoothed_local_slope,
  log2_rna_per_atac = log2(slope$smoothed_local_slope),
  log2_atac_per_rna = -log2(slope$smoothed_local_slope))
turning_ratio <- merge(turning[, c("window_id", "region_id", "pt")],
  ratio[, c("window_id", "log2_rna_per_atac", "log2_atac_per_rna")],
  by = "window_id", sort = FALSE)
turning_ratio <- turning_ratio[match(names(region_colors),
  as.character(turning_ratio$region_id)), ]
turning_ratio$region_id <- factor(turning_ratio$region_id,
  levels = names(region_colors))
stopifnot(all(is.finite(ratio$log2_rna_per_atac)),
  max(abs(ratio$log2_rna_per_atac + ratio$log2_atac_per_rna)) < 1e-12)

y_limit <- ceiling(max(abs(ratio$log2_rna_per_atac)) * 10) / 10
x_range <- range(ratio$pt_median)
x_note <- x_range[[1L]] + 0.015 * diff(x_range)
ratio_plot <- ggplot(ratio, aes(pt_median, log2_rna_per_atac)) +
  geom_rect(data = turning_bands,
    aes(xmin = pt_min, xmax = pt_max, fill = region_id),
    ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.18) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.7,
    color = "#444444") +
  geom_line(color = "#4B1D82", linewidth = 1.05) +
  geom_point(data = turning_ratio,
    aes(pt, log2_rna_per_atac, color = region_id), inherit.aes = FALSE,
    size = 2.5) +
  geom_text(data = turning_ratio,
    aes(pt, log2_rna_per_atac, label = region_id, color = region_id),
    inherit.aes = FALSE, vjust = -0.9, fontface = "bold", size = 3.3,
    show.legend = FALSE) +
  annotate("text", x = x_note, y = 0.91 * y_limit, label = "RNA faster",
    hjust = 0, color = "#8E0152", fontface = "bold", size = 3.5) +
  annotate("text", x = x_note, y = -0.91 * y_limit, label = "ATAC faster",
    hjust = 0, color = "#2166AC", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = region_colors, guide = "none") +
  scale_color_manual(values = region_colors, guide = "none") +
  scale_y_continuous(limits = c(-y_limit, y_limit),
    breaks = seq(-1, 1, by = 0.5),
    sec.axis = sec_axis(~ -.,
      name = expression(log[2](Delta*ATAC/Delta*RNA)))) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 8.5, hjust = 0),
    axis.title.y.right = element_text(color = "#2166AC"),
    axis.title.y.left = element_text(color = "#8E0152")) +
  labs(title = "Local RNA–ATAC progression balance along Bdev pseudotime",
    subtitle = "Explicit log2 ratio; T1–T4 are the existing fixed turning windows",
    x = "Bdev pseudotime (window median pt)",
    y = expression(log[2](Delta*RNA/Delta*ATAC)),
    caption = paste(
      "Zero is matched progression. +1 means RNA changes 2× faster; −1 means",
      "ATAC changes 2× faster. The right axis is the exact inverse, not a separate fit."))

ggsave(file.path(stage_dir, "rna_atac_log2_ratio_with_inverse_axis.png"),
  ratio_plot, width = 9.2, height = 5.8, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "rna_atac_log2_ratio_with_inverse_axis.pdf"),
  ratio_plot, width = 9.2, height = 5.8, device = cairo_pdf, bg = "white")
write_tsv(ratio, file.path(stage_dir, "log2_modality_ratio_windows.tsv"))
write_tsv(turning_ratio,
  file.path(stage_dir, "log2_modality_ratio_turning_points.tsv"))

proposal <- c("# Proposed gene-pattern screen", "",
  "## Target", "",
  "Use the saved per-window `log2(delta RNA / delta ATAC)` curve as the target.",
  "The inverse target is exactly its negative: `log2(delta ATAC / delta RNA)`.",
  "Therefore one fitted score is sufficient: positive coefficients identify the",
  "RNA/ATAC pattern and negative coefficients identify the inverse pattern.", "",
  "## Two complementary definitions", "",
  "1. **State-following genes:** correlate each gene's standardized, smoothed",
  "   30-cell-window `SCT_data` expression profile with the target curve.",
  "2. **Change-following genes:** correlate the centered pseudotime derivative",
  "   of each expression profile with the target. This asks which genes rise or",
  "   fall specifically during modality-imbalanced intervals and is closer to a",
  "   dynamic interpretation.", "",
  "Report both scores. A gene can track the modality state without changing at",
  "the same instant, so the two screens should not be collapsed.", "",
  "## Recommended inference", "",
  "- Filter genes by detection and trajectory variance before testing.",
  "- Use the same 30-cell/5-cell windows and smoothing span for comparability.",
  "- Estimate zero-lag and prespecified lagged associations; a positive lag can",
  "  test whether accessibility imbalance precedes RNA expression.",
  "- Calibrate P values with pseudotime-block permutations or circular shifts,",
  "  not ordinary correlation tests, because overlapping windows are strongly",
  "  autocorrelated.",
  "- Control BH FDR across all tested genes and require an effect-size threshold",
  "  plus stability across full LSI1–50 and LSI2–50 sensitivity curves.",
  "- To limit circularity from deriving the RNA trajectory with the same genes",
  "  being tested, use feature-fold cross-fitting: compute the RNA path from one",
  "  gene fold and score genes held out from that fold, then rotate folds.", "",
  "## Interpretation", "",
  "Positive association = follows RNA-faster intervals. Negative association =",
  "follows the inverse, ATAC-faster intervals. Running a second ATAC/RNA screen",
  "would return identical magnitudes and P values with reversed signs/ranks.")
writeLines(proposal, file.path(stage_dir, "GENE_PATTERN_PROPOSAL.md"))

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  windows = nrow(ratio), target = "log2(smoothed delta-RNA/delta-ATAC)",
  inverse = "negative target; log2(smoothed delta-ATAC/delta-RNA)",
  inverse_max_identity_error = max(abs(
    ratio$log2_rna_per_atac + ratio$log2_atac_per_rna)),
  y_limit = y_limit, turning_regions = nrow(turning), bundle_mutated = FALSE,
  capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Preliminary explicit log2 modality-ratio figure", "",
  "The saved robust-LOESS-smoothed local RNA/ATAC slope is displayed as",
  "`log2(delta RNA / delta ATAC)`. Zero denotes matched local progression;",
  "positive values are RNA-faster and negative values ATAC-faster. The right",
  "axis displays the exact inverse. T1–T4 are unchanged fixed references.", "",
  "See `GENE_PATTERN_PROPOSAL.md` for a proposed direct/inverse gene screen."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed explicit log2 modality-ratio report at", output_dir, "\n")
