#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L)
  stop("Usage: 11_trajectory_geneset_region_z_histogram.R <report_root>")
report_root <- normalizePath(args[[1L]], mustWork = TRUE)
input_path <- file.path(report_root, "all_significant_positive_regions.tsv.gz")
stopifnot(file.exists(input_path))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

input_connection <- gzfile(input_path, "rt")
regions <- as.data.table(read.delim(input_connection, check.names = FALSE,
  stringsAsFactors = FALSE))
close(input_connection)
stopifnot(nrow(regions) > 0L, all(is.finite(regions$peak_z)))
labels <- c(cytokines = "Cytokines", bp = "GO BP", mf = "GO MF", cc = "GO CC")
regions[, group_label := factor(labels[analysis_group], levels = unname(labels))]
summary <- regions[, .(n_regions = .N, minimum_z = min(peak_z),
  q25_z = quantile(peak_z, .25), median_z = median(peak_z),
  mean_z = mean(peak_z), q75_z = quantile(peak_z, .75),
  maximum_z = max(peak_z)), by = .(analysis_group, group_label)]
setorder(summary, group_label)

annotation <- summary[, .(group_label,
  label = paste0("n = ", format(n_regions, big.mark = ","),
    "\nmedian = ", sprintf("%.2f", median_z)), x = 49, y = Inf)]
colors <- c("Cytokines" = "#B91C1C", "GO BP" = "#2563EB",
  "GO MF" = "#7C3AED", "GO CC" = "#047857")
p <- ggplot(regions, aes(peak_z, fill = group_label)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "white", linewidth = .2) +
  geom_vline(data = summary, aes(xintercept = median_z),
    linewidth = .6, linetype = "22", color = "#252525") +
  geom_text(data = annotation, aes(x = x, y = y, label = label),
    inherit.aes = FALSE, hjust = 1, vjust = 1.25, size = 3.2) +
  facet_wrap(~group_label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = colors, guide = "none") +
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, 10),
    expand = expansion(mult = c(0, .01))) +
  labs(title = "Peak z-scores of significant trajectory regions",
    subtitle = "Final q < 0.01, prominence-filtered and merged regions; dashed line is the group median",
    x = "Peak z-score", y = "Number of regions") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F3F4F6", color = "#D1D5DB"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"))

output_dir <- file.path(report_root, "summary_figures")
dir.create(output_dir, showWarnings = FALSE)
ggsave(file.path(output_dir, "significant_region_peak_z_histogram.png"), p,
  width = 10, height = 7.5, dpi = 220, bg = "white")
ggsave(file.path(output_dir, "significant_region_peak_z_histogram.pdf"), p,
  width = 10, height = 7.5, device = cairo_pdf, bg = "white")
fwrite(summary, file.path(output_dir, "significant_region_peak_z_summary.tsv"),
  sep = "\t", quote = FALSE, na = "NA")
