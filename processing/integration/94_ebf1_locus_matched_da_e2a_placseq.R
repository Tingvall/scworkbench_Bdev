#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 94_ebf1_locus_matched_da_e2a_placseq.R",
  "SOURCE_TCF3_REPORT MATCHED_DA_REPORT OUTPUT_DIR SCRIPT_VERSION"
), call. = FALSE)

source_report <- normalizePath(args[[1L]], mustWork = TRUE)
matched_report <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
script_version <- args[[4L]]

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
  library(digest)
})

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

chromosome <- "chr11"
ebf1_tss <- 44617317L
locus_start <- ebf1_tss - 2000000L
locus_end <- ebf1_tss + 2000000L

coverage_path <- file.path(source_report,
  "condition_mean_full_region_1kb_coverage.tsv")
tcf3_peak_path <- file.path(source_report, "eight_tcf3_peak_loci.tsv")
link_path <- file.path(source_report,
  "ebf1_promoter_placseq_links_in_evaluated_region.tsv")
da_path <- file.path(matched_report,
  "all_targets_depth_complexity_matched_ebf1_pm2mb_da.tsv")
for (path in c(coverage_path, tcf3_peak_path, link_path, da_path))
  stopifnot(file.exists(path))

coverage <- fread(coverage_path)[midpoint_mb >= locus_start / 1e6 &
  midpoint_mb <= locus_end / 1e6]
tcf3_peaks <- fread(tcf3_peak_path)
links <- fread(link_path)[distal_end >= locus_start & distal_start <= locus_end]
da <- fread(da_path)
stopifnot(nrow(da) == 502L, uniqueN(da$feature_id) == 251L,
  da[target == "Flt3l", sum(matched_significant_bh05)] == 5L,
  da[target == "Il7", sum(matched_significant_bh05)] == 0L)

da[, `:=`(
  peak_midpoint_mb = (peak_start_1based + peak_end_1based) / 2e6,
  target_label = factor(target, levels = c("Flt3l", "Il7"),
    labels = c("Flt3L matched: WT - KO", "IL-7 matched: WT - KO")),
  significance = fifelse(matched_significant_bh05 &
    matched_direction == "WT_higher", "WT higher (BH < 0.05)",
    fifelse(matched_significant_bh05 & matched_direction == "KO_higher",
      "KO higher (BH < 0.05)", "Not significant")))]
significant <- da[matched_significant_bh05 == TRUE]
significant[, label := fifelse(
  peak_start_1based <= ebf1_tss & peak_end_1based >= ebf1_tss,
  "Ebf1 TSS", sprintf("%+.2f Mb", midpoint_relative_to_tss_bp / 1e6))]
setorder(significant, peak_start_1based)
write_tsv(significant, file.path(stage_dir,
  "five_remaining_depth_complexity_matched_peaks.tsv"))

condition_order <- c("control", "TCF3_depleted")
coverage[, condition := factor(condition, levels = condition_order)]
condition_labels <- c(control = "Control mean (n=3)",
  TCF3_depleted = "E2A/TCF3 depleted mean (n=3)")
condition_colors <- c(control = "#D55E00", TCF3_depleted = "#3264A8")
significance_colors <- c(
  "WT higher (BH < 0.05)" = "#C84545",
  "KO higher (BH < 0.05)" = "#3B77A8",
  "Not significant" = "grey78")
tcf3_marks <- tcf3_peaks[peak_midpoint >= locus_start &
  peak_midpoint <= locus_end, .(x = peak_midpoint / 1e6,
    label = paste0("P", peak_index))]

highlight_rects <- significant[, .(
  xmin = peak_start_1based / 1e6,
  xmax = peak_end_1based / 1e6,
  significance)]
ymax <- max(coverage$smoothed_coverage, na.rm = TRUE)
coverage_plot <- ggplot(coverage,
  aes(midpoint_mb, smoothed_coverage, color = condition,
    fill = condition)) +
  geom_rect(data = highlight_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
      fill = significance), inherit.aes = FALSE, alpha = 0.14,
    color = NA) +
  geom_vline(xintercept = ebf1_tss / 1e6, color = "#222222",
    linewidth = 0.45, linetype = "dashed") +
  geom_vline(data = tcf3_marks, aes(xintercept = x), inherit.aes = FALSE,
    color = "#8B1A1A", linewidth = 0.25, alpha = 0.7) +
  geom_area(alpha = 0.60, color = NA, na.rm = TRUE) +
  geom_line(linewidth = 0.52, na.rm = TRUE) +
  geom_text(data = tcf3_marks, aes(x, ymax * 1.01, label = label),
    inherit.aes = FALSE, size = 2.4, color = "#8B1A1A", vjust = 0) +
  facet_grid(condition ~ ., labeller = as_labeller(condition_labels)) +
  scale_color_manual(values = condition_colors, guide = "none") +
  scale_fill_manual(values = c(condition_colors, significance_colors),
    guide = "none") +
  scale_x_continuous(limits = c(locus_start, locus_end) / 1e6,
    expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, ymax * 1.08), clip = "off") +
  labs(x = NULL,
    y = "Mean fragment coverage CPM\n(1-kb bins; 5-kb rolling mean)",
    title = "Ebf1 locus after depth/complexity-matched WT/KO analysis",
    subtitle = paste(
      "E2A/TCF3 coverage, five retained Flt3L DA peaks, and",
      "Ebf1-promoter PLAC-seq; mm10")) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
    strip.text.y = element_text(face = "bold"),
    plot.margin = margin(8, 8, 2, 8))

da_plot <- ggplot(da, aes(peak_midpoint_mb,
  matched_detection_difference_wt_minus_ko, color = significance)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_vline(xintercept = ebf1_tss / 1e6, color = "#222222",
    linewidth = 0.45, linetype = "dashed") +
  geom_segment(aes(xend = peak_midpoint_mb, y = 0,
    yend = matched_detection_difference_wt_minus_ko),
    linewidth = 0.55, alpha = 0.85) +
  geom_point(size = 1.5, alpha = 0.95) +
  geom_text(data = significant,
    aes(x = peak_midpoint_mb,
      y = matched_detection_difference_wt_minus_ko,
      label = label), inherit.aes = FALSE,
    color = "#222222", size = 2.6, vjust = -0.7,
    check_overlap = TRUE) +
  facet_grid(target_label ~ .) +
  scale_color_manual(values = significance_colors,
    breaks = names(significance_colors), name = NULL) +
  scale_x_continuous(limits = c(locus_start, locus_end) / 1e6,
    expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Accessible-cell\nfraction: WT - KO",
    subtitle = paste0(
      "Within-bin ATAC depth/complexity matching: Flt3L 5/251 BH-significant; ",
      "IL-7 0/251")) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
    strip.text.y = element_text(face = "bold"),
    plot.margin = margin(1, 8, 2, 8))

links[, `:=`(
  promoter_midpoint = (promoter_start + promoter_end) / 2,
  distal_midpoint = (distal_start + distal_end) / 2,
  minus_log10_q = -log10(pmax(`Q-Value_Bias`, 1e-300)))]
arc_rows <- rbindlist(lapply(seq_len(nrow(links)), function(i) {
  x1 <- links$promoter_midpoint[[i]] / 1e6
  x2 <- links$distal_midpoint[[i]] / 1e6
  t <- seq(0, 1, length.out = 80L)
  data.table(interaction_id = links$interaction_id[[i]],
    x = x1 + (x2 - x1) * t,
    y = -sin(pi * t) * sqrt(abs(x2 - x1)),
    minus_log10_q = links$minus_log10_q[[i]])
}))
plac_plot <- ggplot(arc_rows,
  aes(x, y, group = interaction_id, color = minus_log10_q)) +
  geom_path(linewidth = 0.52, alpha = 0.72) +
  geom_vline(xintercept = ebf1_tss / 1e6, color = "#222222",
    linewidth = 0.45, linetype = "dashed") +
  geom_vline(data = significant,
    aes(xintercept = peak_midpoint_mb), inherit.aes = FALSE,
    color = "grey55", linewidth = 0.25, linetype = "dotted") +
  scale_color_viridis_c(option = "C", direction = -1,
    name = expression(-log[10](PLAC~q))) +
  scale_x_continuous(limits = c(locus_start, locus_end) / 1e6,
    expand = expansion(mult = c(0, 0))) +
  labs(x = "chr11 position (Mb; mm10)", y = NULL,
    subtitle = paste0(nrow(links),
      " 230238 H3K4me3 PLAC-seq Ebf1-promoter links")) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank(), legend.position = "right",
    plot.margin = margin(0, 8, 8, 8))

figure <- coverage_plot / da_plot / plac_plot +
  plot_layout(heights = c(2.25, 1.8, 1.35))
stem <- "ebf1_pm2mb_e2a_tcf3_matched_wt_ko_da_placseq"
ggsave(file.path(stage_dir, "figures", paste0(stem, ".png")), figure,
  width = 15, height = 10.5, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "figures", paste0(stem, ".pdf")), figure,
  width = 15, height = 10.5, device = cairo_pdf, bg = "white")

inputs <- data.table(
  role = c("source_coverage", "source_tcf3_peaks",
    "source_placseq_links", "matched_wt_ko_da"),
  path = c(coverage_path, tcf3_peak_path, link_path, da_path))
inputs[, `:=`(size_bytes = file.info(path)$size,
  sha256 = vapply(path, sha256_file, character(1L)))]
write_tsv(inputs, file.path(stage_dir, "input_manifest.tsv"))

checks <- list(
  tested_peaks_per_target = 251L,
  retained_significant_peaks = nrow(significant),
  retained_flt3l = significant[target == "Flt3l", .N],
  retained_il7 = significant[target == "Il7", .N],
  tss_peak_retained = significant[peak_start_1based <= ebf1_tss &
    peak_end_1based >= ebf1_tss, .N] == 1L,
  source_coverage_values_changed = FALSE,
  source_placseq_values_changed = FALSE,
  figures_exist = all(file.exists(file.path(stage_dir, "figures",
    paste0(stem, c(".png", ".pdf"))))),
  promotions = 0L)
stopifnot(checks$retained_significant_peaks == 5L,
  checks$retained_flt3l == 5L, checks$retained_il7 == 0L,
  checks$tss_peak_retained, checks$figures_exist,
  checks$promotions == 0L)
write_json(list(status = "complete_unpromoted_report_only",
  script_version = script_version,
  region = list(chromosome = chromosome, start_1based = locus_start,
    end_1based = locus_end, ebf1_tss_1based = ebf1_tss),
  encoding = list(DA_y = "matched accessible-cell fraction WT minus KO",
    WT_higher = "red", KO_higher = "blue", nonsignificant = "gray"),
  checks = checks), file.path(stage_dir, "diagnostics.json"))

writeLines(c(
  "# Ebf1 locus with depth/complexity-matched WT/KO DA", "",
  "This Ebf1 ±2 Mb composite retains the original E2A/TCF3 coverage and Ebf1-promoter PLAC-seq context but replaces the pooled WT/KO tracks with the within-bin depth/complexity-matched sensitivity results.", "",
  paste0("![Matched Ebf1 locus](figures/", stem, ".png)"), "",
  "- Five Flt3L peaks remain BH-significant: one WT-higher Ebf1 TSS peak and four KO-higher distal peaks.",
  "- No IL-7 peaks remain BH-significant.",
  "- Colored translucent bands in the coverage panels locate the five retained peaks; P1-P8 remain the original E2A/TCF3-dependent peak labels.",
  "- Coverage and PLAC-seq values are unchanged from the source report.", "",
  paste0("- [Vector PDF](figures/", stem, ".pdf)"),
  "- [Five retained peaks](five_remaining_depth_complexity_matched_peaks.tsv)",
  "- [Input manifest](input_manifest.tsv)"
), file.path(stage_dir, "REPORT.md"))

capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
if (!file.rename(stage_dir, output_dir))
  stop("Failed to finalize output directory", call. = FALSE)
cat("Completed matched Ebf1 locus figure at", output_dir, "\n")
