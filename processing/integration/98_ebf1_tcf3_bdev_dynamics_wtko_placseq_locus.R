#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(paste(
  "Usage: 98_ebf1_tcf3_bdev_dynamics_wtko_placseq_locus.R BUNDLE",
  "FIVE_REGION_REPORT RATIO_REPORT TCF3_PEAKS PLAC_LINKS OUTPUT_DIR",
  "EXPECTED_FINGERPRINT SCRIPT_VERSION"), call. = FALSE)

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
five_dir <- normalizePath(args[[2L]], mustWork = TRUE)
ratio_dir <- normalizePath(args[[3L]], mustWork = TRUE)
tcf3_path <- normalizePath(args[[4L]], mustWork = TRUE)
plac_path <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
expected_fp <- args[[7L]]
script_version <- args[[8L]]

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
})
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

chromosome <- "chr11"
ebf1_tss <- 44617317L
locus_start <- ebf1_tss - 250000L
locus_end <- ebf1_tss + 2100000L

tcf3 <- fread(tcf3_path)[tcf3_route_pass == TRUE]
setorder(tcf3, peak_start_1based, peak_end_1based)
stopifnot(nrow(tcf3) == 8L, !anyDuplicated(tcf3$bdev_feature_id),
  all(tcf3$tcf3_best_log2fc < -1), all(tcf3$tcf3_best_padj < 0.05))
tcf3[, `:=`(feature_id = bdev_feature_id,
  midpoint_1based = (peak_start_1based + peak_end_1based) / 2,
  midpoint_mb = (peak_start_1based + peak_end_1based) / 2e6,
  locus_label = paste0("P", seq_len(.N)))]
peak_ids <- tcf3$feature_id

assignments <- fread(cmd = paste("zcat", shQuote(file.path(five_dir,
  "bdev_five_region_cell_assignments.tsv.gz"))))
regions <- fread(file.path(five_dir, "bdev_five_region_definition.tsv"))
dynamic <- fread(file.path(five_dir, "ebf1_pm2mb_dynamic_peak_summary.tsv"))
stopifnot(!anyDuplicated(assignments$cell_id), nrow(regions) == 5L,
  all(regions$region_order == 1:5))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fp),
  all(assignments$cell_id %in% wb_cells(bundle, "set:bdev_reference")),
  all(peak_ids %in% wb_feature_ids(bundle, "ATAC_bdev_mm10_peaks")))
counts <- as.matrix(wb_matrix_read_features(bundle, "ATAC_bdev_mm10_peaks",
  peak_ids, assignments$cell_id, missing = "NA"))
stopifnot(identical(rownames(counts), peak_ids),
  identical(colnames(counts), assignments$cell_id), !anyNA(counts))

bdev_profiles <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
  r <- regions[i]
  idx <- which(assignments$region_id == r$region_id)
  data.table(feature_id = peak_ids, region_order = r$region_order,
    region_label = r$region_label, bdev_cells = length(idx),
    accessible_percent = 100 * rowMeans(counts[, idx, drop = FALSE] > 0))
}))
bdev_profiles <- merge(bdev_profiles, tcf3[, .(feature_id, locus_label,
  midpoint_mb)], by = "feature_id", all.x = TRUE)
bdev_profiles[, region_factor := factor(region_label,
  levels = rev(regions$region_label))]

dynamic_class <- tcf3[, .(feature_id, locus_label, midpoint_mb)]
dynamic_class <- merge(dynamic_class,
  dynamic[, .(feature_id, statistically_dynamic, robust_dynamic)],
  by = "feature_id", all.x = TRUE)
dynamic_class[, dynamic_class := fcase(
  robust_dynamic == TRUE, "Robust Bdev dynamic",
  statistically_dynamic == TRUE, "BH-significant Bdev dynamic",
  is.na(statistically_dynamic), "Outside ±2 Mb dynamic screen",
  default = "Not Bdev dynamic")]

ratios <- fread(file.path(ratio_dir,
  "all_tcf3_peaks_normalized_wt_ko_ratios.tsv"))
stopifnot(nrow(ratios) == 16L, setequal(ratios$feature_id, peak_ids))
ratios <- merge(ratios, tcf3[, .(feature_id, locus_label, midpoint_mb)],
  by = "feature_id", all.x = TRUE)
ratios[, target_label := factor(knockout_target,
  levels = c("Flt3l", "Il7"),
  labels = c("Flt3L WT / FLT3LKO", "IL-7 WT / IL7KO"))]
ratio_cap <- 2.25
ratios[, display_log2_ratio := pmax(-ratio_cap,
  pmin(ratio_cap, normalized_log2_wt_over_ko))]

links <- fread(plac_path)[distal_end >= locus_start & distal_start <= locus_end]
links[, `:=`(promoter_midpoint = (promoter_start + promoter_end) / 2,
  distal_midpoint = (distal_start + distal_end) / 2,
  minus_log10_q = -log10(pmax(`Q-Value_Bias`, 1e-300)))]
arc_rows <- rbindlist(lapply(seq_len(nrow(links)), function(i) {
  x1 <- links$promoter_midpoint[[i]] / 1e6
  x2 <- links$distal_midpoint[[i]] / 1e6
  z <- seq(0, 1, length.out = 80L)
  data.table(interaction_id = links$interaction_id[[i]],
    x = x1 + (x2 - x1) * z,
    y = -sin(pi * z) * sqrt(abs(x2 - x1)),
    minus_log10_q = links$minus_log10_q[[i]])
}))

x_limits <- c(locus_start, locus_end) / 1e6
tss_mb <- ebf1_tss / 1e6
common_theme <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
    plot.margin = margin(2, 8, 2, 8))

p_bdev <- ggplot(bdev_profiles,
  aes(midpoint_mb, region_factor, size = accessible_percent,
    fill = accessible_percent)) +
  geom_vline(xintercept = tss_mb, linetype = "dashed",
    linewidth = 0.45, color = "grey25") +
  geom_point(shape = 21, color = "grey25", stroke = 0.25) +
  geom_text(data = tcf3[, .(midpoint_mb, locus_label)],
    aes(midpoint_mb, y = 5.72, label = locus_label), inherit.aes = FALSE,
    size = 2.7, fontface = "bold", vjust = 0) +
  annotate("text", x = tss_mb, y = 5.72, label = "Ebf1 TSS",
    size = 2.7, hjust = -0.08, vjust = 0) +
  scale_size_area(max_size = 9, limits = c(0, 30),
    breaks = c(1, 5, 10, 20, 30), name = "Accessible Bdev cells") +
  scale_fill_viridis_c(option = "C", limits = c(0, 30),
    oob = scales::squish, labels = scales::label_percent(scale = 1),
    name = "Accessible Bdev cells") +
  scale_x_continuous(limits = x_limits, expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0.5, 5.9), clip = "off") +
  labs(x = NULL, y = NULL,
    title = "Ebf1 locus: Bdev dynamics, WT/KO balance, TCF3 dependence, and promoter contacts",
    subtitle = "Eight significant E2A/TCF3-depletion peaks; mm10") +
  common_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "right")

p_ratio <- ggplot(ratios,
  aes(midpoint_mb, target_label, fill = display_log2_ratio,
    shape = ratio_status)) +
  geom_vline(xintercept = tss_mb, linetype = "dashed",
    linewidth = 0.45, color = "grey25") +
  geom_point(size = 3.4, color = "grey20", stroke = 0.35) +
  scale_fill_gradient2(low = "#3264A8", mid = "white", high = "#C43C39",
    midpoint = 0, limits = c(-ratio_cap, ratio_cap),
    name = "log2 normalized\nWT / KO",
    breaks = c(-2, -1, 0, 1, 2)) +
  scale_shape_manual(values = c(finite = 21, WT_only = 24,
    KO_only = 25, neither_detected = 4), name = "Support") +
  scale_x_continuous(limits = x_limits, expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = NULL,
    subtitle = paste("Detection-normalized WT/KO ratios; triangles indicate",
      "detection in only one genotype")) +
  common_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "right")

p_tcf3 <- ggplot(tcf3,
  aes(midpoint_mb, tcf3_best_log2fc)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_vline(xintercept = tss_mb, linetype = "dashed",
    linewidth = 0.45, color = "grey25") +
  geom_segment(aes(xend = midpoint_mb, y = 0, yend = tcf3_best_log2fc),
    linewidth = 0.8, color = "grey60") +
  geom_point(data = merge(tcf3, dynamic_class[, .(feature_id,
    dynamic_class)], by = "feature_id"),
    aes(color = dynamic_class), size = 3) +
  geom_text(aes(label = locus_label), size = 2.7, vjust = 1.7) +
  scale_color_manual(values = c("Robust Bdev dynamic" = "#7A3E9D",
    "BH-significant Bdev dynamic" = "#E69F00",
    "Not Bdev dynamic" = "grey45",
    "Outside ±2 Mb dynamic screen" = "grey75"), name = NULL) +
  scale_x_continuous(limits = x_limits, expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "TCF3 depletion\nlog2FC",
    subtitle = "All eight peaks lose accessibility after E2A/TCF3 depletion") +
  common_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "right")

p_plac <- ggplot(arc_rows,
  aes(x, y, group = interaction_id, color = minus_log10_q)) +
  geom_path(linewidth = 0.55, alpha = 0.75) +
  geom_vline(xintercept = tss_mb, linetype = "dashed",
    linewidth = 0.45, color = "grey25") +
  scale_color_viridis_c(option = "C", direction = -1,
    name = expression(-log[10](PLAC~q))) +
  scale_x_continuous(limits = x_limits, expand = expansion(mult = c(0, 0))) +
  labs(x = "chr11 position (Mb; mm10)", y = NULL,
    subtitle = paste0(nrow(links),
      " H3K4me3 PLAC-seq links anchored at the Ebf1 promoter")) +
  common_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank(), legend.position = "right")

figure <- p_bdev / p_ratio / p_tcf3 / p_plac +
  plot_layout(heights = c(2.4, 1.25, 1.35, 1.5), guides = "keep")
stem <- "ebf1_tcf3_bdev_dynamics_normalized_wtko_placseq_locus"
ggsave(file.path(stage_dir, "figures", paste0(stem, ".png")), figure,
  width = 15, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "figures", paste0(stem, ".pdf")), figure,
  width = 15, height = 12, device = cairo_pdf, bg = "white")

summary_table <- merge(tcf3[, .(feature_id, locus_label,
  peak_start_1based, peak_end_1based, midpoint_relative_to_tss_bp,
  tcf3_log2fc_depleted_vs_control = tcf3_best_log2fc,
  tcf3_bh_q = tcf3_best_padj, placseq_support = placseq_route_pass)],
  dynamic_class[, .(feature_id, dynamic_class)], by = "feature_id")
ratio_export <- dcast(ratios, feature_id ~ knockout_target,
  value.var = c("normalized_wt_over_ko_ratio", "ratio_status"))
summary_table <- merge(summary_table, ratio_export, by = "feature_id")
write_tsv(summary_table, file.path(stage_dir,
  "eight_peak_integrated_locus_summary.tsv"))
write_tsv(bdev_profiles, file.path(stage_dir,
  "eight_peak_five_region_bdev_accessibility.tsv"))
write_tsv(links, file.path(stage_dir, "placseq_links_in_displayed_locus.tsv"))

write_json(list(status = "complete_unpromoted_report_only",
  script_version = script_version, bundle_fingerprint = fingerprint_before,
  genome = "GRCm38/mm10", chromosome = chromosome,
  locus = list(start_1based = locus_start, end_1based = locus_end,
    ebf1_tss_1based = ebf1_tss), peaks = 8L,
  bdev_states = as.character(regions$region_label),
  ratio_definition = "mean[10000*I(accessible)/cell detected peaks] WT divided by KO",
  tcf3_definition = "original E2A/TCF3-depleted versus control log2FC",
  placseq_links = nrow(links), bundle_mutated = FALSE, promotion = FALSE),
  file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Integrated Ebf1 locus summary", "",
  "This figure aligns all eight significant E2A/TCF3-depletion peaks across a focused Ebf1 interval that extends to +2.1 Mb so P8 is retained.", "",
  "- Bdev bubble tracks show the accessible-cell percentage in the five ordered trajectory regions.",
  "- WT/KO tracks show the requested detected-peak-normalized ratio separately for Flt3L and IL-7; upward and downward triangles indicate WT-only and KO-only detection.",
  "- The TCF3 track shows the original depletion-versus-control log2FC. Purple marks the robust Bdev-dynamic P7 element; orange marks statistically dynamic P3.",
  "- The bottom track shows Ebf1-promoter H3K4me3 PLAC-seq arcs.", "",
  paste0("![Integrated Ebf1 locus](figures/", stem, ".png)"), "",
  paste0("- [Vector PDF](figures/", stem, ".pdf)"),
  "- [Integrated eight-peak table](eight_peak_integrated_locus_summary.tsv)",
  "- [Five-region Bdev values](eight_peak_five_region_bdev_accessibility.tsv)",
  "- [Displayed PLAC-seq links](placseq_links_in_displayed_locus.tsv)"),
  file.path(stage_dir, "REPORT.md"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed integrated Ebf1 locus summary at", output_dir, "\n")
