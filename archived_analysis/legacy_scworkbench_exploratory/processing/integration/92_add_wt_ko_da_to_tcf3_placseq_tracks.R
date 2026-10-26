#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 92_add_wt_ko_da_to_tcf3_placseq_tracks.R",
  "SOURCE_TCF3_REPORT SHARED_SUPPORT_DA_TABLE OUTPUT_DIR SCRIPT_VERSION"
), call. = FALSE)

source_report <- normalizePath(args[[1L]], mustWork = TRUE)
da_path <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
script_version <- args[[4L]]

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
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
region_start <- ebf1_tss - 3000000L
region_end <- ebf1_tss + 3000000L
focus_start <- 44000000L
focus_end <- 46000000L

coverage_path <- file.path(source_report,
  "condition_mean_full_region_1kb_coverage.tsv")
peak_path <- file.path(source_report, "eight_tcf3_peak_loci.tsv")
link_path <- file.path(source_report,
  "ebf1_promoter_placseq_links_in_evaluated_region.tsv")
stopifnot(file.exists(coverage_path), file.exists(peak_path),
  file.exists(link_path))

coverage <- fread(coverage_path)
peaks <- fread(peak_path)
links <- fread(link_path)
da <- fread(da_path)

required_da <- c("feature_id", "pct_a", "pct_b", "p_val_adj",
  "chromosome", "peak_start_1based", "peak_end_1based", "target",
  "direction", "significant_bh05")
stopifnot(all(required_da %in% names(da)), nrow(da) == 502L,
  uniqueN(da$feature_id) == 251L,
  setequal(unique(da$target), c("Flt3l", "Il7")),
  all(da$chromosome == chromosome))

da[, `:=`(
  peak_midpoint_mb = (peak_start_1based + peak_end_1based) / 2e6,
  accessible_fraction_difference = pct_a - pct_b,
  target_label = factor(target, levels = c("Flt3l", "Il7"),
    labels = c("Flt3L: WT - KO", "IL-7: WT - KO")),
  significance = fifelse(significant_bh05 & direction == "WT_higher",
    "WT higher (BH < 0.05)",
    fifelse(significant_bh05 & direction == "KO_higher",
      "KO higher (BH < 0.05)", "Not significant")))]
stopifnot(da[target == "Flt3l", sum(significant_bh05)] == 86L,
  da[target == "Il7", sum(significant_bh05)] == 26L)
setorder(da, target, peak_start_1based)
write_tsv(da, file.path(stage_dir,
  "wt_ko_da_tracks_ebf1_pm2mb.tsv"))

condition_order <- c("control", "TCF3_depleted")
coverage[, condition := factor(condition, levels = condition_order)]
condition_labels <- c(control = "Control mean (n=3)",
  TCF3_depleted = "E2A/TCF3 depleted mean (n=3)")
condition_colors <- c(control = "#D55E00", TCF3_depleted = "#3264A8")
significance_colors <- c(
  "WT higher (BH < 0.05)" = "#C84545",
  "KO higher (BH < 0.05)" = "#3B77A8",
  "Not significant" = "grey76")
peak_marks <- peaks[, .(x = peak_midpoint / 1e6,
  label = paste0("P", peak_index))]

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

coverage_panel <- function(xmin, xmax, title) {
  d <- coverage[midpoint_mb >= xmin / 1e6 & midpoint_mb <= xmax / 1e6]
  marks <- peak_marks[x >= xmin / 1e6 & x <= xmax / 1e6]
  ymax <- max(d$smoothed_coverage, na.rm = TRUE)
  ggplot(d, aes(midpoint_mb, smoothed_coverage, color = condition,
    fill = condition)) +
    geom_vline(xintercept = ebf1_tss / 1e6, color = "#222222",
      linewidth = 0.4, linetype = "dashed") +
    geom_vline(data = marks, aes(xintercept = x), inherit.aes = FALSE,
      color = "#8B1A1A", linewidth = 0.25, alpha = 0.72) +
    geom_area(alpha = 0.62, color = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.52, na.rm = TRUE) +
    geom_text(data = marks, aes(x, ymax * 1.01, label = label),
      inherit.aes = FALSE, size = 2.5, color = "#8B1A1A", vjust = 0) +
    facet_grid(condition ~ ., labeller = as_labeller(condition_labels)) +
    scale_color_manual(values = condition_colors, guide = "none") +
    scale_fill_manual(values = condition_colors, guide = "none") +
    scale_x_continuous(limits = c(xmin, xmax) / 1e6,
      expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = c(0, ymax * 1.08), clip = "off") +
    labs(x = NULL,
      y = "Mean fragment coverage CPM\n(1-kb bins; 5-kb rolling mean)",
      title = title,
      subtitle = paste(
        "Acute E2A/TCF3 degradation; equal-weight mean of three",
        "replicate-level CPM tracks per condition; mm10")) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
      strip.text.y = element_text(face = "bold"),
      plot.margin = margin(8, 8, 2, 8))
}

da_panel <- function(xmin, xmax) {
  d <- da[peak_end_1based >= xmin & peak_start_1based <= xmax]
  ggplot(d, aes(peak_midpoint_mb, accessible_fraction_difference,
    color = significance)) +
    geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
    geom_vline(xintercept = ebf1_tss / 1e6, color = "#222222",
      linewidth = 0.4, linetype = "dashed") +
    geom_segment(aes(xend = peak_midpoint_mb, y = 0,
      yend = accessible_fraction_difference), linewidth = 0.55,
      alpha = 0.85) +
    geom_point(size = 1.5, alpha = 0.95) +
    facet_grid(target_label ~ .) +
    scale_color_manual(values = significance_colors,
      breaks = names(significance_colors), name = NULL) +
    scale_x_continuous(limits = c(xmin, xmax) / 1e6,
      expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "Accessible-cell\nfraction: WT - KO",
      subtitle = paste0(
        "Shared-support DA: colored peaks pass BH < 0.05; gray peaks do not (",
        "Flt3L 86/251; IL-7 26/251)")) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
      strip.text.y = element_text(face = "bold"),
      plot.margin = margin(1, 8, 2, 8))
}

plac_panel <- function(xmin, xmax) {
  keep_links <- links[distal_end >= xmin & distal_start <= xmax]
  d <- arc_rows[interaction_id %in% keep_links$interaction_id]
  ggplot(d, aes(x, y, group = interaction_id, color = minus_log10_q)) +
    geom_path(linewidth = 0.52, alpha = 0.72) +
    geom_vline(xintercept = ebf1_tss / 1e6, color = "#222222",
      linewidth = 0.4, linetype = "dashed") +
    scale_color_viridis_c(option = "C", direction = -1,
      name = expression(-log[10](PLAC~q))) +
    scale_x_continuous(limits = c(xmin, xmax) / 1e6,
      expand = expansion(mult = c(0, 0))) +
    labs(x = "chr11 position (Mb; mm10)", y = NULL,
      subtitle = paste0(nrow(keep_links),
        " 230238 H3K4me3 PLAC-seq Ebf1-promoter links")) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      panel.grid = element_blank(), legend.position = "right",
      plot.margin = margin(0, 8, 8, 8))
}

save_composite <- function(xmin, xmax, stem, title) {
  p <- coverage_panel(xmin, xmax, title) /
    da_panel(xmin, xmax) /
    plac_panel(xmin, xmax) +
    plot_layout(heights = c(2.25, 1.75, 1.35))
  ggsave(file.path(stage_dir, "figures", paste0(stem, ".png")), p,
    width = 15, height = 10.2, dpi = 260, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(stem, ".pdf")), p,
    width = 15, height = 10.2, device = cairo_pdf, bg = "white")
}

full_stem <- "e2a_tcf3_coverage_wt_ko_da_placseq_ebf1_pm3mb"
focus_stem <- "e2a_tcf3_coverage_wt_ko_da_placseq_chr11_44_46mb"
save_composite(region_start, region_end, full_stem,
  "E2A/TCF3 + WT/KO DA + PLAC-seq: Ebf1 ±3 Mb")
save_composite(focus_start, focus_end, focus_stem,
  "E2A/TCF3 + WT/KO DA + PLAC-seq: chr11:44-46 Mb")

inputs <- data.table(
  role = c("source_coverage", "source_tcf3_peaks", "source_placseq_links",
    "shared_support_wt_ko_da"),
  path = c(coverage_path, peak_path, link_path, da_path))
inputs[, `:=`(size_bytes = file.info(path)$size,
  sha256 = vapply(path, sha256_file, character(1L)))]
write_tsv(inputs, file.path(stage_dir, "input_manifest.tsv"))

checks <- list(
  da_rows = nrow(da),
  unique_da_peaks = uniqueN(da$feature_id),
  flt3l_significant = da[target == "Flt3l", sum(significant_bh05)],
  il7_significant = da[target == "Il7", sum(significant_bh05)],
  source_tcf3_peaks = nrow(peaks),
  placseq_links = nrow(links),
  figures_exist = all(file.exists(file.path(stage_dir, "figures",
    paste0(rep(c(full_stem, focus_stem), each = 2L),
      rep(c(".png", ".pdf"), times = 2L))))),
  source_bundles_mutated = FALSE,
  promotions = 0L)
stopifnot(checks$da_rows == 502L, checks$unique_da_peaks == 251L,
  checks$flt3l_significant == 86L, checks$il7_significant == 26L,
  checks$source_tcf3_peaks == 8L, checks$placseq_links == 28L,
  checks$figures_exist, !checks$source_bundles_mutated,
  checks$promotions == 0L)
write_json(list(status = "complete_unpromoted_report_only",
  script_version = script_version, checks = checks),
  file.path(stage_dir, "diagnostics.json"))

writeLines(c(
  "# E2A/TCF3 degradation, shared-support WT/KO DA, and Ebf1 PLAC-seq", "",
  "The existing E2A/TCF3 degradation coverage and 230238 Ebf1-promoter PLAC-seq display now includes aligned Flt3L and IL-7 WT-versus-KO differential-accessibility tracks.", "",
  paste0("![Full region](figures/", full_stem, ".png)"), "",
  paste0("![Focused region](figures/", focus_stem, ".png)"), "",
  "## Differential-accessibility encoding", "",
  "- Track height is accessible-cell fraction in WT minus KO.",
  "- Red marks BH-significant WT-higher peaks; blue marks BH-significant KO-higher peaks; gray marks nonsignificant tested peaks.",
  "- Each cytokine test contains 251 Bdev peaks within Ebf1 TSS ±2 Mb and is restricted to consecutive pseudotime bins containing at least 10 WT and 10 KO cells per bin.",
  "- BH adjustment was performed independently within each cytokine comparison.",
  "- The E2A/TCF3 coverage and PLAC-seq data are unchanged from the v4 source report.", "",
  "Because each WT/KO genotype has one biological sample, the DA tracks show descriptive cell-level differences rather than replicate-level inference.", "",
  "## Downloads", "",
  paste0("- [Full-region PDF](figures/", full_stem, ".pdf)"),
  paste0("- [Focused PDF](figures/", focus_stem, ".pdf)"),
  "- [Plotted WT/KO DA values](wt_ko_da_tracks_ebf1_pm2mb.tsv)",
  "- [Input manifest](input_manifest.tsv)"
), file.path(stage_dir, "REPORT.md"))

capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
if (!file.rename(stage_dir, output_dir))
  stop("Failed to finalize output directory", call. = FALSE)
cat("Completed E2A/TCF3 + WT/KO DA + PLAC-seq report at", output_dir, "\n")
