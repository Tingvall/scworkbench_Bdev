#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 97_tcf3_peaks_detection_normalized_wt_ko_ratio.R BUNDLE",
  "SHARED_SUPPORT_REPORT COMPLEXITY_REPORT TCF3_TABLE OUTPUT_DIR",
  "EXPECTED_FINGERPRINT SCRIPT_VERSION"), call. = FALSE)

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
support_dir <- normalizePath(args[[2L]], mustWork = TRUE)
complexity_dir <- normalizePath(args[[3L]], mustWork = TRUE)
tcf3_path <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
expected_fp <- args[[6L]]
script_version <- args[[7L]]

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
})
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

assignments <- fread(cmd = paste("zcat", shQuote(file.path(support_dir,
  "shared_region_cell_assignments.tsv.gz"))))
complexity <- fread(cmd = paste("zcat", shQuote(file.path(complexity_dir,
  "cell_depth_complexity_and_match_status.tsv.gz"))))
stopifnot(!anyDuplicated(assignments$cell_id),
  !anyDuplicated(complexity$cell_id),
  all(c("cell_id", "detected_peaks") %in% names(complexity)))
cells <- merge(assignments, complexity[, .(cell_id, detected_peaks)],
  by = "cell_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(cells) == nrow(assignments), !anyNA(cells$detected_peaks),
  all(cells$detected_peaks > 0),
  cells[knockout_target == "Flt3l", range(bin)] == c(14L, 25L),
  cells[knockout_target == "Il7", range(bin)] == c(22L, 30L))
setorder(cells, knockout_target, genotype, mapped_pt, cell_id)

tcf3 <- fread(tcf3_path)[tcf3_route_pass == TRUE]
setorder(tcf3, peak_start_1based, peak_end_1based)
stopifnot(nrow(tcf3) == 8L, !anyDuplicated(tcf3$bdev_feature_id),
  all(tcf3$tcf3_best_log2fc < -1), all(tcf3$tcf3_best_padj < 0.05))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fp),
  all(cells$cell_id %in% wb_cells(bundle, "set:ko_query")),
  all(tcf3$bdev_feature_id %in%
    wb_feature_ids(bundle, "ATAC_bdev_mm10_peaks")))

peak_ids <- as.character(tcf3$bdev_feature_id)
counts <- as.matrix(wb_matrix_read_features(bundle, "ATAC_bdev_mm10_peaks",
  peak_ids, cells$cell_id, missing = "NA"))
stopifnot(identical(rownames(counts), peak_ids),
  identical(colnames(counts), cells$cell_id), !anyNA(counts))

long <- rbindlist(lapply(seq_along(peak_ids), function(i) data.table(
  cell_id = cells$cell_id, feature_id = peak_ids[[i]],
  accessible = counts[i, ] > 0,
  detected_peaks = cells$detected_peaks,
  normalized_detection_per_10k =
    10000 * as.numeric(counts[i, ] > 0) / cells$detected_peaks,
  knockout_target = cells$knockout_target, genotype = cells$genotype,
  bin = cells$bin)))

group_summary <- long[, .(
  cells = .N,
  accessible_cells = sum(accessible),
  accessible_percent = 100 * mean(accessible),
  mean_detected_peaks_per_cell = mean(detected_peaks),
  median_detected_peaks_per_cell = as.numeric(median(detected_peaks)),
  mean_detection_normalized_per_10k = mean(normalized_detection_per_10k)),
  by = .(knockout_target, feature_id, genotype)]

wide <- dcast(group_summary,
  knockout_target + feature_id ~ genotype,
  value.var = c("cells", "accessible_cells", "accessible_percent",
    "mean_detected_peaks_per_cell", "median_detected_peaks_per_cell",
    "mean_detection_normalized_per_10k"))
stopifnot(!anyNA(wide$mean_detection_normalized_per_10k_WT),
  !anyNA(wide$mean_detection_normalized_per_10k_KO))
wide[, `:=`(
  raw_accessible_percent_wt_over_ko = accessible_percent_WT /
    accessible_percent_KO,
  normalized_wt_over_ko_ratio = mean_detection_normalized_per_10k_WT /
    mean_detection_normalized_per_10k_KO,
  normalized_log2_wt_over_ko = log2(mean_detection_normalized_per_10k_WT /
    mean_detection_normalized_per_10k_KO))]
wide[, ratio_status := fcase(
  mean_detection_normalized_per_10k_WT == 0 &
    mean_detection_normalized_per_10k_KO == 0, "neither_detected",
  mean_detection_normalized_per_10k_KO == 0, "WT_only",
  mean_detection_normalized_per_10k_WT == 0, "KO_only",
  default = "finite")]
wide[, normalized_direction := fifelse(normalized_wt_over_ko_ratio > 1,
  "WT_higher", fifelse(normalized_wt_over_ko_ratio < 1, "KO_higher", "equal"))]
wide <- merge(wide, tcf3[, .(feature_id = bdev_feature_id,
  peak_start_1based, peak_end_1based, midpoint_relative_to_tss_bp,
  tcf3_log2fc_depleted_vs_control = tcf3_best_log2fc,
  tcf3_bh_q = tcf3_best_padj, placseq_support = placseq_route_pass)],
  by = "feature_id", all.x = TRUE)
setorder(wide, knockout_target, peak_start_1based)

write_tsv(group_summary, file.path(stage_dir,
  "all_group_normalized_detection_values.tsv"))
write_tsv(wide, file.path(stage_dir,
  "all_tcf3_peaks_normalized_wt_ko_ratios.tsv"))
fwrite(long[, .(cell_id, knockout_target, genotype, bin, feature_id,
  accessible, detected_peaks, normalized_detection_per_10k)],
  file.path(stage_dir, "cell_level_normalized_peak_detection.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip")

plot_dt <- copy(wide)
plot_dt[, peak_label := sprintf("%s (%+.2f Mb)", feature_id,
  midpoint_relative_to_tss_bp / 1e6)]
peak_levels <- unique(plot_dt[order(peak_start_1based), peak_label])
plot_dt[, peak_label := factor(peak_label, levels = rev(peak_levels))]
plot_dt[, target_label := factor(knockout_target,
  levels = c("Flt3l", "Il7"), labels = c("Flt3L WT / FLT3LKO",
    "IL-7 WT / IL7KO"))]
finite_limit <- max(2.25,
  max(abs(plot_dt[is.finite(normalized_log2_wt_over_ko),
    normalized_log2_wt_over_ko])))
plot_dt[, display_log2_ratio := pmax(-finite_limit,
  pmin(finite_limit, normalized_log2_wt_over_ko))]
p <- ggplot(plot_dt, aes(normalized_log2_wt_over_ko, peak_label,
  color = normalized_direction)) +
  geom_vline(xintercept = 0, color = "grey45", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = display_log2_ratio,
    yend = peak_label), linewidth = 0.7, color = "grey70") +
  geom_point(aes(x = display_log2_ratio, shape = ratio_status), size = 3) +
  facet_wrap(~target_label, nrow = 1) +
  scale_color_manual(values = c(WT_higher = "#C43C39",
    KO_higher = "#3264A8", equal = "grey45"), name = NULL) +
  scale_shape_manual(values = c(finite = 16, WT_only = 24,
    KO_only = 25, neither_detected = 4), name = "Ratio support") +
  scale_x_continuous(limits = c(-finite_limit, finite_limit),
    oob = scales::squish) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    axis.title.y = element_blank()) +
  labs(x = "log2(WT / KO) of detection-normalized accessibility",
    title = "All eight TCF3-dependent peaks in shared-support WT/KO regions",
    subtitle = paste("Per-cell score = 10,000 × accessible indicator /",
      "number of detected Bdev-aligned peaks; group means are ratioed"))
ggsave(file.path(stage_dir, "figures",
  "all_tcf3_peaks_normalized_wt_ko_ratio.png"), p,
  width = 12, height = 6.5, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "all_tcf3_peaks_normalized_wt_ko_ratio.pdf"), p,
  width = 12, height = 6.5, device = cairo_pdf, bg = "white")

target_summary <- wide[, .(
  peaks = .N, WT_higher = sum(normalized_direction == "WT_higher"),
  KO_higher = sum(normalized_direction == "KO_higher"),
  finite_ratios = sum(ratio_status == "finite"),
  WT_only = sum(ratio_status == "WT_only"),
  KO_only = sum(ratio_status == "KO_only"),
  median_finite_wt_over_ko = median(
    normalized_wt_over_ko_ratio[ratio_status == "finite"]),
  minimum_wt_over_ko = min(normalized_wt_over_ko_ratio),
  maximum_wt_over_ko = max(normalized_wt_over_ko_ratio)),
  by = knockout_target]
write_tsv(target_summary, file.path(stage_dir, "target_summary.tsv"))
cell_counts <- cells[, .N, by = .(knockout_target, genotype)]
write_json(list(status = "complete_unpromoted_report_only",
  script_version = script_version, bundle_fingerprint = fingerprint_before,
  peaks = length(peak_ids),
  cells = split(cell_counts, cell_counts$knockout_target),
  regions = list(Flt3l = "bins 14-25", Il7 = "bins 22-30"),
  cell_score = "10000 * I(peak count > 0) / detected peaks among all 97,788 Bdev-aligned peaks",
  comparison = "ratio of WT and KO group mean cell scores; no pseudocount",
  target_summary = split(target_summary, target_summary$knockout_target),
  inference = "descriptive; one biological sample per genotype",
  bundle_mutated = FALSE, promotion = FALSE),
  file.path(stage_dir, "diagnostics.json"))

summary_lines <- apply(target_summary, 1L, function(x) paste0("- `",
  x[["knockout_target"]], "`: ", x[["WT_higher"]],
  " WT-higher and ", x[["KO_higher"]],
  " KO-higher peaks; median finite WT/KO = ",
  sprintf("%.3f", as.numeric(x[["median_finite_wt_over_ko"]])),
  "; WT-only = ", x[["WT_only"]], "; KO-only = ", x[["KO_only"]], "."))
writeLines(c("# Detection-normalized WT/KO ratios for TCF3-dependent peaks", "",
  "All eight peaks passing the existing E2A/TCF3-depletion differential-accessibility criterion were evaluated in the same shared-support mapped-cell regions used for the Ebf1 analysis: Flt3L bins 14-25 and IL-7 bins 22-30.", "",
  "For each cell and peak, the score is `10,000 × I(peak count > 0) / number of detected peaks across all 97,788 Bdev-aligned peaks`. The reported WT/KO value is the WT group mean divided by the KO group mean, without a pseudocount. The factor 10,000 improves readability and cancels in the ratio.", "",
  summary_lines, "",
  "These ratios are descriptive because there is one biological sample per genotype.", "",
  "- [Complete ratio table](all_tcf3_peaks_normalized_wt_ko_ratios.tsv)",
  "- [Group summaries](all_group_normalized_detection_values.tsv)",
  "- [Cell-level values](cell_level_normalized_peak_detection.tsv.gz)",
  "- [Ratio figure](figures/all_tcf3_peaks_normalized_wt_ko_ratio.png)"),
  file.path(stage_dir, "REPORT.md"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed TCF3 normalized WT/KO ratio report at", output_dir, "\n")
