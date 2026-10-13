#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(paste(
  "Usage: 36_tcf3_ebf1_coverage_placseq_tracks.R SAMPLE_MANIFEST BAM_DIR",
  "QC_DIR THREE_WAY_REPORT OUTPUT_DIR SCRIPT_VERSION EBF1_TSS_1BASED SEED"
), call. = FALSE)

sample_manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
bam_dir <- normalizePath(args[[2L]], mustWork = TRUE)
qc_dir <- normalizePath(args[[3L]], mustWork = TRUE)
three_way_report <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
script_version <- args[[6L]]
ebf1_tss <- as.integer(args[[7L]])
seed <- as.integer(args[[8L]])
stopifnot(ebf1_tss == 44617317L, is.finite(seed))

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)
dir.create(file.path(stage_dir, "tracks"), recursive = TRUE,
  showWarnings = FALSE)
set.seed(seed)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(Rsamtools)
  library(GenomicRanges)
  library(IRanges)
  library(jsonlite)
  library(digest)
})

log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    paste(..., collapse = " "), "\n")
  flush.console()
}
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
roll_mean <- function(x, k) as.numeric(stats::filter(x, rep(1 / k, k),
  sides = 2, circular = FALSE))

chromosome <- "chr11"
region_flank_bp <- 3000000L
region_start <- ebf1_tss - region_flank_bp
region_end <- ebf1_tss + region_flank_bp
region_length <- region_end - region_start + 1L
fine_bin_bp <- 50L
full_bin_bp <- 1000L
full_smoothing_bins <- 5L
zoom_flank_bp <- 10000L
zoom_smoothing_bins <- 5L
max_insert_bp <- 2000L

manifest <- fread(sample_manifest_path)
required_manifest <- c("sample_id", "condition", "treatment_time_h",
  "geo_accession", "sra_run")
stopifnot(all(required_manifest %in% names(manifest)), nrow(manifest) == 6L,
  !anyDuplicated(manifest$sample_id),
  identical(manifest[, .N, by = condition][order(condition)]$N, c(3L, 3L)))
condition_order <- c("control", "TCF3_depleted")
stopifnot(setequal(manifest$condition, condition_order))
manifest[, `:=`(
  bam_path = file.path(bam_dir, paste0(sample_id, ".filtered.bam")),
  bai_path = file.path(bam_dir, paste0(sample_id, ".filtered.bam.bai")),
  flagstat_path = file.path(qc_dir, paste0(sample_id, ".flagstat.txt")))]
stopifnot(all(file.exists(manifest$bam_path)), all(file.exists(manifest$bai_path)),
  all(file.exists(manifest$flagstat_path)))

parse_read1_count <- function(path) {
  line <- grep("^[0-9]+ \\+ [0-9]+ read1$", readLines(path), value = TRUE)
  stopifnot(length(line) == 1L)
  as.numeric(sub(" .*", "", line))
}
manifest[, normalization_fragments := vapply(flagstat_path,
  parse_read1_count, numeric(1L))]
stopifnot(all(manifest$normalization_fragments > 1e6))

master_path <- file.path(three_way_report,
  "all_259_bdev_ebf1_window_peaks_three_way_tests.tsv")
links_path <- file.path(three_way_report,
  "route3_all_ebf1_promoter_placseq_interactions.tsv")
stopifnot(file.exists(master_path), file.exists(links_path))
master <- fread(master_path)
peaks <- master[tcf3_route_pass == TRUE]
setorder(peaks, peak_start_1based, peak_end_1based)
stopifnot(nrow(peaks) == 8L, all(peaks$chromosome == chromosome),
  all(peaks$tcf3_best_log2fc < -1), all(peaks$tcf3_best_padj < 0.05))
peaks[, `:=`(
  peak_index = seq_len(.N),
  peak_midpoint = (peak_start_1based + peak_end_1based) / 2,
  zoom_start = pmax(region_start, floor((peak_start_1based +
    peak_end_1based) / 2) - zoom_flank_bp),
  zoom_end = pmin(region_end, floor((peak_start_1based +
    peak_end_1based) / 2) + zoom_flank_bp))]
peaks[, panel_label := sprintf(
  "P%d  %s\nTCF3 log2FC %.2f; padj %.2g%s", peak_index,
  bdev_feature_id, tcf3_best_log2fc, tcf3_best_padj,
  fifelse(placseq_route_pass, "; PLAC-seq+", ""))]
write_tsv(peaks, file.path(stage_dir, "eight_tcf3_peak_loci.tsv"))

links <- fread(links_path)
required_links <- c("interaction_id", "promoter_start", "promoter_end",
  "distal_chr", "distal_start", "distal_end", "Q-Value_Bias", "cc")
stopifnot(all(required_links %in% names(links)), nrow(links) == 28L)
links <- links[distal_chr == chromosome & distal_end >= region_start &
  distal_start <= region_end]
links[, `:=`(
  promoter_midpoint = (promoter_start + promoter_end) / 2,
  distal_midpoint = (distal_start + distal_end) / 2,
  minus_log10_q = -log10(pmax(`Q-Value_Bias`, 1e-300)))]
write_tsv(links, file.path(stage_dir,
  "ebf1_promoter_placseq_links_in_evaluated_region.tsv"))

# Derive equal-width regional coverage from one record per paired fragment.
# The first-mate record and absolute TLEN reconstruct the fragment span. The
# same first-mate count in the filtered BAM is the replicate CPM denominator.
fine_n <- ceiling(region_length / fine_bin_bp)
fine_index <- ((seq_len(region_length) - 1L) %/% fine_bin_bp) + 1L
fine_widths <- tabulate(fine_index, nbins = fine_n)
fine_starts <- region_start + (seq_len(fine_n) - 1L) * fine_bin_bp
fine_ends <- pmin(region_end, fine_starts + fine_bin_bp - 1L)
sample_tracks <- vector("list", nrow(manifest))
sample_audit <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  sample <- manifest[i]
  log_event("Reading regional fragments for", sample$sample_id)
  query <- GRanges(chromosome, IRanges(max(1L, region_start -
    max_insert_bp), region_end + max_insert_bp))
  param <- ScanBamParam(which = query,
    what = c("pos", "mpos", "isize", "flag"))
  records <- scanBam(sample$bam_path, param = param)[[1L]]
  keep <- bitwAnd(as.integer(records$flag), 64L) != 0L &
    is.finite(records$pos) & is.finite(records$mpos) &
    is.finite(records$isize) & abs(records$isize) > 0L &
    abs(records$isize) <= max_insert_bp
  fragment_start <- pmin(as.integer(records$pos[keep]),
    as.integer(records$mpos[keep]))
  fragment_end <- fragment_start + abs(as.integer(records$isize[keep])) - 1L
  overlap <- fragment_end >= region_start & fragment_start <= region_end
  fragment_start <- pmax(fragment_start[overlap], region_start)
  fragment_end <- pmin(fragment_end[overlap], region_end)
  stopifnot(length(fragment_start) > 0L,
    length(fragment_start) == length(fragment_end),
    all(fragment_start <= fragment_end))

  start_idx <- fragment_start - region_start + 1L
  end_idx <- fragment_end - region_start + 1L
  delta <- as.numeric(tabulate(start_idx, nbins = region_length + 1L)) -
    as.numeric(tabulate(end_idx + 1L, nbins = region_length + 1L))
  base_coverage <- cumsum(delta)[seq_len(region_length)] /
    sample$normalization_fragments * 1e6
  fine_value <- as.numeric(rowsum(matrix(base_coverage, ncol = 1L),
    fine_index, reorder = FALSE)[, 1L]) / fine_widths
  sample_tracks[[i]] <- data.table(
    chromosome = chromosome,
    start_1based = fine_starts,
    end_1based = fine_ends,
    sample_id = sample$sample_id,
    condition = sample$condition,
    treatment_time_h = sample$treatment_time_h,
    normalized_fragment_coverage = fine_value)
  sample_audit[[i]] <- data.table(
    sample_id = sample$sample_id,
    condition = sample$condition,
    treatment_time_h = sample$treatment_time_h,
    normalization_fragments = sample$normalization_fragments,
    regional_first_mate_records = sum(keep),
    regional_overlapping_fragments = length(fragment_start),
    maximum_regional_fragment_length = max(fragment_end - fragment_start + 1L),
    mean_regional_fragment_length = mean(fragment_end - fragment_start + 1L))
  rm(records, fragment_start, fragment_end, delta, base_coverage)
  invisible(gc())
}

sample_tracks <- rbindlist(sample_tracks)
sample_audit <- rbindlist(sample_audit)
stopifnot(nrow(sample_tracks) == 6L * fine_n,
  sample_tracks[, uniqueN(sample_id), by = .(condition, start_1based)][, all(V1 == 3L)])
write_tsv(sample_audit, file.path(stage_dir,
  "replicate_cpm_normalization_and_fragment_audit.tsv"))
fwrite(sample_tracks, file.path(stage_dir,
  "replicate_normalized_50bp_coverage.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")

group_tracks <- sample_tracks[, .(
  normalized_fragment_coverage = mean(normalized_fragment_coverage)),
  by = .(chromosome, start_1based, end_1based, condition)]
group_tracks[, condition := factor(condition, levels = condition_order)]
stopifnot(nrow(group_tracks) == 2L * fine_n)
fwrite(group_tracks, file.path(stage_dir,
  "condition_mean_normalized_50bp_coverage.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")
for (condition_value in condition_order) {
  out <- group_tracks[condition == condition_value,
    .(chromosome, start_0based = start_1based - 1L,
      end_0based_exclusive = end_1based,
      normalized_fragment_coverage)]
  fwrite(out, file.path(stage_dir, "tracks", paste0(condition_value,
    "_mean_fragment_CPM_50bp.bedGraph.gz")), sep = "\t", col.names = FALSE,
    compress = "gzip")
}

# Full ±3 Mb plot: aggregate 50-bp values to 1-kb means and display a centered
# 5-kb rolling mean. The condition tracks remain equal-weight replicate means.
group_tracks[, full_bin := (start_1based - region_start) %/% full_bin_bp + 1L]
full_tracks <- group_tracks[, .(
  start_1based = min(start_1based), end_1based = max(end_1based),
  normalized_fragment_coverage = weighted.mean(
    normalized_fragment_coverage, end_1based - start_1based + 1L)),
  by = .(condition, full_bin)]
setorder(full_tracks, condition, full_bin)
full_tracks[, smoothed_coverage := roll_mean(normalized_fragment_coverage,
  full_smoothing_bins), by = condition]
full_tracks[, midpoint_mb := ((start_1based + end_1based) / 2) / 1e6]
write_tsv(full_tracks, file.path(stage_dir,
  "condition_mean_full_region_1kb_coverage.tsv"))

condition_labels <- c(control = "Control mean (n=3)",
  TCF3_depleted = "TCF3 depleted mean (n=3)")
condition_colors <- c(control = "#D55E00", TCF3_depleted = "#3264A8")
peak_marks <- peaks[, .(x = peak_midpoint / 1e6,
  label = paste0("P", peak_index))]
y_top <- max(full_tracks$smoothed_coverage, na.rm = TRUE)
coverage_plot <- ggplot(full_tracks,
  aes(midpoint_mb, smoothed_coverage, color = condition,
    fill = condition)) +
  geom_vline(xintercept = ebf1_tss / 1e6, color = "#222222",
    linewidth = 0.45, linetype = "dashed") +
  geom_vline(data = peak_marks, aes(xintercept = x), inherit.aes = FALSE,
    color = "#8B1A1A", linewidth = 0.28, alpha = 0.8) +
  geom_area(alpha = 0.28, position = "identity", color = NA,
    na.rm = TRUE) +
  geom_line(linewidth = 0.48, na.rm = TRUE) +
  geom_text(data = peak_marks, aes(x, y_top * 1.01, label = label),
    inherit.aes = FALSE, size = 2.7, color = "#8B1A1A", vjust = 0) +
  scale_color_manual(values = condition_colors, labels = condition_labels,
    breaks = condition_order, name = NULL) +
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_x_continuous(limits = c(region_start, region_end) / 1e6,
    expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, y_top * 1.08), clip = "off") +
  labs(x = NULL,
    y = "Mean fragment coverage CPM\n(1-kb bins; 5-kb rolling mean)",
    title = "TCF3-dependent chromatin accessibility across the Ebf1 ±3 Mb region",
    subtitle = "Equal-weight mean of three replicate-level CPM tracks per condition; mm10") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 3, 8))

arc_rows <- rbindlist(lapply(seq_len(nrow(links)), function(i) {
  x1 <- links$promoter_midpoint[[i]] / 1e6
  x2 <- links$distal_midpoint[[i]] / 1e6
  t <- seq(0, 1, length.out = 80L)
  span <- abs(x2 - x1)
  data.table(interaction_id = links$interaction_id[[i]],
    x = x1 + (x2 - x1) * t,
    y = -sin(pi * t) * sqrt(span),
    minus_log10_q = links$minus_log10_q[[i]])
}))
link_plot <- ggplot(arc_rows, aes(x, y, group = interaction_id,
  color = minus_log10_q)) +
  geom_path(linewidth = 0.52, alpha = 0.72) +
  geom_vline(xintercept = ebf1_tss / 1e6, color = "#222222",
    linewidth = 0.45, linetype = "dashed") +
  scale_color_viridis_c(option = "C", direction = -1,
    name = expression(-log[10](PLAC~q))) +
  scale_x_continuous(limits = c(region_start, region_end) / 1e6,
    expand = expansion(mult = c(0, 0))) +
  labs(x = "chr11 position (Mb; mm10)", y = NULL,
    subtitle = paste0(nrow(links),
      " 230238 H3K4me3 PLAC-seq Ebf1-promoter links")) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank(), legend.position = "right",
    plot.margin = margin(0, 8, 8, 8))

full_figure <- coverage_plot / link_plot +
  plot_layout(heights = c(1.6, 1.35))
ggsave(file.path(stage_dir, "figures",
  "tcf3_mean_coverage_ebf1_pm3mb_with_230238_links.png"), full_figure,
  width = 15, height = 6.2, dpi = 260, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "tcf3_mean_coverage_ebf1_pm3mb_with_230238_links.pdf"), full_figure,
  width = 15, height = 6.2, device = cairo_pdf, bg = "white")

# Focused view requested over chr11:44-46 Mb, retaining the same normalization,
# binning, smoothing, colors, and PLAC-seq link encoding as the ±3 Mb view.
focus_start <- 44000000L
focus_end <- 46000000L
focus_tracks <- full_tracks[midpoint_mb >= focus_start / 1e6 &
  midpoint_mb <= focus_end / 1e6]
focus_peaks <- peak_marks[x >= focus_start / 1e6 & x <= focus_end / 1e6]
focus_links <- links[distal_end >= focus_start & distal_start <= focus_end]
focus_arc_rows <- arc_rows[interaction_id %in% focus_links$interaction_id]
focus_distal_rects <- focus_links[, .(
  xmin = pmax(distal_start, focus_start) / 1e6,
  xmax = pmin(distal_end, focus_end) / 1e6)]
focus_y_top <- max(focus_tracks$smoothed_coverage, na.rm = TRUE)
focus_coverage_plot <- ggplot(focus_tracks,
  aes(midpoint_mb, smoothed_coverage, color = condition,
    fill = condition)) +
  geom_vline(xintercept = ebf1_tss / 1e6, color = "#222222",
    linewidth = 0.45, linetype = "dashed") +
  geom_vline(data = focus_peaks, aes(xintercept = x), inherit.aes = FALSE,
    color = "#8B1A1A", linewidth = 0.3, alpha = 0.8) +
  geom_area(alpha = 0.28, position = "identity", color = NA,
    na.rm = TRUE) +
  geom_line(linewidth = 0.52, na.rm = TRUE) +
  geom_text(data = focus_peaks,
    aes(x, focus_y_top * 1.01, label = label), inherit.aes = FALSE,
    size = 2.8, color = "#8B1A1A", vjust = 0) +
  scale_color_manual(values = condition_colors, labels = condition_labels,
    breaks = condition_order, name = NULL) +
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_x_continuous(limits = c(focus_start, focus_end) / 1e6,
    expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, focus_y_top * 1.08), clip = "off") +
  labs(x = NULL,
    y = "Mean fragment coverage CPM\n(1-kb bins; 5-kb rolling mean)",
    title = "TCF3-dependent chromatin accessibility across chr11:44-46 Mb",
    subtitle = "Equal-weight mean of three replicate-level CPM tracks per condition; mm10") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 3, 8))

focus_link_plot <- ggplot(focus_arc_rows,
  aes(x, y, group = interaction_id, color = minus_log10_q)) +
  geom_path(linewidth = 0.58, alpha = 0.75) +
  geom_vline(xintercept = ebf1_tss / 1e6, color = "#222222",
    linewidth = 0.45, linetype = "dashed") +
  scale_color_viridis_c(option = "C", direction = -1,
    name = expression(-log[10](PLAC~q))) +
  scale_x_continuous(limits = c(focus_start, focus_end) / 1e6,
    expand = expansion(mult = c(0, 0))) +
  labs(x = "chr11 position (Mb; mm10)", y = NULL,
    subtitle = paste0(nrow(focus_links),
      " 230238 H3K4me3 PLAC-seq Ebf1-promoter links in the focused region")) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank(), legend.position = "right",
    plot.margin = margin(0, 8, 8, 8))

focus_figure <- focus_coverage_plot / focus_link_plot +
  plot_layout(heights = c(1.6, 1.35))
ggsave(file.path(stage_dir, "figures",
  "tcf3_mean_coverage_chr11_44_46mb_with_230238_links.png"), focus_figure,
  width = 15, height = 6.2, dpi = 260, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "tcf3_mean_coverage_chr11_44_46mb_with_230238_links.pdf"), focus_figure,
  width = 15, height = 6.2, device = cairo_pdf, bg = "white")

# Additional layouts place control and depleted means in separate stacked
# tracks with identical y limits, followed by the downward-facing link track.
separate_full_coverage <- ggplot(full_tracks,
  aes(midpoint_mb, smoothed_coverage, color = condition,
    fill = condition)) +
  geom_vline(xintercept = ebf1_tss / 1e6, color = "#222222",
    linewidth = 0.4, linetype = "dashed") +
  geom_vline(data = peak_marks, aes(xintercept = x), inherit.aes = FALSE,
    color = "#8B1A1A", linewidth = 0.25, alpha = 0.75) +
  geom_area(alpha = 0.62, color = NA, na.rm = TRUE) +
  geom_line(linewidth = 0.52, na.rm = TRUE) +
  geom_text(data = peak_marks, aes(x, y_top * 1.01, label = label),
    inherit.aes = FALSE, size = 2.5, color = "#8B1A1A", vjust = 0) +
  facet_grid(condition ~ ., labeller = as_labeller(condition_labels)) +
  scale_color_manual(values = condition_colors, guide = "none") +
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_x_continuous(limits = c(region_start, region_end) / 1e6,
    expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, y_top * 1.08), clip = "off") +
  labs(x = NULL,
    y = "Mean fragment coverage CPM\n(1-kb bins; 5-kb rolling mean)",
    title = "Separate control and TCF3-depleted coverage across Ebf1 ±3 Mb",
    subtitle = "Equal-weight mean of three replicate-level CPM tracks per condition; shared y-axis; mm10") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
    strip.text.y = element_text(face = "bold"),
    plot.margin = margin(8, 8, 3, 8))

separate_full_figure <- separate_full_coverage / link_plot +
  plot_layout(heights = c(2.4, 1.35))
ggsave(file.path(stage_dir, "figures",
  "tcf3_separate_coverage_ebf1_pm3mb_with_230238_links.png"),
  separate_full_figure, width = 15, height = 7.2, dpi = 260, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "tcf3_separate_coverage_ebf1_pm3mb_with_230238_links.pdf"),
  separate_full_figure, width = 15, height = 7.2, device = cairo_pdf,
  bg = "white")

separate_focus_coverage <- ggplot(focus_tracks,
  aes(midpoint_mb, smoothed_coverage, color = condition,
    fill = condition)) +
  geom_vline(xintercept = ebf1_tss / 1e6, color = "#222222",
    linewidth = 0.4, linetype = "dashed") +
  geom_vline(data = focus_peaks, aes(xintercept = x), inherit.aes = FALSE,
    color = "#8B1A1A", linewidth = 0.27, alpha = 0.75) +
  geom_area(alpha = 0.62, color = NA, na.rm = TRUE) +
  geom_line(linewidth = 0.56, na.rm = TRUE) +
  geom_text(data = focus_peaks,
    aes(x, focus_y_top * 1.01, label = label), inherit.aes = FALSE,
    size = 2.7, color = "#8B1A1A", vjust = 0) +
  facet_grid(condition ~ ., labeller = as_labeller(condition_labels)) +
  scale_color_manual(values = condition_colors, guide = "none") +
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_x_continuous(limits = c(focus_start, focus_end) / 1e6,
    expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, focus_y_top * 1.08), clip = "off") +
  labs(x = NULL,
    y = "Mean fragment coverage CPM\n(1-kb bins; 5-kb rolling mean)",
    title = "Separate control and TCF3-depleted coverage across chr11:44-46 Mb",
    subtitle = "Equal-weight mean of three replicate-level CPM tracks per condition; shared y-axis; mm10") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
    strip.text.y = element_text(face = "bold"),
    plot.margin = margin(8, 8, 3, 8))

separate_focus_figure <- separate_focus_coverage / focus_link_plot +
  plot_layout(heights = c(2.4, 1.35))
ggsave(file.path(stage_dir, "figures",
  "tcf3_separate_coverage_chr11_44_46mb_with_230238_links.png"),
  separate_focus_figure, width = 15, height = 7.2, dpi = 260, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "tcf3_separate_coverage_chr11_44_46mb_with_230238_links.pdf"),
  separate_focus_figure, width = 15, height = 7.2, device = cairo_pdf,
  bg = "white")

# Peak-centered views use the native 50-bp regional bins and a 250-bp rolling
# mean. Both conditions share a y scale within a peak, while peaks are free.
zoom_list <- vector("list", nrow(peaks))
rectangles <- peaks[, .(panel_label, xmin = -(
  peak_midpoint - peak_start_1based) / 1000,
  xmax = (peak_end_1based - peak_midpoint) / 1000)]
for (i in seq_len(nrow(peaks))) {
  peak <- peaks[i]
  z <- copy(group_tracks[start_1based <= peak$zoom_end &
    end_1based >= peak$zoom_start])
  setorder(z, condition, start_1based)
  z[, smoothed_coverage := roll_mean(normalized_fragment_coverage,
    zoom_smoothing_bins), by = condition]
  z[, `:=`(
    peak_index = peak$peak_index,
    panel_label = peak$panel_label,
    relative_kb = ((start_1based + end_1based) / 2 -
      peak$peak_midpoint) / 1000)]
  zoom_list[[i]] <- z
}
zoom_tracks <- rbindlist(zoom_list)
zoom_tracks[, panel_label := factor(panel_label, levels = peaks$panel_label)]
rectangles[, panel_label := factor(panel_label, levels = peaks$panel_label)]
fwrite(zoom_tracks, file.path(stage_dir,
  "eight_peak_zoom_condition_mean_50bp_coverage.tsv.gz"), sep = "\t",
  na = "NA", compress = "gzip")

zoom_plot <- ggplot(zoom_tracks,
  aes(relative_kb, smoothed_coverage, color = condition)) +
  geom_rect(data = rectangles,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#BDBDBD", alpha = 0.42) +
  geom_line(linewidth = 0.58, na.rm = TRUE) +
  facet_wrap(~panel_label, ncol = 2, scales = "free_y") +
  scale_color_manual(values = condition_colors, labels = condition_labels,
    breaks = condition_order, name = NULL) +
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10),
    limits = c(-10, 10), expand = expansion(mult = c(0, 0))) +
  labs(x = "Position relative to Bdev peak midpoint (kb; mm10)",
    y = "Mean fragment coverage CPM\n(50-bp bins; 250-bp rolling mean)",
    title = "Normalized TCF3 ATAC coverage at the eight passing peaks",
    subtitle = paste(
      "Gray shading marks the Bdev peak; PLAC-seq+ uses the ±20-kb",
      "distal-anchor overlap rule")) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8, face = "bold", lineheight = 0.95),
    axis.title = element_text(size = 9))
ggsave(file.path(stage_dir, "figures",
  "tcf3_mean_coverage_eight_passing_peaks_zoom.png"), zoom_plot,
  width = 13, height = 14, dpi = 260, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "tcf3_mean_coverage_eight_passing_peaks_zoom.pdf"), zoom_plot,
  width = 13, height = 14, device = cairo_pdf, bg = "white")

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
inputs <- rbind(
  data.table(role = c("sample_manifest", "three_way_master",
    "placseq_promoter_links", "analysis_script"),
    path = c(sample_manifest_path, master_path, links_path, script_path)),
  manifest[, .(role = paste0("filtered_bam_", sample_id), path = bam_path)],
  manifest[, .(role = paste0("filtered_bam_index_", sample_id), path = bai_path)],
  manifest[, .(role = paste0("flagstat_", sample_id), path = flagstat_path)])
inputs[, `:=`(size_bytes = file.info(path)$size,
  sha256 = vapply(path, sha256_file, character(1L)))]
write_tsv(inputs, file.path(stage_dir, "input_manifest.tsv"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

checks <- list(
  six_samples = nrow(manifest) == 6L,
  three_replicates_per_condition = all(manifest[, .N,
    by = condition]$N == 3L),
  equal_weight_replicate_mean = TRUE,
  eight_tcf3_passing_peaks = nrow(peaks) == 8L,
  four_tcf3_placseq_peaks = sum(peaks$placseq_route_pass) == 4L,
  placseq_promoter_links_in_region = nrow(links),
  full_track_grid_complete = nrow(group_tracks) == 2L * fine_n,
  all_regional_fragments_present = all(
    sample_audit$regional_overlapping_fragments > 0),
  figures_exist = all(file.exists(file.path(stage_dir, "figures", c(
    "tcf3_mean_coverage_ebf1_pm3mb_with_230238_links.png",
    "tcf3_mean_coverage_chr11_44_46mb_with_230238_links.png",
    "tcf3_separate_coverage_ebf1_pm3mb_with_230238_links.png",
    "tcf3_separate_coverage_chr11_44_46mb_with_230238_links.png",
    "tcf3_mean_coverage_eight_passing_peaks_zoom.png")))),
  outside_scworkbench_approved = TRUE,
  promotions = 0L)
stopifnot(all(unlist(checks[c("six_samples", "three_replicates_per_condition",
  "equal_weight_replicate_mean", "eight_tcf3_passing_peaks",
  "four_tcf3_placseq_peaks", "full_track_grid_complete",
  "all_regional_fragments_present", "figures_exist",
  "outside_scworkbench_approved")])), checks$promotions == 0L)
write_json(list(status = "passed", checks = checks),
  file.path(stage_dir, "checks.json"))

writeLines(c(
  "# TCF3 ATAC coverage across the Ebf1 enhancer region", "",
  "This report averages three replicate-level normalized tracks for control ",
    "and three for TCF3 depletion. The depleted mean contains the two 6-hour ",
    "and one 8-hour libraries specified in the source manifest.", "",
  "![Full region](figures/tcf3_mean_coverage_ebf1_pm3mb_with_230238_links.png)", "",
  "## Full region with separate condition tracks", "",
  "![Separate full-region tracks](figures/tcf3_separate_coverage_ebf1_pm3mb_with_230238_links.png)", "",
  "## Focused chr11:44-46 Mb view", "",
  "![Focused region](figures/tcf3_mean_coverage_chr11_44_46mb_with_230238_links.png)", "",
  "### Focused view with separate condition tracks", "",
  "![Separate focused tracks](figures/tcf3_separate_coverage_chr11_44_46mb_with_230238_links.png)", "",
  paste0("The focused panel contains ", nrow(focus_links),
    " Ebf1-promoter PLAC-seq links whose distal anchor intersects chr11:44-46 Mb."), "",
  "## Eight passing peaks", "",
  "![Eight peaks](figures/tcf3_mean_coverage_eight_passing_peaks_zoom.png)", "",
  "## Method", "",
  "- Genome and interval: mm10 chr11:41,617,317-47,617,317 (Ebf1 TSS ±3 Mb).",
  "- Fragment definition: the first-mate BAM record and absolute template length reconstruct one paired-end fragment span; inserts are limited to the 2-kb alignment maximum.",
  "- Normalization: fragment-span coverage divided by each filtered BAM's first-mate fragment count and multiplied by one million. Three normalized replicates are averaged with equal weight within each condition.",
  "- Full-region display: 1-kb means with a centered 5-kb rolling mean. The lower panel contains the 230238 H3K4me3 PLAC-seq links whose promoter anchor contains the Ebf1 TSS.",
  "- Zoom display: ±10 kb around each of the eight passing Bdev peaks, using 50-bp means and a centered 250-bp rolling mean.",
  "- This explicitly approved report-only calculation is outside scWorkbench because no registered coverage-track operation is available. No bundle was modified or promoted.", "",
  "## Downloads", "",
  "- [Full-region vector PDF](figures/tcf3_mean_coverage_ebf1_pm3mb_with_230238_links.pdf)",
  "- [Full-region separate-track PDF](figures/tcf3_separate_coverage_ebf1_pm3mb_with_230238_links.pdf)",
  "- [Focused chr11:44-46 Mb vector PDF](figures/tcf3_mean_coverage_chr11_44_46mb_with_230238_links.pdf)",
  "- [Focused separate-track PDF](figures/tcf3_separate_coverage_chr11_44_46mb_with_230238_links.pdf)",
  "- [Eight-peak vector PDF](figures/tcf3_mean_coverage_eight_passing_peaks_zoom.pdf)",
  "- [Condition mean 50-bp values](condition_mean_normalized_50bp_coverage.tsv.gz)",
  "- [Replicate-level 50-bp values](replicate_normalized_50bp_coverage.tsv.gz)",
  "- [Control mean bedGraph](tracks/control_mean_fragment_CPM_50bp.bedGraph.gz)",
  "- [TCF3-depleted mean bedGraph](tracks/TCF3_depleted_mean_fragment_CPM_50bp.bedGraph.gz)",
  "- [Normalization and fragment audit](replicate_cpm_normalization_and_fragment_audit.tsv)"
), file.path(stage_dir, "REPORT.md"))

html <- c(
  "<!doctype html><html><head><meta charset='utf-8'>",
  "<title>TCF3 ATAC coverage across Ebf1</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1500px;margin:2rem auto;padding:0 1rem;line-height:1.45}img{max-width:100%;height:auto}code{background:#f2f2f2;padding:.1rem .25rem}</style></head><body>",
  "<h1>TCF3 ATAC coverage across the Ebf1 enhancer region</h1>",
  "<p>Each condition is the equal-weight mean of three replicate-level fragment-CPM tracks. The depleted group contains two 6-hour and one 8-hour libraries.</p>",
  "<h2>Ebf1 ±3 Mb and 230238 PLAC-seq links</h2>",
  "<img src='figures/tcf3_mean_coverage_ebf1_pm3mb_with_230238_links.png'>",
  "<h2>Full region: separate control and depleted tracks</h2>",
  "<img src='figures/tcf3_separate_coverage_ebf1_pm3mb_with_230238_links.png'>",
  "<h2>Focused chr11:44-46 Mb view</h2>",
  paste0("<p>", nrow(focus_links), " Ebf1-promoter links have distal anchors intersecting this interval.</p>"),
  "<img src='figures/tcf3_mean_coverage_chr11_44_46mb_with_230238_links.png'>",
  "<h2>Focused view: separate control and depleted tracks</h2>",
  "<img src='figures/tcf3_separate_coverage_chr11_44_46mb_with_230238_links.png'>",
  "<h2>Eight passing peaks</h2>",
  "<img src='figures/tcf3_mean_coverage_eight_passing_peaks_zoom.png'>",
  "<h2>Downloads</h2><ul>",
  "<li><a href='figures/tcf3_mean_coverage_ebf1_pm3mb_with_230238_links.pdf'>Full-region PDF</a></li>",
  "<li><a href='figures/tcf3_separate_coverage_ebf1_pm3mb_with_230238_links.pdf'>Full-region separate-track PDF</a></li>",
  "<li><a href='figures/tcf3_mean_coverage_chr11_44_46mb_with_230238_links.pdf'>Focused 44-46 Mb PDF</a></li>",
  "<li><a href='figures/tcf3_separate_coverage_chr11_44_46mb_with_230238_links.pdf'>Focused separate-track PDF</a></li>",
  "<li><a href='figures/tcf3_mean_coverage_eight_passing_peaks_zoom.pdf'>Eight-peak PDF</a></li>",
  "<li><a href='condition_mean_normalized_50bp_coverage.tsv.gz'>Condition mean values</a></li>",
  "<li><a href='replicate_normalized_50bp_coverage.tsv.gz'>Replicate values</a></li>",
  "<li><a href='tracks/control_mean_fragment_CPM_50bp.bedGraph.gz'>Control bedGraph</a></li>",
  "<li><a href='tracks/TCF3_depleted_mean_fragment_CPM_50bp.bedGraph.gz'>TCF3-depleted bedGraph</a></li>",
  "</ul></body></html>")
writeLines(html, file.path(stage_dir, "REPORT.html"))

write_json(list(
  schema_version = "tcf3_ebf1_coverage_placseq_tracks.v4",
  status = "complete_unpromoted_report_only",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script_version = script_version,
  seed = seed,
  genome = "GRCm38/mm10",
  region = list(chromosome = chromosome, start_1based = region_start,
    end_1based = region_end, ebf1_tss_1based = ebf1_tss),
  normalization = "per_replicate_fragment_span_CPM_then_equal_weight_condition_mean",
  sample_counts = stats::setNames(manifest[, .N, by = condition]$N,
    manifest[, .N, by = condition]$condition),
  full_region = list(bin_width_bp = full_bin_bp,
    smoothing_width_bp = full_bin_bp * full_smoothing_bins),
  focused_region = list(chromosome = chromosome,
    start_1based = focus_start, end_1based = focus_end,
    placseq_links = nrow(focus_links)),
  zoom = list(peaks = nrow(peaks), flank_bp = zoom_flank_bp,
    bin_width_bp = fine_bin_bp,
    smoothing_width_bp = fine_bin_bp * zoom_smoothing_bins),
  placseq_links = nrow(links),
  implementation = "user_approved_outside_scworkbench_report_only",
  promotions = 0L,
  checks = checks
), file.path(stage_dir, "completion.json"))

if (!file.rename(stage_dir, output_dir))
  stop("Failed to finalize output directory")
log_event("Completed TCF3 Ebf1 coverage and PLAC-seq track report at",
  output_dir)
