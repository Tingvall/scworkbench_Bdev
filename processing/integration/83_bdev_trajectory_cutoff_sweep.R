#!/usr/bin/env Rscript

# Sensitivity analysis for progressively less strict frozen-Bdev trajectory
# support cutoffs. Report-only; no model, bundle, or mapped pseudotime changes.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 83_bdev_trajectory_cutoff_sweep.R BASE_REPORT ORIGINAL_REGIONS",
  "BDEV_BUNDLE QUERY_PROJECTION OUTPUT_DIR VERSION"), call. = FALSE)
base_dir <- normalizePath(args[[1L]], mustWork = TRUE)
region_path <- normalizePath(args[[2L]], mustWork = TRUE)
bdev_bundle_path <- normalizePath(args[[3L]], mustWork = TRUE)
projection_path <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
script_version <- args[[6L]]
if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output already exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")

query_path <- file.path(base_dir,
  "query_trajectory_support_and_acceptance.tsv.gz")
cv_path <- file.path(base_dir,
  "bdev_leave_one_sample_out_trajectory_support.tsv.gz")
query <- fread(cmd = paste("zcat", shQuote(query_path)))
cv <- fread(cmd = paste("zcat", shQuote(cv_path)))
base_density <- fread(file.path(base_dir,
  "query_bdev_normalized_density_by_policy.tsv"))
regions <- fread(region_path)[target %in% c("Flt3l", "Il7")]
stopifnot(nrow(query) == 23252L, sum(query$current_accepted) == 22471L,
  nrow(cv) == 26723L)

primary_cutoff <- 0.623237373726582
cutoffs <- c(primary_cutoff, .60, .55, .50, .45, .40, .35, .30, .25,
  .20, .10, 0)
cutoff_labels <- sprintf("%.3f", cutoffs)
cutoff_labels[[1L]] <- "0.623 (primary)"

bdev_grid <- unique(base_density[, .(bin, bin_left, bin_right, bin_midpoint,
  bdev_n_cells, bdev_fraction, bdev_rolling_fraction)])
setorder(bdev_grid, bin)
stopifnot(nrow(bdev_grid) == 100L)
breaks <- c(bdev_grid$bin_left[[1L]], bdev_grid$bin_right)
assign_bin <- function(pt) pmin(100L, pmax(1L,
  findInterval(pt, breaks, all.inside = TRUE)))
current_totals <- query[current_accepted & is.finite(mapped_pt),
  .(fixed_denominator = .N), by = .(knockout_target, genotype)]

acceptance <- rbindlist(lapply(seq_along(cutoffs), function(i) {
  cutoff_i <- cutoffs[[i]]
  query[, .(n_currently_mapped = sum(current_accepted),
    n_accepted = sum(current_accepted & trajectory_support >= cutoff_i),
    accepted_fraction = sum(current_accepted &
      trajectory_support >= cutoff_i) / sum(current_accepted)),
    by = .(sample_id, knockout_target, genotype)][,
      `:=`(cutoff = cutoff_i, cutoff_label = cutoff_labels[[i]])]
}))

density <- rbindlist(lapply(seq_along(cutoffs), function(i) {
  cutoff_i <- cutoffs[[i]]
  rbindlist(lapply(c("Flt3l", "Il7"), function(target_i) {
    rbindlist(lapply(c("WT", "KO"), function(genotype_i) {
      cells <- query[current_accepted & trajectory_support >= cutoff_i &
        knockout_target == target_i & genotype == genotype_i &
        is.finite(mapped_pt)]
      counts <- cells[, .(n_cells = .N), by = .(bin = assign_bin(mapped_pt))]
      ans <- merge(copy(bdev_grid), counts, by = "bin", all.x = TRUE)
      ans[is.na(n_cells), n_cells := 0L]
      den <- current_totals[knockout_target == target_i &
        genotype == genotype_i, fixed_denominator]
      ans[, fixed_fraction := n_cells / den]
      ans[, rolling_fraction := frollmean(fixed_fraction, 3L,
        align = "center", fill = NA_real_)]
      ans[, enrichment_ratio := rolling_fraction / bdev_rolling_fraction]
      ans[!is.finite(enrichment_ratio), enrichment_ratio := NA_real_]
      ans[, `:=`(cutoff = cutoff_i, cutoff_label = cutoff_labels[[i]],
        target = target_i, genotype = genotype_i,
        total_accepted = nrow(cells), fixed_denominator = den)]
      ans
    }))
  }))
}))

wt_above <- density[genotype == "WT" & is.finite(enrichment_ratio) &
  enrichment_ratio > 1]
setorder(wt_above, cutoff, target, bin)
wt_above[, run := cumsum(seq_len(.N) == 1L |
  cutoff != shift(cutoff, fill = cutoff[[1L]]) |
  target != shift(target, fill = target[[1L]]) |
  bin != shift(bin, fill = bin[[1L]] - 1L) + 1L)]
peaks <- wt_above[, .(first_bin = min(bin), last_bin = max(bin),
  n_bins = .N, pt_left = min(bin_left), pt_right = max(bin_right),
  maximum_wt_ratio = max(enrichment_ratio)),
  by = .(cutoff, cutoff_label, target, run)]
setorder(peaks, -cutoff, target, pt_left)
peaks[, peak_number := seq_len(.N), by = .(cutoff, target)]

peak_summary <- density[genotype == "WT", .(
  maximum_wt_ratio = max(enrichment_ratio, na.rm = TRUE),
  bins_above_one = sum(enrichment_ratio > 1, na.rm = TRUE)),
  by = .(cutoff, cutoff_label, target)]
peak_summary[, any_enrichment_peak := bins_above_one > 0]

region_retention <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
  reg <- regions[i]
  rbindlist(lapply(seq_along(cutoffs), function(j) {
    d <- query[current_accepted & knockout_target == reg$target &
      is.finite(mapped_pt) & mapped_pt >= reg$pt_left &
      mapped_pt < reg$pt_right]
    d[, .(n_original_region_cells = .N,
      n_accepted = sum(trajectory_support >= cutoffs[[j]]),
      accepted_fraction = mean(trajectory_support >= cutoffs[[j]])),
      by = genotype][, `:=`(region_id = reg$region_id, target = reg$target,
        cutoff = cutoffs[[j]], cutoff_label = cutoff_labels[[j]])]
  }))
}))

write_tsv(acceptance, file.path(stage_dir, "cell_retention_by_cutoff.tsv"))
write_tsv(region_retention, file.path(stage_dir,
  "original_region_retention_by_cutoff.tsv"))
write_tsv(peak_summary, file.path(stage_dir,
  "maximum_wt_enrichment_by_cutoff.tsv"))
write_tsv(peaks[, !"run"], file.path(stage_dir,
  "wt_enrichment_intervals_by_cutoff.tsv"))
write_tsv(density, file.path(stage_dir, "density_curves_by_cutoff.tsv.gz"))

genotype_cols <- c("WT" = "#0072B2", "KO" = "#D55E00")
target_cols <- c("Flt3l" = "#E69F00", "Il7" = "#0072B2")
base_theme <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
  strip.text = element_text(face = "bold"), legend.position = "bottom")
p_cells <- ggplot(acceptance,
    aes(cutoff, accepted_fraction, colour = sample_id)) +
  geom_line(linewidth = .8) + geom_point(size = 1.8) +
  geom_vline(xintercept = primary_cutoff, linetype = "dashed") +
  scale_x_reverse(breaks = c(.623, .5, .4, .3, .2, .1, 0)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Trajectory-support cutoff (less strict →)",
    y = "Currently mapped cells retained", colour = "Sample",
    title = "A. Cell recovery as the cutoff is relaxed") + base_theme
p_peak <- ggplot(peak_summary,
    aes(cutoff, maximum_wt_ratio, colour = target)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line(linewidth = .8) + geom_point(size = 1.8) +
  geom_vline(xintercept = primary_cutoff, linetype = "dashed") +
  scale_x_reverse(breaks = c(.623, .5, .4, .3, .2, .1, 0)) +
  scale_colour_manual(values = target_cols) +
  labs(x = "Trajectory-support cutoff (less strict →)",
    y = "Maximum WT/Bdev enrichment ratio", colour = "Target",
    title = "B. Strongest apparent WT enrichment") + base_theme
p_region <- ggplot(region_retention,
    aes(cutoff, accepted_fraction, colour = genotype)) +
  geom_line(linewidth = .8) + geom_point(size = 1.6) +
  geom_vline(xintercept = primary_cutoff, linetype = "dashed") +
  facet_wrap(~region_id, nrow = 1) +
  scale_x_reverse(breaks = c(.623, .5, .4, .3, .2, .1, 0)) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = genotype_cols) +
  labs(x = "Trajectory-support cutoff (less strict →)",
    y = "Original-region cells retained", colour = "Genotype",
    title = "C. Recovery of cells in the original mapped peaks") + base_theme
sweep_figure <- (p_cells | p_peak) / p_region +
  plot_annotation(title = "Sensitivity to progressively less strict Bdev trajectory gates",
    subtitle = "Dashed line: primary sqrt-on cutoff 0.623; density normalization remains fixed to all currently mapped cells")
ggsave(file.path(stage_dir, "figures", "cutoff_sweep_summary.pdf"),
  sweep_figure, width = 13, height = 8.5, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "cutoff_sweep_summary.png"),
  sweep_figure, width = 13, height = 8.5, dpi = 220, bg = "white")

representative <- c(primary_cutoff, .55, .50, .40)
curve_plot <- density[cutoff %in% representative]
curve_plot[, cutoff_display := factor(sprintf("cutoff %.3f", cutoff),
  levels = sprintf("cutoff %.3f", representative))]
p_curves <- ggplot(curve_plot,
    aes(bin_midpoint, enrichment_ratio, colour = genotype)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line(linewidth = .8, na.rm = TRUE) +
  facet_grid(target ~ cutoff_display, scales = "free_y") +
  scale_colour_manual(values = genotype_cols) +
  labs(x = "Mapped Bdev pseudotime",
    y = "Retained query / Bdev rolling fraction", colour = "Genotype",
    title = "Enrichment curves at representative permissive cutoffs",
    subtitle = "All panels retain the original currently mapped sample total as denominator") +
  base_theme
ggsave(file.path(stage_dir, "figures", "representative_cutoff_curves.pdf"),
  p_curves, width = 14, height = 7, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "representative_cutoff_curves.png"),
  p_curves, width = 14, height = 7, dpi = 220, bg = "white")

# Render every tested cutoff in a readable three-page PDF (four cutoffs/page),
# including the unfiltered cutoff 0 where the rejected early peaks return.
make_curve_page <- function(page_cutoffs, page_number) {
  d <- density[cutoff %in% page_cutoffs]
  d[, cutoff_display := factor(sprintf("cutoff %.3f", cutoff),
    levels = sprintf("cutoff %.3f", page_cutoffs))]
  ggplot(d, aes(bin_midpoint, enrichment_ratio, colour = genotype)) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_line(linewidth = .8, na.rm = TRUE) +
    facet_grid(target ~ cutoff_display, scales = "free_y") +
    scale_colour_manual(values = genotype_cols) +
    labs(x = "Mapped Bdev pseudotime",
      y = "Retained query / Bdev rolling fraction", colour = "Genotype",
      title = sprintf("Enrichment curves across all cutoffs — page %d of 3",
        page_number),
      subtitle = "Fixed denominator: all currently mapped cells in each target/genotype") +
    base_theme
}
curve_pages <- lapply(seq_len(3L), function(i) {
  idx <- ((i - 1L) * 4L + 1L):(i * 4L)
  make_curve_page(cutoffs[idx], i)
})
pdf_path <- file.path(stage_dir, "figures", "all_cutoff_enrichment_curves.pdf")
cairo_pdf(pdf_path, width = 14, height = 7)
for (p in curve_pages) print(p)
dev.off()
for (i in seq_along(curve_pages)) {
  ggsave(file.path(stage_dir, "figures", sprintf(
    "all_cutoff_enrichment_curves_page%d.png", i)), curve_pages[[i]],
    width = 14, height = 7, dpi = 220, bg = "white")
}

# Display the held-out Bdev trajectory-support values in the native Bdev UMAP.
# These are cross-sample prediction scores, not in-sample neighbor fractions.
bdev <- wb_open_bundle(bdev_bundle_path, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
bdev_ids <- as.character(wb_cells(bdev, "all"))
bdev_md <- as.data.table(wb_meta(bdev, cells = bdev_ids))
bdev_umap <- as.data.table(wb_space(bdev, "umap.atac", cells = bdev_ids))
stopifnot(identical(as.character(bdev_md$cell_id), bdev_ids),
  identical(as.character(bdev_umap$cell_id), bdev_ids),
  "peaks_snn_res.1.2" %in% names(bdev_md))
cv[, source_cell_id := sub("^bdev::", "", cell_id)]
bdev_plot <- merge(bdev_md[, .(source_cell_id = cell_id,
  cluster = as.character(`peaks_snn_res.1.2`), pt)],
  bdev_umap[, .(source_cell_id = cell_id, x, y)], by = "source_cell_id")
bdev_plot <- merge(bdev_plot, cv[, .(source_cell_id, held_out_sample,
  trajectory_truth, trajectory_support)], by = "source_cell_id")
stopifnot(nrow(bdev_plot) == 26723L, !anyDuplicated(bdev_plot$source_cell_id),
  identical(is.finite(bdev_plot$pt), bdev_plot$trajectory_truth))
bdev_plot[, cluster18 := cluster == "18"]

cluster_summary <- bdev_plot[, .(n_cells = .N,
  n_on_trajectory = sum(trajectory_truth),
  on_trajectory_fraction = mean(trajectory_truth),
  support_mean = mean(trajectory_support),
  support_median = median(trajectory_support),
  support_q25 = quantile(trajectory_support, .25),
  support_q75 = quantile(trajectory_support, .75),
  fraction_at_least_0623 = mean(trajectory_support >= primary_cutoff)),
  by = cluster]
cluster_summary[, cluster_number := as.integer(cluster)]
setorder(cluster_summary, cluster_number)
write_tsv(cluster_summary, file.path(stage_dir,
  "bdev_trajectory_support_by_cluster.tsv"))

score_scale <- scale_colour_viridis_c(option = "magma", limits = c(0, 1),
  breaks = c(0, .25, .5, .75, 1), name = "Trajectory\nsupport")
p_bdev_score <- ggplot(bdev_plot[order(trajectory_support)],
    aes(x, y, colour = trajectory_support)) +
  geom_point(size = .22, alpha = .75) + score_scale + coord_equal() +
  labs(x = "Bdev UMAP 1", y = "Bdev UMAP 2",
    title = "A. All Bdev cells",
    subtitle = "Held-out trajectory-support score") +
  base_theme + theme(plot.title = element_text(size = 11),
    plot.subtitle = element_text(size = 9))
p_t_cluster <- ggplot() +
  geom_point(data = bdev_plot[cluster18 == FALSE], aes(x, y), colour = "grey85",
    size = .18) +
  geom_point(data = bdev_plot[cluster18 == TRUE],
    aes(x, y, colour = trajectory_support), size = .65, alpha = .9) +
  score_scale + coord_equal() +
  labs(x = "Bdev UMAP 1", y = "Bdev UMAP 2",
    title = "B. T-lineage cluster 18 only",
    subtitle = sprintf("n=%d; finite-pt=%d; median support=%.3f",
      cluster_summary[cluster == "18", n_cells],
      cluster_summary[cluster == "18", n_on_trajectory],
      cluster_summary[cluster == "18", support_median])) + base_theme +
  theme(plot.title = element_text(size = 11),
    plot.subtitle = element_text(size = 9))
cluster_order <- cluster_summary[order(support_median), cluster]
bdev_plot[, cluster_ordered := factor(cluster, levels = cluster_order)]
p_cluster <- ggplot(bdev_plot,
    aes(cluster_ordered, trajectory_support,
      fill = cluster18)) +
  geom_boxplot(outlier.shape = NA, linewidth = .3, width = .7) +
  geom_hline(yintercept = primary_cutoff, linetype = "dashed",
    colour = "#CC79A7") +
  scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = "#D55E00"),
    labels = c("Other clusters", "T-lineage cluster 18")) +
  coord_flip() +
  labs(x = "Bdev cluster (ordered by median support)",
    y = "Held-out trajectory-support score", fill = NULL,
    title = "C. Score distribution by Bdev cluster",
    subtitle = "Dashed line: primary cutoff 0.623") + base_theme
bdev_figure <- (p_bdev_score | p_t_cluster) / p_cluster +
  plot_annotation(title = "Where the Bdev-derived trajectory-support score lies in the reference") +
  plot_layout(heights = c(1, 1.15))
ggsave(file.path(stage_dir, "figures", "bdev_trajectory_support_umap.pdf"),
  bdev_figure, width = 13, height = 11, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "bdev_trajectory_support_umap.png"),
  bdev_figure, width = 13, height = 11, dpi = 220, bg = "white")

display_cutoffs <- c(.50, .40, .30, .20, .10)
cutoff_umap <- rbindlist(lapply(display_cutoffs, function(cutoff_i) {
  ans <- copy(bdev_plot)
  ans[, `:=`(cutoff = cutoff_i,
    included = trajectory_support >= cutoff_i)]
  ans
}))
cutoff_counts <- cutoff_umap[, .(n_cells = .N,
  n_included = sum(included), included_fraction = mean(included),
  n_cluster18_included = sum(included & cluster18)), by = cutoff]
cutoff_counts[, panel_label := sprintf(
  "cutoff %.1f\nincluded %s/%s (%.1f%%); T18 %s/447",
  cutoff, format(n_included, big.mark = ","), format(n_cells, big.mark = ","),
  100 * included_fraction, format(n_cluster18_included, big.mark = ","))]
cutoff_umap <- cutoff_counts[, .(cutoff, panel_label)][cutoff_umap,
  on = "cutoff"]
cutoff_umap[, panel_label := factor(panel_label,
  levels = cutoff_counts[match(display_cutoffs, cutoff), panel_label])]
cutoff_umap[, inclusion_status := fifelse(included, "Included", "Excluded")]
setorder(cutoff_umap, included)
p_cutoff_umap <- ggplot(cutoff_umap,
    aes(x, y, colour = inclusion_status)) +
  geom_point(size = .16, alpha = .68) +
  facet_wrap(~panel_label, nrow = 1) + coord_equal() +
  scale_colour_manual(values = c("Excluded" = "grey82",
    "Included" = "#0072B2")) +
  labs(x = "Bdev UMAP 1", y = "Bdev UMAP 2", colour = "Bdev cell",
    title = "Bdev cells included by increasingly permissive trajectory-support cutoffs",
    subtitle = "Scores are from leave-one-original-sample-out Bdev calibration; T18 denotes T-lineage cluster 18") +
  base_theme
ggsave(file.path(stage_dir, "figures", "bdev_cells_included_by_cutoff.pdf"),
  p_cutoff_umap, width = 17, height = 5, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "bdev_cells_included_by_cutoff.png"),
  p_cutoff_umap, width = 17, height = 5, dpi = 220, bg = "white")
write_tsv(cutoff_counts, file.path(stage_dir,
  "bdev_cell_inclusion_counts_by_cutoff.tsv"))

# Overlay the thresholded WT/KO query projection on the corresponding Bdev
# accepted-cell background. Current general-support mapping remains a required
# prerequisite; each displayed cutoff is then applied to trajectory support.
projection <- fread(cmd = paste("zcat", shQuote(projection_path)))[,
  .(cell_id, projected_x, projected_y)]
query_xy <- merge(query, projection, by = "cell_id", all.x = TRUE)
stopifnot(nrow(query_xy) == nrow(query),
  sum(query_xy$current_accepted & is.finite(query_xy$projected_x)) == 22471L)
query_cutoff <- rbindlist(lapply(display_cutoffs, function(cutoff_i) {
  ans <- copy(query_xy[current_accepted & is.finite(projected_x) &
    is.finite(projected_y)])
  ans[, `:=`(cutoff = cutoff_i,
    included = trajectory_support >= cutoff_i)]
  ans
}))
query_counts <- query_cutoff[, .(n_currently_mapped = .N,
  n_included = sum(included), included_fraction = mean(included)),
  by = .(cutoff, sample_id, genotype, knockout_target)]
query_counts[, count_label := sprintf("query %s/%s (%.1f%%)",
  format(n_included, big.mark = ","),
  format(n_currently_mapped, big.mark = ","), 100 * included_fraction)]
query_cutoff[, cutoff_display := factor(sprintf("cutoff %.1f", cutoff),
  levels = sprintf("cutoff %.1f", display_cutoffs))]
query_cutoff[, sample_display := factor(sample_id,
  levels = c("FLT3LKO_BM_WT", "FLT3LKO_BM_KO", "IL7KO_BM_WT",
    "IL7KO_BM_KO"))]
query_cutoff[, query_status := fifelse(included,
  as.character(genotype), "Rejected query")]
setorder(query_cutoff, included)

bdev_cutoff_background <- cutoff_umap[included == TRUE,
  .(x, y, cutoff_display = factor(sprintf("cutoff %.1f", cutoff),
    levels = sprintf("cutoff %.1f", display_cutoffs)))]
annotation <- copy(query_counts)
annotation[, `:=`(
  cutoff_display = factor(sprintf("cutoff %.1f", cutoff),
    levels = sprintf("cutoff %.1f", display_cutoffs)),
  sample_display = factor(sample_id,
    levels = c("FLT3LKO_BM_WT", "FLT3LKO_BM_KO", "IL7KO_BM_WT",
      "IL7KO_BM_KO")),
  x = min(bdev_plot$x), y = max(bdev_plot$y))]
p_query_cutoff <- ggplot() +
  geom_point(data = bdev_plot, aes(x, y), colour = "grey93", size = .06) +
  geom_point(data = bdev_cutoff_background, aes(x, y),
    colour = "#B8D8CF", size = .07, alpha = .7) +
  geom_point(data = query_cutoff,
    aes(projected_x, projected_y, colour = query_status),
    size = .18, alpha = .62) +
  geom_text(data = annotation, aes(x, y, label = count_label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.5,
    colour = "black") +
  facet_grid(sample_display ~ cutoff_display) + coord_equal() +
  scale_colour_manual(values = c("Rejected query" = "grey68",
    "WT" = "#0072B2", "KO" = "#D55E00")) +
  labs(x = "Bdev UMAP 1", y = "Bdev UMAP 2", colour = "Query status",
    title = "Mapped WT/KO cells retained at each Bdev trajectory-support cutoff",
    subtitle = "Pale green: Bdev cells passing the cutoff; grey query: mapped but rejected by the cutoff") +
  base_theme + theme(axis.text = element_text(size = 7),
    strip.text.y = element_text(size = 8))
ggsave(file.path(stage_dir, "figures",
  "mapped_wt_ko_by_bdev_cutoff.pdf"), p_query_cutoff,
  width = 17, height = 14, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures",
  "mapped_wt_ko_by_bdev_cutoff.png"), p_query_cutoff,
  width = 17, height = 14, dpi = 220, bg = "white")
write_tsv(query_counts, file.path(stage_dir,
  "mapped_wt_ko_inclusion_counts_by_cutoff.tsv"))

report <- c("# Bdev trajectory-support cutoff sweep", "",
  "This report relaxes the primary sample-aware cutoff from 0.623 to 0 while leaving the frozen Bdev model, query projection, mapped pseudotime, and source bundles unchanged.", "",
  "Query densities use the original currently mapped cell count for each target/genotype as a fixed denominator. Consequently, lowering the cutoff can add density but cannot inflate a retained subset merely through renormalization.", "",
  "## Outputs", "",
  "- `cell_retention_by_cutoff.tsv`",
  "- `original_region_retention_by_cutoff.tsv`",
  "- `maximum_wt_enrichment_by_cutoff.tsv`",
  "- `wt_enrichment_intervals_by_cutoff.tsv`",
  "- `bdev_trajectory_support_by_cluster.tsv`",
  "- `bdev_cell_inclusion_counts_by_cutoff.tsv`",
  "- `mapped_wt_ko_inclusion_counts_by_cutoff.tsv`",
  "- `figures/cutoff_sweep_summary.pdf`",
  "- `figures/representative_cutoff_curves.pdf`",
  "- `figures/all_cutoff_enrichment_curves.pdf`",
  "- `figures/bdev_trajectory_support_umap.pdf`", "",
  "- `figures/bdev_cells_included_by_cutoff.pdf`",
  "- `figures/mapped_wt_ko_by_bdev_cutoff.pdf`", "",
  sprintf("Script version: `%s`.", script_version))
writeLines(report, file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete_unpromoted_report_only",
  script_version = script_version, primary_cutoff = primary_cutoff,
  tested_cutoffs = cutoffs, parent_report = base_dir,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  file.path(stage_dir, "receipt.json"), pretty = TRUE, auto_unbox = TRUE)

stopifnot(file.exists(file.path(stage_dir, "figures",
  "cutoff_sweep_summary.pdf")), file.exists(pdf_path),
  file.exists(file.path(stage_dir, "figures",
    "bdev_trajectory_support_umap.pdf")),
  file.exists(file.path(stage_dir, "figures",
    "bdev_cells_included_by_cutoff.pdf")),
  file.exists(file.path(stage_dir, "figures",
    "mapped_wt_ko_by_bdev_cutoff.pdf")))
if (!file.rename(stage_dir, output_dir))
  stop("Failed to finalize output directory", call. = FALSE)
cat("Completed cutoff sweep: ", output_dir, "\n", sep = "")
