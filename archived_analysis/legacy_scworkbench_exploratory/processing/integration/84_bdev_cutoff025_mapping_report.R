#!/usr/bin/env Rscript

# Focused report for the user-selected Bdev trajectory-support cutoff 0.25.
# Report-only: the frozen Bdev model, query projection, bundles, and preserved
# Bdev pseudotime are read unchanged.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 84_bdev_cutoff025_mapping_report.R BASE_REPORT BDEV_BUNDLE",
  "QUERY_PROJECTION BDEV_SIGNAL_REGIONS OUTPUT_DIR VERSION"), call. = FALSE)
base_dir <- normalizePath(args[[1L]], mustWork = TRUE)
bdev_bundle_path <- normalizePath(args[[2L]], mustWork = TRUE)
projection_path <- normalizePath(args[[3L]], mustWork = TRUE)
signal_region_path <- normalizePath(args[[4L]], mustWork = TRUE)
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
  library(grid)
})
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
cutoff <- 0.25
n_bins <- 100L
n_background_draws <- 1000L
seed <- 1729L
set.seed(seed)

cv_path <- file.path(base_dir,
  "bdev_leave_one_sample_out_trajectory_support.tsv.gz")
query_path <- file.path(base_dir,
  "query_trajectory_support_and_acceptance.tsv.gz")
cv <- fread(cmd = paste("zcat", shQuote(cv_path)))
query <- fread(cmd = paste("zcat", shQuote(query_path)))
base_density <- fread(file.path(base_dir,
  "query_bdev_normalized_density_by_policy.tsv"))
projection <- fread(cmd = paste("zcat", shQuote(projection_path)))[,
  .(cell_id, projected_x, projected_y)]
signal_regions <- fread(signal_region_path)[feature_set %in% c(
  "parse_cytokine_rest_only_refinement_all_top100_v1.flt3l.up",
  "parse_cytokine_rest_only_refinement_all_top100_v1.il_7.up")]
signal_regions[, target := fifelse(grepl("flt3l.up$", feature_set),
  "Flt3l", "Il7")]
signal_regions <- signal_regions[, .(target, feature_set, region_id,
  region_left, region_right, peak_pseudotime, p_region_fwer, q_dynamic)]
stopifnot(nrow(signal_regions[target == "Flt3l"]) == 2L,
  nrow(signal_regions[target == "Il7"]) == 1L)
stopifnot(nrow(cv) == 26723L, nrow(query) == 23252L)

bdev <- wb_open_bundle(bdev_bundle_path, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
bdev_ids <- as.character(wb_cells(bdev, "all"))
bdev_md <- as.data.table(wb_meta(bdev, cells = bdev_ids))
bdev_umap <- as.data.table(wb_space(bdev, "umap.atac", cells = bdev_ids))
stopifnot(identical(as.character(bdev_md$cell_id), bdev_ids),
  identical(as.character(bdev_umap$cell_id), bdev_ids))
cv[, source_cell_id := sub("^bdev::", "", cell_id)]
bdev_plot <- merge(bdev_md[, .(source_cell_id = cell_id,
  cluster = as.character(`peaks_snn_res.1.2`), pt)],
  bdev_umap[, .(source_cell_id = cell_id, x, y)], by = "source_cell_id")
bdev_plot <- merge(bdev_plot, cv[, .(source_cell_id, held_out_sample,
  trajectory_truth, trajectory_support)], by = "source_cell_id")
stopifnot(nrow(bdev_plot) == 26723L,
  identical(is.finite(bdev_plot$pt), bdev_plot$trajectory_truth))
bdev_plot[, `:=`(cluster18 = cluster == "18",
  included = trajectory_support >= cutoff)]

t_background <- bdev_plot[cluster18 == TRUE, .(
  n_cells = .N, n_defined_pt = sum(trajectory_truth),
  support_q75 = quantile(trajectory_support, .75),
  support_q90 = quantile(trajectory_support, .90),
  support_q95 = quantile(trajectory_support, .95),
  n_included = sum(included), included_fraction = mean(included))]
calibration <- data.table(cutoff = cutoff,
  sensitivity = bdev_plot[trajectory_truth == TRUE, mean(included)],
  specificity = bdev_plot[trajectory_truth == FALSE, mean(!included)],
  balanced_accuracy = (bdev_plot[trajectory_truth == TRUE, mean(included)] +
    bdev_plot[trajectory_truth == FALSE, mean(!included)]) / 2,
  t_cluster_q90 = t_background$support_q90,
  margin_above_t_cluster_q90 = cutoff - t_background$support_q90,
  t_cluster_excluded_fraction = 1 - t_background$included_fraction)

query <- merge(query, projection, by = "cell_id", all.x = TRUE)
query[, accepted := current_accepted & trajectory_support >= cutoff]
acceptance <- query[, .(n_query_cells = .N,
  n_currently_mapped = sum(current_accepted), n_accepted = sum(accepted),
  accepted_fraction_currently_mapped = sum(accepted) / sum(current_accepted)),
  by = .(sample_id, knockout_target, genotype)]

bdev_grid <- unique(base_density[, .(bin, bin_left, bin_right, bin_midpoint,
  bdev_n_cells, bdev_fraction, bdev_rolling_fraction)])
setorder(bdev_grid, bin)
breaks <- c(bdev_grid$bin_left[[1L]], bdev_grid$bin_right)
assign_bin <- function(pt) pmin(n_bins, pmax(1L,
  findInterval(pt, breaks, all.inside = TRUE)))
current_totals <- query[current_accepted & is.finite(mapped_pt),
  .(fixed_denominator = .N), by = .(knockout_target, genotype)]
density <- rbindlist(lapply(c("Flt3l", "Il7"), function(target_i) {
  rbindlist(lapply(c("WT", "KO"), function(genotype_i) {
    cells <- query[accepted & knockout_target == target_i &
      genotype == genotype_i & is.finite(mapped_pt)]
    counts <- cells[, .(n_cells = .N), by = .(bin = assign_bin(mapped_pt))]
    ans <- merge(copy(bdev_grid), counts, by = "bin", all.x = TRUE)
    ans[is.na(n_cells), n_cells := 0L]
    den <- current_totals[knockout_target == target_i &
      genotype == genotype_i, fixed_denominator]
    ans[, fixed_fraction := n_cells / den]
    ans[, conditional_fraction := n_cells / nrow(cells)]
    ans[, rolling_fraction := frollmean(fixed_fraction, 3L,
      align = "center", fill = NA_real_)]
    ans[, conditional_rolling_fraction := frollmean(conditional_fraction,
      3L, align = "center", fill = NA_real_)]
    ans[, enrichment_ratio := rolling_fraction / bdev_rolling_fraction]
    ans[, conditional_enrichment_ratio := conditional_rolling_fraction /
      bdev_rolling_fraction]
    ans[!is.finite(enrichment_ratio), enrichment_ratio := NA_real_]
    ans[!is.finite(conditional_enrichment_ratio),
      conditional_enrichment_ratio := NA_real_]
    ans[, `:=`(target = target_i, genotype = genotype_i,
      total_accepted = nrow(cells), fixed_denominator = den)]
    ans
  }))
}))

# Define genotype-specific enrichment thresholds from Bdev background rather
# than using a fixed ratio. For each target and genotype, draw that group's
# retained sample size from the Bdev pseudotime distribution 1,000 times and
# take the bin-wise 95th percentile of the resulting conditional
# Bdev-normalized rolling ratios. WT remains the primary peak definition; KO
# above-threshold bins and intervals are retained for future use.
background_threshold <- rbindlist(lapply(c("Flt3l", "Il7"), function(target_i) {
  rbindlist(lapply(c("WT", "KO"), function(genotype_i) {
    n_retained <- density[target == target_i & genotype == genotype_i,
      unique(total_accepted)]
    stopifnot(length(n_retained) == 1L, n_retained > 0L)
    sims <- replicate(n_background_draws,
      as.numeric(rmultinom(1L, n_retained, bdev_grid$bdev_fraction)))
    sim_fraction <- sims / n_retained
    sim_rolling <- apply(sim_fraction, 2L, function(z)
      frollmean(z, 3L, align = "center", fill = NA_real_))
    sim_ratio <- sim_rolling / bdev_grid$bdev_rolling_fraction
    data.table(target = target_i, genotype = genotype_i,
      bin = bdev_grid$bin, bin_midpoint = bdev_grid$bin_midpoint,
      background_ratio_median = apply(sim_ratio, 1L, median, na.rm = TRUE),
      enrichment_threshold_q95 = apply(sim_ratio, 1L, quantile,
        probs = .95, na.rm = TRUE), retained_n = n_retained)
  }))
}))
background_threshold[!is.finite(enrichment_threshold_q95),
  enrichment_threshold_q95 := NA_real_]

density <- background_threshold[, .(target, genotype, bin,
  enrichment_threshold_q95)][density, on = c("target", "genotype", "bin")]
above <- density[
  is.finite(conditional_enrichment_ratio) &
  is.finite(enrichment_threshold_q95) &
  conditional_enrichment_ratio > enrichment_threshold_q95]
setorder(above, target, genotype, bin)
above[, run := cumsum(bin != shift(bin, fill = first(bin) - 1L) + 1L),
  by = .(target, genotype)]
regions <- above[, .(first_bin = min(bin), last_bin = max(bin),
  n_bins = .N, pt_left = min(bin_left), pt_right = max(bin_right),
  peak_pt = bin_midpoint[which.max(conditional_enrichment_ratio -
    enrichment_threshold_q95)],
  maximum_ratio = max(conditional_enrichment_ratio),
  maximum_excess_over_background = max(conditional_enrichment_ratio -
    enrichment_threshold_q95)), by = .(target, genotype, run)]
setorder(regions, target, genotype, pt_left)
regions[, region_number := seq_len(.N), by = .(target, genotype)]
wt_regions <- copy(regions[genotype == "WT"])
wt_regions[, `:=`(peak_number = region_number,
  maximum_wt_ratio = maximum_ratio)]
ko_regions <- copy(regions[genotype == "KO"])
ko_regions[, ko_region_number := region_number]

write_tsv(calibration, file.path(stage_dir, "cutoff025_calibration.tsv"))
write_tsv(t_background, file.path(stage_dir,
  "t_cluster18_background_summary.tsv"))
write_tsv(acceptance, file.path(stage_dir, "query_acceptance_cutoff025.tsv"))
write_tsv(density, file.path(stage_dir, "query_density_cutoff025.tsv.gz"))
write_tsv(wt_regions[, !c("run", "region_number")], file.path(stage_dir,
  "wt_enrichment_regions_cutoff025.tsv"))
write_tsv(ko_regions[, !c("run", "region_number")], file.path(stage_dir,
  "ko_above_background_regions_cutoff025.tsv"))
write_tsv(regions[, !"run"], file.path(stage_dir,
  "genotype_above_background_regions_cutoff025.tsv"))
write_tsv(above[, .(target, genotype, bin, bin_left, bin_right,
  bin_midpoint, conditional_enrichment_ratio, enrichment_threshold_q95,
  excess_over_background = conditional_enrichment_ratio -
    enrichment_threshold_q95)], file.path(stage_dir,
  "genotype_above_background_pt_bins_cutoff025.tsv"))
write_tsv(above[genotype == "KO", .(target, genotype, bin, bin_left,
  bin_right, bin_midpoint, conditional_enrichment_ratio,
  enrichment_threshold_q95, excess_over_background =
    conditional_enrichment_ratio - enrichment_threshold_q95)],
  file.path(stage_dir, "ko_above_background_pt_bins_cutoff025.tsv"))
write_tsv(background_threshold, file.path(stage_dir,
  "bdev_background_enrichment_threshold_cutoff025.tsv"))
write_tsv(signal_regions, file.path(stage_dir,
  "bdev_significant_cytokine_signal_regions.tsv"))

genotype_cols <- c("WT" = "#0072B2", "KO" = "#D55E00")
base_theme <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
  strip.text = element_text(face = "bold"), legend.position = "bottom")

p_score <- ggplot(bdev_plot[order(trajectory_support)],
    aes(x, y, colour = trajectory_support)) +
  geom_point(size = .22, alpha = .75) + coord_equal() +
  scale_colour_viridis_c(option = "magma", limits = c(0, 1),
    breaks = c(0, .25, .5, .75, 1), name = "Trajectory\nsupport") +
  labs(x = "Bdev UMAP 1", y = "Bdev UMAP 2",
    title = "A. Held-out support score for every Bdev cell") + base_theme
p_status <- ggplot(bdev_plot[order(included)],
    aes(x, y, colour = fifelse(included, "Included", "Below 0.25"))) +
  geom_point(size = .22, alpha = .72) + coord_equal() +
  scale_colour_manual(values = c("Below 0.25" = "grey80",
    "Included" = "#009E73"), name = "Cutoff status") +
  labs(x = "Bdev UMAP 1", y = "Bdev UMAP 2",
    title = "B. Bdev cells retained at cutoff 0.25") + base_theme
score_density <- rbindlist(list(
  bdev_plot[, .(trajectory_support,
    group = fifelse(trajectory_truth, "Defined pt", "Missing pt"))],
  bdev_plot[cluster18 == TRUE, .(trajectory_support,
    group = "T-lineage cluster 18")]))
p_background <- ggplot(score_density,
    aes(trajectory_support, colour = group, fill = group)) +
  geom_density(alpha = .16, linewidth = .75) +
  geom_vline(xintercept = t_background$support_q90,
    linetype = "dotted", colour = "#D55E00", linewidth = .8) +
  geom_vline(xintercept = cutoff, colour = "black", linewidth = .85) +
  scale_colour_manual(values = c("Defined pt" = "#0072B2",
    "Missing pt" = "grey45", "T-lineage cluster 18" = "#D55E00")) +
  scale_fill_manual(values = c("Defined pt" = "#56B4E9",
    "Missing pt" = "grey65", "T-lineage cluster 18" = "#E69F00")) +
  labs(x = "Held-out trajectory-support score", y = "Density",
    colour = NULL, fill = NULL,
    title = "C. Reference support-score context",
    subtitle = sprintf("Orange dotted: T18 q90 %.3f; black: selected cutoff 0.25",
      t_background$support_q90)) + base_theme
bdev_figure <- (p_score | p_status) / p_background +
  plot_annotation(title = "Bdev trajectory support and selected cutoff 0.25") +
  plot_layout(heights = c(1, .8))
ggsave(file.path(stage_dir, "figures", "bdev_support_and_threshold.pdf"),
  bdev_figure, width = 13, height = 10, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "bdev_support_and_threshold.png"),
  bdev_figure, width = 13, height = 10, dpi = 220, bg = "white")

query_plot <- query[current_accepted & is.finite(projected_x) &
  is.finite(projected_y)]
query_plot[, status := fifelse(accepted, as.character(genotype),
  "Rejected by 0.25")]
query_plot[, sample_display := factor(sample_id,
  levels = c("FLT3LKO_BM_WT", "FLT3LKO_BM_KO", "IL7KO_BM_WT",
    "IL7KO_BM_KO"))]
setorder(query_plot, accepted)
annotation <- copy(acceptance)
annotation[, `:=`(sample_display = factor(sample_id,
  levels = c("FLT3LKO_BM_WT", "FLT3LKO_BM_KO", "IL7KO_BM_WT",
    "IL7KO_BM_KO")), x = min(bdev_plot$x), y = max(bdev_plot$y),
  label = sprintf("retained %s/%s (%.1f%%)",
    format(n_accepted, big.mark = ","),
    format(n_currently_mapped, big.mark = ","),
    100 * accepted_fraction_currently_mapped))]
p_query <- ggplot() +
  geom_point(data = bdev_plot[!is.na(x)], aes(x, y), colour = "grey91",
    size = .08) +
  geom_point(data = bdev_plot[included == TRUE], aes(x, y),
    colour = "#C4DDD6", size = .09) +
  geom_point(data = query_plot,
    aes(projected_x, projected_y, colour = status), size = .25, alpha = .65) +
  geom_text(data = annotation, aes(x, y, label = label), hjust = 0,
    vjust = 1, size = 3) + facet_wrap(~sample_display, ncol = 2) +
  coord_equal() + scale_colour_manual(values = c(
    "Rejected by 0.25" = "grey65", genotype_cols)) +
  labs(x = "Bdev UMAP 1", y = "Bdev UMAP 2", colour = "Query status",
    title = "Mapped WT/KO cells retained at trajectory-support cutoff 0.25",
    subtitle = "Pale green background: Bdev cells passing the same cutoff") +
  base_theme
ggsave(file.path(stage_dir, "figures", "retained_wt_ko_umap_cutoff025.pdf"),
  p_query, width = 11, height = 10, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "retained_wt_ko_umap_cutoff025.png"),
  p_query, width = 11, height = 10, dpi = 220, bg = "white")

p_density <- ggplot(density,
    aes(bin_midpoint, conditional_enrichment_ratio, colour = genotype)) +
  geom_rect(data = signal_regions,
    aes(xmin = region_left, xmax = region_right, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#CC79A7", alpha = .16) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line(data = background_threshold[genotype == "WT"],
    aes(bin_midpoint, enrichment_threshold_q95), inherit.aes = FALSE,
    colour = "black", linetype = "dashed", linewidth = .7) +
  geom_line(linewidth = .85, na.rm = TRUE) + facet_wrap(~target, ncol = 2,
    scales = "free_y") + scale_colour_manual(values = genotype_cols) +
  labs(x = "Mapped Bdev pseudotime",
    y = "Conditional query / Bdev rolling fraction", colour = "Genotype",
    title = "A. WT/KO enrichment along pseudotime at cutoff 0.25",
    subtitle = "Purple: significant Bdev residual-UCell regions; black dashed: WT-specific 95% Bdev background (KO intervals stored separately)") + base_theme
trajectory_figure <- p_density +
  plot_annotation(title = "Mapped enrichment with significant Bdev cytokine-signal regions")
ggsave(file.path(stage_dir, "figures",
  "wt_ko_pt_enrichment_bdev_signal_regions_cutoff025.pdf"),
  trajectory_figure, width = 13, height = 6, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures",
  "wt_ko_pt_enrichment_bdev_signal_regions_cutoff025.png"),
  trajectory_figure, width = 13, height = 6, dpi = 220, bg = "white")

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
acceptance_lines <- acceptance[order(knockout_target, genotype), sprintf(
  "%s %s: %s/%s (%s)", knockout_target, genotype,
  format(n_accepted, big.mark = ","),
  format(n_currently_mapped, big.mark = ","),
  fmt_pct(accepted_fraction_currently_mapped))]
region_lines <- wt_regions[, sprintf("%s peak %d: pt %.2f to %.2f, max %.2fx",
  target, peak_number, pt_left, pt_right, maximum_wt_ratio)]
ko_region_lines <- ko_regions[, sprintf(
  "%s KO region %d: pt %.2f to %.2f, peak pt %.2f, max %.2fx",
  target, ko_region_number, pt_left, pt_right, peak_pt, maximum_ratio)]
signal_region_lines <- signal_regions[order(target, region_left), sprintf(
  "%s Bdev signal: pt %.2f to %.2f", target, region_left, region_right)]
report_text <- paste(c(
  "Bdev trajectory-support cutoff 0.25 mapping report", "",
  "Approach",
  "The frozen Bdev-only PeakVI representation and existing query projection were reused without refitting. Each Bdev score is from leave-one-original-sample-out neighbor mapping. Query cells must first pass the existing general-support gate, then have trajectory support ≥0.25. Bdev pseudotime values and missing entries are unchanged.", "",
  "Support-cutoff context",
  sprintf("T-lineage cluster 18 has %d cells, %d with defined pt, q90 support %.3f, and q95 %.3f. The user-selected cutoff 0.25 is %.3f above q90 and excludes %s of this cluster.",
    t_background$n_cells, t_background$n_defined_pt,
    t_background$support_q90, t_background$support_q95,
    cutoff - t_background$support_q90,
    fmt_pct(1 - t_background$included_fraction)),
  sprintf("Across all held-out Bdev cells, sensitivity is %s and specificity is %s. Thus 0.25 is deliberately permissive and is not the ROC-optimal cutoff.",
    fmt_pct(calibration$sensitivity), fmt_pct(calibration$specificity)), "",
  "Retained mapped query cells", acceptance_lines, "",
  "Pseudotime-region enrichment threshold",
  sprintf("For each target and genotype, %d multinomial samples matching that group's retained cell count were drawn from the Bdev pseudotime distribution. The bin-wise 95th percentile of the simulated conditional Bdev-normalized rolling density is the genotype-specific background threshold. Consecutive WT bins above the WT threshold define the primary enriched regions.", n_background_draws), "",
  "WT enrichment intervals", region_lines, "",
  "KO above-background intervals retained for later use", ko_region_lines, "",
  "Significant Bdev cytokine-signal regions highlighted in the curves",
  signal_region_lines, "",
  "Interpretation",
  "This threshold recovers more mapped cells while still suppressing most of the low-support T-lineage island. Both Flt3L and IL-7 WT enrichment intervals are present at this permissive setting. Purple bands show significant cytokine-signal regions derived only from the default Bdev residual-UCell permutation analysis; no WT/KO signal score was calculated or plotted.", "",
  sprintf("Script version: %s", script_version)), collapse = "\n")
wrapped_report <- paste(vapply(strsplit(report_text, "\n", fixed = TRUE)[[1L]],
  function(x) paste(strwrap(x, width = 105), collapse = "\n"), character(1)),
  collapse = "\n")

report_pdf <- file.path(stage_dir, "CUTOFF025_MAPPING_REPORT.pdf")
cairo_pdf(report_pdf, width = 13, height = 10)
grid.newpage()
grid.text(wrapped_report, x = unit(.06, "npc"), y = unit(.96, "npc"),
  just = c("left", "top"), gp = gpar(fontsize = 9.5, lineheight = 1.12))
print(bdev_figure)
print(p_query)
print(trajectory_figure)
dev.off()

report_md <- c("# Bdev trajectory-support cutoff 0.25 mapping report", "",
  "The frozen Bdev model and mapped pseudotime were unchanged. The selected 0.25 cutoff is user-chosen and report-only.", "",
  "## Support-cutoff context", "",
  sprintf("Cluster 18 q90 support is %.3f; cutoff 0.25 is %.3f higher and excludes %s of cluster 18.",
    t_background$support_q90, cutoff - t_background$support_q90,
    fmt_pct(1 - t_background$included_fraction)),
  sprintf("Held-out Bdev sensitivity is %s and specificity is %s.",
    fmt_pct(calibration$sensitivity), fmt_pct(calibration$specificity)), "",
  "## Pseudotime-region threshold", "",
  sprintf("WT-enriched regions exceed the bin-wise 95th-percentile background from %d Bdev pseudotime resamples matching the retained WT cell count. KO regions and individual pseudotime bins exceeding a separately sample-size-matched KO threshold are retained as future-use artifacts.",
    n_background_draws), "",
  "## Output", "", "- `CUTOFF025_MAPPING_REPORT.pdf`", "",
  sprintf("Script version: `%s`.", script_version))
writeLines(report_md, file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete_unpromoted_report_only",
  script_version = script_version, selected_cutoff = cutoff,
  threshold_basis = "user_selected_trajectory_support_cutoff",
  enrichment_region_threshold = "genotype_specific_binwise_q95_from_1000_Bdev_background_multinomial_draws_matching_retained_n_primary_WT_regions",
  ko_above_background_pt_retained = TRUE,
  bdev_signal_region_source = signal_region_path,
  seed = seed,
  t_cluster18_q90 = t_background$support_q90,
  source_bundles_preserved = TRUE, completed_at = format(Sys.time(),
    "%Y-%m-%dT%H:%M:%S%z")), file.path(stage_dir, "receipt.json"),
  pretty = TRUE, auto_unbox = TRUE)

stopifnot(file.exists(report_pdf), file.info(report_pdf)$size > 0,
  file.exists(file.path(stage_dir, "figures",
    "bdev_support_and_threshold.pdf")),
  file.exists(file.path(stage_dir, "figures",
    "retained_wt_ko_umap_cutoff025.pdf")),
  file.exists(file.path(stage_dir, "figures",
    "wt_ko_pt_enrichment_bdev_signal_regions_cutoff025.pdf")))
if (!file.rename(stage_dir, output_dir))
  stop("Failed to finalize output directory", call. = FALSE)
cat("Completed cutoff-0.25 mapping report: ", output_dir, "\n", sep = "")
