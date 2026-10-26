#!/usr/bin/env Rscript

# Define supported WT-enrichment runs and visualize cell support using centered
# three-bin rolling windows. This is a read-only report over the accepted
# cutoff-0.25 mapping output; it does not modify a bundle or mapping result.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 85_wt_enrichment_rolling_window_cell_support.R",
  "<cutoff025_report> <output_dir> <minimum_cells>",
  "<rolling_width> <minimum_run_bins>"), call. = FALSE)
source_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
minimum_cells <- as.integer(args[[3L]])
rolling_width <- as.integer(args[[4L]])
minimum_run_bins <- as.integer(args[[5L]])
stopifnot(minimum_cells >= 1L, rolling_width >= 3L,
  rolling_width %% 2L == 1L, minimum_run_bins >= rolling_width)
if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output already exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

region_path <- file.path(source_dir, "wt_enrichment_regions_cutoff025.tsv")
density_path <- file.path(source_dir, "query_density_cutoff025.tsv.gz")
stopifnot(file.exists(region_path), file.exists(density_path))
regions <- fread(region_path)
density <- fread(cmd = paste("zcat", shQuote(density_path)))
required_region <- c("target", "first_bin", "last_bin", "peak_number")
required_density <- c("target", "genotype", "bin", "bin_left",
  "bin_right", "bin_midpoint", "n_cells")
stopifnot(all(required_region %in% names(regions)),
  all(required_density %in% names(density)),
  setequal(unique(density$genotype), c("WT", "KO")),
  !anyDuplicated(density[, .(target, genotype, bin)]))

counts <- dcast(density[, .(target, genotype, bin, bin_left, bin_right,
    bin_midpoint, n_cells)],
  target + bin + bin_left + bin_right + bin_midpoint ~ genotype,
  value.var = "n_cells", fill = 0L)
setnames(counts, c("WT", "KO"), c("wt_cells", "ko_cells"))
setorder(counts, target, bin)

source_rows <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
  region <- regions[i]
  ans <- counts[target == region$target & bin >= region$first_bin &
    bin <= region$last_bin]
  stopifnot(nrow(ans) == region$last_bin - region$first_bin + 1L)
  ans[, `:=`(wt_peak = region$peak_number,
    source_region_pt_left = region$pt_left,
    source_region_pt_right = region$pt_right)]
  ans
}))
setorder(source_rows, target, wt_peak, bin)
source_rows[, `:=`(wt_pass = wt_cells >= minimum_cells,
  ko_pass = ko_cells >= minimum_cells)]

# Identify maximal consecutive runs of source bins that all pass the WT rule.
passing <- source_rows[wt_pass == TRUE]
passing[, run_local := cumsum(bin != shift(bin, fill = first(bin) - 1L) + 1L),
  by = .(target, wt_peak)]
run_table <- passing[, .(first_bin = min(bin), last_bin = max(bin),
    n_source_bins = .N, source_pt_left = min(bin_left),
    source_pt_right = max(bin_right), wt_cells = sum(wt_cells),
    ko_cells = sum(ko_cells), ko_failed_source_bins = sum(!ko_pass)),
  by = .(target, wt_peak, run_local)]
setorder(run_table, target, wt_peak, first_bin)
run_table[, retained_run := n_source_bins >= minimum_run_bins]
run_table[, run_number := seq_len(.N), by = .(target, wt_peak)]
run_table[, region_id := sprintf("%s_wt_peak%d_run%d", tolower(target),
  wt_peak, run_number)]

source_rows <- merge(source_rows,
  passing[, .(target, wt_peak, bin, run_local)],
  by = c("target", "wt_peak", "bin"), all.x = TRUE, sort = FALSE)
source_rows <- merge(source_rows,
  run_table[, .(target, wt_peak, run_local, region_id, retained_run)],
  by = c("target", "wt_peak", "run_local"), all.x = TRUE, sort = FALSE)
source_rows[is.na(retained_run), retained_run := FALSE]
setorder(source_rows, target, wt_peak, bin)

half_width <- (rolling_width - 1L) %/% 2L
retained <- run_table[retained_run == TRUE]
rolling <- rbindlist(lapply(seq_len(nrow(retained)), function(i) {
  run <- retained[i]
  bins <- seq.int(run$first_bin, run$last_bin)
  centers <- bins[(1L + half_width):(length(bins) - half_width)]
  rbindlist(lapply(centers, function(center_bin) {
    members <- source_rows[target == run$target & wt_peak == run$wt_peak &
      bin >= center_bin - half_width & bin <= center_bin + half_width]
    stopifnot(nrow(members) == rolling_width, all(members$wt_pass))
    data.table(target = run$target, wt_peak = run$wt_peak,
      region_id = run$region_id, center_bin = center_bin,
      window_first_bin = min(members$bin),
      window_last_bin = max(members$bin),
      window_pt_left = min(members$bin_left),
      window_pt_right = max(members$bin_right),
      center_pt = members[bin == center_bin, bin_midpoint],
      wt_cells = sum(members$wt_cells),
      ko_cells = sum(members$ko_cells),
      center_bin_wt_cells = members[bin == center_bin, wt_cells],
      center_bin_ko_cells = members[bin == center_bin, ko_cells],
      ko_failed_source_bins_in_window = sum(!members$ko_pass))
  }))
}))
stopifnot(nrow(rolling) > 0L,
  all(rolling$wt_cells >= rolling_width * minimum_cells),
  !anyDuplicated(rolling[, .(region_id, center_bin)]))

# KO support is a property of the raw middle bin, not the overlapping rolling
# count. Exclude regions with no KO-supported middle bin at all. For retained
# regions, draw one boundary between the last supported middle bin and the
# first subsequent unsupported middle bin; do not mark individual windows.
rolling[, center_bin_ko_pass := center_bin_ko_cells >= minimum_cells]
ko_region_support <- rolling[, .(
  ko_supported_center_bins = sum(center_bin_ko_pass),
  ko_unsupported_center_bins = sum(!center_bin_ko_pass)),
  by = .(target, wt_peak, region_id)]
ko_region_support[, retained_for_visualization :=
  ko_supported_center_bins > 0L]
excluded_no_ko <- ko_region_support[retained_for_visualization == FALSE]
rolling <- rolling[region_id %in%
  ko_region_support[retained_for_visualization == TRUE, region_id]]
stopifnot(nrow(rolling) > 0L)

support_boundaries <- rolling[order(center_bin), {
  supported <- which(center_bin_ko_pass)
  last_supported <- max(supported)
  if (last_supported < .N &&
      all(!center_bin_ko_pass[seq.int(last_supported + 1L, .N)])) {
    list(last_supported_center_bin = center_bin[[last_supported]],
      first_unsupported_center_bin = center_bin[[last_supported + 1L]],
      ko_support_boundary_pt = mean(center_pt[c(last_supported,
        last_supported + 1L)]))
  } else {
    list(last_supported_center_bin = center_bin[[last_supported]],
      first_unsupported_center_bin = NA_integer_,
      ko_support_boundary_pt = NA_real_)
  }
}, by = .(target, wt_peak, region_id)]

included_regions <- rolling[, .(
  source_first_bin = min(window_first_bin),
  source_last_bin = max(window_last_bin),
  source_pt_left = min(window_pt_left),
  source_pt_right = max(window_pt_right),
  first_center_bin = min(center_bin),
  last_center_bin = max(center_bin),
  center_pt_left = min(center_pt),
  center_pt_right = max(center_pt),
  n_source_bins = max(window_last_bin) - min(window_first_bin) + 1L,
  n_rolling_windows = .N,
  middle_bins_ko_below_10 = sum(!center_bin_ko_pass),
  middle_bins_ko_at_least_10 = sum(center_bin_ko_pass)),
  by = .(target, wt_peak, region_id)]
included_regions <- merge(included_regions, support_boundaries,
  by = c("target", "wt_peak", "region_id"), all.x = TRUE, sort = FALSE)
setorder(included_regions, target, wt_peak)

plot_data <- melt(rolling,
  id.vars = c("target", "wt_peak", "region_id", "center_bin", "center_pt",
    "center_bin_ko_pass"), measure.vars = c("wt_cells", "ko_cells"),
  variable.name = "group", value.name = "n_cells")
plot_data[, group := factor(group, levels = c("wt_cells", "ko_cells"),
  labels = c("WT", "KO"))]
plot_data[, region_label := sprintf("%s WT peak %d", target, wt_peak)]
region_levels <- unique(plot_data[order(target, wt_peak), region_label])
plot_data[, region_label := factor(region_label, levels = region_levels)]
boundary_plot <- merge(support_boundaries,
  unique(plot_data[, .(target, wt_peak, region_id, region_label)]),
  by = c("target", "wt_peak", "region_id"), all.x = TRUE, sort = FALSE)

p <- ggplot(plot_data, aes(center_pt, n_cells, colour = group,
    group = group)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2) +
  geom_vline(data = boundary_plot[is.finite(ko_support_boundary_pt)],
    aes(xintercept = ko_support_boundary_pt), inherit.aes = FALSE,
    linetype = "dashed", colour = "#D55E00", linewidth = 0.8) +
  facet_wrap(~region_label, scales = "free_x", ncol = 1) +
  scale_colour_manual(values = c(WT = "#0072B2", KO = "#D55E00")) +
  scale_x_continuous(breaks = function(x) pretty(x, n = 6)) +
  theme_bw(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")) +
  labs(title = "Cell support in WT-enriched pseudotime regions",
    subtitle = paste0("Centered ", rolling_width,
      "-bin windows; retained source bins each contain at least ",
      minimum_cells, " WT cells"), x = "Middle-bin Bdev pseudotime",
    y = "Cells in centered three-bin window", colour = NULL,
    caption = paste0("Orange vertical line: boundary after the last middle source bin with at least ",
      minimum_cells, " KO cells; middle bins to its right are below threshold."))
ggsave(file.path(stage_dir, "figures",
  "wt_enrichment_rolling_3bin_cell_support.png"), p,
  width = 10, height = max(6, 2.4 * length(region_levels)), dpi = 300,
  bg = "white")
ggsave(file.path(stage_dir, "figures",
  "wt_enrichment_rolling_3bin_cell_support.pdf"), p,
  width = 10, height = max(6, 2.4 * length(region_levels)),
  device = cairo_pdf, bg = "white")

write_tsv(source_rows, file.path(stage_dir, "source_bin_cell_counts.tsv"))
write_tsv(rolling, file.path(stage_dir,
  "rolling_3bin_window_cell_counts.tsv"))
write_tsv(included_regions, file.path(stage_dir, "included_regions.tsv"))
write_tsv(excluded_no_ko, file.path(stage_dir,
  "excluded_regions_without_ko_supported_middle_bin.tsv"))
write_tsv(support_boundaries, file.path(stage_dir,
  "ko_middle_bin_support_boundaries.tsv"))
write_tsv(run_table[retained_run == FALSE], file.path(stage_dir,
  "excluded_short_wt_runs.tsv"))
write_tsv(run_table, file.path(stage_dir, "all_wt_supported_runs.tsv"))

checks <- list(
  source_density_rows_unique = !anyDuplicated(
    density[, .(target, genotype, bin)]),
  source_regions_exactly_subset = nrow(source_rows) ==
    sum(regions$last_bin - regions$first_bin + 1L),
  retained_source_bins_all_wt_ge_minimum = all(
    source_rows[retained_run == TRUE]$wt_pass),
  retained_regions_have_ko_supported_middle_bin = all(
    rolling[, any(center_bin_ko_pass), by = region_id]$V1),
  rolling_windows_exact_width = all(
    rolling$window_last_bin - rolling$window_first_bin + 1L == rolling_width),
  rolling_windows_centered = all(
    rolling$center_bin - rolling$window_first_bin == half_width &
    rolling$window_last_bin - rolling$center_bin == half_width),
  rolling_wt_counts_at_least_width_times_minimum = all(
    rolling$wt_cells >= rolling_width * minimum_cells))
stopifnot(all(unlist(checks)))

diagnostics <- list(status = "complete_unpromoted_report_only",
  source_report = source_dir, source_files = list(
    regions = region_path, density = density_path),
  rules = list(minimum_cells_per_source_wt_bin = minimum_cells,
    retain_maximal_consecutive_wt_passing_runs = TRUE,
    minimum_source_bins_per_retained_run = minimum_run_bins,
    rolling_width_bins = rolling_width,
    rolling_alignment = "centered_on_middle_source_bin",
    rolling_window_cell_count = "sum_of_cells_in_three_source_bins",
    ko_support_assessment = "raw_cell_count_in_middle_source_bin",
    exclude_region_when_no_ko_supported_middle_bin = TRUE,
    ko_highlight = paste0("one vertical boundary after the last middle bin with at least ",
      minimum_cells, " KO cells")),
  retained_regions = nrow(included_regions),
  retained_rolling_windows = nrow(rolling),
  middle_bins_below_ko_threshold = sum(!rolling$center_bin_ko_pass),
  checks = checks, bundle_mutated = FALSE, promotion = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))

region_lines <- vapply(seq_len(nrow(included_regions)), function(i) {
  x <- included_regions[i]
  sprintf("- `%s`: source bins %d-%d (pt %.3f to %.3f); displayed centers %d-%d (pt %.3f to %.3f); %d windows; %d middle bins below the KO threshold.",
    x$region_id, x$source_first_bin, x$source_last_bin,
    x$source_pt_left, x$source_pt_right, x$first_center_bin,
    x$last_center_bin, x$center_pt_left, x$center_pt_right,
    x$n_rolling_windows, x$middle_bins_ko_below_10)
}, character(1L))
writeLines(c("# WT-enrichment centered rolling-window cell support", "",
  paste0("A source bin is retained when it belongs to a maximal consecutive run in a WT-defined enrichment peak and contains at least **",
    minimum_cells, " WT cells**. A run must contain at least **",
    minimum_run_bins, " bins**, the minimum needed for a centered ",
    rolling_width, "-bin window."), "",
  paste0("Each displayed point is calculated from the union of the middle source bin and its ",
    half_width, " immediate bin on each side. Its x-coordinate is the middle-bin pseudotime. KO support is assessed from the unsmoothed middle source bin, not the rolling total. A region is excluded when none of its middle bins has at least ",
    minimum_cells, " KO cells. In each retained region, one orange vertical line is drawn between the last supported middle bin and the first subsequent unsupported middle bin; individual bins are not marked."), "",
  "## Included regions", "", region_lines, "",
  "![Rolling-window cell support](figures/wt_enrichment_rolling_3bin_cell_support.png)", "",
  "This is descriptive cell support. Each target/genotype is represented by one biological sample, and overlapping rolling windows are not independent."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed rolling-window support report at", output_dir, "\n")
