#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 64_cytokine_ucell_permutation.R",
  "<native_ucell_report> <output_dir> <seed> <permutations>"))
source_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
base_seed <- as.integer(args[[3L]])
n_permutations <- as.integer(args[[4L]])
stopifnot(is.finite(base_seed), identical(n_permutations, 1000L))
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(matrixStats)
})

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA")
write_tsv_gz <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA", compress = "gzip")
short_label <- function(x) sub(
  "^parse_cytokine_rest_only_refinement_all_top100_v1\\.", "", x)
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ...,
    "\n")
  flush.console()
}
read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}

completion <- read_json(file.path(source_dir, "completion.json"),
  simplifyVector = TRUE)
stopifnot(identical(completion$status, "complete"),
  completion$passing_signatures == 38L,
  completion$defined_pt_cells == 12642L,
  completion$missing_pt_cells == 14081L,
  completion$max_rank == 1500L)

pt_table <- fread(file.path(source_dir, "source_pseudotime.tsv"))
setnames(pt_table, c("cell_id", "pt"))
stopifnot(nrow(pt_table) == 12642L, !anyDuplicated(pt_table$cell_id),
  all(is.finite(pt_table$pt)))
setorder(pt_table, cell_id)
pt <- pt_table$pt

log_event("Reading retained native scWorkbench UCell scores")
scores_long <- read_tsv_gz(file.path(source_dir, "ucell_cell_scores.tsv.gz"))[
  , .(key, cell_id, value)]
scores_long <- scores_long[cell_id %chin% pt_table$cell_id]
stopifnot(uniqueN(scores_long$key) == 38L,
  uniqueN(scores_long$cell_id) == 12642L,
  nrow(scores_long) == 38L * 12642L,
  !anyDuplicated(scores_long, by = c("key", "cell_id")),
  all(is.finite(scores_long$value)), all(scores_long$value >= 0))
score_wide <- dcast(scores_long, cell_id ~ key, value.var = "value")
setorder(score_wide, cell_id)
stopifnot(identical(score_wide$cell_id, pt_table$cell_id))
feature_sets <- sort(setdiff(names(score_wide), "cell_id"))
score_matrix <- t(as.matrix(score_wide[, ..feature_sets]))
rownames(score_matrix) <- feature_sets
rm(scores_long, score_wide); invisible(gc())

width_fraction <- 0.05
step_fraction <- 0.01
min_cells <- 20L
q_cutoff <- 0.01
merge_gap_fraction <- width_fraction
pt_limits <- range(pt)
pt_span <- diff(pt_limits)
centers <- pt_limits[[1L]] + seq(0, 1, by = step_fraction) * pt_span
window_width <- width_fraction * pt_span
windows <- data.table(grid_window_index = seq_along(centers),
  window_center = centers,
  window_left = pmax(pt_limits[[1L]], centers - window_width / 2),
  window_right = pmin(pt_limits[[2L]], centers + window_width / 2))
membership_index <- lapply(seq_len(nrow(windows)), function(i)
  which(pt >= windows$window_left[[i]] & pt <= windows$window_right[[i]]))
windows[, n_cells := lengths(membership_index)]
keep <- windows$n_cells >= min_cells
windows <- windows[keep]
membership_index <- membership_index[keep]
windows[, window_index := seq_len(.N)]
stopifnot(nrow(windows) == 97L, min(windows$n_cells) >= min_cells)

# Validate reconstructed medians against the durable native scWorkbench report.
native <- fread(file.path(source_dir, "sliding_window_summary.tsv"))[
  abs(width_fraction - 0.05) < 1e-12 & eligible_window == TRUE,
  .(feature_set, grid_window_index = window_index,
    native_window_center = pseudotime, native_n_cells = n_cells,
    native_median = score_median)]
stopifnot(nrow(native) == 38L * 97L,
  setequal(native$feature_set, feature_sets))

# A plateau contributes one peak at its midpoint, rather than one candidate per
# tied window. Endpoints may be peaks when they exceed their sole neighbour.
local_maxima_plateau <- function(values) {
  rr <- rle(values)
  ends <- cumsum(rr$lengths)
  starts <- ends - rr$lengths + 1L
  left <- c(-Inf, head(rr$values, -1L))
  right <- c(tail(rr$values, -1L), -Inf)
  keep <- rr$values >= left & rr$values >= right &
    (rr$values > left | rr$values > right)
  as.integer(floor((starts[keep] + ends[keep]) / 2))
}

analyze_one <- function(ii) {
  fs <- feature_sets[[ii]]
  scores <- as.numeric(score_matrix[ii, ])
  baseline <- median(scores)
  observed <- vapply(membership_index,
    function(idx) median(scores[idx]), numeric(1))
  effect <- observed - baseline
  peak_index <- local_maxima_plateau(observed)
  peak_index <- peak_index[effect[peak_index] > .Machine$double.eps]
  curve <- data.table(analysis_group = "cytokines", feature_set = fs,
    windows, ucell_median = observed, trajectory_median = baseline,
    ucell_effect = effect)
  if (!length(peak_index)) return(list(curve = curve, candidates = NULL,
    permutation = data.table(analysis_group = "cytokines", feature_set = fs,
      candidate_peaks = 0L, permutations_used = 0L,
      null_max_mean = NA_real_, null_max_sd = NA_real_,
      null_max_q95 = NA_real_, null_max_q99 = NA_real_,
      stop_reason = "no_positive_local_peak")))

  null_max <- numeric(n_permutations)
  batch_size <- 100L
  set.seed(base_seed + ii - 1L)
  for (start in seq.int(1L, n_permutations, by = batch_size)) {
    take <- min(batch_size, n_permutations - start + 1L)
    permuted <- replicate(take, sample(scores, length(scores), replace = FALSE))
    if (is.null(dim(permuted))) permuted <- matrix(permuted, ncol = 1L)
    batch_max <- rep(-Inf, take)
    for (idx in membership_index) {
      window_median <- matrixStats::colMedians(permuted[idx, , drop = FALSE])
      batch_max <- pmax(batch_max, window_median - baseline)
    }
    null_max[start:(start + take - 1L)] <- batch_max
  }
  null_sd <- sd(null_max)
  candidates <- data.table(analysis_group = "cytokines", feature_set = fs,
    peak_window = peak_index,
    peak_pseudotime = windows$window_center[peak_index],
    region_left = windows$window_left[peak_index],
    region_right = windows$window_right[peak_index],
    n_cells = windows$n_cells[peak_index],
    peak_ucell_median = observed[peak_index],
    trajectory_median = baseline,
    peak_ucell_effect = effect[peak_index],
    null_max_mean = mean(null_max), null_max_sd = null_sd,
    permutation_z = if (is.finite(null_sd) && null_sd > 0)
      (effect[peak_index] - mean(null_max)) / null_sd else NA_real_,
    exceedances = vapply(effect[peak_index],
      function(x) sum(null_max >= x), integer(1)),
    permutations_used = n_permutations)
  candidates[, p_region := (1 + exceedances) / (1 + permutations_used)]
  list(curve = curve, candidates = candidates,
    permutation = data.table(analysis_group = "cytokines", feature_set = fs,
      candidate_peaks = length(peak_index), permutations_used = n_permutations,
      null_max_mean = mean(null_max), null_max_sd = null_sd,
      null_max_q95 = unname(quantile(null_max, .95)),
      null_max_q99 = unname(quantile(null_max, .99)),
      stop_reason = "complete"))
}

workers_requested <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", "8")))
if (!is.finite(workers_requested) || workers_requested < 1L)
  stop("N_WORKERS must be a positive integer")
workers <- min(workers_requested, parallel::detectCores(), length(feature_sets))
log_event("Starting", n_permutations, "median-UCell permutations for",
  length(feature_sets), "cytokine sets with", workers, "workers")
results <- parallel::mclapply(seq_along(feature_sets), analyze_one,
  mc.cores = workers, mc.preschedule = TRUE)
if (any(vapply(results, inherits, logical(1), "try-error")))
  stop("At least one cytokine permutation worker failed")
curves <- rbindlist(lapply(results, `[[`, "curve"))
candidate_rows <- Filter(Negate(is.null), lapply(results, `[[`, "candidates"))
candidates <- if (length(candidate_rows)) rbindlist(candidate_rows) else
  data.table()
permutation_summary <- rbindlist(lapply(results, `[[`, "permutation"))

native_check <- merge(curves[, .(feature_set, grid_window_index,
  window_center, n_cells, ucell_median)], native,
  by = c("feature_set", "grid_window_index"), all = TRUE)
stopifnot(nrow(native_check) == 38L * 97L, !anyNA(native_check),
  max(abs(native_check$window_center - native_check$native_window_center)) < 1e-10,
  identical(native_check$n_cells, native_check$native_n_cells),
  max(abs(native_check$ucell_median - native_check$native_median)) < 1e-12)

if (nrow(candidates)) {
  candidates[, q_region := p.adjust(p_region, method = "BH")]
  candidates[, significant := q_region < q_cutoff]
  setorder(candidates, q_region, p_region, -peak_ucell_effect,
    feature_set, peak_pseudotime)
  candidates[, peak_id := sprintf("ucell_cytokine_peak_%04d", seq_len(.N))]
  significant <- candidates[significant == TRUE]
} else significant <- data.table()

merge_gap_pt <- merge_gap_fraction * pt_span
merge_peaks <- function(x) {
  x <- copy(x[order(region_left, region_right)])
  gid <- integer(nrow(x)); gid[[1L]] <- 1L
  if (nrow(x) > 1L) for (i in 2:nrow(x)) {
    gap <- x$region_left[[i]] - x$region_right[[i - 1L]]
    gid[[i]] <- gid[[i - 1L]] + as.integer(gap >= merge_gap_pt)
  }
  rbindlist(lapply(split(seq_len(nrow(x)), gid), function(idx) {
    g <- x[idx]
    representative <- order(g$q_region, g$p_region,
      -g$peak_ucell_effect)[[1L]]
    out <- copy(g[representative])
    out[, `:=`(region_left = min(g$region_left),
      region_right = max(g$region_right),
      n_constituent_peaks = nrow(g),
      constituent_peak_pts = paste(signif(g$peak_pseudotime, 7), collapse = ";"),
      constituent_peak_ids = paste(g$peak_id, collapse = ";"))]
    out
  }))
}
if (nrow(significant)) {
  final_regions <- rbindlist(lapply(split(seq_len(nrow(significant)),
    significant$feature_set), function(idx) merge_peaks(significant[idx])))
  setorder(final_regions, q_region, p_region, -peak_ucell_effect,
    feature_set, region_left)
  final_regions[, region_id := sprintf("ucell_cytokine_region_%04d",
    seq_len(.N))]
  significant_sets <- final_regions[, .(
    short_name = short_label(feature_set),
    n_regions = .N,
    best_q = min(q_region),
    best_p = min(p_region),
    strongest_permutation_z = max(permutation_z, na.rm = TRUE),
    strongest_ucell_effect = max(peak_ucell_effect),
    peak_pseudotimes = paste(signif(peak_pseudotime[order(region_left)], 6),
      collapse = ";")), by = feature_set]
  setorder(significant_sets, best_q, best_p, -strongest_permutation_z,
    feature_set)
} else {
  final_regions <- data.table()
  significant_sets <- data.table()
}

write_tsv_gz(curves, file.path(stage_dir, "window_curves.tsv.gz"))
write_tsv_gz(candidates, file.path(stage_dir, "all_candidate_peaks.tsv.gz"))
write_tsv(significant,
  file.path(stage_dir, "significant_individual_peaks.tsv"))
write_tsv(final_regions,
  file.path(stage_dir, "significant_positive_regions.tsv"))
write_tsv(significant_sets,
  file.path(stage_dir, "significant_gene_sets.tsv"))
write_tsv(permutation_summary,
  file.path(stage_dir, "permutation_summary.tsv"))

if (nrow(final_regions)) {
  plot_order <- rev(significant_sets$feature_set)
  plot_regions <- copy(final_regions)
  plot_regions[, y := match(feature_set, plot_order)]
  q_floor <- 1 / (n_permutations + 1)
  p <- ggplot(plot_regions) +
    geom_segment(aes(x = region_left, xend = region_right, y = y, yend = y,
      color = -log10(pmax(q_region, q_floor))), linewidth = 5,
      lineend = "butt") +
    geom_point(aes(peak_pseudotime, y), shape = 21, size = 2.1,
      fill = "white", color = "#7F0000", stroke = .7) +
    scale_color_gradient(low = "#FCAE91", high = "#99000D",
      name = expression(-log[10](q))) +
    scale_y_continuous(breaks = seq_along(plot_order),
      labels = short_label(plot_order), expand = expansion(add = .6)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = "Significant positive cytokine UCell regions",
      subtitle = paste0("1,000 pseudotime-label permutations; BH q < 0.01; ",
        "no UCell magnitude or prominence cutoff"),
      x = "Existing Bdev pseudotime (pt)", y = "Cytokine gene set") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"))
  plot_height <- max(5, 1.8 + .25 * nrow(significant_sets))
  ggsave(file.path(stage_dir, "significant_ucell_regions.png"), p,
    width = 11, height = plot_height, dpi = 220, bg = "white")
  ggsave(file.path(stage_dir, "significant_ucell_regions.pdf"), p,
    width = 11, height = plot_height, device = cairo_pdf, bg = "white")
}

fmt_set <- function(i) paste0("| ", significant_sets$short_name[[i]], " | ",
  significant_sets$n_regions[[i]], " | ",
  signif(significant_sets$best_p[[i]], 4), " | ",
  signif(significant_sets$best_q[[i]], 4), " | ",
  signif(significant_sets$strongest_permutation_z[[i]], 4), " | ",
  significant_sets$peak_pseudotimes[[i]], " |")
set_rows <- if (nrow(significant_sets))
  vapply(seq_len(nrow(significant_sets)), fmt_set, character(1)) else
  "| None | 0 | NA | NA | NA | NA |"
report <- c("# Cytokine UCell permutation analysis", "",
  paste0("- Significant gene sets: ", nrow(significant_sets), " / 38"),
  paste0("- Significant individual peaks: ", nrow(significant)),
  paste0("- Final merged regions: ", nrow(final_regions)),
  "- Score: native scWorkbench UCell on RNA counts (`max_rank=1500`, no smoothing).",
  "- Statistic: positive deviation of the 5%-window median UCell score from the full-trajectory median.",
  "- Null: shuffle the same cells' UCell scores across fixed pseudotime positions; compare each peak with the maximum positive window deviation in each permutation.",
  "- Windows: 5% width, 1% step, minimum 20 cells; 97 retained windows.",
  "- Inference: 1,000 permutations and BH correction across all cytokine candidate peaks; q < 0.01.",
  "- Multiple peaks are retained. Significant peak windows merge only when the uncovered gap is strictly less than one 5% window.",
  "- No raw-UCell magnitude, z-score, relative-prominence, or gene-contribution cutoff was applied.", "",
  "## Significant gene sets", "",
  "| Gene set | Regions | Best p | Best q | Strongest permutation z | Peak pseudotime(s) |",
  "|---|---:|---:|---:|---:|---|", set_rows)
writeLines(report, file.path(stage_dir, "REPORT.md"))

write_json(list(status = "complete",
  analysis = "standalone native-UCell trajectory permutation",
  tested_gene_sets = 38L,
  significant_gene_sets = nrow(significant_sets),
  significant_individual_peaks = nrow(significant),
  significant_positive_regions = nrow(final_regions),
  score_source = "native scWorkbench UCell RNA_counts max_rank=1500",
  window_statistic = "median UCell",
  test_statistic = "window median minus full-trajectory median",
  permutation_null = "maximum positive window deviation across trajectory",
  permutations = n_permutations, multiple_testing =
    "BH across all cytokine candidate peaks", q_cutoff = q_cutoff,
  magnitude_or_prominence_filter = FALSE,
  window_width_fraction = width_fraction,
  window_step_fraction = step_fraction,
  minimum_window_cells = min_cells,
  retained_windows = nrow(windows),
  merge_gap_fraction_strictly_less_than = merge_gap_fraction,
  defined_pt_cells = nrow(pt_table),
  missing_pt_cells_preserved = completion$missing_pt_cells,
  native_window_medians_reproduced = TRUE,
  random_seed = base_seed, workers = workers,
  source_ucell_report = source_dir, promoted = FALSE),
  file.path(stage_dir, "completion.json"), pretty = TRUE, auto_unbox = TRUE,
  digits = NA)

stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed cytokine UCell permutation report at", output_dir)
