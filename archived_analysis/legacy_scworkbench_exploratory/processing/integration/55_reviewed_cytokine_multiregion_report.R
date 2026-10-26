#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(7L, 10L, 11L)) stop(paste(
  "Usage: 55_reviewed_cytokine_multiregion_report.R",
  "<bundle_dir> <inventory.tsv> <membership.tsv> <validated_report_dir>",
  "<output_dir> <n_permutations> <seed>",
  "[both|positive_only] [cluster_alpha] [region_q_cutoff]",
  "[cluster|peak_window]"))

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
inventory_path <- normalizePath(args[[2L]], mustWork = TRUE)
membership_path <- normalizePath(args[[3L]], mustWork = TRUE)
validated_dir <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
n_permutations <- as.integer(args[[6L]])
permutation_seed <- as.integer(args[[7L]])
direction_mode <- if (length(args) >= 10L) args[[8L]] else "both"
cluster_alpha <- if (length(args) >= 10L) as.numeric(args[[9L]]) else 0.05
region_q_cutoff <- if (length(args) >= 10L) as.numeric(args[[10L]]) else 0.05
inference_mode <- if (length(args) == 11L) args[[11L]] else "cluster"
stopifnot(n_permutations > 0L, is.finite(permutation_seed),
  direction_mode %in% c("both", "positive_only"),
  inference_mode %in% c("cluster", "peak_window"),
  is.finite(cluster_alpha), cluster_alpha > 0, cluster_alpha < 0.5,
  is.finite(region_q_cutoff), region_q_cutoff > 0, region_q_cutoff < 1)
if (inference_mode == "peak_window" && direction_mode != "positive_only")
  stop("peak_window inference currently requires positive_only mode")
tested_directions <- if (direction_mode == "positive_only") "positive" else
  c("positive", "negative")
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)
figure_dir <- file.path(stage_dir, "figures")
dir.create(figure_dir)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
  library(splines)
})

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt"); on.exit(close(con), add = TRUE)
  write.table(as.data.frame(x), con, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ..., "\n")
  flush.console()
}
as_flag <- function(x) isTRUE(x) || identical(tolower(as.character(x)), "true")
slugify <- function(x) {
  y <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("(^_+|_+$)", "", y)
}

inventory <- fread(inventory_path)
members <- fread(membership_path)
stopifnot(nrow(inventory) == 38L, !anyDuplicated(inventory$signature_id),
  setequal(unique(members$signature_id), inventory$signature_id))
validated <- read_json(file.path(validated_dir, "completion.json"),
  simplifyVector = TRUE)
stopifnot(identical(validated$status, "complete"),
  identical(as.integer(validated$tested_sets), 38L),
  as_flag(validated$native_validation$passed),
  identical(as.integer(validated$n_permutations), n_permutations))
validated_curve_connection <- gzfile(file.path(validated_dir,
  "native_trajectory_curves.tsv.gz"), "rt")
validated_curves <- as.data.table(read.delim(validated_curve_connection,
  check.names = FALSE, stringsAsFactors = FALSE))
close(validated_curve_connection)

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before,
  as.character(validated$bundle_fingerprint_after)))
cells <- as.character(wb_cells(bundle, "all"))
metadata <- as.data.table(wb_meta(bundle, cells = cells))
umap <- as.data.table(wb_space(bundle, "umap.atac", cells = cells))
stopifnot(length(cells) == 26723L,
  identical(as.character(metadata$cell_id), cells), "pt" %in% names(metadata),
  "nCount_RNA" %in% names(metadata),
  identical(as.character(umap$cell_id), cells),
  all(is.finite(umap$x)), all(is.finite(umap$y)),
  all(is.finite(metadata$nCount_RNA)), all(metadata$nCount_RNA > 0))
pt_all <- as.numeric(metadata$pt)
pt_defined <- is.finite(pt_all)
pt <- pt_all[pt_defined]
pt_cells <- cells[pt_defined]
stopifnot(length(pt) == 12642L, sum(!pt_defined) == 14081L)

# Match the earlier Bdev conditional RNA-depth sensitivity analysis: estimate
# the depth coefficient after allowing a flexible pseudotime trend, then remove
# only that depth component. The biological pseudotime component and residual
# are retained. Scaling is learned on the trajectory cells and applied to all
# cells so adjusted per-cell scores can also be saved and plotted.
depth_spline_df <- 8L
log_rna_depth_all <- log1p(as.numeric(metadata$nCount_RNA))
depth_center <- mean(log_rna_depth_all[pt_defined])
depth_scale <- sd(log_rna_depth_all[pt_defined])
stopifnot(is.finite(depth_scale), depth_scale > 0)
depth_z_all <- (log_rna_depth_all - depth_center) / depth_scale
depth_z <- depth_z_all[pt_defined]
depth_design <- cbind(intercept = 1, ns_pt = ns(pt, df = depth_spline_df),
  log_RNA_depth_z = depth_z)
depth_qr <- qr(depth_design)
pt_only_qr <- qr(depth_design[, -ncol(depth_design), drop = FALSE])
partial_depth <- qr.resid(pt_only_qr, depth_z)

# Fixed analysis grid: 5%-wide windows, centers every 1% of the unchanged pt
# range. Edge windows are clipped, just as in the existing Bdev summaries.
width_fraction <- 0.05
step_fraction <- 0.01
min_cells <- 30L
minimum_relative_prominence <- 0.10
pt_limits <- range(pt)
pt_span <- diff(pt_limits)
centers <- pt_limits[[1L]] + seq(0, 1, by = step_fraction) * pt_span
window_width <- width_fraction * pt_span
windows <- data.table(window_index = seq_along(centers),
  window_center = centers,
  window_left = pmax(pt_limits[[1L]], centers - window_width / 2),
  window_right = pmin(pt_limits[[2L]], centers + window_width / 2))
windows[, edge_truncated := (window_right - window_left) < window_width - 1e-12]
membership_index <- lapply(seq_len(nrow(windows)), function(i)
  which(pt >= windows$window_left[[i]] & pt <= windows$window_right[[i]]))
windows[, n_cells := lengths(membership_index)]
windows[, eligible := n_cells >= min_cells]
stopifnot(any(windows$eligible))
membership_index <- membership_index[windows$eligible]
windows <- windows[eligible == TRUE]
windows[, grid_window_index := window_index]
windows[, window_index := seq_len(.N)]
stopifnot(all(vapply(membership_index,
  function(i) length(i) >= min_cells, logical(1))))
weight_matrix <- matrix(0, nrow(windows), length(pt))
for (i in seq_len(nrow(windows)))
  weight_matrix[i, membership_index[[i]]] <- 1 / length(membership_index[[i]])

cluster_threshold <- qnorm(1 - cluster_alpha)
cluster_runs <- function(z, direction) {
  active <- if (direction == "positive") z >= cluster_threshold else
    z <= -cluster_threshold
  active[!is.finite(active)] <- FALSE
  rr <- rle(active)
  ends <- cumsum(rr$lengths)
  starts <- ends - rr$lengths + 1L
  Map(function(a, b) seq.int(a, b), starts[rr$values], ends[rr$values])
}
cluster_mass <- function(z, idx, direction) {
  if (direction == "positive") sum(z[idx] - cluster_threshold) else
    sum(-z[idx] - cluster_threshold)
}
max_cluster_mass <- function(z) {
  masses <- unlist(lapply(tested_directions, function(direction) {
    runs <- cluster_runs(z, direction)
    vapply(runs, function(i) cluster_mass(z, i, direction), numeric(1))
  }), use.names = FALSE)
  if (length(masses)) max(masses) else 0
}
local_extrema <- function(z, idx, direction) {
  if (length(idx) == 1L) return(idx)
  s <- if (direction == "positive") z else -z
  keep <- vapply(idx, function(i) {
    left_ok <- i == idx[[1L]] || s[[i]] >= s[[i - 1L]]
    right_ok <- i == idx[[length(idx)]] || s[[i]] >= s[[i + 1L]]
    left_ok && right_ok
  }, logical(1))
  idx[keep]
}
# Standard topographic prominence on the complete ordered score curve. Peak
# discovery does not use this quantity: it is attached for the later, separate
# effect-size filtering stage. Searching across sparse-grid gaps prevents a
# segment endpoint on a continuing slope from acquiring artificial prominence.
peak_prominence <- function(values, index) {
  peak <- values[[index]]
  left_higher <- if (index > 1L)
    which(values[seq.int(index - 1L, 1L)] > peak) else integer()
  left_boundary <- if (length(left_higher)) index - left_higher[[1L]] else 1L
  right_higher <- if (index < length(values))
    which(values[seq.int(index + 1L, length(values))] > peak) else integer()
  right_boundary <- if (length(right_higher)) index + right_higher[[1L]] else
    length(values)
  left_base <- min(values[seq.int(left_boundary, index)])
  right_base <- min(values[seq.int(index, right_boundary)])
  peak - max(left_base, right_base)
}

raw_score_matrix <- matrix(NA_real_, nrow = 38L, ncol = length(cells),
  dimnames = list(inventory$signature_id, cells))
score_matrix <- matrix(NA_real_, nrow = 38L, ncol = length(cells),
  dimnames = list(inventory$signature_id, cells))
window_rows <- vector("list", 38L)
region_rows <- list()
depth_diagnostic_rows <- vector("list", 38L)
validation_differences <- numeric(38L)
region_cursor <- 0L

log_event("Computing regional permutation tests for 38 reviewed cytokine signatures")
for (ii in seq_len(nrow(inventory))) {
  feature_set <- inventory$signature_id[[ii]]
  genes <- members[signature_id == feature_set, feature_id]
  expression <- as.matrix(wb_matrix_read_features(bundle,
    "RNA_log1p_cp10k", genes, cells, missing = "NA"))
  raw_scores <- colMeans(expression, na.rm = TRUE)
  stopifnot(length(raw_scores) == length(cells), all(is.finite(raw_scores)))
  raw_score_matrix[feature_set, ] <- raw_scores
  raw_pt_scores <- raw_scores[pt_defined]
  depth_coefficients <- qr.coef(depth_qr, raw_pt_scores)
  depth_beta <- unname(depth_coefficients[[length(depth_coefficients)]])
  stopifnot(is.finite(depth_beta))
  scores <- raw_scores - depth_beta * depth_z_all
  pt_scores <- scores[pt_defined]
  score_matrix[feature_set, ] <- scores
  raw_partial_score <- qr.resid(pt_only_qr, raw_pt_scores)
  adjusted_partial_score <- qr.resid(pt_only_qr, pt_scores)
  depth_diagnostic_rows[[ii]] <- data.table(feature_set = feature_set,
    cytokine = inventory$cytokine[[ii]],
    signature_direction = inventory$direction[[ii]],
    depth_beta_per_sd_log1p_nCount_RNA = depth_beta,
    conditional_depth_effect_score_sd = abs(depth_beta) / sd(raw_pt_scores),
    marginal_depth_correlation_raw = cor(raw_pt_scores,
      log_rna_depth_all[pt_defined]),
    marginal_depth_correlation_adjusted = cor(pt_scores,
      log_rna_depth_all[pt_defined]),
    partial_depth_correlation_raw = cor(raw_partial_score, partial_depth),
    partial_depth_correlation_adjusted = cor(adjusted_partial_score,
      partial_depth))
  observed <- as.numeric(weight_matrix %*% pt_scores)
  raw_observed <- as.numeric(weight_matrix %*% raw_pt_scores)

  # Confirm that score construction reproduces the already native-validated
  # 30-bin trajectory before adding the finer regional inference layer.
  bins30 <- as.integer(cut(pt, breaks = 30L, include.lowest = TRUE,
    labels = FALSE))
  check_curve <- vapply(seq_len(30L), function(b)
    mean(raw_pt_scores[bins30 == b]), numeric(1))
  target_id <- feature_set
  old_curve <- validated_curves[validated_curves$feature_set == target_id][order(bin)]
  stopifnot(nrow(old_curve) >= 5L,
    all(old_curve$bin %in% seq_len(30L)))
  validation_differences[[ii]] <- max(abs(
    check_curve[old_curve$bin] - old_curve$value))

  set.seed(permutation_seed + ii - 1L)
  null_curves <- matrix(NA_real_, nrow(windows), n_permutations)
  batch_size <- 100L
  for (start in seq.int(1L, n_permutations, by = batch_size)) {
    take <- min(batch_size, n_permutations - start + 1L)
    permuted <- replicate(take, sample(pt_scores, length(pt_scores),
      replace = FALSE))
    if (is.null(dim(permuted))) permuted <- matrix(permuted, ncol = 1L)
    null_curves[, start:(start + take - 1L)] <- weight_matrix %*% permuted
  }
  null_mean <- rowMeans(null_curves)
  null_sd <- apply(null_curves, 1L, sd)
  stopifnot(all(is.finite(null_sd)), all(null_sd > 0))
  observed_z <- (observed - null_mean) / null_sd
  null_z <- sweep(sweep(null_curves, 1L, null_mean, "-"), 1L, null_sd, "/")
  null_max <- if (inference_mode == "peak_window")
    apply(null_z, 2L, max) else apply(null_z, 2L, max_cluster_mass)

  window_rows[[ii]] <- data.table(feature_set = feature_set,
    cytokine = inventory$cytokine[[ii]],
    signature_direction = inventory$direction[[ii]], windows,
    score_mean = observed, raw_score_mean = raw_observed,
    null_mean = null_mean, null_sd = null_sd,
    z_score = observed_z,
    positive_cluster_forming = observed_z >= cluster_threshold,
    negative_cluster_forming = observed_z <= -cluster_threshold)

  if (inference_mode == "peak_window") {
    # Peak locations describe the enrichment curve itself, not its standardized
    # test statistic. Work within contiguous portions of the original 1% grid
    # so unsupported gaps are never treated as adjacent windows.
    segment_id <- cumsum(c(TRUE, diff(windows$grid_window_index) != 1L))
    segments <- split(seq_along(observed), segment_id)
    peak_indices <- unlist(lapply(segments, function(idx)
      local_extrema(observed, idx, "positive")), use.names = FALSE)
    peak_indices <- peak_indices[observed_z[peak_indices] >= cluster_threshold]
    runs_by_direction <- list(positive = lapply(peak_indices, function(i) i))
  } else {
    runs_by_direction <- setNames(lapply(tested_directions,
      function(direction) cluster_runs(observed_z, direction)), tested_directions)
  }
  for (regional_change in names(runs_by_direction)) {
    for (idx in runs_by_direction[[regional_change]]) {
      region_cursor <- region_cursor + 1L
      mass <- if (inference_mode == "peak_window") observed_z[[idx]] else
        cluster_mass(observed_z, idx, regional_change)
      peak_index <- if (regional_change == "positive")
        idx[[which.max(observed_z[idx])]] else idx[[which.min(observed_z[idx])]]
      extrema <- local_extrema(observed_z, idx, regional_change)
      if (inference_mode == "peak_window") {
        left_boundary <- windows$window_left[[peak_index]]
        right_boundary <- windows$window_right[[peak_index]]
      } else {
        left_boundary <- if (min(idx) > 1L) mean(windows$window_center[
          c(min(idx) - 1L, min(idx))]) else max(pt_limits[[1L]],
          windows$window_center[[min(idx)]] - step_fraction * pt_span / 2)
        right_boundary <- if (max(idx) < nrow(windows)) mean(windows$window_center[
          c(max(idx), max(idx) + 1L)]) else min(pt_limits[[2L]],
          windows$window_center[[max(idx)]] + step_fraction * pt_span / 2)
      }
      region_rows[[region_cursor]] <- data.table(
        feature_set = feature_set,
        cytokine = inventory$cytokine[[ii]],
        signature_direction = inventory$direction[[ii]],
        regional_change = regional_change,
        first_window = min(idx), last_window = max(idx),
        n_windows = length(idx),
        region_left = left_boundary,
        region_right = right_boundary,
        first_center = min(windows$window_center[idx]),
        last_center = max(windows$window_center[idx]),
        peak_window = peak_index,
        peak_pseudotime = windows$window_center[[peak_index]],
        peak_score = observed[[peak_index]],
        raw_peak_score = raw_observed[[peak_index]],
        peak_prominence = if (inference_mode == "peak_window")
          peak_prominence(observed, peak_index) else NA_real_,
        relative_peak_prominence = if (inference_mode == "peak_window")
          peak_prominence(observed, peak_index) / diff(range(observed)) else
          NA_real_,
        peak_excess = observed[[peak_index]] - null_mean[[peak_index]],
        peak_z = observed_z[[peak_index]],
        cluster_mass = mass,
        n_descriptive_local_extrema = length(extrema),
        descriptive_local_extrema_pt = paste(signif(
          windows$window_center[extrema], 7), collapse = ";"),
        p_region = (1 + sum(null_max >= mass)) / (1 + n_permutations))
    }
  }
  if (ii %% 5L == 0L || ii == 38L)
    log_event("Regional tests completed", ii, "of 38")
}
stopifnot(max(validation_differences) <= 1e-12)

window_statistics <- rbindlist(window_rows)
regions <- if (length(region_rows)) rbindlist(region_rows) else data.table()
if (nrow(regions)) {
  regions[, q_region := p.adjust(p_region, method = "BH")]
  setorder(regions, q_region, p_region, feature_set, region_left)
  regions[, region_id := sprintf("region_%03d", seq_len(.N))]
  regions[, significant := q_region < region_q_cutoff]
  regions[, passes_prominence := if (inference_mode == "peak_window")
    is.finite(relative_peak_prominence) &
      relative_peak_prominence >= minimum_relative_prominence else TRUE]
  regions[, retained_for_reporting := significant & passes_prominence]
} else {
  stop("No supra-threshold candidate regions were found")
}
significant_individual_peaks <- regions[significant == TRUE]
significant_constituent_peaks <- regions[retained_for_reporting == TRUE]
merge_gap_fraction <- width_fraction
merge_gap_pt <- merge_gap_fraction * pt_span
merge_significant_peaks <- function(x) {
  x <- copy(x[order(region_left, region_right)])
  merge_group <- integer(nrow(x)); merge_group[[1L]] <- 1L
  if (nrow(x) > 1L) for (i in 2:nrow(x)) {
    gap <- x$region_left[[i]] - x$region_right[[i - 1L]]
    merge_group[[i]] <- merge_group[[i - 1L]] + as.integer(gap >= merge_gap_pt)
  }
  groups <- split(seq_len(nrow(x)), merge_group)
  rbindlist(lapply(groups, function(idx) {
    g <- x[idx]
    representative <- which.max(g$peak_score)
    out <- copy(g[representative])
    out[, `:=`(
      first_window = min(g$first_window), last_window = max(g$last_window),
      n_windows = max(g$last_window) - min(g$first_window) + 1L,
      region_left = min(g$region_left), region_right = max(g$region_right),
      first_center = min(g$first_center), last_center = max(g$last_center),
      p_region = max(g$p_region), q_region = max(g$q_region),
      n_descriptive_local_extrema = nrow(g),
      descriptive_local_extrema_pt = paste(signif(g$peak_pseudotime, 7),
        collapse = ";"),
      n_constituent_peaks = nrow(g),
      constituent_peak_pts = paste(signif(g$peak_pseudotime, 7), collapse = ";"),
      constituent_region_ids = paste(g$region_id, collapse = ";"))]
    out
  }))
}
if (inference_mode == "peak_window" && nrow(significant_constituent_peaks)) {
  merge_key <- paste(significant_constituent_peaks$feature_set,
    significant_constituent_peaks$regional_change, sep = "\r")
  significant_regions <- rbindlist(lapply(
    split(seq_len(nrow(significant_constituent_peaks)), merge_key), function(idx)
      merge_significant_peaks(significant_constituent_peaks[idx])))
  setorder(significant_regions, q_region, p_region, feature_set, region_left)
  significant_regions[, region_id := sprintf("merged_region_%03d", seq_len(.N))]
} else {
  significant_regions <- copy(significant_constituent_peaks)
  significant_regions[, `:=`(n_constituent_peaks = 1L,
    constituent_peak_pts = as.character(signif(peak_pseudotime, 7)),
    constituent_region_ids = region_id)]
}
positive_regions <- significant_regions[regional_change == "positive"]
negative_regions <- significant_regions[regional_change == "negative"]

# Assign significant region IDs back to windows for auditability.
window_statistics[, `:=`(significant_positive_region = NA_character_,
  significant_negative_region = NA_character_)]
for (i in seq_len(nrow(significant_regions))) {
  r <- significant_regions[i]
  column <- if (r$regional_change == "positive")
    "significant_positive_region" else "significant_negative_region"
  rows <- window_statistics$feature_set == r$feature_set &
    window_statistics$window_index >= r$first_window &
    window_statistics$window_index <= r$last_window
  set(window_statistics, which(rows), column, r$region_id)
}

write_tsv_gz(window_statistics,
  file.path(stage_dir, "window_statistics.tsv.gz"))
write_tsv(regions, file.path(stage_dir, "all_candidate_regions.tsv"))
write_tsv(significant_individual_peaks,
  file.path(stage_dir, "significant_individual_peaks_unfiltered.tsv"))
write_tsv(significant_constituent_peaks,
  file.path(stage_dir, "prominence_retained_constituent_peaks.tsv"))
depth_diagnostics <- rbindlist(depth_diagnostic_rows)
write_tsv(depth_diagnostics,
  file.path(stage_dir, "rna_depth_adjustment_diagnostics.tsv"))
write_tsv(significant_regions,
  file.path(stage_dir, "significant_regions.tsv"))
write_tsv(positive_regions,
  file.path(stage_dir, "significant_positive_regions.tsv"))
write_tsv(negative_regions,
  file.path(stage_dir, "significant_negative_regions.tsv"))
saveRDS(list(cell_id = cells, pt = pt_all, raw_scores = raw_score_matrix,
  adjusted_scores = score_matrix, scores = score_matrix,
  matrix = "RNA_log1p_cp10k", score_method = "mean",
  depth_adjustment = paste0("subtract conditional log1p(nCount_RNA) ",
    "coefficient after ns(pt,df=", depth_spline_df, ")")),
  file.path(stage_dir, "cell_scores.rds"), compress = "xz")

base_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"))
sets_with_positive_peak <- unique(positive_regions$feature_set)
figure_index <- list()

for (ii in seq_along(sets_with_positive_peak)) {
  feature_set <- sets_with_positive_peak[[ii]]
  target_id <- feature_set
  inv <- inventory[inventory$signature_id == target_id]
  pos <- positive_regions[positive_regions$feature_set == target_id][order(region_left)]
  neg <- negative_regions[negative_regions$feature_set == target_id][order(region_left)]
  curve <- window_statistics[window_statistics$feature_set == target_id]
  positive_mask <- rep(FALSE, length(cells))
  for (j in seq_len(nrow(pos)))
    positive_mask <- positive_mask |
      (is.finite(pt_all) & pt_all >= pos$region_left[[j]] &
        pt_all <= pos$region_right[[j]])
  negative_mask <- rep(FALSE, length(cells))
  for (j in seq_len(nrow(neg)))
    negative_mask <- negative_mask |
      (is.finite(pt_all) & pt_all >= neg$region_left[[j]] &
        pt_all <= neg$region_right[[j]])
  stopifnot(!any(positive_mask & negative_mask))
  score <- score_matrix[feature_set, ]
  score_limits <- as.numeric(quantile(score, c(.01, .99), names = FALSE))
  if (diff(score_limits) == 0)
    score_limits <- range(score) + c(-1, 1) * .Machine$double.eps^0.5
  scaled_score <- pmin(1, pmax(0,
    (score - score_limits[[1L]]) / diff(score_limits)))
  red_palette <- grDevices::colorRampPalette(c("#FECACA", "#B91C1C"))(256L)
  # Lower scores represent stronger depletion in negative regions, so the
  # blue scale is darkest at its low-score end.
  blue_palette <- grDevices::colorRampPalette(c("#1D4ED8", "#BFDBFE"))(256L)
  color_index <- pmax(1L, pmin(256L, 1L + floor(scaled_score * 255L)))
  d <- data.table(x = umap$x, y = umap$y, score = score,
    scaled_score = scaled_score,
    status = fifelse(!is.finite(pt_all), "not_in_trajectory",
      fifelse(positive_mask, "significant_positive",
        fifelse(negative_mask, "significant_negative", "not_significant"))))
  d[, point_color := fifelse(status == "not_in_trajectory", "#CDD2D8",
    fifelse(status == "not_significant", "#AEB5BD",
      fifelse(status == "significant_positive", red_palette[color_index],
        blue_palette[color_index])))]

  umap_plot <- ggplot() +
    geom_point(data = d[status == "not_in_trajectory"], aes(x, y),
      color = "#CDD2D8", size = .17, alpha = .90, stroke = 0) +
    geom_point(data = d[status == "not_significant"], aes(x, y),
      color = "#AEB5BD", size = .19, alpha = .92, stroke = 0) +
    geom_point(data = d[status == "significant_negative"][order(-scaled_score)],
      aes(x, y, color = point_color), size = .48, alpha = .96, stroke = 0) +
    geom_point(data = d[status == "significant_positive"][order(scaled_score)],
      aes(x, y, color = point_color), size = .48, alpha = .96, stroke = 0) +
    scale_color_identity() +
    coord_equal() + theme_void(base_size = 10) +
    theme(legend.position = "none", plot.title = element_text(face = "bold")) +
    labs(title = paste(inv$cytokine, toupper(inv$direction),
        "depth-adjusted score on Bdev UMAP"),
      subtitle = paste(if (direction_mode == "positive_only")
          "Red: significant positive enrichment" else
          "Red: increase; blue: decrease",
        "Very light grey: pt missing; light grey: not significant", sep = "\n"),
      caption = if (direction_mode == "positive_only")
        "Darker red = higher depth-adjusted gene-set score." else
        "Dark red = higher adjusted score; dark blue = lower adjusted score.")

  trajectory_plot <- ggplot(curve, aes(window_center, score_mean)) +
    geom_rect(data = pos, aes(xmin = region_left, xmax = region_right,
      ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
      fill = "#D81B60", alpha = .13) +
    geom_rect(data = neg, aes(xmin = region_left, xmax = region_right,
      ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
      fill = "#1976D2", alpha = .12) +
    geom_line(aes(y = raw_score_mean), color = "#8B929A", linewidth = .65,
      linetype = "22") +
    geom_line(color = "#252525", linewidth = .8) +
    geom_point(data = pos, aes(x = peak_pseudotime, y = peak_score),
      inherit.aes = FALSE, color = "#D81B60", size = 2) +
    geom_point(data = neg, aes(x = peak_pseudotime, y = peak_score),
      inherit.aes = FALSE, color = "#1976D2", size = 2) +
    base_theme +
    labs(title = "Regional trajectory test",
      subtitle = if (direction_mode == "positive_only")
        "Black: depth-adjusted; dashed grey: raw; pink: retained region" else
        "Pink: significant positive; blue: significant negative",
      x = "Existing Bdev pseudotime (pt)",
      y = "Mean gene-set expression score",
      caption = paste0("5%-wide windows, 1% steps; q < ", region_q_cutoff,
        "; relative prominence >= ", minimum_relative_prominence, "."))

  composite <- (umap_plot | trajectory_plot) +
    plot_annotation(title = paste(inv$cytokine, toupper(inv$direction),
      "multi-region enrichment"), tag_levels = "A")
  filename <- paste0(sprintf("%02d", ii), "_", slugify(inv$cytokine), "_",
    slugify(inv$direction), "_significant_regions.png")
  ggsave(file.path(figure_dir, filename), composite, width = 12.5,
    height = 5.6, dpi = 220, bg = "white")
  figure_index[[ii]] <- data.table(feature_set = feature_set,
    cytokine = inv$cytokine, signature_direction = inv$direction,
    n_significant_positive_regions = nrow(pos),
    n_significant_negative_regions = nrow(neg),
    positive_cells = sum(positive_mask), negative_cells = sum(negative_mask),
    nonsignificant_trajectory_cells = sum(pt_defined & !positive_mask & !negative_mask),
    missing_pt_cells = sum(!pt_defined), total_cells = length(cells),
    figure = file.path("figures", filename))
}
figure_index <- if (length(figure_index)) rbindlist(figure_index) else data.table()
write_tsv(figure_index, file.path(stage_dir, "umap_figure_index.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
completion <- list(status = "complete",
  collection = "reviewed cytokine signatures", tested_sets = 38L,
  sets_with_significant_positive_regions = length(sets_with_positive_peak),
  significant_positive_regions = nrow(positive_regions),
  significant_negative_regions = nrow(negative_regions),
  significant_individual_peaks_before_prominence =
    nrow(significant_individual_peaks),
  prominence_retained_constituent_peaks = nrow(significant_constituent_peaks),
  n_permutations = n_permutations, permutation_seed = permutation_seed,
  matrix = "RNA_log1p_cp10k", score_method = "depth-adjusted mean",
  depth_adjustment = list(covariate = "scaled log1p(nCount_RNA)",
    conditional_model = paste0("score ~ ns(pt, df=", depth_spline_df,
      ") + scaled log1p(nCount_RNA)"),
    removed_component = "RNA-depth coefficient times scaled log-depth only",
    retained_components = "pseudotime spline plus residual",
    median_absolute_conditional_effect_score_sd = median(
      depth_diagnostics$conditional_depth_effect_score_sd),
    maximum_absolute_conditional_effect_score_sd = max(
      depth_diagnostics$conditional_depth_effect_score_sd)),
  pseudotime_source = "exact Bdev metadata pt",
  defined_pt_cells = length(pt), missing_pt_cells = sum(!pt_defined),
  windows = list(width_fraction = width_fraction,
    step_fraction = step_fraction, n_windows = nrow(windows),
    min_cells = min_cells, edge_policy = "clip to observed pt range"),
  inference = list(cluster_forming_one_sided_alpha = cluster_alpha,
    cluster_forming_z = cluster_threshold,
    inference_mode = inference_mode,
    test_statistic = if (inference_mode == "peak_window")
      paste("local maxima located on the depth-adjusted mean gene-set score; peak z tested",
        "against the permutation maximum z across windows") else
      "sum of z excess beyond the forming threshold",
    tested_directions = tested_directions,
    within_signature_control = if (inference_mode == "peak_window")
      "permutation maximum positive z across all pseudotime windows" else
      paste("permutation maximum cluster mass across",
        paste(tested_directions, collapse = " and "), "direction(s)"),
    collection_adjustment = "BH across all candidate regions in all 38 signatures",
    region_q_cutoff = region_q_cutoff,
    downstream_prominence_filter = if (inference_mode == "peak_window") list(
      relative_prominence_minimum = minimum_relative_prominence,
      denominator = "full eligible-window range of the adjusted score curve",
      order = "applied after peak discovery, permutation testing, and BH adjustment; before merging") else "none",
    reporting_merge = if (inference_mode == "peak_window") list(
      rule = "merge adjacent significant peak windows when uncovered gap is strictly less than one window",
      gap_fraction_strictly_less_than = merge_gap_fraction,
      gap_pt_units_strictly_less_than = merge_gap_pt,
      implication_on_1pct_grid = "0% through 4% gaps merge; 5% gaps do not",
      representative_peak = "retained constituent peak with highest depth-adjusted mean gene-set score",
      constituent_peak_statistics_preserved = TRUE) else "none",
    displayed_region_boundary = if (inference_mode == "peak_window")
      paste("the union of actual 5%-wide scoring windows centered on retained",
        "depth-adjusted-score peaks; unsupported grid gaps are separate during",
        "initial peak discovery") else
      paste("midpoints between adjacent tested window centers; keeps clusters",
        "separated despite overlapping window footprints")),
  umap_encoding = list(all_cells_visible = TRUE,
    missing_pt = "visible very light grey",
    trajectory_not_significant = "light grey, visibly darker than missing pt",
    significant_positive = "red; darker with higher gene-set score",
    significant_negative = if (direction_mode == "positive_only")
      "not tested" else "blue; darker with lower gene-set score"),
  validation = list(source_report = validated_dir,
    prior_native_validation_passed = TRUE,
    max_30_bin_score_difference = max(validation_differences)),
  promoted = FALSE, bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = fingerprint_after, bundle_unchanged = TRUE)
write_json(completion, file.path(stage_dir, "completion.json"))

region_table <- function(x) {
  if (!nrow(x)) return("No significant regions.")
  c("| Cytokine | Signature | Change | Region pt | Representative peak pt | Peaks merged | Relative prominence | Peak z | p | q |",
    "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
    vapply(seq_len(nrow(x)), function(i) paste0("| ", x$cytokine[[i]], " | ",
      x$signature_direction[[i]], " | ", x$regional_change[[i]], " | ",
      signif(x$region_left[[i]], 5), "-", signif(x$region_right[[i]], 5), " | ",
      signif(x$peak_pseudotime[[i]], 5), " | ", x$n_constituent_peaks[[i]],
      " | ", signif(x$relative_peak_prominence[[i]], 4),
      " | ", signif(x$peak_z[[i]], 4),
      " | ", signif(x$p_region[[i]], 4), " | ", signif(x$q_region[[i]], 4),
      " |"), character(1)))
}
figure_lines <- if (nrow(figure_index)) vapply(seq_len(nrow(figure_index)),
  function(i) paste0("- [", figure_index$cytokine[[i]], " ",
    toupper(figure_index$signature_direction[[i]]), "](",
    figure_index$figure[[i]], ") — ",
    figure_index$n_significant_positive_regions[[i]], " positive region(s), ",
    figure_index$n_significant_negative_regions[[i]], " negative region(s); ",
    figure_index$total_cells[[i]], " cells shown."), character(1)) else
  "No gene set had a significant positive region."
report <- c("# Reviewed cytokine multi-region trajectory analysis", "",
  paste0("The 38 reviewed cytokine gene sets were tested with ",
    n_permutations, " permutations each in ", direction_mode, " mode. ",
    "Inference mode: ", inference_mode, ". ",
    nrow(positive_regions),
    " positive and ", nrow(negative_regions),
    " negative final reporting regions remained after significance, prominence,",
    " and merging."), "",
  "## Interpretation", "",
  paste("A positive region means the depth-adjusted gene-set score is higher than expected",
    "under exchangeability along pseudotime; a negative region means it is lower."),
  if (inference_mode == "peak_window") paste(
    "Analysis is staged explicitly. First, every eligible local maximum on the",
    "depth-adjusted enrichment curve is discovered and permutation-tested; BH",
    "adjustment is then applied across all tested peaks and all 38 gene sets.") else paste(
    "Separated significant clusters are reported independently, allowing",
    "more than one enrichment stage for a gene set."),
  if (inference_mode == "peak_window") paste(
    "Second, q-significant individual peaks must have relative topographic",
    "prominence of at least 10% of that gene set's complete adjusted trajectory",
    "range. This is a downstream effect-size filter, not part of discovery or",
    "multiple-testing correction, and the number of retained peaks is not capped.") else paste(
    "Each region spans consecutive windows above the cluster-forming threshold."),
  if (inference_mode == "peak_window") paste(
    "Third, retained peak windows are merged for reporting when their uncovered",
    "gap is strictly less than one window (5% of the trajectory): 0% through 4%",
    "gaps merge, while 5% gaps do not. Unfiltered q-significant peaks and retained",
    "constituent peaks are preserved in separate audit tables.")
    else NULL,
  paste("Before peak analysis, each cell-level mean log1p(CP10K) score is modeled",
    "as ns(pt, df=8) plus scaled log1p(nCount_RNA). Only the fitted RNA-depth",
    "component is subtracted, preserving the pseudotime spline and residual.",
    "Figures show the adjusted curve in black and the raw curve in dashed grey."),
  paste("The UMAP encoding is region-level. All cells are shown: missing-pt cells",
    "are visible very light grey, trajectory cells outside significant regions",
    "are a noticeably darker light grey, and significant positive regions are red.",
    if (direction_mode == "positive_only") "Negative enrichment was not tested."
      else "Significant negative regions are blue.",
    "This is not a per-cell significance claim."), "",
  "## Significant regions", "", region_table(significant_regions), "",
  "## UMAP and trajectory figures", "", figure_lines, "",
  "## Statistical specification", "",
  paste0("Windows span 5% of the observed pt range and advance by 1% (",
    nrow(windows), " windows; minimum 30 cells). A one-sided z = ",
    signif(cluster_threshold, 5),
    if (direction_mode == "positive_only")
      " threshold defines eligible positive peak candidates."
    else " threshold forms positive and negative clusters."),
  if (inference_mode == "peak_window") paste(
    "For each permutation, the largest positive window z provides the null statistic.",
    "Observed local-peak p-values were BH-adjusted together across all 38 gene sets.")
    else paste("For each permutation, the largest cluster mass across tested directions",
      "provides the null statistic. Candidate-region p-values were BH-adjusted",
      "together across all 38 gene sets."),
  if (inference_mode == "peak_window") paste(
    "Reported region bounds are the union of retained 5%-wide peak windows after",
    "the separate 10% prominence filter and strict less-than-5% gap merge.")
    else paste("Reported region bounds lie halfway between adjacent tested window centers,",
      "so distinct clusters remain disjoint even though the underlying 5% windows overlap."),
  paste("The unadjusted score calculation reproduced the earlier native-validated",
    "30-bin curves with maximum absolute difference",
    signif(max(validation_differences), 4), "."))
writeLines(report, file.path(stage_dir, "REPORT.md"))

wb_close_bundle(bundle)
file.rename(stage_dir, output_dir)
log_event("Completed cytokine multi-region report at", output_dir)
