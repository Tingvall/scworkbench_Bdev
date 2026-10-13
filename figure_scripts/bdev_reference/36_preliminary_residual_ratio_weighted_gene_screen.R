#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 36_preliminary_residual_ratio_weighted_gene_screen.R",
  "<bdev_bundle> <window_report> <ratio_depth_report> <output_dir>",
  "<expected_bundle_fingerprint> <batch_size>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
ratio_depth_report <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_fingerprint <- args[[5L]]
batch_size <- as.integer(args[[6L]])
stopifnot(batch_size >= 16L)
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(Matrix)
  library(matrixStats)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
  library(splines)
})
read_tsv <- function(path) read.delim(path, check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = "NA")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con))
  write.table(x, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
row_cosine <- function(x, target, weights) {
  weights <- weights / sum(weights)
  numerator <- rowSums(sweep(x, 2L, weights * target, "*"))
  denominator <- sqrt(rowSums(sweep(x^2, 2L, weights, "*")) *
    sum(weights * target^2))
  answer <- numerator / denominator
  answer[!is.finite(answer)] <- 0
  answer
}
row_weighted_correlation <- function(x, target, weights) {
  weights <- weights / sum(weights)
  x_center <- rowSums(sweep(x, 2L, weights, "*"))
  y_center <- sum(weights * target)
  row_cosine(x - x_center, target - y_center, weights)
}
row_scale_invariant_direction <- function(x, target, weights) {
  gene_scale <- rowQuantiles(abs(x), probs = 0.75, drop = TRUE)
  gene_scale <- pmax(gene_scale, .Machine$double.eps)
  target_scale <- as.numeric(quantile(abs(target), 0.75,
    names = FALSE, type = 7))
  stopifnot(is.finite(target_scale), target_scale > 0)
  gene_direction <- tanh(x / gene_scale)
  target_direction <- tanh(target / target_scale)
  row_cosine(gene_direction, target_direction, weights)
}
row_pearson <- function(x, y) {
  yc <- y - mean(y)
  xc <- x - rowMeans(x)
  answer <- rowSums(sweep(xc, 2L, yc, "*")) /
    sqrt(rowSums(xc^2) * sum(yc^2))
  answer[!is.finite(answer)] <- 0
  answer
}
robust_unit <- function(x) {
  scale <- as.numeric(quantile(x, 0.95, names = FALSE, type = 7))
  stopifnot(is.finite(scale), scale > 0)
  pmin(x / scale, 1)
}

trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
windows <- read_tsv(file.path(window_report, "window_distance_and_slope.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
ratio <- read_tsv(file.path(ratio_depth_report,
  "ratio_total_signal_fitted_and_residual_windows.tsv"))
stopifnot(nrow(trajectory) == 12642L, nrow(windows) == 2524L,
  nrow(ratio) == nrow(windows), identical(windows$window_id, ratio$window_id),
  all(table(membership$window_id) == 30L),
  all(is.finite(ratio$full_residual)))

cell_index <- match(membership$cell_id, trajectory$cell_id)
window_index <- match(membership$window_id, windows$window_id)
stopifnot(!anyNA(cell_index), !anyNA(window_index))
window_design <- sparseMatrix(i = cell_index, j = window_index,
  x = rep(1 / 30, nrow(membership)),
  dims = c(nrow(trajectory), nrow(windows)))
membership_split <- split(cell_index, factor(membership$window_id,
  levels = windows$window_id))
smoothing_fraction <- 0.03
smoothing_half_width <- floor(nrow(windows) * smoothing_fraction / 2)
smoothing_rows <- unlist(lapply(seq_len(nrow(windows)), function(j)
  seq.int(max(1L, j - smoothing_half_width),
    min(nrow(windows), j + smoothing_half_width))), use.names = FALSE)
smoothing_columns <- rep(seq_len(nrow(windows)),
  vapply(seq_len(nrow(windows)), function(j)
    min(nrow(windows), j + smoothing_half_width) -
      max(1L, j - smoothing_half_width) + 1L, numeric(1)))
smoothing_counts <- tabulate(smoothing_columns, nbins = nrow(windows))
smoothing_design <- sparseMatrix(i = smoothing_rows, j = smoothing_columns,
  x = 1 / smoothing_counts[smoothing_columns],
  dims = c(nrow(windows), nrow(windows)))

# Apply exactly the same smoother to the corrected ratio and every gene. This
# prevents high-frequency target structure that no equally smoothed gene could
# reproduce from distorting the ranking. Association is then scale invariant.
target_raw <- ratio$full_residual
target <- as.numeric(t(target_raw) %*% smoothing_design)
target_derivative <- c(target[[2L]] - target[[1L]],
  (target[3:length(target)] - target[1:(length(target) - 2L)]) / 2,
  target[[length(target)]] - target[[length(target) - 1L]])
target_curvature <- c(target_derivative[[2L]] - target_derivative[[1L]],
  (target_derivative[3:length(target_derivative)] -
    target_derivative[1:(length(target_derivative) - 2L)]) / 2,
  target_derivative[[length(target_derivative)]] -
    target_derivative[[length(target_derivative) - 1L]])
magnitude_component <- robust_unit(abs(target))
derivative_component <- robust_unit(abs(target_derivative))
curvature_component <- robust_unit(abs(target_curvature))
interest <- 0.40 * magnitude_component + 0.40 * derivative_component +
  0.20 * curvature_component
baseline_weight <- 0.25
window_weight <- baseline_weight + (1 - baseline_weight) * interest
window_weight <- window_weight / sum(window_weight)
# Eight segments contain equal cumulative target-interest weight. The score
# compares normalized, clipped slopes, making it primarily an up/down pattern
# score rather than a magnitude score.
n_segments <- 8L
cumulative_interest <- cumsum(window_weight) - window_weight / 2
segment <- pmin(n_segments, floor(cumulative_interest * n_segments) + 1L)
stopifnot(identical(as.integer(sort(unique(segment))), seq_len(n_segments)))
weight_table <- data.frame(window_id = windows$window_id,
  window_index = windows$window_index, pt = windows$pt_median,
  corrected_ratio_raw = target_raw, corrected_ratio_matched_smooth = target,
  ratio_derivative = target_derivative,
  ratio_curvature = target_curvature,
  magnitude_component = magnitude_component,
  derivative_component = derivative_component,
  curvature_component = curvature_component,
  continuous_interest = interest, normalized_window_weight = window_weight,
  interest_balanced_segment = segment)
write_tsv(weight_table, file.path(stage_dir,
  "corrected_ratio_continuous_interest_weights.tsv"))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))
metadata <- wb_meta(bundle, cells = trajectory$cell_id,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt", "nCount_RNA"))
stopifnot(identical(as.character(metadata$cell_id), trajectory$cell_id),
  all(is.finite(metadata$nCount_RNA) & metadata$nCount_RNA > 0))
log_depth <- log1p(as.numeric(metadata$nCount_RNA))
depth_z <- as.numeric(scale(log_depth))
depth_design <- cbind(intercept = 1, ns_pt = ns(trajectory$pt, df = 8),
  log_RNA_depth_z = depth_z)
depth_solver <- solve(crossprod(depth_design))
features <- wb_feature_ids(bundle, "SCT_data")
stopifnot(length(features) > 10000L, !anyDuplicated(features))

minimum_detection_fraction <- 0.01
global_shape_weight <- 0.30
segment_shape_weight <- 0.70
segment_consistency_penalty <- 0.75
amplitude_half_saturation <- 0.25
results_list <- list()
batch_starts <- seq.int(1L, length(features), by = batch_size)
for (batch_number in seq_along(batch_starts)) {
  first <- batch_starts[[batch_number]]
  last <- min(first + batch_size - 1L, length(features))
  batch_features <- features[first:last]
  message(sprintf("Batch %d/%d: features %d-%d", batch_number,
    length(batch_starts), first, last))
  expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data",
    batch_features, trajectory$cell_id, missing = "NA"))
  stopifnot(identical(rownames(expression), batch_features),
    identical(colnames(expression), trajectory$cell_id),
    all(is.finite(expression)))
  detection_fraction <- rowMeans(expression > 0)
  keep <- detection_fraction >= minimum_detection_fraction
  if (!any(keep)) next
  ids <- batch_features[keep]
  expression <- expression[keep, , drop = FALSE]

  coefficients <- depth_solver %*% crossprod(depth_design, t(expression))
  depth_beta <- as.numeric(coefficients[nrow(coefficients), ])
  adjusted <- expression - tcrossprod(depth_beta, depth_z)
  depth_cor_before <- row_pearson(expression, log_depth)
  depth_cor_after <- row_pearson(adjusted, log_depth)
  conditional_depth_effect_sd <- abs(depth_beta) / pmax(rowSds(expression),
    .Machine$double.eps)

  profiles <- as.matrix((adjusted %*% window_design) %*% smoothing_design)
  derivatives <- (profiles[, 3:nrow(windows), drop = FALSE] -
    profiles[, 1:(nrow(windows) - 2L), drop = FALSE]) / 2
  derivative_indices <- 2:(nrow(windows) - 1L)
  target_derivative_center <- target_derivative[derivative_indices]
  derivative_weight <- window_weight[derivative_indices]
  global_profile_correlation <- row_weighted_correlation(profiles, target,
    window_weight)
  global_direction_similarity <- row_scale_invariant_direction(derivatives,
    target_derivative_center, derivative_weight)
  segment_scores <- do.call(cbind, lapply(seq_len(n_segments), function(k) {
    ix <- which(segment[derivative_indices] == k)
    row_scale_invariant_direction(derivatives[, ix, drop = FALSE],
      target_derivative_center[ix], derivative_weight[ix])
  }))
  segment_mean <- rowMeans(segment_scores)
  segment_sd <- rowSds(segment_scores)
  raw_direct <- global_shape_weight * global_direction_similarity +
    segment_shape_weight *
    (segment_mean - segment_consistency_penalty * segment_sd)
  raw_inverse <- -global_shape_weight * global_direction_similarity +
    segment_shape_weight *
    (-segment_mean - segment_consistency_penalty * segment_sd)

  profile_quantiles <- rowQuantiles(profiles, probs = c(0.05, 0.95),
    drop = FALSE)
  dynamic_range <- profile_quantiles[, 2L] - profile_quantiles[, 1L]
  raw_q25 <- raw_q75 <- raw_sd <- matrix(NA_real_, nrow = nrow(adjusted),
    ncol = nrow(windows))
  for (j in seq_len(nrow(windows))) {
    ix <- membership_split[[j]]
    q <- rowQuantiles(adjusted[, ix, drop = FALSE],
      probs = c(0.25, 0.75), drop = FALSE)
    raw_q25[, j] <- q[, 1L]
    raw_q75[, j] <- q[, 2L]
    raw_sd[, j] <- rowSds(adjusted[, ix, drop = FALSE])
  }
  smoothed_q25 <- as.matrix(raw_q25 %*% smoothing_design)
  smoothed_q75 <- as.matrix(raw_q75 %*% smoothing_design)
  median_iqr <- rowMedians(smoothed_q75 - smoothed_q25)
  median_sd <- rowMedians(raw_sd)
  conservative_spread <- pmax(median_iqr, 1.349 * median_sd)
  precision <- dynamic_range / (dynamic_range + conservative_spread)
  precision[!is.finite(precision)] <- 0
  amplitude <- dynamic_range / (dynamic_range + amplitude_half_saturation)
  amplitude[!is.finite(amplitude)] <- 0
  evidence_weight <- precision * amplitude
  result <- data.frame(feature_id = ids,
    detection_fraction = detection_fraction[keep], depth_beta = depth_beta,
    conditional_depth_effect_sd = conditional_depth_effect_sd,
    correlation_total_RNA_before = depth_cor_before,
    correlation_total_RNA_after = depth_cor_after,
    adjusted_dynamic_range_p95_p05 = dynamic_range,
    adjusted_median_smoothed_window_iqr = median_iqr,
    adjusted_median_window_sd = median_sd,
    precision_reliability = precision, amplitude_weight = amplitude,
    evidence_weight = evidence_weight,
    global_shape_correlation_direct = global_profile_correlation,
    global_direction_similarity_direct = global_direction_similarity,
    interest_segment_direction_mean_direct = segment_mean,
    interest_segment_direction_sd = segment_sd,
    raw_shape_score_direct = raw_direct,
    raw_shape_score_inverse = raw_inverse,
    evidence_weighted_score_direct = raw_direct * evidence_weight,
    evidence_weighted_score_inverse = raw_inverse * evidence_weight,
    stringsAsFactors = FALSE)
  for (k in seq_len(n_segments))
    result[[paste0("segment_", k, "_direction_similarity_direct")]] <-
      segment_scores[, k]
  results_list[[length(results_list) + 1L]] <- result
  rm(expression, adjusted, profiles, derivatives, raw_q25, raw_q75, raw_sd,
    smoothed_q25, smoothed_q75, segment_scores)
  if (batch_number %% 4L == 0L) gc(verbose = FALSE)
}

results <- do.call(rbind, results_list)
stopifnot(nrow(results) > 10000L, !anyDuplicated(results$feature_id),
  all(is.finite(results$evidence_weighted_score_direct)),
  all(is.finite(results$evidence_weighted_score_inverse)))
results$direct_rank <- rank(-results$evidence_weighted_score_direct,
  ties.method = "min")
results$inverse_rank <- rank(-results$evidence_weighted_score_inverse,
  ties.method = "min")
results <- results[order(results$feature_id), ]
write_tsv_gz(results, file.path(stage_dir,
  "residual_ratio_weighted_gene_scores.tsv.gz"))
top_direct <- head(results[order(results$direct_rank, results$feature_id), ],
  100L)
top_inverse <- head(results[order(results$inverse_rank, results$feature_id), ],
  100L)
write_tsv(top_direct, file.path(stage_dir,
  "residual_ratio_weighted_top100_direct.tsv"))
write_tsv(top_inverse, file.path(stage_dir,
  "residual_ratio_weighted_top100_inverse.tsv"))

# Re-read only displayed genes and repeat the identical depth adjustment.
top_n <- 20L
direct20 <- head(top_direct, top_n)
inverse20 <- head(top_inverse, top_n)
plot_features <- unique(c(direct20$feature_id, inverse20$feature_id))
plot_expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data",
  plot_features, trajectory$cell_id, missing = "NA"))
plot_coef <- depth_solver %*% crossprod(depth_design, t(plot_expression))
plot_beta <- as.numeric(plot_coef[nrow(plot_coef), ])
plot_adjusted <- plot_expression - tcrossprod(plot_beta, depth_z)
plot_mean <- as.matrix((plot_adjusted %*% window_design) %*% smoothing_design)
plot_q25 <- plot_q75 <- matrix(NA_real_, nrow = length(plot_features),
  ncol = nrow(windows))
for (j in seq_len(nrow(windows))) {
  q <- rowQuantiles(plot_adjusted[, membership_split[[j]], drop = FALSE],
    probs = c(0.25, 0.75), drop = FALSE)
  plot_q25[, j] <- q[, 1L]
  plot_q75[, j] <- q[, 2L]
}
plot_q25 <- as.matrix(plot_q25 %*% smoothing_design)
plot_q75 <- as.matrix(plot_q75 %*% smoothing_design)
rownames(plot_mean) <- rownames(plot_q25) <- rownames(plot_q75) <- plot_features
plot_center <- rowMeans(plot_mean)
plot_sd <- rowSds(plot_mean)
plot_mean_z <- (plot_mean - plot_center) / plot_sd
plot_q25_z <- (plot_q25 - plot_center) / plot_sd
plot_q75_z <- (plot_q75 - plot_center) / plot_sd
target_z <- as.numeric(scale(target))

make_plot_data <- function(top, direction) {
  rows <- lapply(seq_len(nrow(top)), function(i) {
    gene <- top$feature_id[[i]]
    data.frame(pt = windows$pt_median,
      panel = sprintf("%02d %s | score %.3f | depth beta %.3f",
        i, gene, top[[paste0("evidence_weighted_score_", direction)]][[i]],
        top$depth_beta[[i]]),
      mean_z = plot_mean_z[gene, ], q25_z = plot_q25_z[gene, ],
      q75_z = plot_q75_z[gene, ],
      target_z = if (direction == "direct") target_z else -target_z)
  })
  out <- do.call(rbind, rows)
  out$panel <- factor(out$panel, levels = unique(out$panel))
  out
}
make_top_plot <- function(top, direction) {
  color <- if (direction == "direct") "#C51B5A" else "#2166AC"
  label <- if (direction == "direct") "direct" else "inverse"
  ggplot(make_plot_data(top, direction), aes(pt)) +
    geom_hline(yintercept = 0, color = "#999999", linetype = "dotted",
      linewidth = 0.25) +
    geom_ribbon(aes(ymin = q25_z, ymax = q75_z), fill = color,
      alpha = 0.18, color = NA) +
    geom_line(aes(y = mean_z), color = color, linewidth = 0.8) +
    geom_line(aes(y = target_z), color = "#111111", linewidth = 0.55) +
    facet_wrap(~ panel, ncol = 4L) +
    coord_cartesian(ylim = c(-3.5, 3.5)) +
    theme_classic(base_size = 9) +
    theme(strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 7.2),
      panel.spacing = grid::unit(5, "mm"),
      plot.title = element_text(face = "bold", size = 15)) +
    labs(title = paste("Top 20", label,
        "depth-adjusted corrected-ratio associations"),
      subtitle = paste("Complete pseudotime is scored; black is the corrected",
        "ratio target; ribbon is adjusted-expression IQR"),
      x = "Bdev pseudotime", y = "Standardized profile")
}
direct_plot <- make_top_plot(direct20, "direct")
inverse_plot <- make_top_plot(inverse20, "inverse")
ggsave(file.path(stage_dir, "residual_ratio_weighted_top20_direct.png"),
  direct_plot, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "residual_ratio_weighted_top20_direct.pdf"),
  direct_plot, width = 14, height = 12, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "residual_ratio_weighted_top20_inverse.png"),
  inverse_plot, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "residual_ratio_weighted_top20_inverse.pdf"),
  inverse_plot, width = 14, height = 12, device = cairo_pdf, bg = "white")

target_data <- data.frame(pt = windows$pt_median, ratio = target,
  weight_relative = window_weight / min(window_weight),
  interest = interest)
target_panel <- ggplot(target_data, aes(pt, ratio)) +
  geom_hline(yintercept = 0, color = "#888888", linetype = "dotted") +
  geom_line(color = "#008837", linewidth = 0.9) +
  theme_classic(base_size = 11) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), plot.title = element_text(face = "bold")) +
  labs(title = "Corrected RNA/ATAC ratio", y = "Residual log2\nRNA / ATAC")
weight_panel <- ggplot(target_data, aes(pt, weight_relative)) +
  geom_area(fill = "#E6AB02", alpha = 0.45, position = "identity") +
  geom_line(color = "#A6761D", linewidth = 0.7) +
  geom_hline(yintercept = 1, color = "#777777", linetype = "dotted") +
  theme_classic(base_size = 11) +
  labs(title = "Continuous target-derived analysis weight",
    subtitle = "Every window contributes; magnitude, slope, and curvature increase weight",
    x = "Bdev pseudotime", y = "Weight relative\nto minimum")
weight_figure <- target_panel / weight_panel +
  plot_annotation(title = "Full-pseudotime weighting from the corrected ratio",
    theme = theme(plot.title = element_text(face = "bold", size = 15)))
ggsave(file.path(stage_dir, "corrected_ratio_continuous_weights.png"),
  weight_figure, width = 11, height = 7, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "corrected_ratio_continuous_weights.pdf"),
  weight_figure, width = 11, height = 7, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  target = "full-path residual log2 RNA/ATAC ratio",
  target_external_regions_used = FALSE, full_pseudotime_scored = TRUE,
  interest_formula = "0.40*|ratio| + 0.40*|first derivative| + 0.20*|second derivative| after componentwise p95 scaling",
  baseline_weight = baseline_weight, maximum_to_minimum_weight =
    max(window_weight) / min(window_weight), interest_balanced_segments = n_segments,
  gene_depth_adjustment = "SCT_data minus conditional log1p(nCount_RNA) coefficient after ns(pt,df=8)",
  target_smoothing = "same 3% trajectory smoother as every gene",
  score = "30% global plus 70% segment-balanced normalized clipped-slope similarity",
  slope_normalization = "each gene and target divided by its own p75 absolute slope then tanh transformed",
  segment_evidence = paste(segment_consistency_penalty,
    "times between-segment correlation SD penalty"),
  tested_genes = nrow(results), direct_top20 = as.list(direct20$feature_id),
  inverse_top20 = as.list(inverse20$feature_id),
  feature_fold_cross_fitting = FALSE,
  inference = "descriptive depth-adjusted screen; cross-fitting remains required for unbiased inference",
  bundle_mutated = FALSE, capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
report <- c("# Corrected-ratio full-pseudotime weighted gene screen", "",
  paste0("All ", nrow(results),
    " detected SCT_data genes were conditionally adjusted for total RNA and ranked."),
  "No BB or T-region definitions were used.", "",
  paste("Every pseudotime window contributes. Continuous weights are derived",
    "only from the corrected ratio's magnitude, first derivative, and second",
    "derivative, with a 0.25 weight floor."),
  paste("The corrected ratio and every gene receive the identical 3% trajectory",
    "smoother. Slopes are divided by each curve's own p75 absolute slope and",
    "tanh transformed, so up/down agreement dominates while magnitude may differ."),
  paste("The trajectory is divided into eight equal target-interest-mass",
    "segments. Seventy percent of the shape score is segment-balanced and",
    "subtracts 0.75 between-segment SD, smoothly penalizing genes whose",
    "agreement is confined to one feature without imposing a cutoff."), "",
  paste("Top direct genes:", paste(direct20$feature_id, collapse = ", ")), "",
  paste("Top inverse genes:", paste(inverse20$feature_id, collapse = ", ")), "",
  paste("This is a descriptive depth-adjusted screen. The ratio still contains",
    "RNA features used to construct the RNA PCA; feature-fold cross-fitting is",
    "required before interpreting ranks as unbiased evidence."))
writeLines(report, file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed residual-ratio weighted gene screen at", output_dir, "\n")
