#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 30_preliminary_t1_t3_continuous_gene_screen.R",
  "<bdev_bundle> <window_report> <smoothed_report> <turn_pattern_report>",
  "<output_dir> <expected_bundle_fingerprint> <batch_size>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
smoothed_report <- normalizePath(args[[3L]], mustWork = TRUE)
turn_pattern_report <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
expected_fingerprint <- args[[6L]]
batch_size <- as.integer(args[[7L]])
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
  library(jsonlite)
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
row_cosine <- function(x, target, weights = rep(1 / length(target),
    length(target))) {
  numerator <- rowSums(sweep(x, 2L, weights * target, "*"))
  denominator <- sqrt(rowSums(sweep(x^2, 2L, weights, "*")) *
    sum(weights * target^2))
  answer <- numerator / denominator
  answer[!is.finite(answer)] <- 0
  answer
}

trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
windows <- read_tsv(file.path(smoothed_report, "smoothed_slope_windows.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
event_definition <- read_tsv(file.path(turn_pattern_report,
  "turn_event_definition.tsv"))
event_definition <- event_definition[match(paste0("T", 1:3),
  event_definition$region_id), ]
stopifnot(nrow(trajectory) == 12642L, nrow(windows) == 2524L,
  !anyDuplicated(trajectory$cell_id), !anyDuplicated(windows$window_id),
  all(windows$smoothed_local_slope > 0), all(table(membership$window_id) == 30L),
  identical(event_definition$region_id, paste0("T", 1:3)))

cell_index <- match(membership$cell_id, trajectory$cell_id)
window_index <- match(membership$window_id, windows$window_id)
stopifnot(!anyNA(cell_index), !anyNA(window_index))
window_design <- sparseMatrix(i = cell_index, j = window_index,
  x = rep(1 / 30, nrow(membership)),
  dims = c(nrow(trajectory), nrow(windows)))
membership_split <- split(cell_index, factor(membership$window_id,
  levels = windows$window_id))

smoothing_fraction <- 0.08
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
stopifnot(all(abs(colSums(window_design) - 1) < 1e-12),
  all(abs(colSums(smoothing_design) - 1) < 1e-12))

# T4 is deliberately excluded. The analysis interval starts at the left T1
# flank and ends at the right T3 flank, retaining all intervening trajectory.
interval_indices <- seq.int(event_definition$left_window_index[[1L]],
  event_definition$right_window_index[[3L]])
needed_iqr_indices <- seq.int(min(interval_indices) - smoothing_half_width,
  max(interval_indices) + smoothing_half_width)
stopifnot(min(needed_iqr_indices) >= 1L,
  max(needed_iqr_indices) <= nrow(windows))
interval_smoothing_design <- smoothing_design[needed_iqr_indices,
  interval_indices, drop = FALSE]

target <- log2(windows$smoothed_local_slope)
target_derivative <- (target[3:length(target)] -
  target[1:(length(target) - 2L)]) / 2
derivative_center_indices <- interval_indices[
  interval_indices > 1L & interval_indices < nrow(windows)]
derivative_column_indices <- derivative_center_indices - 1L
target_derivative_interval <- target_derivative[derivative_column_indices]
derivative_weights <- abs(target_derivative_interval)
derivative_weights <- derivative_weights / sum(derivative_weights)
event_indices <- event_definition$window_index
left_indices <- event_definition$left_window_index
right_indices <- event_definition$right_window_index
target_event_vector <- as.vector(rbind(target[event_indices] -
  target[left_indices], target[right_indices] - target[event_indices]))
stopifnot(length(target_event_vector) == 6L)

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))
features <- wb_feature_ids(bundle, "SCT_data")
stopifnot(length(features) > 10000L, !anyDuplicated(features))

minimum_detection_fraction <- 0.01
derivative_weight <- 0.60
event_vector_weight <- 0.40
amplitude_half_saturation <- 0.50
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
    identical(colnames(expression), trajectory$cell_id), all(is.finite(expression)))
  detection_fraction <- rowMeans(expression > 0)
  keep <- detection_fraction >= minimum_detection_fraction
  if (!any(keep)) next
  ids <- batch_features[keep]
  kept_expression <- expression[keep, , drop = FALSE]
  profiles <- as.matrix((kept_expression %*% window_design) %*%
    smoothing_design)
  rownames(profiles) <- ids
  interval_profiles <- profiles[, interval_indices, drop = FALSE]
  profile_quantiles <- rowQuantiles(interval_profiles,
    probs = c(0.05, 0.95), drop = FALSE)
  dynamic_range <- profile_quantiles[, 2L] - profile_quantiles[, 1L]
  derivatives <- (profiles[, 3:nrow(windows), drop = FALSE] -
    profiles[, 1:(nrow(windows) - 2L), drop = FALSE]) / 2
  derivative_cosine <- row_cosine(
    derivatives[, derivative_column_indices, drop = FALSE],
    target_derivative_interval, derivative_weights)
  gene_event_vector <- do.call(cbind, lapply(seq_len(3L), function(k)
    cbind(profiles[, event_indices[[k]]] - profiles[, left_indices[[k]]],
      profiles[, right_indices[[k]]] - profiles[, event_indices[[k]]])))
  event_vector_cosine <- row_cosine(gene_event_vector, target_event_vector)
  raw_shape_score <- derivative_weight * derivative_cosine +
    event_vector_weight * event_vector_cosine

  # Exact cell-level IQR is calculated for each source 30-cell window needed
  # to smooth the T1-T3 interval. SD is retained as a safeguard against a
  # zero IQR in sparse, zero-inflated genes.
  raw_q25 <- matrix(NA_real_, nrow = nrow(kept_expression),
    ncol = length(needed_iqr_indices))
  raw_q75 <- raw_q25
  raw_sd <- raw_q25
  for (j in seq_along(needed_iqr_indices)) {
    ix <- membership_split[[needed_iqr_indices[[j]]]]
    q <- rowQuantiles(kept_expression[, ix, drop = FALSE],
      probs = c(0.25, 0.75), drop = FALSE)
    raw_q25[, j] <- q[, 1L]
    raw_q75[, j] <- q[, 2L]
    raw_sd[, j] <- rowSds(kept_expression[, ix, drop = FALSE])
  }
  smoothed_q25 <- as.matrix(raw_q25 %*% interval_smoothing_design)
  smoothed_q75 <- as.matrix(raw_q75 %*% interval_smoothing_design)
  median_iqr <- rowMedians(smoothed_q75 - smoothed_q25)
  interval_raw_columns <- match(interval_indices, needed_iqr_indices)
  median_window_sd <- rowMedians(raw_sd[, interval_raw_columns, drop = FALSE])
  conservative_cell_spread <- pmax(median_iqr, 1.349 * median_window_sd)
  precision_reliability <- dynamic_range /
    (dynamic_range + conservative_cell_spread)
  precision_reliability[!is.finite(precision_reliability)] <- 0
  amplitude_weight <- dynamic_range /
    (dynamic_range + amplitude_half_saturation)
  amplitude_weight[!is.finite(amplitude_weight)] <- 0
  evidence_weight <- precision_reliability * amplitude_weight
  direct_score <- raw_shape_score * evidence_weight
  inverse_score <- -raw_shape_score * evidence_weight
  results_list[[length(results_list) + 1L]] <- data.frame(
    feature_id = ids, detection_fraction = detection_fraction[keep],
    t1_t3_dynamic_range_p95_p05 = dynamic_range,
    t1_t3_median_smoothed_window_iqr = median_iqr,
    t1_t3_median_window_sd = median_window_sd,
    t1_t3_precision_reliability = precision_reliability,
    t1_t3_amplitude_weight = amplitude_weight,
    t1_t3_evidence_weight = evidence_weight,
    t1_t3_derivative_cosine_direct = derivative_cosine,
    t1_t3_six_flank_vector_cosine_direct = event_vector_cosine,
    t1_t3_raw_shape_score_direct = raw_shape_score,
    t1_t3_evidence_weighted_score_direct = direct_score,
    t1_t3_raw_shape_score_inverse = -raw_shape_score,
    t1_t3_evidence_weighted_score_inverse = inverse_score,
    stringsAsFactors = FALSE)
  rm(expression, kept_expression, profiles, interval_profiles, derivatives,
    raw_q25, raw_q75, raw_sd, smoothed_q25, smoothed_q75)
  if (batch_number %% 5L == 0L) gc(verbose = FALSE)
}

results <- do.call(rbind, results_list)
stopifnot(nrow(results) > 10000L, !anyDuplicated(results$feature_id),
  all(is.finite(results$t1_t3_raw_shape_score_direct)),
  all(is.finite(results$t1_t3_evidence_weighted_score_direct)))
results$direct_shape_rank <- rank(-results$t1_t3_raw_shape_score_direct,
  ties.method = "min")
results$inverse_shape_rank <- rank(-results$t1_t3_raw_shape_score_inverse,
  ties.method = "min")
results$direct_evidence_rank <- rank(
  -results$t1_t3_evidence_weighted_score_direct, ties.method = "min")
results$inverse_evidence_rank <- rank(
  -results$t1_t3_evidence_weighted_score_inverse, ties.method = "min")
results <- results[order(results$feature_id), ]
write_tsv_gz(results, file.path(stage_dir,
  "t1_t3_continuous_gene_scores.tsv.gz"))

top_direct <- head(results[order(results$direct_evidence_rank,
  results$feature_id), ], 100L)
top_inverse <- head(results[order(results$inverse_evidence_rank,
  results$feature_id), ], 100L)
write_tsv(top_direct, file.path(stage_dir,
  "t1_t3_continuous_top100_direct.tsv"))
write_tsv(top_inverse, file.path(stage_dir,
  "t1_t3_continuous_top100_inverse.tsv"))

genes_of_prior_interest <- c("Nipbl", "Mau2", "Wapl", "Ctcf", "Rad21",
  "Smc1a", "Smc3", "Stag1", "Stag2", "Pds5a", "Pds5b", "Baz2b",
  "Smarca4", "Arid1b", "Bcl7a", "Ezh2", "Eed", "Suz12")
interest_table <- results[match(genes_of_prior_interest, results$feature_id), ]
stopifnot(!anyNA(interest_table$feature_id))
write_tsv(interest_table, file.path(stage_dir,
  "prior_interest_genes_continuous_scores.tsv"))

# Plot the unbiased top 20 from each direction with exact IQR. These genes are
# chosen solely by the genome-wide evidence rank above.
top_n <- 20L
direct20 <- head(top_direct, top_n)
inverse20 <- head(top_inverse, top_n)
plot_features <- unique(c(direct20$feature_id, inverse20$feature_id))
plot_expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data",
  plot_features, trajectory$cell_id, missing = "NA"))
plot_raw_mean <- as.matrix(plot_expression %*% window_design)
plot_profiles <- as.matrix(plot_raw_mean %*% smoothing_design)
plot_q25 <- matrix(NA_real_, nrow = length(plot_features),
  ncol = length(needed_iqr_indices))
plot_q75 <- plot_q25
for (j in seq_along(needed_iqr_indices)) {
  ix <- membership_split[[needed_iqr_indices[[j]]]]
  q <- rowQuantiles(plot_expression[, ix, drop = FALSE],
    probs = c(0.25, 0.75), drop = FALSE)
  plot_q25[, j] <- q[, 1L]
  plot_q75[, j] <- q[, 2L]
}
plot_mean <- plot_profiles[, interval_indices, drop = FALSE]
plot_smoothed_q25 <- as.matrix(plot_q25 %*% interval_smoothing_design)
plot_smoothed_q75 <- as.matrix(plot_q75 %*% interval_smoothing_design)
rownames(plot_mean) <- rownames(plot_smoothed_q25) <-
  rownames(plot_smoothed_q75) <- plot_features
plot_center <- rowMeans(plot_mean)
plot_sd <- sqrt(rowSums((plot_mean - plot_center)^2) /
  (ncol(plot_mean) - 1L))
plot_mean_z <- (plot_mean - plot_center) / plot_sd
plot_q25_z <- (plot_smoothed_q25 - plot_center) / plot_sd
plot_q75_z <- (plot_smoothed_q75 - plot_center) / plot_sd
direct_target_z <- as.numeric(scale(target[interval_indices]))
inverse_target_z <- -direct_target_z

make_plot_data <- function(top, direction) {
  rows <- lapply(seq_len(nrow(top)), function(i) {
    gene <- top$feature_id[[i]]
    data.frame(pt = windows$pt_median[interval_indices],
      gene_panel = sprintf("%02d  %s | shape %.2f | evidence %.2f", i, gene,
        top[[paste0("t1_t3_raw_shape_score_", direction)]][[i]],
        top[[paste0("t1_t3_evidence_weighted_score_", direction)]][[i]]),
      mean_z = plot_mean_z[gene, ], q25_z = plot_q25_z[gene, ],
      q75_z = plot_q75_z[gene, ],
      target_z = if (direction == "direct") direct_target_z else inverse_target_z,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$gene_panel <- factor(out$gene_panel, levels = unique(out$gene_panel))
  out
}
make_figure <- function(data, direction) {
  color <- if (direction == "direct") "#C51B5A" else "#2166AC"
  label <- if (direction == "direct") "direct RNA/ATAC" else
    "inverse ATAC/RNA"
  ggplot(data, aes(pt)) +
    geom_hline(yintercept = 0, color = "#999999", linetype = "dotted",
      linewidth = 0.25) +
    geom_ribbon(aes(ymin = q25_z, ymax = q75_z), fill = color,
      alpha = 0.18, color = NA) +
    geom_line(aes(y = mean_z), color = color, linewidth = 0.8) +
    geom_line(aes(y = target_z), color = "#111111", linewidth = 0.55) +
    geom_vline(data = event_definition, aes(xintercept = pt),
      inherit.aes = FALSE, color = "#555555", linetype = "dashed",
      linewidth = 0.25) +
    facet_wrap(~ gene_panel, ncol = 4L) +
    coord_cartesian(ylim = c(-3.5, 3.5)) +
    theme_classic(base_size = 9) +
    theme(strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 8.0),
      panel.spacing = grid::unit(5, "mm"),
      plot.title = element_text(face = "bold", size = 15)) +
    labs(title = paste("Unbiased top 20 continuous T1-T3", label, "hits"),
      subtitle = paste(
        "Colored line/ribbon: standardized gene mean and IQR; black: target;",
        "T4 excluded; ranking includes continuous amplitude/precision weight"),
      x = "Bdev pseudotime", y = "Standardized profile")
}
direct_plot <- make_figure(make_plot_data(direct20, "direct"), "direct")
inverse_plot <- make_figure(make_plot_data(inverse20, "inverse"), "inverse")
ggsave(file.path(stage_dir, "t1_t3_continuous_top20_direct_with_iqr.png"),
  direct_plot, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "t1_t3_continuous_top20_direct_with_iqr.pdf"),
  direct_plot, width = 14, height = 12, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "t1_t3_continuous_top20_inverse_with_iqr.png"),
  inverse_plot, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "t1_t3_continuous_top20_inverse_with_iqr.pdf"),
  inverse_plot, width = 14, height = 12, device = cairo_pdf, bg = "white")

nipbl <- results[results$feature_id == "Nipbl", ]
scatter_plot <- ggplot(results, aes(t1_t3_raw_shape_score_direct,
    t1_t3_evidence_weighted_score_direct)) +
  geom_point(color = "#BDBDBD", alpha = 0.45, size = 0.65) +
  geom_point(data = nipbl, color = "#7A0177", size = 2.8) +
  geom_text(data = nipbl, aes(label = paste0("Nipbl: rank ",
    direct_evidence_rank)), color = "#7A0177", fontface = "bold",
    nudge_y = 0.025, size = 3.5) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Continuous T1-T3 shape versus amplitude/IQR-weighted evidence",
    subtitle = "All detected genes were scored identically; T4 was excluded",
    x = "Raw direct shape score",
    y = "Evidence-weighted direct score")
ggsave(file.path(stage_dir, "t1_t3_shape_vs_evidence_nipbl.png"),
  scatter_plot, width = 7.5, height = 6.5, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "t1_t3_shape_vs_evidence_nipbl.pdf"),
  scatter_plot, width = 7.5, height = 6.5, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  expression_matrix = "SCT_data", trajectory_cells = nrow(trajectory),
  windows = nrow(windows), window_size = 30L, stride = 5L,
  excluded_turn = "T4", included_turns = as.list(paste0("T", 1:3)),
  interval_window_indices = c(min(interval_indices), max(interval_indices)),
  smoothing_fraction = smoothing_fraction,
  raw_shape_score = paste0(derivative_weight,
    " * target-change-weighted derivative cosine + ", event_vector_weight,
    " * six-flank-vector cosine"),
  precision_reliability = paste(
    "dynamic range / (dynamic range + conservative cell spread);",
    "spread is max(median smoothed IQR, 1.349 * median window SD)"),
  amplitude_weight = paste0("dynamic range / (dynamic range + ",
    amplitude_half_saturation, ")"),
  evidence_score = "raw shape score * precision reliability * amplitude weight",
  hard_turn_or_flank_cutoffs = FALSE,
  tested_genes = nrow(results), nipbl = as.list(nipbl[1, ]),
  direct_top20 = as.list(direct20$feature_id),
  inverse_top20 = as.list(inverse20$feature_id),
  bundle_mutated = FALSE, capsule_registered = FALSE,
  inference = "unbiased genome-wide descriptive ranking; no significance cutoff")
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
report <- c("# Continuous T1-T3 gene-pattern screen", "",
  paste0("All ", nrow(results),
    " detected genes were ranked without gene-class preselection. T4 was excluded."),
  "", "## Method", "",
  paste("The raw shape score is 60% continuous derivative cosine across the",
    "entire T1-T3 interval and 40% cosine of the six T1-T3 flank changes.",
    "No flank or complete-turn pass count is used."),
  paste("The ranking score smoothly downweights low-amplitude or imprecise",
    "profiles. Precision uses the gene's p95-p05 dynamic range relative to a",
    "conservative within-window spread derived from both exact IQR and SD.",
    "There is no hard expression-amplitude threshold."), "",
  "Top direct genes:", paste(direct20$feature_id, collapse = ", "), "",
  "Top inverse genes:", paste(inverse20$feature_id, collapse = ", "), "",
  paste0("Nipbl direct evidence rank: ", nipbl$direct_evidence_rank,
    " of ", nrow(results), "; raw shape score: ",
    sprintf("%.3f", nipbl$t1_t3_raw_shape_score_direct),
    "; evidence-weighted score: ",
    sprintf("%.3f", nipbl$t1_t3_evidence_weighted_score_direct), "."), "",
  "This is an unbiased descriptive ranking, not a differential-expression or",
  "significance test. Exact scores for every tested gene are exported.")
writeLines(report, file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed continuous T1-T3 gene screen at", output_dir, "\n")
