#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 27_preliminary_modality_turn_pattern_screen.R",
  "<bdev_bundle> <window_report> <smoothed_report> <turning_report>",
  "<output_dir> <expected_bundle_fingerprint> <batch_size>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
smoothed_report <- normalizePath(args[[3L]], mustWork = TRUE)
turning_report <- normalizePath(args[[4L]], mustWork = TRUE)
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

trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
windows <- read_tsv(file.path(smoothed_report, "smoothed_slope_windows.tsv"))
turns <- read_tsv(file.path(turning_report, "selected_turning_regions.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
turns <- turns[match(paste0("T", 1:4), turns$region_id), ]
stopifnot(nrow(trajectory) == 12642L, nrow(windows) == 2524L,
  !anyDuplicated(trajectory$cell_id), !anyDuplicated(windows$window_id),
  all(windows$smoothed_local_slope > 0), all(table(membership$window_id) == 30L),
  identical(turns$region_id, paste0("T", 1:4)))

cell_index <- match(membership$cell_id, trajectory$cell_id)
window_index <- match(membership$window_id, windows$window_id)
stopifnot(!anyNA(cell_index), !anyNA(window_index))
window_design <- sparseMatrix(i = cell_index, j = window_index,
  x = rep(1 / 30, nrow(membership)),
  dims = c(nrow(trajectory), nrow(windows)))
membership_split <- split(cell_index, factor(membership$window_id,
  levels = windows$window_id))

# Keep the fixed smoother used by the preceding screen, but score only the
# direction and amplitude of changes, especially around the four accepted turns.
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

target <- log2(windows$smoothed_local_slope)
target_derivative <- (target[3:length(target)] -
  target[1:(length(target) - 2L)]) / 2
strong_change_quantile <- 0.70
strong_change_threshold <- unname(quantile(abs(target_derivative),
  strong_change_quantile, names = FALSE))
strong_mask <- abs(target_derivative) >= strong_change_threshold &
  abs(target_derivative) > 0
strong_weights <- abs(target_derivative[strong_mask])
strong_weights <- strong_weights / sum(strong_weights)
target_strong <- target_derivative[strong_mask]

event_flank_windows <- 100L
event_index <- as.integer(turns$window_index)
stopifnot(all(event_index - event_flank_windows >= 1L),
  all(event_index + event_flank_windows <= nrow(windows)))
target_pre <- target[event_index] - target[event_index - event_flank_windows]
target_post <- target[event_index + event_flank_windows] - target[event_index]
target_pre_sign <- sign(target_pre)
target_post_sign <- sign(target_post)
stopifnot(all(target_pre_sign != 0), all(target_post_sign != 0))
event_definition <- data.frame(region_id = turns$region_id,
  window_index = event_index, pt = turns$pt,
  flank_windows = event_flank_windows,
  left_window_index = event_index - event_flank_windows,
  right_window_index = event_index + event_flank_windows,
  target_pre_change = target_pre, target_post_change = target_post,
  direct_pre_direction = ifelse(target_pre_sign > 0, "up", "down"),
  direct_post_direction = ifelse(target_post_sign > 0, "up", "down"),
  stringsAsFactors = FALSE)
write_tsv(event_definition, file.path(stage_dir, "turn_event_definition.tsv"))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))
features <- wb_feature_ids(bundle, "SCT_data")
stopifnot(length(features) > 10000L, !anyDuplicated(features))

minimum_detection_fraction <- 0.01
minimum_dynamic_range <- 0.50
minimum_local_change_fraction <- 0.05
minimum_weighted_sign_concordance <- 0.65
minimum_weighted_cosine <- 0.20
minimum_matched_flanks <- 6L
minimum_complete_turns <- 2L
results_list <- list()
events_list <- list()
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
  profiles <- as.matrix((expression[keep, , drop = FALSE] %*%
    window_design) %*% smoothing_design)
  rownames(profiles) <- ids
  profile_quantiles <- apply(profiles, 1L, quantile,
    probs = c(0.05, 0.95), names = FALSE, type = 7)
  dynamic_range <- profile_quantiles[2L, ] - profile_quantiles[1L, ]
  derivatives <- (profiles[, 3:nrow(windows), drop = FALSE] -
    profiles[, 1:(nrow(windows) - 2L), drop = FALSE]) / 2
  gene_strong <- derivatives[, strong_mask, drop = FALSE]
  target_sign_matrix <- matrix(sign(target_strong), nrow = nrow(gene_strong),
    ncol = length(target_strong), byrow = TRUE)
  direct_sign_concordance <- as.numeric((sign(gene_strong) ==
    target_sign_matrix) %*% strong_weights)
  inverse_sign_concordance <- as.numeric((sign(gene_strong) ==
    -target_sign_matrix) %*% strong_weights)
  weighted_numerator <- rowSums(sweep(gene_strong, 2L,
    strong_weights * target_strong, "*"))
  weighted_denominator <- sqrt(rowSums(sweep(gene_strong^2, 2L,
    strong_weights, "*")) * sum(strong_weights * target_strong^2))
  direct_cosine <- weighted_numerator / weighted_denominator
  direct_cosine[!is.finite(direct_cosine)] <- 0
  inverse_cosine <- -direct_cosine

  gene_pre <- vapply(event_index, function(i)
    profiles[, i] - profiles[, i - event_flank_windows], numeric(nrow(profiles)))
  gene_post <- vapply(event_index, function(i)
    profiles[, i + event_flank_windows] - profiles[, i], numeric(nrow(profiles)))
  local_minimum <- minimum_local_change_fraction * dynamic_range
  direct_pre_match <- sweep(sign(gene_pre), 2L, target_pre_sign, "==") &
    abs(gene_pre) >= local_minimum
  direct_post_match <- sweep(sign(gene_post), 2L, target_post_sign, "==") &
    abs(gene_post) >= local_minimum
  inverse_pre_match <- sweep(sign(gene_pre), 2L, -target_pre_sign, "==") &
    abs(gene_pre) >= local_minimum
  inverse_post_match <- sweep(sign(gene_post), 2L, -target_post_sign, "==") &
    abs(gene_post) >= local_minimum
  direct_matched_flanks <- rowSums(direct_pre_match) + rowSums(direct_post_match)
  inverse_matched_flanks <- rowSums(inverse_pre_match) + rowSums(inverse_post_match)
  direct_complete_turns <- rowSums(direct_pre_match & direct_post_match)
  inverse_complete_turns <- rowSums(inverse_pre_match & inverse_post_match)
  direct_score <- 0.45 * direct_cosine +
    0.35 * (2 * direct_sign_concordance - 1) +
    0.20 * (direct_matched_flanks / 4 - 1)
  inverse_score <- 0.45 * inverse_cosine +
    0.35 * (2 * inverse_sign_concordance - 1) +
    0.20 * (inverse_matched_flanks / 4 - 1)
  direct_selected <- dynamic_range >= minimum_dynamic_range &
    direct_sign_concordance >= minimum_weighted_sign_concordance &
    direct_cosine >= minimum_weighted_cosine &
    direct_matched_flanks >= minimum_matched_flanks &
    direct_complete_turns >= minimum_complete_turns
  inverse_selected <- dynamic_range >= minimum_dynamic_range &
    inverse_sign_concordance >= minimum_weighted_sign_concordance &
    inverse_cosine >= minimum_weighted_cosine &
    inverse_matched_flanks >= minimum_matched_flanks &
    inverse_complete_turns >= minimum_complete_turns
  results_list[[length(results_list) + 1L]] <- data.frame(
    feature_id = ids, detection_fraction = detection_fraction[keep],
    expression_dynamic_range_p95_p05 = dynamic_range,
    direct_weighted_sign_concordance = direct_sign_concordance,
    direct_weighted_derivative_cosine = direct_cosine,
    direct_matched_flanks_of_8 = direct_matched_flanks,
    direct_complete_turns_of_4 = direct_complete_turns,
    direct_composite_score = direct_score, direct_selected = direct_selected,
    inverse_weighted_sign_concordance = inverse_sign_concordance,
    inverse_weighted_derivative_cosine = inverse_cosine,
    inverse_matched_flanks_of_8 = inverse_matched_flanks,
    inverse_complete_turns_of_4 = inverse_complete_turns,
    inverse_composite_score = inverse_score, inverse_selected = inverse_selected,
    stringsAsFactors = FALSE)
  events_list[[length(events_list) + 1L]] <- do.call(rbind,
    lapply(seq_len(4L), function(k) data.frame(feature_id = ids,
      region_id = turns$region_id[[k]], gene_pre_change = gene_pre[, k],
      gene_post_change = gene_post[, k], local_minimum_change = local_minimum,
      direct_pre_match = direct_pre_match[, k],
      direct_post_match = direct_post_match[, k],
      direct_complete_turn = direct_pre_match[, k] & direct_post_match[, k],
      inverse_pre_match = inverse_pre_match[, k],
      inverse_post_match = inverse_post_match[, k],
      inverse_complete_turn = inverse_pre_match[, k] & inverse_post_match[, k],
      stringsAsFactors = FALSE)))
  rm(expression, profiles, derivatives, gene_strong)
  if (batch_number %% 10L == 0L) gc(verbose = FALSE)
}

results <- do.call(rbind, results_list)
events <- do.call(rbind, events_list)
stopifnot(nrow(results) > 10000L, !anyDuplicated(results$feature_id),
  nrow(events) == 4L * nrow(results),
  !any(results$direct_selected & results$inverse_selected))
results <- results[order(results$feature_id), ]
write_tsv_gz(results, file.path(stage_dir, "modality_turn_pattern_scores.tsv.gz"))
write_tsv_gz(events, file.path(stage_dir, "gene_turn_event_matches.tsv.gz"))

rank_direction <- function(direction, n = 100L) {
  selected <- results[[paste0(direction, "_selected")]]
  score <- results[[paste0(direction, "_composite_score")]]
  amplitude_eligible <- results$expression_dynamic_range_p95_p05 >=
    minimum_dynamic_range
  ordering <- order(!amplitude_eligible, !selected, -score,
    -results$expression_dynamic_range_p95_p05, results$feature_id)
  head(results[ordering, ], n)
}
top_direct <- rank_direction("direct")
top_inverse <- rank_direction("inverse")
top_direct$rank <- seq_len(nrow(top_direct))
top_inverse$rank <- seq_len(nrow(top_inverse))
write_tsv(top_direct, file.path(stage_dir, "direct_turn_pattern_top100.tsv"))
write_tsv(top_inverse, file.path(stage_dir, "inverse_turn_pattern_top100.tsv"))
write_tsv(results[results$direct_selected, ],
  file.path(stage_dir, "direct_turn_pattern_selected.tsv"))
write_tsv(results[results$inverse_selected, ],
  file.path(stage_dir, "inverse_turn_pattern_selected.tsv"))

sensitivity <- do.call(rbind, lapply(c(0.25, 0.50, 1.00), function(cutoff) {
  direct <- results$expression_dynamic_range_p95_p05 >= cutoff &
    results$direct_weighted_sign_concordance >= minimum_weighted_sign_concordance &
    results$direct_weighted_derivative_cosine >= minimum_weighted_cosine &
    results$direct_matched_flanks_of_8 >= minimum_matched_flanks &
    results$direct_complete_turns_of_4 >= minimum_complete_turns
  inverse <- results$expression_dynamic_range_p95_p05 >= cutoff &
    results$inverse_weighted_sign_concordance >= minimum_weighted_sign_concordance &
    results$inverse_weighted_derivative_cosine >= minimum_weighted_cosine &
    results$inverse_matched_flanks_of_8 >= minimum_matched_flanks &
    results$inverse_complete_turns_of_4 >= minimum_complete_turns
  data.frame(dynamic_range_cutoff = cutoff,
    direct_selected = sum(direct), inverse_selected = sum(inverse))
}))
write_tsv(sensitivity, file.path(stage_dir, "dynamic_range_sensitivity.tsv"))

# Re-read only the displayed genes to add within-window cell-level IQR.
top_n <- 20L
direct20 <- head(top_direct, top_n)
inverse20 <- head(top_inverse, top_n)
plot_features <- unique(c(direct20$feature_id, inverse20$feature_id))
plot_expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data",
  plot_features, trajectory$cell_id, missing = "NA"))
raw_mean <- as.matrix(plot_expression %*% window_design)
raw_q25 <- vapply(membership_split, function(ix)
  apply(plot_expression[, ix, drop = FALSE], 1L, quantile, probs = 0.25,
    names = FALSE, type = 7), numeric(length(plot_features)))
raw_q75 <- vapply(membership_split, function(ix)
  apply(plot_expression[, ix, drop = FALSE], 1L, quantile, probs = 0.75,
    names = FALSE, type = 7), numeric(length(plot_features)))
rownames(raw_q25) <- rownames(raw_q75) <- plot_features
smoothed_mean <- as.matrix(raw_mean %*% smoothing_design)
smoothed_q25 <- as.matrix(raw_q25 %*% smoothing_design)
smoothed_q75 <- as.matrix(raw_q75 %*% smoothing_design)
rownames(smoothed_mean) <- rownames(smoothed_q25) <-
  rownames(smoothed_q75) <- plot_features

make_plot_data <- function(top, direction) {
  rows <- lapply(seq_len(nrow(top)), function(i) {
    gene <- top$feature_id[[i]]
    data.frame(pt = windows$pt_median,
      gene_panel = sprintf("%02d  %s  (score = %.2f; turns = %d/4)", i, gene,
        top[[paste0(direction, "_composite_score")]][[i]],
        top[[paste0(direction, "_complete_turns_of_4")]][[i]]),
      mean = smoothed_mean[gene, ], q25 = smoothed_q25[gene, ],
      q75 = smoothed_q75[gene, ], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$gene_panel <- factor(out$gene_panel, levels = unique(out$gene_panel))
  out
}
make_expression_figure <- function(data, direction) {
  direction_label <- if (direction == "direct") "RNA/ATAC direct" else
    "ATAC/RNA inverse"
  color <- if (direction == "direct") "#C51B5A" else "#2166AC"
  ggplot(data, aes(pt)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color, alpha = 0.20,
      color = NA) +
    geom_line(aes(y = mean), color = color, linewidth = 0.8) +
    geom_vline(data = turns, aes(xintercept = pt), inherit.aes = FALSE,
      color = "#333333", linetype = "dashed", linewidth = 0.25) +
    facet_wrap(~ gene_panel, ncol = 4L, scales = "free_y") +
    theme_classic(base_size = 9) +
    theme(strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 8.2),
      panel.spacing = grid::unit(5, "mm"),
      plot.title = element_text(face = "bold", size = 15)) +
    labs(title = paste("Top 20 event-pattern genes:", direction_label),
      subtitle = paste(
        "Genes are ranked after amplitude gating; ribbon is within-window IQR;",
        "dashed lines are fixed T1-T4"),
      x = "Bdev pseudotime", y = "Smoothed SCT_data expression")
}
direct_plot <- make_expression_figure(make_plot_data(direct20, "direct"),
  "direct")
inverse_plot <- make_expression_figure(make_plot_data(inverse20, "inverse"),
  "inverse")
ggsave(file.path(stage_dir, "direct_turn_pattern_top20_with_iqr.png"),
  direct_plot, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "direct_turn_pattern_top20_with_iqr.pdf"),
  direct_plot, width = 14, height = 12, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "inverse_turn_pattern_top20_with_iqr.png"),
  inverse_plot, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "inverse_turn_pattern_top20_with_iqr.pdf"),
  inverse_plot, width = 14, height = 12, device = cairo_pdf, bg = "white")

counts <- list(features_in_sct_data = length(features), tested_features = nrow(results),
  direct_selected = sum(results$direct_selected),
  inverse_selected = sum(results$inverse_selected))
diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  expression_matrix = "SCT_data", trajectory_cells = nrow(trajectory),
  windows = nrow(windows), window_size = 30L, stride = 5L,
  smoothing_fraction = smoothing_fraction,
  smoothing_half_width_windows = smoothing_half_width,
  target = "log2(smoothed delta-RNA/delta-ATAC)",
  strong_change_mask = paste0("top ", 100 * (1 - strong_change_quantile),
    "% of absolute target derivatives"),
  event_flank_windows = event_flank_windows,
  minimum_detection_fraction = minimum_detection_fraction,
  minimum_dynamic_range_p95_p05 = minimum_dynamic_range,
  minimum_local_change_fraction_of_dynamic_range = minimum_local_change_fraction,
  minimum_weighted_sign_concordance = minimum_weighted_sign_concordance,
  minimum_weighted_derivative_cosine = minimum_weighted_cosine,
  minimum_matched_flanks_of_8 = minimum_matched_flanks,
  minimum_complete_turns_of_4 = minimum_complete_turns,
  composite_score = paste(
    "0.45 weighted derivative cosine + 0.35 oriented sign agreement +",
    "0.20 T1-T4 flank agreement"),
  counts = counts, direct_top20 = as.list(direct20$feature_id),
  inverse_top20 = as.list(inverse20$feature_id),
  bundle_mutated = FALSE, capsule_registered = FALSE,
  inference = "descriptive ranking and threshold screen; no p-values")
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
report <- c("# Event-focused modality-pattern gene screen", "",
  paste0("Tested ", nrow(results), " genes; selected ", counts$direct_selected,
    " direct RNA/ATAC-pattern genes and ", counts$inverse_selected,
    " inverse ATAC/RNA-pattern genes."), "",
  "## Why this replaces the state-wide ranking", "",
  paste("The previous whole-trajectory state correlation could reward a nearly",
    "flat gene for agreeing across long quiet portions. This screen requires a",
    "p95-p05 expression range of at least 0.50 SCT_data, emphasizes the 30%",
    "strongest target changes, and explicitly evaluates both sides of T1-T4."), "",
  "A selected gene must match at least six of eight event flanks, reproduce at",
  "least two of four complete turns, have weighted direction concordance >= 0.65,",
  "and have weighted derivative cosine >= 0.20. A local event change smaller than",
  "5% of that gene's p95-p05 range cannot count as a match.", "",
  "Top direct genes:", paste(direct20$feature_id, collapse = ", "), "",
  "Top inverse genes:", paste(inverse20$feature_id, collapse = ", "), "",
  "## Interpretation", "",
  paste("Direct means the gene rises and falls in the same orientation as",
    "log2(delta RNA/delta ATAC); inverse means the same criteria were applied",
    "against its exact sign reversal, log2(delta ATAC/delta RNA). The complete",
    "per-gene/per-turn decisions and amplitude sensitivity are exported."),
  "This is a deterministic descriptive screen, not a significance test.")
writeLines(report, file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed event-focused modality-pattern screen at", output_dir, "\n")
