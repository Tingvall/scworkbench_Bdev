#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 25_preliminary_modality_gene_pattern_screen.R",
  "<bdev_bundle> <window_report> <smoothed_report> <output_dir>",
  "<expected_bundle_fingerprint> <batch_size> <n_shifts>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
smoothed_report <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_fingerprint <- args[[5L]]
batch_size <- as.integer(args[[6L]])
n_shifts_requested <- as.integer(args[[7L]])
stopifnot(batch_size >= 16L, n_shifts_requested >= 50L)
set.seed(20261009L)
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
row_standardize <- function(x) {
  mu <- rowMeans(x)
  centered <- x - mu
  sd <- sqrt(rowSums(centered^2) / (ncol(x) - 1L))
  list(z = centered / sd, sd = sd)
}
make_shift_matrix <- function(target, requested) {
  n <- length(target)
  offsets <- unique(as.integer(round(seq(ceiling(0.05 * n),
    floor(0.95 * n), length.out = requested))))
  shifted <- vapply(offsets, function(k) c(target[(k + 1L):n],
    target[seq_len(k)]), numeric(n))
  list(matrix = shifted, offsets = offsets)
}
pooled_tail_p <- function(score, null_abs) {
  score <- as.numeric(score)
  null_abs <- as.numeric(null_abs)
  sorted <- sort(null_abs)
  (length(sorted) - findInterval(abs(score), sorted) + 1) /
    (length(sorted) + 1)
}

trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
windows <- read_tsv(file.path(smoothed_report, "smoothed_slope_windows.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
stopifnot(nrow(trajectory) == 12642L, !anyDuplicated(trajectory$cell_id),
  nrow(windows) == 2524L, !anyDuplicated(windows$window_id),
  all(windows$smoothed_local_slope > 0),
  all(table(membership$window_id) == 30L),
  all(membership$cell_id %in% trajectory$cell_id))

cell_index <- match(membership$cell_id, trajectory$cell_id)
window_index <- match(membership$window_id, windows$window_id)
stopifnot(!anyNA(cell_index), !anyNA(window_index))
window_design <- sparseMatrix(i = cell_index, j = window_index,
  x = rep(1 / 30, nrow(membership)),
  dims = c(nrow(trajectory), nrow(windows)))
stopifnot(all(abs(colSums(window_design) - 1) < 1e-12))

# Apply one fixed, linear smoother to every gene before either screen. Its
# total neighborhood is approximately 8% of the ordered windows, matching the
# accepted display span without gene-specific fitting or outcome-dependent
# smoothing. Edge columns use all available neighboring windows and reweight
# to sum to one.
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
stopifnot(all(abs(colSums(smoothing_design) - 1) < 1e-12))

target <- log2(windows$smoothed_local_slope)
target_z <- as.numeric(scale(target))
target_change_z <- as.numeric(scale(target[2:(length(target) - 1L)]))
state_shifts <- make_shift_matrix(target_z, n_shifts_requested)
change_shifts <- make_shift_matrix(target_change_z, n_shifts_requested)
n_shifts <- min(ncol(state_shifts$matrix), ncol(change_shifts$matrix))
state_shift_matrix <- state_shifts$matrix[, seq_len(n_shifts), drop = FALSE]
change_shift_matrix <- change_shifts$matrix[, seq_len(n_shifts), drop = FALSE]

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))
features <- wb_feature_ids(bundle, "SCT_data")
stopifnot(length(features) > 10000L, !anyDuplicated(features))

minimum_detection_fraction <- 0.01
minimum_profile_sd <- 1e-8
minimum_abs_correlation <- 0.30
results_list <- list()
state_null_list <- list()
change_null_list <- list()
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
  cell_mean <- rowMeans(expression)
  cell_sd <- sqrt(pmax(rowMeans(expression^2) - cell_mean^2, 0))
  keep_initial <- detection_fraction >= minimum_detection_fraction &
    cell_sd > minimum_profile_sd
  if (!any(keep_initial)) next
  kept_expression <- expression[keep_initial, , drop = FALSE]
  raw_profiles <- as.matrix(kept_expression %*% window_design)
  profiles <- as.matrix(raw_profiles %*% smoothing_design)
  state_standard <- row_standardize(profiles)
  # Use ordered window index as the derivative coordinate, matching the
  # turning-point workflow. Source pt contains legitimate ties, so dividing by
  # pt differences would create zero denominators.
  derivatives <- (profiles[, 3:nrow(windows), drop = FALSE] -
    profiles[, 1:(nrow(windows) - 2L), drop = FALSE]) / 2
  change_standard <- row_standardize(derivatives)
  keep_profile <- is.finite(state_standard$sd) &
    state_standard$sd > minimum_profile_sd &
    is.finite(change_standard$sd) & change_standard$sd > minimum_profile_sd
  if (!any(keep_profile)) next
  profiles_z <- state_standard$z[keep_profile, , drop = FALSE]
  derivatives_z <- change_standard$z[keep_profile, , drop = FALSE]
  state_score <- as.numeric(profiles_z %*% target_z /
    (length(target_z) - 1L))
  change_score <- as.numeric(derivatives_z %*% target_change_z /
    (length(target_change_z) - 1L))
  state_null <- profiles_z %*% state_shift_matrix /
    (length(target_z) - 1L)
  change_null <- derivatives_z %*% change_shift_matrix /
    (length(target_change_z) - 1L)
  initial_indices <- which(keep_initial)[keep_profile]
  result_features <- as.character(batch_features[initial_indices])
  stopifnot(length(result_features) == nrow(profiles_z),
    !anyDuplicated(result_features))
  results_list[[length(results_list) + 1L]] <- data.frame(
    feature_id = result_features,
    detection_fraction = detection_fraction[initial_indices],
    cell_mean_sct_data = cell_mean[initial_indices],
    cell_sd_sct_data = cell_sd[initial_indices],
    window_profile_sd = state_standard$sd[keep_profile],
    derivative_profile_sd = change_standard$sd[keep_profile],
    state_correlation = state_score,
    state_own_shift_p = (1 + rowSums(abs(state_null) >= abs(state_score))) /
      (n_shifts + 1),
    change_correlation = change_score,
    change_own_shift_p = (1 + rowSums(abs(change_null) >= abs(change_score))) /
      (n_shifts + 1), stringsAsFactors = FALSE)
  state_null_list[[length(state_null_list) + 1L]] <- as.numeric(abs(state_null))
  change_null_list[[length(change_null_list) + 1L]] <-
    as.numeric(abs(change_null))
  rm(expression, kept_expression, raw_profiles, profiles, derivatives, profiles_z,
    derivatives_z, state_null, change_null)
  if (batch_number %% 10L == 0L) gc(verbose = FALSE)
}

results <- as.data.frame(do.call(rbind, results_list), stringsAsFactors = FALSE,
  check.names = FALSE)
numeric_result_columns <- setdiff(names(results), "feature_id")
results[numeric_result_columns] <- lapply(results[numeric_result_columns],
  function(x) as.numeric(unlist(x, use.names = FALSE)))
feature_id_vector <- as.character(unlist(results[["feature_id"]],
  use.names = FALSE))
state_score_vector <- as.double(unlist(results[["state_correlation"]],
  use.names = FALSE))
change_score_vector <- as.double(unlist(results[["change_correlation"]],
  use.names = FALSE))
stopifnot(length(feature_id_vector) == nrow(results),
  length(state_score_vector) == nrow(results),
  length(change_score_vector) == nrow(results))
results[["feature_id"]] <- feature_id_vector
results[["state_correlation"]] <- state_score_vector
results[["change_correlation"]] <- change_score_vector
state_null_abs <- unlist(state_null_list, use.names = FALSE)
change_null_abs <- unlist(change_null_list, use.names = FALSE)
stopifnot(nrow(results) > 1000L, length(state_null_abs) == nrow(results) * n_shifts,
  length(change_null_abs) == nrow(results) * n_shifts,
  !anyDuplicated(results$feature_id), all(is.finite(results$state_correlation)),
  all(is.finite(results$change_correlation)))

results$state_pooled_null_p <- pooled_tail_p(state_score_vector,
  state_null_abs)
results$state_p_adj_bh <- p.adjust(results$state_pooled_null_p, method = "BH")
results$state_direction <- ifelse(state_score_vector >= 0,
  "RNA_faster_pattern", "ATAC_faster_inverse_pattern")
results$state_selected <- results$state_p_adj_bh < 0.05 &
  abs(state_score_vector) >= minimum_abs_correlation
results$change_pooled_null_p <- pooled_tail_p(change_score_vector,
  change_null_abs)
results$change_p_adj_bh <- p.adjust(results$change_pooled_null_p, method = "BH")
results$change_direction <- ifelse(change_score_vector >= 0,
  "RNA_faster_pattern", "ATAC_faster_inverse_pattern")
results$change_selected <- results$change_p_adj_bh < 0.05 &
  abs(change_score_vector) >= minimum_abs_correlation
results <- results[order(feature_id_vector), ]
write_tsv_gz(results, file.path(stage_dir,
  "modality_gene_pattern_scores.tsv.gz"))

top_table <- function(score_column, direction, n = 100L) {
  score <- results[[score_column]]
  ordering <- if (direction == "direct") order(-score, results$feature_id) else
    order(score, results$feature_id)
  results[head(ordering, n), ]
}
top_state_direct <- top_table("state_correlation", "direct")
top_state_inverse <- top_table("state_correlation", "inverse")
top_change_direct <- top_table("change_correlation", "direct")
top_change_inverse <- top_table("change_correlation", "inverse")
write_tsv(top_state_direct,
  file.path(stage_dir, "state_top100_rna_faster_pattern.tsv"))
write_tsv(top_state_inverse,
  file.path(stage_dir, "state_top100_atac_faster_inverse_pattern.tsv"))
write_tsv(top_change_direct,
  file.path(stage_dir, "change_top100_rna_faster_pattern.tsv"))
write_tsv(top_change_inverse,
  file.path(stage_dir, "change_top100_atac_faster_inverse_pattern.tsv"))
write_tsv(results[results$state_selected | results$change_selected, ],
  file.path(stage_dir, "selected_gene_pattern_hits.tsv"))

top_n_consensus <- 20L
membership_rows <- rbind(
  data.frame(approach = "State", direction = "RNA-faster",
    rank = seq_len(top_n_consensus),
    feature_id = top_state_direct$feature_id[seq_len(top_n_consensus)]),
  data.frame(approach = "State", direction = "ATAC-faster inverse",
    rank = seq_len(top_n_consensus),
    feature_id = top_state_inverse$feature_id[seq_len(top_n_consensus)]),
  data.frame(approach = "Change", direction = "RNA-faster",
    rank = seq_len(top_n_consensus),
    feature_id = top_change_direct$feature_id[seq_len(top_n_consensus)]),
  data.frame(approach = "Change", direction = "ATAC-faster inverse",
    rank = seq_len(top_n_consensus),
    feature_id = top_change_inverse$feature_id[seq_len(top_n_consensus)]))
write_tsv(membership_rows,
  file.path(stage_dir, "top20_consensus_membership.tsv"))
consensus_features <- unique(membership_rows$feature_id)
consensus_expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data",
  consensus_features, trajectory$cell_id, missing = "NA"))
consensus_profiles <- as.matrix(
  (consensus_expression %*% window_design) %*% smoothing_design)
rownames(consensus_profiles) <- consensus_features
consensus_state_z <- row_standardize(consensus_profiles)$z
consensus_derivatives <- (consensus_profiles[, 3:nrow(windows), drop = FALSE] -
  consensus_profiles[, 1:(nrow(windows) - 2L), drop = FALSE]) / 2
consensus_change_z <- row_standardize(consensus_derivatives)$z

consensus_plot_rows <- lapply(seq_len(4L), function(i) {
  group <- unique(membership_rows[, c("approach", "direction")])[i, ]
  ids <- membership_rows$feature_id[membership_rows$approach == group$approach &
    membership_rows$direction == group$direction]
  sign <- if (group$direction == "RNA-faster") 1 else -1
  if (group$approach == "State") {
    consensus <- colMeans(consensus_state_z[ids, , drop = FALSE]) * sign
    pt <- windows$pt_median
    target_values <- target_z
  } else {
    consensus <- colMeans(consensus_change_z[ids, , drop = FALSE]) * sign
    pt <- windows$pt_median[2:(nrow(windows) - 1L)]
    target_values <- target_change_z
  }
  data.frame(approach = group$approach, direction = group$direction,
    pt = pt, consensus = consensus, target = target_values,
    stringsAsFactors = FALSE)
})
consensus_plot_data <- do.call(rbind, consensus_plot_rows)
consensus_plot_data$panel <- factor(paste(consensus_plot_data$approach,
  consensus_plot_data$direction, sep = " — "), levels = c(
    "State — RNA-faster", "State — ATAC-faster inverse",
    "Change — RNA-faster", "Change — ATAC-faster inverse"))
consensus_plot <- ggplot(consensus_plot_data, aes(pt)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "#666666") +
  geom_line(aes(y = target, color = "Modality target"), linewidth = 0.8) +
  geom_line(aes(y = consensus, color = "Oriented top-20 consensus"),
    linewidth = 0.85) +
  facet_wrap(~ panel, ncol = 2L, scales = "free_y") +
  scale_color_manual(values = c("Modality target" = "#222222",
    "Oriented top-20 consensus" = "#C51B5A")) +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(),
    strip.text = element_text(face = "bold"), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Top gene-pattern consensus profiles",
    subtitle = paste(
      "Inverse profiles are sign-flipped only in this display so agreement",
      "with the RNA/ATAC target is directly visible"),
    x = "Bdev pseudotime", y = "Standardized profile", color = NULL)
ggsave(file.path(stage_dir, "top20_gene_pattern_consensus.png"),
  consensus_plot, width = 11, height = 7.5, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "top20_gene_pattern_consensus.pdf"),
  consensus_plot, width = 11, height = 7.5, device = cairo_pdf, bg = "white")

results$selection_class <- "Not selected"
results$selection_class[results$state_selected & !results$change_selected] <-
  "State only"
results$selection_class[!results$state_selected & results$change_selected] <-
  "Change only"
results$selection_class[results$state_selected & results$change_selected] <-
  "Both"
scatter_plot <- ggplot(results, aes(state_correlation, change_correlation,
    color = selection_class)) +
  geom_hline(yintercept = 0, color = "#AAAAAA", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "#AAAAAA", linewidth = 0.4) +
  geom_point(size = 0.65, alpha = 0.65, stroke = 0) +
  scale_color_manual(values = c("Not selected" = "#BDBDBD",
    "State only" = "#0072B2", "Change only" = "#D55E00",
    "Both" = "#8E0152")) +
  theme_classic(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", size = 14),
    legend.position = "right") +
  labs(title = "State-following and change-following gene scores",
    subtitle = paste0("Selection: pooled-null BH < 0.05 and |correlation| ≥ ",
      minimum_abs_correlation),
    x = "State correlation", y = "Change correlation", color = NULL)
ggsave(file.path(stage_dir, "state_vs_change_gene_scores.png"),
  scatter_plot, width = 7.5, height = 6.8, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "state_vs_change_gene_scores.pdf"),
  scatter_plot, width = 7.5, height = 6.8, device = cairo_pdf, bg = "white")

counts <- list(features_in_sct_data = length(features),
  tested_features = nrow(results),
  state_selected = sum(results$state_selected),
  state_direct = sum(results$state_selected & results$state_correlation > 0),
  state_inverse = sum(results$state_selected & results$state_correlation < 0),
  change_selected = sum(results$change_selected),
  change_direct = sum(results$change_selected & results$change_correlation > 0),
  change_inverse = sum(results$change_selected & results$change_correlation < 0),
  selected_by_either = sum(results$state_selected | results$change_selected),
  selected_by_both = sum(results$state_selected & results$change_selected))
diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  seed = 20261009L, expression_matrix = "SCT_data",
  trajectory_cells = nrow(trajectory), windows = nrow(windows),
  window_size = 30L, stride = 5L,
  expression_smoother = "centered moving mean over ordered windows",
  smoothing_fraction = smoothing_fraction,
  smoothing_half_width_windows = smoothing_half_width,
  target = "log2(smoothed delta-RNA/delta-ATAC)",
  state_definition = "Pearson correlation of window-mean expression and target",
  change_definition = paste(
    "Pearson correlation of centered expression derivative over ordered",
    "window index and target; window-index derivative avoids tied-pt denominators"),
  inverse_definition = "negative score; no separate inverse fit",
  minimum_detection_fraction = minimum_detection_fraction,
  minimum_abs_correlation = minimum_abs_correlation,
  circular_shift_offsets_state = as.list(state_shifts$offsets[seq_len(n_shifts)]),
  circular_shift_offsets_change = as.list(change_shifts$offsets[seq_len(n_shifts)]),
  shifts_per_gene = n_shifts,
  p_value = paste(
    "two-sided tail against pooled absolute circular-shift null;",
    "each retained gene contributes equally"),
  multiple_testing = "BH across tested genes separately per approach",
  counts = counts, bundle_mutated = FALSE, capsule_registered = FALSE,
  limitation = paste(
    "descriptive screen; RNA latent trajectory was derived from the same",
    "expression modality, and feature-fold cross-fitting was not performed"))
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))

top_names <- function(x, column, decreasing) paste(head(x$feature_id[
  order(x[[column]], decreasing = decreasing)], 10L), collapse = ", ")
report <- c("# Preliminary modality-pattern gene screen", "",
  paste0("Tested ", nrow(results), " of ", length(features),
    " SCT_data genes after detection and variance filtering."), "",
  "## Results", "",
  paste0("- State-following selected: ", counts$state_selected, " (",
    counts$state_direct, " RNA-faster; ", counts$state_inverse,
    " ATAC-faster inverse)."),
  paste0("- Change-following selected: ", counts$change_selected, " (",
    counts$change_direct, " RNA-faster; ", counts$change_inverse,
    " ATAC-faster inverse)."),
  paste0("- Selected by either approach: ", counts$selected_by_either,
    "; selected by both: ", counts$selected_by_both, "."), "",
  "Top state RNA-faster genes:",
  paste(head(top_state_direct$feature_id, 10L), collapse = ", "), "",
  "Top state inverse genes:",
  paste(head(top_state_inverse$feature_id, 10L), collapse = ", "), "",
  "Top change RNA-faster genes:",
  paste(head(top_change_direct$feature_id, 10L), collapse = ", "), "",
  "Top change inverse genes:",
  paste(head(top_change_inverse$feature_id, 10L), collapse = ", "), "",
  "## Method and caution", "",
  paste("Expression was averaged in the accepted overlapping 30-cell windows",
    "at 5-cell stride, then passed through the same fixed centered 8%-window",
    "smoother for every gene. State scores correlate smoothed expression levels with the saved",
    "log2 RNA/ATAC target; change scores correlate centered expression",
    "derivatives over ordered window index with that target. Window index avoids",
    "zero denominators from legitimate ties in the preserved source pt."),
  paste("Two-sided empirical P values use pooled circular-shift nulls that",
    "preserve curve autocorrelation. BH correction is separate by approach;",
    paste0("selection additionally requires |correlation| ≥ ",
      minimum_abs_correlation, ".")),
  paste("Negative scores are the ATAC/RNA inverse result. A separate inverse",
    "fit would be algebraically redundant."),
  paste("This is descriptive and preliminary. Because the RNA latent path was",
    "derived from the same expression modality, feature-fold cross-fitting",
    "would be required for a stricter non-circular confirmation."))
writeLines(report, file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed modality gene-pattern screen at", output_dir, "\n")
