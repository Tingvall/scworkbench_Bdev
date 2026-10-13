#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 31_preliminary_turn_local_balanced_gene_screen.R",
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
all_event_definition <- read_tsv(file.path(turn_pattern_report,
  "turn_event_definition.tsv"))
event_definition <- all_event_definition[match(paste0("T", 1:3),
  all_event_definition$region_id), ]
stopifnot(nrow(trajectory) == 12642L, nrow(windows) == 2524L,
  all(table(membership$window_id) == 30L),
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

# Each turn is scored only inside its own fixed 100-window pre/post
# neighborhood. Intervening pseudotime outside these neighborhoods never
# enters a similarity score, and T4 is absent by design.
turn_intervals <- lapply(seq_len(3L), function(k)
  seq.int(event_definition$left_window_index[[k]],
    event_definition$right_window_index[[k]]))
scored_indices <- sort(unique(unlist(turn_intervals)))
needed_iqr_indices <- seq.int(min(scored_indices) - smoothing_half_width,
  max(scored_indices) + smoothing_half_width)
scored_smoothing_design <- smoothing_design[needed_iqr_indices,
  scored_indices, drop = FALSE]
scored_column_lookup <- match(scored_indices, needed_iqr_indices)
stopifnot(!anyNA(scored_column_lookup))

target <- log2(windows$smoothed_local_slope)
target_derivative <- (target[3:length(target)] -
  target[1:(length(target) - 2L)]) / 2
flank_definition <- do.call(rbind, lapply(seq_len(3L), function(k) rbind(
  data.frame(turn_index = k, region_id = paste0("T", k), flank = "pre",
    start = event_definition$left_window_index[[k]],
    end = event_definition$window_index[[k]]),
  data.frame(turn_index = k, region_id = paste0("T", k), flank = "post",
    start = event_definition$window_index[[k]],
    end = event_definition$right_window_index[[k]]))))
flank_target <- lapply(seq_len(nrow(flank_definition)), function(j) {
  interval <- seq.int(flank_definition$start[[j]], flank_definition$end[[j]])
  centers <- interval[interval > 1L & interval < nrow(windows)]
  derivative <- target_derivative[centers - 1L]
  weights <- abs(derivative)
  weights <- weights / sum(weights)
  endpoint_change <- target[flank_definition$end[[j]]] -
    target[flank_definition$start[[j]]]
  list(interval = interval, centers = centers, derivative = derivative,
    weights = weights, endpoint_change = endpoint_change)
})

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))
features <- wb_feature_ids(bundle, "SCT_data")
stopifnot(length(features) > 10000L, !anyDuplicated(features))

minimum_detection_fraction <- 0.01
flank_derivative_weight <- 0.70
flank_endpoint_weight <- 0.30
local_amplitude_half_saturation <- 0.25
within_turn_penalty_weight <- 0.50
consistency_penalty_weight <- 0.50
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
  derivatives <- (profiles[, 3:nrow(windows), drop = FALSE] -
    profiles[, 1:(nrow(windows) - 2L), drop = FALSE]) / 2

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
  smoothed_q25 <- as.matrix(raw_q25 %*% scored_smoothing_design)
  smoothed_q75 <- as.matrix(raw_q75 %*% scored_smoothing_design)
  smoothed_iqr <- smoothed_q75 - smoothed_q25

  flank_results <- lapply(seq_len(nrow(flank_definition)), function(j) {
    info <- flank_target[[j]]
    scored_columns <- match(info$interval, scored_indices)
    raw_columns <- match(info$interval, needed_iqr_indices)
    local_profiles <- profiles[, info$interval, drop = FALSE]
    q <- rowQuantiles(local_profiles, probs = c(0.05, 0.95), drop = FALSE)
    dynamic_range <- q[, 2L] - q[, 1L]
    derivative_cosine <- row_cosine(
      derivatives[, info$centers - 1L, drop = FALSE],
      info$derivative, info$weights)
    endpoint_change <- profiles[, flank_definition$end[[j]]] -
      profiles[, flank_definition$start[[j]]]
    oriented_endpoint <- sign(info$endpoint_change) * endpoint_change /
      pmax(dynamic_range, .Machine$double.eps)
    endpoint_score <- tanh(2 * oriented_endpoint)
    raw_shape <- flank_derivative_weight * derivative_cosine +
      flank_endpoint_weight * endpoint_score
    median_iqr <- rowMedians(smoothed_iqr[, scored_columns, drop = FALSE])
    median_sd <- rowMedians(raw_sd[, raw_columns, drop = FALSE])
    conservative_spread <- pmax(median_iqr, 1.349 * median_sd)
    precision <- dynamic_range / (dynamic_range + conservative_spread)
    precision[!is.finite(precision)] <- 0
    amplitude <- dynamic_range /
      (dynamic_range + local_amplitude_half_saturation)
    amplitude[!is.finite(amplitude)] <- 0
    evidence_weight <- precision * amplitude
    data.frame(dynamic_range = dynamic_range, median_iqr = median_iqr,
      precision = precision, amplitude_weight = amplitude,
      derivative_cosine = derivative_cosine,
      endpoint_score_direct = endpoint_score, raw_shape = raw_shape,
      evidence_weight = evidence_weight,
      evidence_score_direct = raw_shape * evidence_weight)
  })
  direct_flanks <- do.call(cbind, lapply(flank_results,
    function(x) x$evidence_score_direct))
  inverse_flanks <- -direct_flanks
  direct_turns <- do.call(cbind, lapply(seq_len(3L), function(k) {
    pair <- direct_flanks[, flank_definition$turn_index == k, drop = FALSE]
    rowMeans(pair) - within_turn_penalty_weight * rowSds(pair)
  }))
  inverse_turns <- do.call(cbind, lapply(seq_len(3L), function(k) {
    pair <- inverse_flanks[, flank_definition$turn_index == k, drop = FALSE]
    rowMeans(pair) - within_turn_penalty_weight * rowSds(pair)
  }))
  direct_mean <- rowMeans(direct_turns)
  direct_sd <- rowSds(direct_turns)
  inverse_mean <- rowMeans(inverse_turns)
  inverse_sd <- rowSds(inverse_turns)
  result <- data.frame(feature_id = ids,
    detection_fraction = detection_fraction[keep], stringsAsFactors = FALSE)
  for (j in seq_len(nrow(flank_definition))) {
    prefix <- paste0(tolower(flank_definition$region_id[[j]]), "_",
      flank_definition$flank[[j]])
    x <- flank_results[[j]]
    result[[paste0(prefix, "_dynamic_range_p95_p05")]] <- x$dynamic_range
    result[[paste0(prefix, "_median_smoothed_window_iqr")]] <- x$median_iqr
    result[[paste0(prefix, "_precision_reliability")]] <- x$precision
    result[[paste0(prefix, "_amplitude_weight")]] <- x$amplitude_weight
    result[[paste0(prefix, "_derivative_cosine_direct")]] <- x$derivative_cosine
    result[[paste0(prefix, "_endpoint_score_direct")]] <-
      x$endpoint_score_direct
    result[[paste0(prefix, "_raw_shape_score_direct")]] <- x$raw_shape
    result[[paste0(prefix, "_evidence_weight")]] <- x$evidence_weight
    result[[paste0(prefix, "_evidence_score_direct")]] <-
      x$evidence_score_direct
    result[[paste0(prefix, "_evidence_score_inverse")]] <-
      -x$evidence_score_direct
  }
  for (k in seq_len(3L)) {
    prefix <- tolower(event_definition$region_id[[k]])
    result[[paste0(prefix, "_balanced_turn_score_direct")]] <-
      direct_turns[, k]
    result[[paste0(prefix, "_balanced_turn_score_inverse")]] <-
      inverse_turns[, k]
  }
  result$mean_local_evidence_direct <- direct_mean
  result$sd_local_evidence_direct <- direct_sd
  result$balanced_local_score_direct <- direct_mean -
    consistency_penalty_weight * direct_sd
  result$mean_local_evidence_inverse <- inverse_mean
  result$sd_local_evidence_inverse <- inverse_sd
  result$balanced_local_score_inverse <- inverse_mean -
    consistency_penalty_weight * inverse_sd
  results_list[[length(results_list) + 1L]] <- result
  rm(expression, kept_expression, profiles, derivatives, raw_q25, raw_q75,
    raw_sd, smoothed_q25, smoothed_q75, smoothed_iqr)
  if (batch_number %% 5L == 0L) gc(verbose = FALSE)
}

results <- do.call(rbind, results_list)
stopifnot(nrow(results) > 10000L, !anyDuplicated(results$feature_id),
  all(is.finite(results$balanced_local_score_direct)),
  all(is.finite(results$balanced_local_score_inverse)))
results$direct_balanced_rank <- rank(-results$balanced_local_score_direct,
  ties.method = "min")
results$inverse_balanced_rank <- rank(-results$balanced_local_score_inverse,
  ties.method = "min")
results <- results[order(results$feature_id), ]
write_tsv_gz(results, file.path(stage_dir,
  "turn_local_balanced_gene_scores.tsv.gz"))
top_direct <- head(results[order(results$direct_balanced_rank,
  results$feature_id), ], 100L)
top_inverse <- head(results[order(results$inverse_balanced_rank,
  results$feature_id), ], 100L)
write_tsv(top_direct, file.path(stage_dir,
  "turn_local_balanced_top100_direct.tsv"))
write_tsv(top_inverse, file.path(stage_dir,
  "turn_local_balanced_top100_inverse.tsv"))

genes_of_prior_interest <- c("Nipbl", "Igll1", "Ebf1", "Lef1", "Smarca4",
  "Baz2b", "Mau2", "Wapl", "Ctcf", "Rad21", "Smc1a", "Smc3", "Stag1",
  "Stag2", "Pds5a", "Pds5b", "Ezh2", "Eed", "Suz12")
interest <- results[match(genes_of_prior_interest, results$feature_id), ]
stopifnot(!anyNA(interest$feature_id))
write_tsv(interest, file.path(stage_dir,
  "prior_interest_genes_turn_local_scores.tsv"))

top_n <- 20L
direct20 <- head(top_direct, top_n)
inverse20 <- head(top_inverse, top_n)
plot_features <- unique(c(direct20$feature_id, inverse20$feature_id))
plot_expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data",
  plot_features, trajectory$cell_id, missing = "NA"))
plot_profiles <- as.matrix((plot_expression %*% window_design) %*%
  smoothing_design)
plot_needed_indices <- seq_len(nrow(windows))
plot_q25 <- matrix(NA_real_, nrow = length(plot_features),
  ncol = length(plot_needed_indices))
plot_q75 <- plot_q25
for (j in seq_along(plot_needed_indices)) {
  ix <- membership_split[[plot_needed_indices[[j]]]]
  q <- rowQuantiles(plot_expression[, ix, drop = FALSE],
    probs = c(0.25, 0.75), drop = FALSE)
  plot_q25[, j] <- q[, 1L]
  plot_q75[, j] <- q[, 2L]
}
plot_mean <- plot_profiles
plot_smoothed_q25 <- as.matrix(plot_q25 %*% smoothing_design)
plot_smoothed_q75 <- as.matrix(plot_q75 %*% smoothing_design)
rownames(plot_mean) <- rownames(plot_smoothed_q25) <-
  rownames(plot_smoothed_q75) <- plot_features
plot_center <- rowMeans(plot_mean)
plot_sd <- sqrt(rowSums((plot_mean - plot_center)^2) /
  (ncol(plot_mean) - 1L))
plot_mean_z <- (plot_mean - plot_center) / plot_sd
plot_q25_z <- (plot_smoothed_q25 - plot_center) / plot_sd
plot_q75_z <- (plot_smoothed_q75 - plot_center) / plot_sd
direct_target_z <- as.numeric(scale(target))
inverse_target_z <- -direct_target_z
highlight_regions <- data.frame(region_id = event_definition$region_id,
  xmin = windows$pt_median[event_definition$left_window_index],
  xmax = windows$pt_median[event_definition$right_window_index])

make_plot_data <- function(top, direction) {
  rows <- lapply(seq_len(nrow(top)), function(i) {
    gene <- top$feature_id[[i]]
    local <- if (direction == "direct") c(
      top$t1_balanced_turn_score_direct[[i]],
      top$t2_balanced_turn_score_direct[[i]],
      top$t3_balanced_turn_score_direct[[i]]) else c(
      top$t1_balanced_turn_score_inverse[[i]],
      top$t2_balanced_turn_score_inverse[[i]],
      top$t3_balanced_turn_score_inverse[[i]])
    data.frame(pt = windows$pt_median,
      gene_panel = sprintf("%02d %s | T1 %.2f T2 %.2f T3 %.2f | bal %.2f",
        i, gene, local[[1L]], local[[2L]], local[[3L]],
        top[[paste0("balanced_local_score_", direction)]][[i]]),
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
    geom_rect(data = highlight_regions,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
        fill = region_id), inherit.aes = FALSE, alpha = 0.10, color = NA) +
    geom_hline(yintercept = 0, color = "#999999", linetype = "dotted",
      linewidth = 0.25) +
    geom_ribbon(aes(ymin = q25_z, ymax = q75_z), fill = color,
      alpha = 0.18, color = NA) +
    geom_line(aes(y = mean_z), color = color, linewidth = 0.8) +
    geom_line(aes(y = target_z), color = "#111111", linewidth = 0.55) +
    geom_vline(data = all_event_definition, aes(xintercept = pt),
      inherit.aes = FALSE, color = "#555555", linetype = "dashed",
      linewidth = 0.25) +
    facet_wrap(~ gene_panel, ncol = 4L) +
    coord_cartesian(ylim = c(-3.5, 3.5)) +
    scale_fill_manual(values = c("T1" = "#E41A1C", "T2" = "#377EB8",
      "T3" = "#4DAF4A"), name = "Scored ±100 windows") +
    theme_classic(base_size = 9) +
    theme(strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 7.5),
      panel.spacing = grid::unit(5, "mm"), legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 15)) +
    labs(title = paste("Top 20 turn-local balanced", label, "hits"),
      subtitle = paste(
        "Entire pseudotime is displayed; shaded ±100-window T1-T3 regions",
        "alone are scored; T4 is unscored context"),
      x = "Bdev pseudotime", y = "Standardized profile")
}
direct_plot <- make_figure(make_plot_data(direct20, "direct"), "direct")
inverse_plot <- make_figure(make_plot_data(inverse20, "inverse"), "inverse")
ggsave(file.path(stage_dir, "turn_local_balanced_top20_direct_with_iqr.png"),
  direct_plot, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "turn_local_balanced_top20_direct_with_iqr.pdf"),
  direct_plot, width = 14, height = 12, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "turn_local_balanced_top20_inverse_with_iqr.png"),
  inverse_plot, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "turn_local_balanced_top20_inverse_with_iqr.pdf"),
  inverse_plot, width = 14, height = 12, device = cairo_pdf, bg = "white")

heatmap_candidates <- unique(c(direct20$feature_id, "Nipbl", "Igll1"))
heatmap_scores <- results[match(heatmap_candidates, results$feature_id), ]
heatmap_rows <- do.call(rbind, lapply(seq_len(nrow(heatmap_scores)), function(i) {
  row <- heatmap_scores[i, ]
  data.frame(feature_id = row$feature_id, turn = factor(paste0("T", 1:3),
    levels = paste0("T", 1:3)), score = c(row$t1_balanced_turn_score_direct,
      row$t2_balanced_turn_score_direct, row$t3_balanced_turn_score_direct))
}))
heatmap_rows$feature_id <- factor(heatmap_rows$feature_id,
  levels = rev(heatmap_candidates))
local_heatmap <- ggplot(heatmap_rows, aes(turn, feature_id, fill = score)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", score)), size = 3.1,
    fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-max(abs(heatmap_rows$score)),
      max(abs(heatmap_rows$score)))) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
    axis.text = element_text(face = "bold"), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Per-turn direct evidence scores for unbiased top hits",
    subtitle = "Nipbl and Igll1 are included for direct comparison",
    fill = "Local score")
ggsave(file.path(stage_dir, "turn_local_direct_score_heatmap.png"),
  local_heatmap, width = 6.5, height = 9, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "turn_local_direct_score_heatmap.pdf"),
  local_heatmap, width = 6.5, height = 9, device = cairo_pdf, bg = "white")

nipbl <- results[results$feature_id == "Nipbl", ]
igll1 <- results[results$feature_id == "Igll1", ]
diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  expression_matrix = "SCT_data", trajectory_cells = nrow(trajectory),
  windows = nrow(windows), window_size = 30L, stride = 5L,
  turns = as.list(paste0("T", 1:3)), excluded_turn = "T4",
  local_neighborhood = "100 windows before through 100 windows after each turn",
  full_pseudotime_or_interturn_similarity = FALSE,
  local_raw_shape = paste0(flank_derivative_weight,
    " * within-flank derivative cosine + ", flank_endpoint_weight,
    " * continuous oriented endpoint score"),
  local_amplitude_half_saturation = local_amplitude_half_saturation,
  local_precision = "dynamic range relative to max(smoothed IQR, 1.349 * window SD)",
  aggregation = paste(
    "within each turn: equal pre/post mean minus", within_turn_penalty_weight,
    "times pre/post SD; across turns: equal-turn mean minus",
    consistency_penalty_weight, "times between-turn SD"),
  plot_scope = "entire pseudotime with scored T1-T3 +/-100-window regions shaded",
  hard_turn_or_flank_cutoffs = FALSE,
  tested_genes = nrow(results), nipbl = as.list(nipbl[1, ]),
  igll1 = as.list(igll1[1, ]), direct_top20 = as.list(direct20$feature_id),
  inverse_top20 = as.list(inverse20$feature_id),
  bundle_mutated = FALSE, capsule_registered = FALSE,
  inference = "unbiased genome-wide descriptive ranking; no significance cutoff")
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
report <- c("# Turn-local balanced T1-T3 gene screen", "",
  paste0("All ", nrow(results),
    " detected genes were ranked without gene-class preselection."), "",
  "## Method", "",
  paste("The six T1-T3 pre/post flanks are scored independently, using only",
    "their own fixed 100-window neighborhoods. Intervening trajectory and T4",
    "are not used. Each flank score combines derivative similarity and a",
    "continuous oriented endpoint score, then is smoothly weighted by local",
    "amplitude and within-window precision."),
  paste("Pre and post flanks receive equal weight within each turn, with a",
    "continuous imbalance penalty. T1, T2, and T3 then receive equal weight",
    "with a second between-turn consistency penalty. No binary cutoff is used.",
    "Figures display the full pseudotime range and shade the scored regions."), "",
  "Top direct genes:", paste(direct20$feature_id, collapse = ", "), "",
  "Top inverse genes:", paste(inverse20$feature_id, collapse = ", "), "",
  paste0("Nipbl direct balanced rank: ", nipbl$direct_balanced_rank,
    "; T1/T2/T3 balanced-turn evidence: ", paste(sprintf("%.3f", c(
      nipbl$t1_balanced_turn_score_direct, nipbl$t2_balanced_turn_score_direct,
      nipbl$t3_balanced_turn_score_direct)), collapse = ", "), "."),
  paste0("Igll1 direct balanced rank: ", igll1$direct_balanced_rank,
    "; T1/T2/T3 balanced-turn evidence: ", paste(sprintf("%.3f", c(
      igll1$t1_balanced_turn_score_direct, igll1$t2_balanced_turn_score_direct,
      igll1$t3_balanced_turn_score_direct)), collapse = ", "), "."), "",
  "This is a descriptive ranking, not a significance test.")
writeLines(report, file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed turn-local balanced gene screen at", output_dir, "\n")
