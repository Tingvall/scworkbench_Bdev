#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 10_all_cytokine_balanced_score_comparison.R",
  "<bdev_bundle> <existing_38_profile_ucell_report>",
  "<existing_driver_report> <output_dir> <expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
ucell_report <- normalizePath(args[[2L]], mustWork = TRUE)
driver_report <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_fingerprint <- args[[5L]]
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(stage_dir, "profile_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(patchwork)
})
source("src/sliding_window_signature_summary.R")

read_tsv <- function(path) read.delim(path, check.names = FALSE,
  na.strings = "NA", stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
slug <- function(x) gsub("(^_+|_+$)", "",
  tolower(gsub("[^A-Za-z0-9]+", "_", x)))
as_flag <- function(x) !is.na(x) & tolower(as.character(x)) %in%
  c("true", "t", "1", "yes")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}

ucell_inputs <- c("completion.json", "passing_signature_inventory.tsv",
  "passing_gene_membership.tsv", "ucell_scores.rds",
  "sliding_window_summary.tsv", "window_support.tsv",
  "window_membership.rds", "trajectory_parameters.json",
  "bb_up_change_points.tsv")
driver_inputs <- c("completion.json", "top5_driver_genes.tsv",
  "top5_driver_mean_expression_curves.tsv", "checks.json")
ucell_paths <- file.path(ucell_report, ucell_inputs)
driver_paths <- file.path(driver_report, driver_inputs)
stopifnot(all(file.exists(ucell_paths)), all(file.exists(driver_paths)))
ucell_completion <- read_json(ucell_paths[[1L]], simplifyVector = TRUE)
driver_completion <- read_json(driver_paths[[1L]], simplifyVector = TRUE)
stopifnot(identical(ucell_completion$status, "complete"),
  identical(ucell_completion$promoted, FALSE),
  identical(as.integer(ucell_completion$passing_signatures), 38L),
  identical(driver_completion$status, "complete"),
  identical(driver_completion$promoted, FALSE),
  identical(as.integer(driver_completion$profiles), 38L),
  identical(driver_completion$driver_summary, "mean"),
  identical(driver_completion$driver_smoothing, "none"),
  identical(driver_completion$driver_plot_scale,
    "per-gene min-max scaling of mean trajectory to [0,1]"))

inventory <- read_tsv(file.path(ucell_report,
  "passing_signature_inventory.tsv"))
members <- read_tsv(file.path(ucell_report, "passing_gene_membership.tsv"))
old_family <- readRDS(file.path(ucell_report, "ucell_scores.rds"))
old_scores <- data.frame(key = as.character(old_family$key),
  cell_id = as.character(old_family$cell_id),
  value = as.numeric(old_family$value), stringsAsFactors = FALSE)
old_trajectory <- read_tsv(file.path(ucell_report,
  "sliding_window_summary.tsv"))
old_support <- read_tsv(file.path(ucell_report, "window_support.tsv"))
old_membership <- readRDS(file.path(ucell_report, "window_membership.rds"))
old_parameters <- read_json(file.path(ucell_report,
  "trajectory_parameters.json"), simplifyVector = TRUE)
change_points <- read_tsv(file.path(ucell_report,
  "bb_up_change_points.tsv"))
change_points$boundary_label <- paste0("BB", seq_len(nrow(change_points)))
driver_table <- read_tsv(file.path(driver_report, "top5_driver_genes.tsv"))
driver_curves <- read_tsv(file.path(driver_report,
  "top5_driver_mean_expression_curves.tsv"))

stopifnot(nrow(inventory) == 38L, all(as_flag(inventory$final_overall_pass)),
  all(inventory$final_size >= 20L), !anyDuplicated(inventory$signature_id),
  setequal(unique(members$signature_id), inventory$signature_id),
  !anyDuplicated(members[c("signature_id", "feature_id")]),
  setequal(unique(old_scores$key), inventory$signature_id),
  setequal(unique(old_trajectory$feature_set), inventory$signature_id),
  setequal(unique(driver_table$feature_set), inventory$signature_id),
  setequal(unique(driver_curves$feature_set), inventory$signature_id),
  "plot_value_internal_0_1" %in% names(driver_curves),
  all(is.finite(driver_curves$plot_value_internal_0_1)),
  all(driver_curves$plot_value_internal_0_1 >= 0 &
    driver_curves$plot_value_internal_0_1 <= 1))
inventory <- inventory[order(inventory$cytokine,
  match(inventory$direction, c("up", "down"))), , drop = FALSE]
members <- members[order(match(members$signature_id,
  inventory$signature_id), members$final_rank), , drop = FALSE]

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint),
  identical(ucell_completion$bdev_bundle_fingerprint, fingerprint_before),
  identical(driver_completion$bundle_fingerprint, fingerprint_before))
cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "pt"))
stopifnot(identical(as.character(metadata$cell_id), cells),
  length(cells) == 26723L, sum(is.finite(metadata$pt)) == 12642L,
  sum(is.na(metadata$pt)) == 14081L)
for (id in inventory$signature_id) {
  local <- old_scores[old_scores$key == id, ]
  stopifnot(identical(local$cell_id, cells), !anyDuplicated(local$cell_id))
}

# User-authorized analysis outside a single built-in scWorkbench score op.
# Each gene is transformed independently: zeros remain zero; positive SCT_data
# values receive their average percentile rank among positive cells. Equal
# averaging then caps every gene's possible contribution at 1 / set size.
genes <- unique(members$feature_id)
stopifnot(all(genes %in% wb_feature_ids(bundle, "SCT_data")))
log_event("Reading SCT_data for", length(genes), "unique signature genes")
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", genes,
  cells, missing = "zero"))
stopifnot(identical(rownames(expression), genes),
  identical(colnames(expression), cells), all(is.finite(expression)),
  all(expression >= 0))

zero_aware_percentile <- function(x) {
  out <- numeric(length(x))
  positive <- x > 0
  if (any(positive)) {
    out[positive] <- rank(x[positive], ties.method = "average") / sum(positive)
  }
  out
}
log_event("Calculating equal-weight zero-aware gene percentiles")
gene_percentiles <- t(apply(expression, 1L, zero_aware_percentile))
dimnames(gene_percentiles) <- dimnames(expression)
gene_qc <- data.frame(feature_id = genes,
  detected_cells = rowSums(expression > 0),
  detection_fraction = rowMeans(expression > 0),
  expression_mean = rowMeans(expression),
  expression_max = apply(expression, 1L, max),
  percentile_mean = rowMeans(gene_percentiles),
  stringsAsFactors = FALSE)
stopifnot(all(gene_percentiles >= 0 & gene_percentiles <= 1),
  all(gene_percentiles[expression == 0] == 0))

balanced_rows <- vector("list", nrow(inventory))
weight_rows <- vector("list", nrow(inventory))
for (i in seq_len(nrow(inventory))) {
  id <- inventory$signature_id[[i]]
  set_genes <- members$feature_id[members$signature_id == id]
  weights <- rep(1 / length(set_genes), length(set_genes))
  values <- as.numeric(crossprod(weights,
    gene_percentiles[set_genes, , drop = FALSE]))
  balanced_rows[[i]] <- data.frame(key = id, cell_id = cells,
    value = values, stringsAsFactors = FALSE)
  weight_rows[[i]] <- data.frame(feature_set = id,
    feature_id = set_genes, weight = weights,
    maximum_single_gene_contribution = weights,
    stringsAsFactors = FALSE)
}
balanced_scores <- do.call(rbind, balanced_rows)
weights <- do.call(rbind, weight_rows)
rownames(balanced_scores) <- rownames(weights) <- NULL
stopifnot(nrow(balanced_scores) == 38L * length(cells),
  all(is.finite(balanced_scores$value)),
  all(balanced_scores$value >= 0 & balanced_scores$value <= 1),
  all(abs(ave(weights$weight, weights$feature_set, FUN = sum) - 1) < 1e-12))
write_tsv(gene_qc, file.path(stage_dir, "gene_expression_qc.tsv"))
write_tsv(weights, file.path(stage_dir, "equal_gene_weights.tsv"))
saveRDS(balanced_scores, file.path(stage_dir,
  "balanced_cell_scores.rds"), compress = "gzip")
write_tsv_gz(balanced_scores, file.path(stage_dir,
  "balanced_cell_scores.tsv.gz"))

pt <- data.frame(cell_id = cells[is.finite(metadata$pt)],
  value = as.numeric(metadata$pt[is.finite(metadata$pt)]),
  stringsAsFactors = FALSE)
log_event("Calculating native sliding-window summaries for balanced scores")
balanced_window <- sliding_window_signature_summary(balanced_scores, pt,
  width_fractions = 0.05, step_fraction = 0.01, min_cells = 20L)
balanced_trajectory <- balanced_window$summaries
balanced_support <- balanced_window$windows
balanced_trajectory$window_midpoint <-
  (balanced_trajectory$window_left + balanced_trajectory$window_right) / 2
balanced_support$window_midpoint <-
  (balanced_support$window_left + balanced_support$window_right) / 2
support_matches <- isTRUE(all.equal(
  balanced_support[, names(old_support), drop = FALSE], old_support,
  tolerance = 1e-12, check.attributes = FALSE
))
stopifnot(identical(balanced_window$membership, old_membership),
  support_matches,
  setequal(unique(balanced_trajectory$feature_set), inventory$signature_id))
write_tsv(balanced_trajectory, file.path(stage_dir,
  "balanced_sliding_window_summary.tsv"))
write_tsv(balanced_support, file.path(stage_dir, "window_support.tsv"))
saveRDS(balanced_window$membership, file.path(stage_dir,
  "window_membership.rds"))
write_json(balanced_window$parameters, file.path(stage_dir,
  "balanced_trajectory_parameters.json"))

# Method 3 is deliberately trajectory-level. Within every existing rolling
# window, calculate the mean expression of each signature gene. Normalize each
# gene's trajectory internally to [0,1], then use the median across genes as a
# robust consensus that no single gene can determine.
log_event("Calculating internally scaled gene-trajectory consensus scores")
window_ids <- as.character(old_support$window_id)
stopifnot(identical(names(old_membership), window_ids))
window_cell_index <- lapply(old_membership, match, table = cells)
stopifnot(!anyNA(unlist(window_cell_index, use.names = FALSE)))
eligible_windows <- as.logical(old_support$eligible_window)
scale_trajectory <- function(x) {
  out <- rep(NA_real_, length(x))
  limits <- range(x[eligible_windows])
  if (diff(limits) == 0) {
    out[eligible_windows] <- 0
  } else {
    out[eligible_windows] <-
      (x[eligible_windows] - limits[[1L]]) / diff(limits)
  }
  out
}
across_genes <- function(x, fun) {
  apply(x, 2L, function(v) {
    if (all(is.na(v))) return(NA_real_)
    fun(v[is.finite(v)])
  })
}
consensus_gene_rows <- vector("list", nrow(inventory))
consensus_summary_rows <- vector("list", nrow(inventory))
for (i in seq_len(nrow(inventory))) {
  id <- inventory$signature_id[[i]]
  set_genes <- members$feature_id[members$signature_id == id]
  mean_expression <- vapply(window_cell_index, function(ix) {
    rowMeans(expression[set_genes, ix, drop = FALSE])
  }, numeric(length(set_genes)))
  rownames(mean_expression) <- set_genes
  colnames(mean_expression) <- window_ids
  internally_scaled <- t(apply(mean_expression, 1L, scale_trajectory))
  dimnames(internally_scaled) <- dimnames(mean_expression)
  consensus_gene_rows[[i]] <- data.frame(feature_set = id,
    feature_id = rep(set_genes, each = length(window_ids)),
    window_id = rep(window_ids, times = length(set_genes)),
    window_midpoint = rep(old_support$window_midpoint,
      times = length(set_genes)),
    eligible_window = rep(eligible_windows, times = length(set_genes)),
    mean_expression = as.vector(t(mean_expression)),
    internally_scaled_0_1 = as.vector(t(internally_scaled)),
    stringsAsFactors = FALSE)
  consensus_summary_rows[[i]] <- data.frame(feature_set = id,
    window_id = window_ids,
    window_midpoint = old_support$window_midpoint,
    eligible_window = eligible_windows,
    n_cells = old_support$n_cells,
    consensus_median = across_genes(internally_scaled, median),
    consensus_q25 = across_genes(internally_scaled,
      function(v) quantile(v, 0.25, names = FALSE, type = 7)),
    consensus_q75 = across_genes(internally_scaled,
      function(v) quantile(v, 0.75, names = FALSE, type = 7)),
    consensus_trimmed_mean = across_genes(internally_scaled,
      function(v) mean(v, trim = 0.1)),
    genes_with_nonzero_window_mean = colSums(mean_expression > 0),
    signature_size = length(set_genes), stringsAsFactors = FALSE)
}
consensus_gene_curves <- do.call(rbind, consensus_gene_rows)
consensus_trajectory <- do.call(rbind, consensus_summary_rows)
rownames(consensus_gene_curves) <- rownames(consensus_trajectory) <- NULL
stopifnot(setequal(unique(consensus_trajectory$feature_set),
    inventory$signature_id),
  all(is.finite(consensus_gene_curves$mean_expression)),
  all(consensus_gene_curves$mean_expression >= 0),
  all(consensus_gene_curves$internally_scaled_0_1[
    consensus_gene_curves$eligible_window] >= 0 &
    consensus_gene_curves$internally_scaled_0_1[
      consensus_gene_curves$eligible_window] <= 1),
  all(is.na(consensus_gene_curves$internally_scaled_0_1[
    !consensus_gene_curves$eligible_window])),
  all(is.finite(consensus_trajectory$consensus_median[
    consensus_trajectory$eligible_window])))
write_tsv(consensus_gene_curves, file.path(stage_dir,
  "trajectory_consensus_gene_curves.tsv"))
write_tsv(consensus_trajectory, file.path(stage_dir,
  "trajectory_consensus_summary.tsv"))

score_summary <- do.call(rbind, lapply(inventory$signature_id, function(id) {
  old <- old_scores$value[old_scores$key == id]
  new <- balanced_scores$value[balanced_scores$key == id]
  data.frame(feature_set = id,
    old_ucell_mean = mean(old), old_ucell_sd = sd(old),
    balanced_mean = mean(new), balanced_sd = sd(new),
    pearson = cor(old, new, method = "pearson"),
    spearman = cor(old, new, method = "spearman"),
    stringsAsFactors = FALSE)
}))
write_tsv(score_summary, file.path(stage_dir, "score_comparison_summary.tsv"))

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) {
    return(NA_real_)
  }
  cor(x[keep], y[keep], method = "spearman")
}
trajectory_comparison <- do.call(rbind, lapply(inventory$signature_id,
  function(id) {
    consensus <- consensus_trajectory[
      consensus_trajectory$feature_set == id &
        consensus_trajectory$eligible_window, ]
    old <- old_trajectory[old_trajectory$feature_set == id &
      old_trajectory$eligible_window, ]
    balanced <- balanced_trajectory[balanced_trajectory$feature_set == id &
      balanced_trajectory$eligible_window, ]
    old <- old[match(consensus$window_id, old$window_id), ]
    balanced <- balanced[match(consensus$window_id, balanced$window_id), ]
    stopifnot(!anyNA(old$score_median), !anyNA(balanced$score_median))
    data.frame(feature_set = id,
      spearman_consensus_vs_ucell = safe_spearman(
        consensus$consensus_median, old$score_median),
      spearman_consensus_vs_balanced = safe_spearman(
        consensus$consensus_median, balanced$score_median),
      consensus_peak_pt = consensus$window_midpoint[
        which.max(consensus$consensus_median)], stringsAsFactors = FALSE)
  }))
stopifnot(all(is.finite(trajectory_comparison$spearman_consensus_vs_ucell) |
    is.na(trajectory_comparison$spearman_consensus_vs_ucell)),
  all(is.finite(trajectory_comparison$spearman_consensus_vs_balanced) |
    is.na(trajectory_comparison$spearman_consensus_vs_balanced)))
write_tsv(trajectory_comparison, file.path(stage_dir,
  "three_method_trajectory_comparison.tsv"))

space <- wb_space(bundle, "umap.atac", cells = cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells))
  coordinates <- data.frame(cell_id = cells, x = as.numeric(space$x),
    y = as.numeric(space$y), stringsAsFactors = FALSE)
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  coordinates <- data.frame(cell_id = cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]), stringsAsFactors = FALSE)
}

map_plot <- function(values, title, subtitle) {
  positive <- values > 0
  cap <- quantile(values[positive], 0.99, names = FALSE, type = 7)
  display <- rep(NA_real_, length(values))
  display[positive] <- pmin(values[positive], cap) / cap
  data <- coordinates
  data$display <- display
  data <- data[order(ifelse(is.na(data$display), -Inf, data$display)), ]
  ggplot(data, aes(x, y)) +
    geom_point(data = data[is.na(data$display), ], color = "#CDD2DA",
      size = 0.11, alpha = 0.7, stroke = 0) +
    geom_point(data = data[!is.na(data$display), ], aes(color = display),
      size = 0.26, alpha = 1, stroke = 0) +
    scale_color_gradientn(colours = c("#0066FF", "#00D9FF", "#00F060",
      "#FFF000", "#FF5A00", "#FF00A8"), limits = c(0, 1),
      oob = scales::squish, name = "Within-panel\nscore") +
    coord_equal() + theme_void(base_size = 9) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
    labs(title = title, subtitle = subtitle,
      caption = "Positive values scaled to within-panel p99; zero is gray.")
}

unsupported <- old_support[!old_support$eligible_window, , drop = FALSE]
unsupported$region_left <- unsupported$window_left
unsupported$region_right <- unsupported$window_right
trajectory_plot <- function(data, title, color, y_label) {
  data <- data[data$eligible_window, , drop = FALSE]
  ggplot(data, aes(window_midpoint, score_median)) +
    geom_rect(data = unsupported,
      aes(xmin = region_left, xmax = region_right, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "#64748B", alpha = 0.13) +
    geom_vline(data = change_points, aes(xintercept = pt_median),
      inherit.aes = FALSE, linewidth = 0.45, linetype = "dotted",
      color = "#111827") +
    geom_ribbon(aes(ymin = score_q25, ymax = score_q75), fill = color,
      alpha = 0.22, colour = NA) +
    geom_line(color = color, linewidth = 0.85) +
    labs(title = title, x = "Existing Bdev pseudotime (pt)", y = y_label,
      caption = "5% windows; 1% steps; median and IQR; no smoothing.") +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"))
}

consensus_trajectory_plot <- function(data) {
  data <- data[data$eligible_window, , drop = FALSE]
  ggplot(data, aes(window_midpoint, consensus_median)) +
    geom_rect(data = unsupported,
      aes(xmin = region_left, xmax = region_right, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "#64748B", alpha = 0.13) +
    geom_vline(data = change_points, aes(xintercept = pt_median),
      inherit.aes = FALSE, linewidth = 0.45, linetype = "dotted",
      color = "#111827") +
    geom_ribbon(aes(ymin = consensus_q25, ymax = consensus_q75),
      fill = "#00897B", alpha = 0.22, colour = NA) +
    geom_line(color = "#00796B", linewidth = 0.9) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "New: gene-trajectory consensus",
      x = "Existing Bdev pseudotime (pt)",
      y = "Consensus (gene-internal 0-1)",
      caption = "Median and IQR across all genes; single genes cannot set the line.") +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"))
}

consensus_method_panel <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.68, label = "Trajectory consensus",
      fontface = "bold", size = 5, color = "#00796B") +
    annotate("text", x = 0.5, y = 0.43,
      label = paste("No per-cell UMAP\n5% windows; all", label,
        "genes\nper-gene mean trajectory -> 0-1\nmedian across genes"),
      size = 3.5, lineheight = 1.25) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_void()
}

gene_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
driver_plot <- function(id) {
  data <- driver_curves[driver_curves$feature_set == id, , drop = FALSE]
  ranks <- driver_table[driver_table$feature_set == id,
    c("feature_id", "driver_rank"), drop = FALSE]
  ranks <- ranks[order(ranks$driver_rank), ]
  genes <- ranks$feature_id
  data$feature_id <- factor(data$feature_id, levels = genes)
  colors <- setNames(gene_palette, genes)
  ends <- do.call(rbind, lapply(genes, function(gene) {
    x <- data[as.character(data$feature_id) == gene, , drop = FALSE]
    x[which.max(x$pseudotime), , drop = FALSE]
  }))
  span <- diff(range(data$pseudotime))
  ggplot(data, aes(pseudotime, plot_value_internal_0_1, color = feature_id,
      group = feature_id)) +
    geom_vline(data = change_points, aes(xintercept = pt_median),
      inherit.aes = FALSE, linewidth = 0.7, linetype = "dashed",
      color = "#374151") +
    geom_text(data = change_points,
      aes(x = pt_median, y = Inf, label = boundary_label),
      inherit.aes = FALSE, vjust = 1.35, fontface = "bold", size = 3,
      color = "#374151") +
    geom_line(linewidth = 0.9) +
    geom_text_repel(data = ends, aes(label = feature_id), direction = "y",
      hjust = 0, nudge_x = 0.02 * span, seed = 1729,
      box.padding = 0.2, point.padding = 0.15, min.segment.length = 0,
      max.overlaps = Inf, show.legend = FALSE) +
    scale_color_manual(values = colors, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.20))) +
    coord_cartesian(clip = "off") +
    labs(title = "Existing top-five member-gene drivers",
      subtitle = paste("Mean-expression dynamic-score rank; each mean curve",
        "scaled internally to 0-1"),
      x = "Existing Bdev pseudotime (pt)",
      y = "Within-gene scaled mean expression",
      caption = paste("Dashed lines mark BB1/BB2 medians.",
        "Raw mean and median SCT_data curves remain in source tables.")) +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(),
      legend.position = "none", plot.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 75, 5.5, 5.5))
}

log_event("Rendering 38 three-method comparison figures")
pdf_path <- file.path(figure_dir,
  "all_38_profiles_three_method_comparison_with_drivers.pdf")
grDevices::pdf(pdf_path, width = 18, height = 13.5, onefile = TRUE)
figure_rows <- vector("list", nrow(inventory))
for (i in seq_len(nrow(inventory))) {
  row <- inventory[i, ]
  id <- row$signature_id
  old_values <- old_scores$value[old_scores$key == id]
  new_values <- balanced_scores$value[balanced_scores$key == id]
  old_curve <- old_trajectory[old_trajectory$feature_set == id, ]
  new_curve <- balanced_trajectory[balanced_trajectory$feature_set == id, ]
  consensus_curve <- consensus_trajectory[
    consensus_trajectory$feature_set == id, ]
  direction_color <- if (row$direction == "up") "#D73027" else "#2166AC"
  label <- paste(row$cytokine, toupper(row$direction), "-",
    row$final_size, "genes")
  top <- map_plot(old_values, "Method 1: UCell UMAP", label) |
    map_plot(new_values, "Method 2: balanced-score UMAP", label) |
    consensus_method_panel(as.integer(row$final_size))
  middle <- trajectory_plot(old_curve, "Old: UCell rolling trajectory",
    direction_color, "UCell score") |
    trajectory_plot(new_curve, "New: balanced rolling trajectory",
      "#7B3294", "Balanced score") |
    consensus_trajectory_plot(consensus_curve)
  combined <- top / middle / driver_plot(id) +
    plot_layout(heights = c(1, 0.9, 0.8)) +
    plot_annotation(title = paste("Three-method comparison:", label),
      subtitle = paste("UCell; per-cell equal-gene balanced score;",
        "robust consensus across internally scaled gene trajectories"),
      theme = theme(plot.title = element_text(size = 17, face = "bold")))
  print(combined)
  file <- paste0(sprintf("%02d_", i), slug(paste(row$cytokine,
    row$direction)), "_three_method_comparison.png")
  ggsave(file.path(figure_dir, file), combined, width = 18, height = 13.5,
    dpi = 220, units = "in", bg = "white")
  comparison <- score_summary[score_summary$feature_set == id, ]
  trajectory_match <- trajectory_comparison[
    trajectory_comparison$feature_set == id, ]
  figure_rows[[i]] <- data.frame(signature_id = id,
    cytokine = row$cytokine, direction = row$direction,
    signature_size = as.integer(row$final_size),
    pearson_old_vs_new = comparison$pearson,
    spearman_old_vs_new = comparison$spearman,
    spearman_consensus_vs_ucell =
      trajectory_match$spearman_consensus_vs_ucell,
    spearman_consensus_vs_balanced =
      trajectory_match$spearman_consensus_vs_balanced,
    consensus_peak_pt = trajectory_match$consensus_peak_pt,
    file = file, stringsAsFactors = FALSE)
}
grDevices::dev.off()
figure_inventory <- do.call(rbind, figure_rows)
write_tsv(figure_inventory, file.path(stage_dir, "figure_inventory.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
checks <- list(source_has_38_passing_profiles = nrow(inventory) == 38L,
  old_ucell_values_unchanged = identical(old_scores$value,
    as.numeric(old_family$value)),
  equal_weights_sum_to_one = all(abs(ave(weights$weight,
    weights$feature_set, FUN = sum) - 1) < 1e-12),
  balanced_scores_finite_and_bounded = all(is.finite(balanced_scores$value)) &&
    all(balanced_scores$value >= 0 & balanced_scores$value <= 1),
  zero_expression_maps_to_zero =
    all(gene_percentiles[expression == 0] == 0),
  driver_plot_values_are_internal_0_1 =
    all(driver_curves$plot_value_internal_0_1 >= 0 &
      driver_curves$plot_value_internal_0_1 <= 1),
  consensus_uses_all_signature_genes =
    nrow(consensus_gene_curves) == nrow(members) * nrow(old_support),
  consensus_gene_curves_bounded =
    all(consensus_gene_curves$internally_scaled_0_1[
      consensus_gene_curves$eligible_window] >= 0 &
      consensus_gene_curves$internally_scaled_0_1[
        consensus_gene_curves$eligible_window] <= 1),
  consensus_summaries_finite = all(is.finite(
    consensus_trajectory$consensus_median[
      consensus_trajectory$eligible_window])),
  source_window_support_reused = support_matches,
  exact_window_membership_reused =
    identical(balanced_window$membership, old_membership),
  exact_source_pt_preserved = identical(pt$value,
    as.numeric(metadata$pt[is.finite(metadata$pt)])),
  source_missing_pt_preserved = sum(is.na(metadata$pt)) == 14081L,
  every_comparison_figure_exists = all(file.exists(file.path(figure_dir,
    figure_inventory$file))), multipage_pdf_exists = file.exists(pdf_path),
  bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after), result_promotions = 0L)
stopifnot(all(unlist(checks[names(checks) != "result_promotions"])),
  checks$result_promotions == 0L)
write_json(list(status = "passed", checks = checks),
  file.path(stage_dir, "checks.json"))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
provenance_inputs <- c(ucell_paths, driver_paths,
  file.path(bundle_dir, "manifest.json"), script_path,
  "figure_scripts/bdev_reference/run_10_all_cytokine_balanced_score_comparison.sh",
  "src/sliding_window_signature_summary.R", "config/bdev_reference.yml")
input_manifest <- data.frame(path = normalizePath(provenance_inputs,
  mustWork = TRUE), sha256 = vapply(provenance_inputs, wb_source_sha256,
  character(1)), stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

html_rows <- vapply(seq_len(nrow(figure_inventory)), function(i) {
  x <- figure_inventory[i, ]
  paste0("<h2>", x$cytokine, " ", toupper(x$direction), "</h2>",
    "<p>Old/new Spearman correlation: ",
    sprintf("%.3f", x$spearman_old_vs_new), "</p>",
    '<img src="profile_figures/', x$file, '" alt="', x$cytokine, " ",
    x$direction, ' three-method score comparison">')
}, character(1))
html <- c('<!doctype html><html><head><meta charset="utf-8">',
  '<title>Three-method cytokine trajectory comparison</title>',
  '<style>body{font:16px system-ui,sans-serif;max-width:1450px;margin:32px auto;padding:0 22px;color:#172033}img{max-width:100%;height:auto;border:1px solid #ddd}p,li{line-height:1.5}code{overflow-wrap:anywhere}a{color:#1261a0}</style></head><body>',
  '<h1>Three-method comparison for 38 cytokine profiles</h1>',
  '<p>Method 1 preserves the original native RNA-count UCell score. Method 2 transforms every signature gene independently: zero expression remains zero and positive <code>SCT_data</code> values become their percentile among positive cells; the per-cell score is the equal-weight mean across genes.</p>',
  '<p>Method 3 is trajectory-level: in every existing 5% pseudotime window it calculates each gene\'s mean <code>SCT_data</code> expression, scales that gene\'s trajectory to 0-1, and takes the median across all signature genes. The ribbon is the gene-wise IQR. A single gene cannot determine this median. Because this method has no per-cell score, no third UMAP is defined.</p>',
  '<p>The five driver identities and ranking are unchanged. For display, each driver mean-expression trajectory is independently min-max scaled to 0-1 so sparse, low-expression genes remain visible; raw mean and median curves remain in source tables.</p>',
  '<p>Dashed BB1/BB2 lines in every driver panel mark the saved median B-cell boundary positions. These custom sensitivity methods were explicitly authorized; matrix access is from the accepted read-only bundle and the exact existing windows and memberships are reused.</p>',
  '<p><a href="profile_figures/all_38_profiles_three_method_comparison_with_drivers.pdf">Download 38-page PDF</a> · <a href="three_method_trajectory_comparison.tsv">Trajectory comparison table</a> · <a href="trajectory_consensus_summary.tsv">Consensus summaries</a> · <a href="trajectory_consensus_gene_curves.tsv">Gene trajectories</a> · <a href="balanced_cell_scores.tsv.gz">Balanced cell scores</a></p>',
  html_rows, '</body></html>')
writeLines(html, file.path(stage_dir, "REPORT.html"))
writeLines(c("# Three-method cytokine trajectory comparison", "",
  "[Open report](REPORT.html)", "",
  "All 38 profiles compare native UCell, the per-cell equal-gene balanced score, and a robust median consensus across internally scaled gene trajectories."),
  file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", promoted = FALSE,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  analysis_scope = paste("user-authorized comparison of native UCell, custom",
    "per-cell balanced scoring, and custom gene-trajectory consensus;",
    "accepted read-only scWorkbench bundle and exact existing windows"),
  source_ucell_report = ucell_report, source_driver_report = driver_report,
  bundle = bundle_dir, bundle_fingerprint = fingerprint_before,
  profiles = 38L, cytokines = length(unique(inventory$cytokine)),
  old_score = "native UCell on RNA_counts; max_rank 1500",
  new_score = paste("equal-weight mean of per-gene zero-aware positive",
    "percentile ranks on SCT_data"),
  zero_policy = "zero expression maps to zero score contribution",
  positive_rank_policy = "average rank among positive cells divided by positive-cell count",
  gene_weight = "1 / signature_size", random_sampling = FALSE,
  trajectory_consensus = list(
    cell_score = FALSE,
    matrix = "SCT_data",
    per_window_gene_summary = "mean",
    per_gene_scaling = "min-max across eligible windows to [0,1]",
    across_gene_center = "median",
    across_gene_band = "IQR",
    single_gene_can_determine_center = FALSE),
  trajectory = list(width_fraction = 0.05, step_fraction = 0.01,
    min_cells = 20L, center = "median", band = "IQR"),
  driver_plot_scale = "per-gene min-max scaling of mean trajectory to [0,1]",
  driver_boundaries = "BB1/BB2 dashed lines at saved bb_up median pseudotimes",
  driver_selection_changed = FALSE,
  defined_pt_cells = 12642L, missing_pt_cells = 14081L,
  result_promotions = 0L,
  score_correlation = list(spearman_min = min(score_summary$spearman),
    spearman_median = median(score_summary$spearman),
    spearman_max = max(score_summary$spearman)),
  multipage_pdf = file.path("profile_figures", basename(pdf_path)),
  figures = figure_inventory$file,
  report_sha256 = wb_source_sha256(file.path(stage_dir, "REPORT.html"))),
  file.path(stage_dir, "completion.json"))

if (!file.rename(stage_dir, output_dir)) {
  stop("Failed to atomically finalize output directory")
}
log_event("Completed three-method comparison at", output_dir)
