#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 11_all_cytokine_neighborhood_coordination_report.R",
  "<bdev_bundle> <existing_38_profile_ucell_report>",
  "<existing_balanced_report> <existing_driver_report>",
  "<output_dir> <expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
ucell_report <- normalizePath(args[[2L]], mustWork = TRUE)
balanced_report <- normalizePath(args[[3L]], mustWork = TRUE)
driver_report <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
expected_fingerprint <- args[[6L]]
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(stage_dir, "profile_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(Matrix)
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
slug <- function(x) gsub("(^_+|_+$)", "", tolower(gsub("[^A-Za-z0-9]+", "_", x)))
as_flag <- function(x) !is.na(x) & tolower(as.character(x)) %in%
  c("true", "t", "1", "yes")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}

ucell_inputs <- c("completion.json", "passing_signature_inventory.tsv",
  "passing_gene_membership.tsv", "ucell_scores.rds",
  "sliding_window_summary.tsv", "window_support.tsv",
  "window_membership.rds", "bb_up_change_points.tsv")
balanced_inputs <- c("completion.json", "balanced_cell_scores.rds",
  "balanced_sliding_window_summary.tsv", "checks.json")
driver_inputs <- c("completion.json", "top5_driver_genes.tsv",
  "top5_driver_mean_expression_curves.tsv", "checks.json")
ucell_paths <- file.path(ucell_report, ucell_inputs)
balanced_paths <- file.path(balanced_report, balanced_inputs)
driver_paths <- file.path(driver_report, driver_inputs)
stopifnot(all(file.exists(ucell_paths)), all(file.exists(balanced_paths)),
  all(file.exists(driver_paths)))

ucell_completion <- read_json(ucell_paths[[1L]], simplifyVector = TRUE)
balanced_completion <- read_json(balanced_paths[[1L]], simplifyVector = TRUE)
driver_completion <- read_json(driver_paths[[1L]], simplifyVector = TRUE)
stopifnot(identical(ucell_completion$status, "complete"),
  identical(ucell_completion$promoted, FALSE),
  identical(as.integer(ucell_completion$passing_signatures), 38L),
  identical(balanced_completion$status, "complete"),
  identical(balanced_completion$promoted, FALSE),
  identical(as.integer(balanced_completion$profiles), 38L),
  identical(driver_completion$status, "complete"),
  identical(driver_completion$promoted, FALSE),
  identical(as.integer(driver_completion$profiles), 38L),
  identical(driver_completion$driver_summary, "mean"),
  identical(driver_completion$driver_plot_scale,
    "per-gene min-max scaling of mean trajectory to [0,1]"))

inventory <- read_tsv(file.path(ucell_report,
  "passing_signature_inventory.tsv"))
members <- read_tsv(file.path(ucell_report, "passing_gene_membership.tsv"))
old_family <- readRDS(file.path(ucell_report, "ucell_scores.rds"))
old_scores <- data.frame(key = as.character(old_family$key),
  cell_id = as.character(old_family$cell_id), value = as.numeric(old_family$value),
  stringsAsFactors = FALSE)
old_trajectory <- read_tsv(file.path(ucell_report,
  "sliding_window_summary.tsv"))
window_support <- read_tsv(file.path(ucell_report, "window_support.tsv"))
window_membership <- readRDS(file.path(ucell_report, "window_membership.rds"))
change_points <- read_tsv(file.path(ucell_report, "bb_up_change_points.tsv"))
change_points$boundary_label <- paste0("BB", seq_len(nrow(change_points)))
balanced_scores <- readRDS(file.path(balanced_report,
  "balanced_cell_scores.rds"))
balanced_trajectory <- read_tsv(file.path(balanced_report,
  "balanced_sliding_window_summary.tsv"))
driver_table <- read_tsv(file.path(driver_report, "top5_driver_genes.tsv"))
driver_curves <- read_tsv(file.path(driver_report,
  "top5_driver_mean_expression_curves.tsv"))

stopifnot(nrow(inventory) == 38L, all(as_flag(inventory$final_overall_pass)),
  all(inventory$final_size >= 20L), !anyDuplicated(inventory$signature_id),
  setequal(unique(members$signature_id), inventory$signature_id),
  !anyDuplicated(members[c("signature_id", "feature_id")]),
  setequal(unique(old_scores$key), inventory$signature_id),
  setequal(unique(balanced_scores$key), inventory$signature_id),
  setequal(unique(driver_table$feature_set), inventory$signature_id),
  "plot_value_internal_0_1" %in% names(driver_curves))
inventory <- inventory[order(inventory$cytokine,
  match(inventory$direction, c("up", "down"))), , drop = FALSE]
members <- members[order(match(members$signature_id,
  inventory$signature_id), members$final_rank), , drop = FALSE]

# Native knn_neighbors requires a durable, non-promoting session. The accepted
# bundle is opened session-capable but every state-changing outcome is forbidden
# and the source fingerprint is checked before and after.
bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint),
  identical(ucell_completion$bdev_bundle_fingerprint, fingerprint_before),
  identical(balanced_completion$bundle_fingerprint, fingerprint_before),
  identical(driver_completion$bundle_fingerprint, fingerprint_before))
cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "pt"))
stopifnot(identical(as.character(metadata$cell_id), cells),
  length(cells) == 26723L, sum(is.finite(metadata$pt)) == 12642L,
  sum(is.na(metadata$pt)) == 14081L)
for (id in inventory$signature_id) {
  stopifnot(identical(old_scores$cell_id[old_scores$key == id], cells),
    identical(balanced_scores$cell_id[balanced_scores$key == id], cells))
}

space_info <- wb_space_info(bundle, "pca")
stopifnot(nrow(space_info) == 1L, as.integer(space_info$dims) == 50L)
session_root <- file.path(dirname(bundle_dir), "sessions",
  "all_cytokine_neighborhood_coordination_v1")
session_id <- "all_cytokine_neighborhood_coordination_v1"
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) {
  wb_session_resume(bundle, workspace)
} else {
  wb_session_create(bundle, workspace_root = session_root,
    session_id = session_id)
}
if (!is.null(session$durable$promotion_receipt)) {
  stop("Neighborhood coordination session was unexpectedly promoted")
}
knn_ids <- names(session$runs)[vapply(session$runs, function(x)
  identical(x$op_name, "knn_neighbors") &&
    startsWith(as.character(x$notes %||% ""), "[rna_pca_knn50] "), logical(1))]
stopifnot(length(knn_ids) <= 1L)
if (!length(knn_ids)) {
  log_event("Calculating native scWorkbench k=50 neighbors in RNA PCA")
  knn_result <- wb_session_run_op(session, bundle, "knn_neighbors",
    inputs = list(space = "pca", cell_set = "all", query_cells = cells,
      candidate_cells = cells, k = 50L, include_query = TRUE),
    notes = paste0("[rna_pca_knn50] Native neighbors for neighborhood ",
      "coordination sensitivity analysis; all 50 PCA dimensions; no promotion."))
  session <- knn_result$session
  knn_run <- knn_result$run
} else {
  log_event("Reusing completed native RNA-PCA k=50 neighbors")
  knn_run <- session$runs[[knn_ids[[1L]]]]
}
neighbors <- knn_run$outputs$neighbors$table
neighbors <- neighbors[order(match(neighbors$query_cell_id, cells),
  neighbors$neighbor_rank), , drop = FALSE]
stopifnot(nrow(neighbors) == length(cells) * 50L,
  identical(unique(neighbors$query_cell_id), cells),
  all(table(neighbors$query_cell_id) == 50L),
  all(neighbors$neighbor_rank == rep(seq_len(50L), times = length(cells))))
neighbor_index <- matrix(match(neighbors$cell_id, cells), nrow = length(cells),
  ncol = 50L, byrow = TRUE)
stopifnot(!anyNA(neighbor_index),
  identical(neighbor_index[, 1L], seq_along(cells)),
  all(neighbors$distance[neighbors$neighbor_rank == 1L] == 0))
write_tsv_gz(neighbors, file.path(stage_dir, "rna_pca_knn50_neighbors.tsv.gz"))
write_tsv(space_info, file.path(stage_dir, "rna_pca_space_info.tsv"))

genes <- unique(members$feature_id)
stopifnot(all(genes %in% wb_feature_ids(bundle, "SCT_data")))
log_event("Reading SCT_data for", length(genes), "unique signature genes")
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", genes,
  cells, missing = "zero"))
stopifnot(identical(rownames(expression), genes),
  identical(colnames(expression), cells), all(is.finite(expression)),
  all(expression >= 0))

zero_aware_percentile <- function(x) {
  out <- numeric(length(x)); positive <- x > 0
  if (any(positive)) {
    out[positive] <- rank(x[positive], ties.method = "average") / sum(positive)
  }
  out
}
k_values <- c(15L, 30L, 50L)
thresholds <- c(0.60, 0.75, 0.90)
variant_inventory <- expand.grid(k = k_values, threshold = thresholds,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
variant_inventory$variant <- sprintf("k%02d_q%02d", variant_inventory$k,
  round(100 * variant_inventory$threshold))
variant_inventory$primary <- variant_inventory$k == 30L &
  variant_inventory$threshold == 0.75
stopifnot(sum(variant_inventory$primary) == 1L)
write_tsv(variant_inventory, file.path(stage_dir,
  "coordination_variant_inventory.tsv"))

n_profiles <- nrow(inventory)
n_variants <- nrow(variant_inventory)
coordination_matrix <- matrix(NA_real_, nrow = n_profiles * n_variants,
  ncol = length(cells))
rownames(coordination_matrix) <- as.vector(t(outer(variant_inventory$variant,
  inventory$signature_id, paste, sep = "::")))
colnames(coordination_matrix) <- cells
primary_intensity <- matrix(NA_real_, nrow = n_profiles, ncol = length(cells),
  dimnames = list(inventory$signature_id, cells))
primary_participation <- primary_intensity
gene_smoothing_rows <- vector("list", length(k_values))

for (ki in seq_along(k_values)) {
  k <- k_values[[ki]]
  log_event("Smoothing signature genes over native RNA-PCA neighborhoods, k =", k)
  adjacency <- sparseMatrix(
    i = rep(seq_along(cells), each = k),
    j = as.vector(t(neighbor_index[, seq_len(k), drop = FALSE])),
    x = rep(1 / k, length(cells) * k),
    dims = c(length(cells), length(cells)))
  smoothed_expression <- as.matrix(expression %*% t(adjacency))
  dimnames(smoothed_expression) <- dimnames(expression)
  gene_percentiles <- t(apply(smoothed_expression, 1L,
    zero_aware_percentile))
  dimnames(gene_percentiles) <- dimnames(expression)
  gene_smoothing_rows[[ki]] <- data.frame(k = k, feature_id = genes,
    positive_cells = rowSums(smoothed_expression > 0),
    positive_fraction = rowMeans(smoothed_expression > 0),
    smoothed_mean = rowMeans(smoothed_expression),
    percentile_mean = rowMeans(gene_percentiles), stringsAsFactors = FALSE)
  stopifnot(all(gene_percentiles >= 0 & gene_percentiles <= 1))
  for (pi in seq_len(n_profiles)) {
    id <- inventory$signature_id[[pi]]
    set_genes <- members$feature_id[members$signature_id == id]
    intensity <- colMeans(gene_percentiles[set_genes, , drop = FALSE])
    for (ti in seq_along(thresholds)) {
      threshold <- thresholds[[ti]]
      participation <- colMeans(
        gene_percentiles[set_genes, , drop = FALSE] >= threshold)
      variant_row <- (which(variant_inventory$k == k &
        variant_inventory$threshold == threshold) - 1L) * n_profiles + pi
      coordination_matrix[variant_row, ] <- intensity * participation
      if (k == 30L && threshold == 0.75) {
        primary_intensity[pi, ] <- intensity
        primary_participation[pi, ] <- participation
      }
    }
  }
  rm(adjacency, smoothed_expression, gene_percentiles)
  invisible(gc())
}
gene_smoothing_qc <- do.call(rbind, gene_smoothing_rows)
write_tsv(gene_smoothing_qc, file.path(stage_dir,
  "gene_neighborhood_smoothing_qc.tsv"))
stopifnot(all(is.finite(coordination_matrix)),
  all(coordination_matrix >= 0 & coordination_matrix <= 1),
  all(is.finite(primary_intensity)), all(primary_intensity >= 0 &
    primary_intensity <= 1),
  all(is.finite(primary_participation)), all(primary_participation >= 0 &
    primary_participation <= 1))
saveRDS(coordination_matrix, file.path(stage_dir,
  "coordination_sensitivity_cell_scores.rds"), compress = "gzip")

primary_variant <- variant_inventory$variant[variant_inventory$primary]
primary_rows <- match(paste(primary_variant, inventory$signature_id, sep = "::"),
  rownames(coordination_matrix))
primary_matrix <- coordination_matrix[primary_rows, , drop = FALSE]
rownames(primary_matrix) <- inventory$signature_id
primary_scores <- do.call(rbind, lapply(seq_len(n_profiles), function(i) {
  data.frame(key = inventory$signature_id[[i]], cell_id = cells,
    value = as.numeric(primary_matrix[i, ]),
    intensity = as.numeric(primary_intensity[i, ]),
    participation = as.numeric(primary_participation[i, ]),
    stringsAsFactors = FALSE)
}))
rownames(primary_scores) <- NULL
saveRDS(primary_scores, file.path(stage_dir,
  "neighborhood_coordination_cell_scores.rds"), compress = "gzip")
write_tsv_gz(primary_scores, file.path(stage_dir,
  "neighborhood_coordination_cell_scores.tsv.gz"))

pt <- data.frame(cell_id = cells[is.finite(metadata$pt)],
  value = as.numeric(metadata$pt[is.finite(metadata$pt)]),
  stringsAsFactors = FALSE)
log_event("Calculating native rolling summaries for primary coordination scores")
coordination_window <- sliding_window_signature_summary(
  primary_scores[, c("key", "cell_id", "value")], pt,
  width_fractions = 0.05, step_fraction = 0.01, min_cells = 20L)
coordination_trajectory <- coordination_window$summaries
coordination_trajectory$window_midpoint <-
  (coordination_trajectory$window_left +
    coordination_trajectory$window_right) / 2
stopifnot(identical(coordination_window$membership, window_membership))
write_tsv(coordination_trajectory, file.path(stage_dir,
  "neighborhood_coordination_sliding_window_summary.tsv"))

# Sensitivity trajectories reuse the exact saved memberships. They are kept
# diagnostic and are not substituted for the native primary summary above.
window_index <- lapply(window_membership, match, table = cells)
sensitivity_rows <- vector("list", nrow(coordination_matrix))
for (ri in seq_len(nrow(coordination_matrix))) {
  pieces <- strsplit(rownames(coordination_matrix)[[ri]], "::", fixed = TRUE)[[1L]]
  values <- coordination_matrix[ri, ]
  sensitivity_rows[[ri]] <- data.frame(variant = pieces[[1L]],
    feature_set = pieces[[2L]], window_id = window_support$window_id,
    window_midpoint = window_support$window_midpoint,
    eligible_window = window_support$eligible_window,
    score_median = vapply(window_index, function(ix) median(values[ix]), numeric(1)),
    score_q25 = vapply(window_index, function(ix)
      quantile(values[ix], 0.25, names = FALSE, type = 7), numeric(1)),
    score_q75 = vapply(window_index, function(ix)
      quantile(values[ix], 0.75, names = FALSE, type = 7), numeric(1)),
    stringsAsFactors = FALSE)
}
sensitivity_trajectory <- do.call(rbind, sensitivity_rows)
write_tsv(sensitivity_trajectory, file.path(stage_dir,
  "coordination_sensitivity_trajectories.tsv"))

score_summary <- do.call(rbind, lapply(seq_len(n_profiles), function(i) {
  id <- inventory$signature_id[[i]]
  old <- old_scores$value[old_scores$key == id]
  balanced <- balanced_scores$value[balanced_scores$key == id]
  coordination <- primary_matrix[i, ]
  data.frame(feature_set = id,
    spearman_coordination_vs_ucell = cor(coordination, old,
      method = "spearman"),
    spearman_coordination_vs_balanced = cor(coordination, balanced,
      method = "spearman"),
    coordination_mean = mean(coordination),
    coordination_sd = sd(coordination),
    intensity_mean = mean(primary_intensity[i, ]),
    participation_mean = mean(primary_participation[i, ]),
    zero_coordination_fraction = mean(coordination == 0),
    stringsAsFactors = FALSE)
}))
write_tsv(score_summary, file.path(stage_dir,
  "per_cell_method_comparison.tsv"))

primary_by_profile <- split(seq_len(nrow(coordination_matrix)),
  sub("^.*::", "", rownames(coordination_matrix)))
sensitivity_summary <- do.call(rbind, lapply(seq_len(nrow(coordination_matrix)),
  function(ri) {
    pieces <- strsplit(rownames(coordination_matrix)[[ri]], "::", fixed = TRUE)[[1L]]
    id <- pieces[[2L]]
    primary <- primary_matrix[match(id, rownames(primary_matrix)), ]
    local <- sensitivity_trajectory[
      sensitivity_trajectory$variant == pieces[[1L]] &
        sensitivity_trajectory$feature_set == id &
        sensitivity_trajectory$eligible_window, ]
    data.frame(variant = pieces[[1L]], feature_set = id,
      spearman_vs_primary_cells = cor(coordination_matrix[ri, ], primary,
        method = "spearman"),
      peak_pt = local$window_midpoint[which.max(local$score_median)],
      max_rolling_median = max(local$score_median), stringsAsFactors = FALSE)
  }))
sensitivity_summary <- merge(sensitivity_summary, variant_inventory,
  by = "variant", all.x = TRUE, sort = FALSE)
write_tsv(sensitivity_summary, file.path(stage_dir,
  "coordination_sensitivity_summary.tsv"))

space <- wb_space(bundle, "umap.atac", cells = cells)
coordinates <- if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells))
  data.frame(cell_id = cells, x = as.numeric(space$x), y = as.numeric(space$y))
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  data.frame(cell_id = cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]))
}

map_plot <- function(values, title, subtitle) {
  positive <- values > 0
  cap <- quantile(values[positive], 0.99, names = FALSE, type = 7)
  display <- rep(NA_real_, length(values))
  display[positive] <- pmin(values[positive], cap) / cap
  data <- coordinates; data$display <- display
  data <- data[order(ifelse(is.na(data$display), -Inf, data$display)), ]
  ggplot(data, aes(x, y)) +
    geom_point(data = data[is.na(data$display), ], color = "#CDD2DA",
      size = 0.1, alpha = 0.7, stroke = 0) +
    geom_point(data = data[!is.na(data$display), ], aes(color = display),
      size = 0.24, alpha = 1, stroke = 0) +
    scale_color_gradientn(colours = c("#0066FF", "#00D9FF", "#00F060",
      "#FFF000", "#FF5A00", "#FF00A8"), limits = c(0, 1),
      oob = scales::squish, name = "Within-panel\nscore") +
    coord_equal() + theme_void(base_size = 9) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
    labs(title = title, subtitle = subtitle,
      caption = "Positive values scaled to within-panel p99; zero is gray.")
}

unsupported <- window_support[!window_support$eligible_window, , drop = FALSE]
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

gene_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
driver_plot <- function(id) {
  data <- driver_curves[driver_curves$feature_set == id, , drop = FALSE]
  ranks <- driver_table[driver_table$feature_set == id,
    c("feature_id", "driver_rank"), drop = FALSE]
  ranks <- ranks[order(ranks$driver_rank), ]; genes <- ranks$feature_id
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
      aes(x = pt_median, y = Inf, label = boundary_label), inherit.aes = FALSE,
      vjust = 1.35, fontface = "bold", size = 3, color = "#374151") +
    geom_line(linewidth = 0.9) +
    geom_text_repel(data = ends, aes(label = feature_id), direction = "y",
      hjust = 0, nudge_x = 0.02 * span, seed = 1729,
      box.padding = 0.2, point.padding = 0.15, min.segment.length = 0,
      max.overlaps = Inf, show.legend = FALSE) +
    scale_color_manual(values = colors, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.20))) +
    coord_cartesian(clip = "off") +
    labs(title = "Existing top-five member-gene drivers",
      subtitle = "Mean-expression dynamic-score rank; each mean curve scaled internally to 0-1",
      x = "Existing Bdev pseudotime (pt)",
      y = "Within-gene scaled mean expression",
      caption = "Dashed lines mark BB1/BB2 medians; driver selection is unchanged.") +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(),
      legend.position = "none", plot.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 75, 5.5, 5.5))
}

log_event("Rendering 38 three-per-cell-method comparison figures")
pdf_path <- file.path(figure_dir,
  "all_38_profiles_ucell_balanced_coordination_with_drivers.pdf")
grDevices::pdf(pdf_path, width = 18, height = 13.5, onefile = TRUE)
figure_rows <- vector("list", n_profiles)
for (i in seq_len(n_profiles)) {
  row <- inventory[i, ]; id <- row$signature_id
  old_values <- old_scores$value[old_scores$key == id]
  balanced_values <- balanced_scores$value[balanced_scores$key == id]
  coordination_values <- primary_matrix[i, ]
  old_curve <- old_trajectory[old_trajectory$feature_set == id, ]
  balanced_curve <- balanced_trajectory[balanced_trajectory$feature_set == id, ]
  coordination_curve <- coordination_trajectory[
    coordination_trajectory$feature_set == id, ]
  label <- paste(row$cytokine, toupper(row$direction), "-", row$final_size,
    "genes")
  direction_color <- if (row$direction == "up") "#D73027" else "#2166AC"
  top <- map_plot(old_values, "Method 1: UCell UMAP", label) |
    map_plot(balanced_values, "Method 2: balanced-score UMAP", label) |
    map_plot(coordination_values, "Method 3: coordination-score UMAP", label)
  middle <- trajectory_plot(old_curve, "UCell rolling trajectory",
      direction_color, "UCell score") |
    trajectory_plot(balanced_curve, "Balanced rolling trajectory",
      "#7B3294", "Balanced score") |
    trajectory_plot(coordination_curve, "Neighborhood coordination trajectory",
      "#00796B", "Intensity x participation")
  combined <- top / middle / driver_plot(id) +
    plot_layout(heights = c(1, 0.9, 0.8)) +
    plot_annotation(title = paste("Three per-cell methods:", label),
      subtitle = paste("Coordination uses native RNA-PCA k=30 neighborhoods;",
        "75th-percentile activity; intensity x participating-gene fraction"),
      theme = theme(plot.title = element_text(size = 17, face = "bold")))
  print(combined)
  file <- paste0(sprintf("%02d_", i), slug(paste(row$cytokine,
    row$direction)), "_ucell_balanced_coordination.png")
  ggsave(file.path(figure_dir, file), combined, width = 18, height = 13.5,
    dpi = 220, units = "in", bg = "white")
  comparison <- score_summary[score_summary$feature_set == id, ]
  figure_rows[[i]] <- data.frame(signature_id = id, cytokine = row$cytokine,
    direction = row$direction, signature_size = as.integer(row$final_size),
    spearman_coordination_vs_ucell = comparison$spearman_coordination_vs_ucell,
    spearman_coordination_vs_balanced =
      comparison$spearman_coordination_vs_balanced,
    file = file, stringsAsFactors = FALSE)
}
grDevices::dev.off()
figure_inventory <- do.call(rbind, figure_rows)
write_tsv(figure_inventory, file.path(stage_dir, "figure_inventory.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
checks <- list(source_has_38_passing_profiles = n_profiles == 38L,
  native_rna_pca_neighbors = identical(knn_run$op_name, "knn_neighbors"),
  fifty_neighbors_per_cell = all(table(neighbors$query_cell_id) == 50L),
  self_is_rank_one = identical(neighbor_index[, 1L], seq_along(cells)),
  all_signature_genes_used = setequal(genes, unique(members$feature_id)),
  coordination_scores_finite_and_bounded = all(is.finite(coordination_matrix)) &&
    all(coordination_matrix >= 0 & coordination_matrix <= 1),
  primary_intensity_finite_and_bounded = all(is.finite(primary_intensity)) &&
    all(primary_intensity >= 0 & primary_intensity <= 1),
  primary_participation_finite_and_bounded =
    all(is.finite(primary_participation)) &&
    all(primary_participation >= 0 & primary_participation <= 1),
  exact_window_membership_reused = identical(coordination_window$membership,
    window_membership),
  exact_source_pt_preserved = identical(pt$value,
    as.numeric(metadata$pt[is.finite(metadata$pt)])),
  source_missing_pt_preserved = sum(is.na(metadata$pt)) == 14081L,
  nine_sensitivity_variants = nrow(variant_inventory) == 9L,
  every_figure_exists = all(file.exists(file.path(figure_dir,
    figure_inventory$file))), multipage_pdf_exists = file.exists(pdf_path),
  bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after), result_promotions = 0L)
stopifnot(all(unlist(checks[names(checks) != "result_promotions"])),
  checks$result_promotions == 0L,
  is.null(session$durable$promotion_receipt))
write_json(list(status = "passed", checks = checks),
  file.path(stage_dir, "checks.json"))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
provenance_inputs <- c(ucell_paths, balanced_paths, driver_paths,
  file.path(bundle_dir, "manifest.json"), script_path,
  "figure_scripts/bdev_reference/run_11_all_cytokine_neighborhood_coordination_report.sh",
  "src/sliding_window_signature_summary.R", "config/bdev_reference.yml")
input_manifest <- data.frame(path = normalizePath(provenance_inputs,
  mustWork = TRUE), sha256 = vapply(provenance_inputs, wb_source_sha256,
  character(1)), stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

html_rows <- vapply(seq_len(nrow(figure_inventory)), function(i) {
  x <- figure_inventory[i, ]
  paste0("<h2>", x$cytokine, " ", toupper(x$direction), "</h2>",
    "<p>Coordination correlation: UCell ",
    sprintf("%.3f", x$spearman_coordination_vs_ucell),
    "; balanced ", sprintf("%.3f", x$spearman_coordination_vs_balanced),
    "</p>", '<img src="profile_figures/', x$file, '" alt="',
    x$cytokine, " ", x$direction, ' three per-cell methods">')
}, character(1))
html <- c('<!doctype html><html><head><meta charset="utf-8">',
  '<title>Neighborhood coordination cytokine scores</title>',
  '<style>body{font:16px system-ui,sans-serif;max-width:1600px;margin:32px auto;padding:0 22px;color:#172033}img{max-width:100%;height:auto;border:1px solid #ddd}p,li{line-height:1.5}code{overflow-wrap:anywhere}a{color:#1261a0}</style></head><body>',
  '<h1>UCell, balanced, and neighborhood-coordination scores</h1>',
  '<p>The primary coordination score uses native scWorkbench neighbors in the existing 50-dimensional RNA-PCA space. For each cell and gene, expression is averaged over the cell plus its 29 nearest neighbors, converted to a zero-aware within-gene percentile, and summarized as <code>mean percentile x fraction of genes at or above percentile 0.75</code>. A single active gene therefore receives a strong participation penalty.</p>',
  '<p>Sensitivity results cover k=15, 30, and 50 with activity thresholds 0.60, 0.75, and 0.90. These are exploratory scores without matched-random-set null calibration. Driver identities are unchanged and BB1/BB2 are labeled in the driver panel.</p>',
  '<p><a href="profile_figures/all_38_profiles_ucell_balanced_coordination_with_drivers.pdf">Download 38-page PDF</a> · <a href="per_cell_method_comparison.tsv">Method comparison</a> · <a href="coordination_sensitivity_summary.tsv">Sensitivity summary</a> · <a href="neighborhood_coordination_cell_scores.tsv.gz">Primary cell scores</a></p>',
  html_rows, '</body></html>')
writeLines(html, file.path(stage_dir, "REPORT.html"))
writeLines(c("# Neighborhood coordination score comparison", "",
  "[Open report](REPORT.html)", "",
  "All 38 profiles compare native UCell, equal-gene balanced, and neighborhood-smoothed coordination scores."),
  file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", promoted = FALSE,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  analysis_scope = paste("user-authorized neighborhood-smoothed coordination",
    "scoring; native scWorkbench kNN and rolling summaries; no promotion"),
  source_ucell_report = ucell_report,
  source_balanced_report = balanced_report,
  source_driver_report = driver_report, bundle = bundle_dir,
  bundle_fingerprint = fingerprint_before, profiles = 38L,
  cytokines = length(unique(inventory$cytokine)),
  neighborhood = list(space = "pca", dimensions = 50L,
    native_operation = "knn_neighbors", primary_k = 30L,
    include_query = TRUE, candidate_cells = "all Bdev cells"),
  score = list(gene_input = "SCT_data neighborhood mean",
    gene_transform = "zero-aware percentile across cells",
    intensity = "equal-gene mean percentile",
    participation = "fraction of genes with percentile >= 0.75",
    combined = "intensity * participation"),
  sensitivity = list(k = k_values, percentile_threshold = thresholds),
  random_matched_null_calibration = FALSE,
  trajectory = list(width_fraction = 0.05, step_fraction = 0.01,
    min_cells = 20L, center = "median", band = "IQR"),
  driver_selection_changed = FALSE,
  driver_boundaries = "BB1/BB2 dashed lines at saved median pseudotimes",
  defined_pt_cells = 12642L, missing_pt_cells = 14081L,
  result_promotions = 0L, native_session_workspace = workspace,
  multipage_pdf = file.path("profile_figures", basename(pdf_path)),
  figures = figure_inventory$file,
  report_sha256 = wb_source_sha256(file.path(stage_dir, "REPORT.html"))),
  file.path(stage_dir, "completion.json"))

if (!file.rename(stage_dir, output_dir)) {
  stop("Failed to atomically finalize output directory")
}
log_event("Completed neighborhood coordination comparison at", output_dir)
