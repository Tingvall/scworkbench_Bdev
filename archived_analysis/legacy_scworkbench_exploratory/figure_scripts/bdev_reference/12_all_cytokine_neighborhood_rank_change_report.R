#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 12_all_cytokine_neighborhood_rank_change_report.R",
  "<bdev_bundle> <existing_38_profile_ucell_report>",
  "<existing_driver_report> <existing_neighborhood_report>",
  "<output_dir> <expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
ucell_report <- normalizePath(args[[2L]], mustWork = TRUE)
driver_report <- normalizePath(args[[3L]], mustWork = TRUE)
neighborhood_report <- normalizePath(args[[4L]], mustWork = TRUE)
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
  "passing_gene_membership.tsv", "window_support.tsv",
  "window_membership.rds", "bb_up_change_points.tsv")
driver_inputs <- c("completion.json", "top5_driver_genes.tsv", "checks.json")
neighborhood_inputs <- c("completion.json", "rna_pca_knn50_neighbors.tsv.gz",
  "rna_pca_space_info.tsv", "checks.json")
ucell_paths <- file.path(ucell_report, ucell_inputs)
driver_paths <- file.path(driver_report, driver_inputs)
neighborhood_paths <- file.path(neighborhood_report, neighborhood_inputs)
stopifnot(all(file.exists(ucell_paths)), all(file.exists(driver_paths)),
  all(file.exists(neighborhood_paths)))

ucell_completion <- read_json(ucell_paths[[1L]], simplifyVector = TRUE)
driver_completion <- read_json(driver_paths[[1L]], simplifyVector = TRUE)
neighborhood_completion <- read_json(neighborhood_paths[[1L]],
  simplifyVector = TRUE)
stopifnot(identical(ucell_completion$status, "complete"),
  identical(ucell_completion$promoted, FALSE),
  identical(as.integer(ucell_completion$passing_signatures), 38L),
  identical(driver_completion$status, "complete"),
  identical(driver_completion$promoted, FALSE),
  identical(as.integer(driver_completion$profiles), 38L),
  identical(neighborhood_completion$status, "complete"),
  identical(neighborhood_completion$promoted, FALSE),
  identical(as.integer(neighborhood_completion$profiles), 38L))

inventory <- read_tsv(file.path(ucell_report,
  "passing_signature_inventory.tsv"))
members <- read_tsv(file.path(ucell_report, "passing_gene_membership.tsv"))
window_support <- read_tsv(file.path(ucell_report, "window_support.tsv"))
window_membership <- readRDS(file.path(ucell_report, "window_membership.rds"))
change_points <- read_tsv(file.path(ucell_report, "bb_up_change_points.tsv"))
change_points$boundary_label <- paste0("BB", seq_len(nrow(change_points)))
driver_table <- read_tsv(file.path(driver_report, "top5_driver_genes.tsv"))

stopifnot(nrow(inventory) == 38L, all(as_flag(inventory$final_overall_pass)),
  all(inventory$final_size >= 20L), !anyDuplicated(inventory$signature_id),
  setequal(unique(members$signature_id), inventory$signature_id),
  !anyDuplicated(members[c("signature_id", "feature_id")]),
  setequal(unique(driver_table$feature_set), inventory$signature_id),
  nrow(driver_table) == 5L * nrow(inventory),
  nrow(window_support) == length(window_membership),
  identical(as.character(window_support$window_id), names(window_membership)))
inventory <- inventory[order(inventory$cytokine,
  match(inventory$direction, c("up", "down"))), , drop = FALSE]
members <- members[order(match(members$signature_id,
  inventory$signature_id), members$final_rank), , drop = FALSE]

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint),
  identical(ucell_completion$bdev_bundle_fingerprint, fingerprint_before),
  identical(driver_completion$bundle_fingerprint, fingerprint_before),
  identical(neighborhood_completion$bundle_fingerprint, fingerprint_before))
cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "pt"))
stopifnot(identical(as.character(metadata$cell_id), cells),
  length(cells) == 26723L, sum(is.finite(metadata$pt)) == 12642L,
  sum(is.na(metadata$pt)) == 14081L)

log_event("Reading and validating saved native RNA-PCA neighbors")
neighbors <- read_tsv(file.path(neighborhood_report,
  "rna_pca_knn50_neighbors.tsv.gz"))
neighbors <- neighbors[order(match(neighbors$query_cell_id, cells),
  neighbors$neighbor_rank), , drop = FALSE]
stopifnot(nrow(neighbors) == length(cells) * 50L,
  identical(unique(neighbors$query_cell_id), cells),
  all(table(neighbors$query_cell_id) == 50L),
  all(neighbors$neighbor_rank == rep(seq_len(50L), times = length(cells))))
neighbor_index <- matrix(match(neighbors$cell_id, cells), nrow = length(cells),
  ncol = 50L, byrow = TRUE)
stopifnot(!anyNA(neighbor_index),
  identical(neighbor_index[, 1L], seq_along(cells)))

genes <- unique(members$feature_id)
stopifnot(all(genes %in% wb_feature_ids(bundle, "SCT_data")))
log_event("Reading SCT_data for", length(genes), "unique signature genes")
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", genes,
  cells, missing = "zero"))
stopifnot(identical(rownames(expression), genes),
  identical(colnames(expression), cells), all(is.finite(expression)),
  all(expression >= 0))

log_event("Calculating k=30 neighborhood means and within-gene percentiles")
k <- 30L
adjacency <- sparseMatrix(
  i = rep(seq_along(cells), each = k),
  j = as.vector(t(neighbor_index[, seq_len(k), drop = FALSE])),
  x = rep(1 / k, length(cells) * k),
  dims = c(length(cells), length(cells)))
smoothed_expression <- as.matrix(expression %*% t(adjacency))
dimnames(smoothed_expression) <- dimnames(expression)
zero_aware_percentile <- function(x) {
  out <- numeric(length(x)); positive <- x > 0
  if (any(positive)) {
    out[positive] <- rank(x[positive], ties.method = "average") / sum(positive)
  }
  out
}
gene_percentiles <- t(apply(smoothed_expression, 1L, zero_aware_percentile))
dimnames(gene_percentiles) <- dimnames(expression)
stopifnot(all(is.finite(gene_percentiles)),
  all(gene_percentiles >= 0 & gene_percentiles <= 1))
gene_qc <- data.frame(feature_id = genes,
  raw_positive_fraction = rowMeans(expression > 0),
  smoothed_positive_fraction = rowMeans(smoothed_expression > 0),
  smoothed_mean = rowMeans(smoothed_expression),
  percentile_mean = rowMeans(gene_percentiles), stringsAsFactors = FALSE)
write_tsv(gene_qc, file.path(stage_dir, "gene_neighborhood_rank_qc.tsv"))
rm(adjacency, expression, smoothed_expression)
invisible(gc())

log_event("Calculating each gene's median rank in the existing 5% windows")
window_index <- lapply(window_membership, match, table = cells)
stopifnot(all(vapply(window_index, function(x) !anyNA(x), logical(1))))
n_windows <- length(window_index)
gene_window_rank <- vapply(window_index, function(ix)
  apply(gene_percentiles[, ix, drop = FALSE], 1L, median),
  numeric(length(genes)))
rownames(gene_window_rank) <- genes
colnames(gene_window_rank) <- window_support$window_id
lag_windows <- 5L
gene_window_change <- matrix(NA_real_, nrow = length(genes),
  ncol = n_windows, dimnames = dimnames(gene_window_rank))
for (wi in seq.int(lag_windows + 1L, n_windows)) {
  if (as_flag(window_support$eligible_window[[wi]]) &&
      as_flag(window_support$eligible_window[[wi - lag_windows]])) {
    gene_window_change[, wi] <- gene_window_rank[, wi] -
      gene_window_rank[, wi - lag_windows]
  }
}
stopifnot(all(gene_window_change[is.finite(gene_window_change)] >= -1),
  all(gene_window_change[is.finite(gene_window_change)] <= 1))

gene_window_table <- do.call(rbind, lapply(seq_along(genes), function(gi)
  data.frame(feature_id = genes[[gi]], window_id = window_support$window_id,
    window_index = window_support$window_index,
    window_midpoint = window_support$window_midpoint,
    median_neighborhood_rank = as.numeric(gene_window_rank[gi, ]),
    rank_change_vs_five_windows_earlier =
      as.numeric(gene_window_change[gi, ]), stringsAsFactors = FALSE)))
write_tsv(gene_window_table, file.path(stage_dir,
  "gene_window_neighborhood_rank_change.tsv"))

log_event("Aggregating gene-level changes within each signature")
summary_rows <- vector("list", nrow(inventory))
driver_rows <- vector("list", nrow(inventory))
for (pi in seq_len(nrow(inventory))) {
  id <- inventory$signature_id[[pi]]
  set_genes <- members$feature_id[members$signature_id == id]
  delta <- gene_window_change[set_genes, , drop = FALSE]
  valid <- colSums(is.finite(delta)) == length(set_genes)
  score_mean <- score_median <- score_q25 <- score_q75 <-
    positive_fraction <- negative_fraction <- strong_positive_fraction <-
    strong_negative_fraction <- rep(NA_real_, n_windows)
  score_mean[valid] <- colMeans(delta[, valid, drop = FALSE])
  score_median[valid] <- apply(delta[, valid, drop = FALSE], 2L, median)
  score_q25[valid] <- apply(delta[, valid, drop = FALSE], 2L, quantile,
    probs = 0.25, names = FALSE, type = 7)
  score_q75[valid] <- apply(delta[, valid, drop = FALSE], 2L, quantile,
    probs = 0.75, names = FALSE, type = 7)
  positive_fraction[valid] <- colMeans(delta[, valid, drop = FALSE] > 0)
  negative_fraction[valid] <- colMeans(delta[, valid, drop = FALSE] < 0)
  strong_positive_fraction[valid] <- colMeans(
    delta[, valid, drop = FALSE] >= 0.05)
  strong_negative_fraction[valid] <- colMeans(
    delta[, valid, drop = FALSE] <= -0.05)
  directional_support_fraction <- ifelse(score_mean > 0,
    strong_positive_fraction, ifelse(score_mean < 0,
      strong_negative_fraction, 0))
  summary_rows[[pi]] <- data.frame(feature_set = id,
    window_id = window_support$window_id,
    window_index = window_support$window_index,
    window_midpoint = window_support$window_midpoint,
    eligible_window = as_flag(window_support$eligible_window),
    score_mean = score_mean, score_median_diagnostic = score_median,
    score_q25 = score_q25,
    score_q75 = score_q75, positive_gene_fraction = positive_fraction,
    negative_gene_fraction = negative_fraction,
    strong_positive_gene_fraction = strong_positive_fraction,
    strong_negative_gene_fraction = strong_negative_fraction,
    directional_support_fraction = directional_support_fraction,
    coordinated_support = directional_support_fraction >= 0.20,
    n_signature_genes = length(set_genes), stringsAsFactors = FALSE)
  drivers <- driver_table[driver_table$feature_set == id, , drop = FALSE]
  drivers <- drivers[order(drivers$driver_rank), , drop = FALSE]
  driver_rows[[pi]] <- do.call(rbind, lapply(seq_len(nrow(drivers)), function(j)
    data.frame(feature_set = id, feature_id = drivers$feature_id[[j]],
      driver_rank = as.integer(drivers$driver_rank[[j]]),
      window_id = window_support$window_id,
      window_index = window_support$window_index,
      window_midpoint = window_support$window_midpoint,
      rank_change_vs_five_windows_earlier = as.numeric(
        gene_window_change[drivers$feature_id[[j]], ]),
      stringsAsFactors = FALSE)))
}
signature_trajectory <- do.call(rbind, summary_rows)
driver_change <- do.call(rbind, driver_rows)
write_tsv(signature_trajectory, file.path(stage_dir,
  "signature_neighborhood_rank_change_trajectory.tsv"))
write_tsv(driver_change, file.path(stage_dir,
  "top5_driver_neighborhood_rank_change.tsv"))

# A cell-specific UMAP value compares that cell's current gene percentiles with
# the per-gene median ranks in the window five steps before its nearest current
# window. The subtraction is performed gene by gene before taking the signature
# median, preserving the intended coordinated-change definition.
valid_window <- is.finite(signature_trajectory$score_mean[
  signature_trajectory$feature_set == inventory$signature_id[[1L]]])
valid_window_indices <- which(valid_window)
valid_midpoints <- window_support$window_midpoint[valid_window_indices]
pt_defined <- which(is.finite(metadata$pt))
nearest_current_window <- rep(NA_integer_, length(cells))
inside <- metadata$pt[pt_defined] >= min(valid_midpoints) &
  metadata$pt[pt_defined] <= max(valid_midpoints)
nearest_local <- vapply(metadata$pt[pt_defined][inside], function(x)
  which.min(abs(valid_midpoints - x)), integer(1))
nearest_current_window[pt_defined[inside]] <-
  valid_window_indices[nearest_local]
scoreable_cells <- which(is.finite(nearest_current_window))
stopifnot(length(scoreable_cells) > 0L,
  all(nearest_current_window[scoreable_cells] > lag_windows))

log_event("Calculating per-cell gene-level change scores for UMAP display")
cell_score_matrix <- matrix(NA_real_, nrow = nrow(inventory),
  ncol = length(cells), dimnames = list(inventory$signature_id, cells))
cell_support_matrix <- cell_score_matrix
for (pi in seq_len(nrow(inventory))) {
  id <- inventory$signature_id[[pi]]
  set_genes <- members$feature_id[members$signature_id == id]
  values <- rep(NA_real_, length(cells))
  for (wi in valid_window_indices) {
    cell_ix <- which(nearest_current_window == wi)
    if (!length(cell_ix)) next
    baseline <- gene_window_rank[set_genes, wi - lag_windows]
    difference <- sweep(gene_percentiles[set_genes, cell_ix, drop = FALSE],
      1L, baseline, "-")
    local_mean <- colMeans(difference)
    positive_support <- colMeans(difference >= 0.05)
    negative_support <- colMeans(difference <= -0.05)
    values[cell_ix] <- local_mean
    cell_support_matrix[pi, cell_ix] <- ifelse(local_mean > 0,
      positive_support, ifelse(local_mean < 0, negative_support, 0))
  }
  cell_score_matrix[pi, ] <- values
}
stopifnot(all(is.finite(cell_score_matrix[, scoreable_cells, drop = FALSE])),
  all(cell_score_matrix[, scoreable_cells, drop = FALSE] >= -1),
  all(cell_score_matrix[, scoreable_cells, drop = FALSE] <= 1),
  all(is.finite(cell_support_matrix[, scoreable_cells, drop = FALSE])),
  all(cell_support_matrix[, scoreable_cells, drop = FALSE] >= 0 &
    cell_support_matrix[, scoreable_cells, drop = FALSE] <= 1),
  all(is.na(cell_score_matrix[, is.na(metadata$pt), drop = FALSE])))
cell_scores <- do.call(rbind, lapply(seq_len(nrow(inventory)), function(i)
  data.frame(key = inventory$signature_id[[i]], cell_id = cells,
    value = as.numeric(cell_score_matrix[i, ]),
    directional_support_fraction = as.numeric(cell_support_matrix[i, ]),
    coordinated_support = as.numeric(cell_support_matrix[i, ]) >= 0.20,
    pt = as.numeric(metadata$pt),
    nearest_window_index = nearest_current_window,
    reference_window_index = nearest_current_window - lag_windows,
    stringsAsFactors = FALSE)))
saveRDS(cell_scores, file.path(stage_dir,
  "neighborhood_rank_change_cell_scores.rds"), compress = "gzip")
write_tsv_gz(cell_scores, file.path(stage_dir,
  "neighborhood_rank_change_cell_scores.tsv.gz"))

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
  finite <- is.finite(values)
  cap <- quantile(abs(values[finite]), 0.99, names = FALSE, type = 7)
  if (!is.finite(cap) || cap <= 0) cap <- 0.01
  data <- coordinates; data$value <- values
  data <- data[order(ifelse(is.finite(data$value), abs(data$value), -Inf)), ]
  ggplot(data, aes(x, y)) +
    geom_point(color = "#B8BEC8", size = 0.12, alpha = 0.65, stroke = 0) +
    geom_point(data = data[is.finite(data$value), ], aes(color = value),
      size = 0.28, alpha = 1, stroke = 0) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-cap, cap), oob = scales::squish,
      name = "Rank\nchange") +
    coord_equal() + theme_void(base_size = 10) +
    theme(panel.background = element_rect(fill = "#E5E7EB", colour = NA),
      legend.position = "bottom", plot.title = element_text(face = "bold")) +
    labs(title = title, subtitle = subtitle,
      caption = "Blue: decrease; white: zero; red: increase. Gray: unavailable pt/lag.")
}

trajectory_plot <- function(data) {
  data <- data[is.finite(data$score_mean), , drop = FALSE]
  cap <- quantile(abs(c(data$score_q25, data$score_q75)), 0.99,
    names = FALSE, type = 7)
  if (!is.finite(cap) || cap <= 0) cap <- 0.01
  ggplot(data, aes(window_midpoint, score_mean)) +
    geom_hline(yintercept = 0, linewidth = 0.45, color = "#6B7280") +
    geom_vline(data = change_points, aes(xintercept = pt_median),
      inherit.aes = FALSE, linewidth = 0.65, linetype = "dashed",
      color = "#374151") +
    geom_ribbon(aes(ymin = score_q25, ymax = score_q75),
      fill = "#9CA3AF", alpha = 0.22) +
    geom_line(color = "#374151", linewidth = 0.55) +
    geom_point(aes(color = score_mean), size = 1.35) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-cap, cap), oob = scales::squish,
      name = "Mean rank\nchange") +
    labs(title = "Coordinated neighborhood-rank change along trajectory",
      subtitle = "Mean across genes; ribbon is gene-wise IQR; support is reported separately",
      x = "Existing Bdev pseudotime (pt)",
      y = "Change versus five windows earlier",
      caption = "Coordination support requires >=20% of genes to change >=0.05 in the score direction.") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "right",
      plot.title = element_text(face = "bold"))
}

gene_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
driver_plot <- function(data) {
  data <- data[is.finite(data$rank_change_vs_five_windows_earlier), , drop = FALSE]
  ranks <- unique(data[c("feature_id", "driver_rank")])
  ranks <- ranks[order(ranks$driver_rank), , drop = FALSE]
  driver_genes <- ranks$feature_id
  data$feature_id <- factor(data$feature_id, levels = driver_genes)
  colors <- setNames(gene_palette, driver_genes)
  ends <- do.call(rbind, lapply(driver_genes, function(gene) {
    x <- data[as.character(data$feature_id) == gene, , drop = FALSE]
    x[which.max(x$window_midpoint), , drop = FALSE]
  }))
  span <- diff(range(data$window_midpoint))
  ggplot(data, aes(window_midpoint,
      rank_change_vs_five_windows_earlier, color = feature_id,
      group = feature_id)) +
    geom_hline(yintercept = 0, linewidth = 0.45, color = "#6B7280") +
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
    labs(title = "Top-five driver-gene rank changes",
      subtitle = "Each line is that gene's k=30 neighborhood-rank change",
      x = "Existing Bdev pseudotime (pt)",
      y = "Gene rank change versus five windows earlier",
      caption = "Driver identities/ranks are unchanged; dashed lines mark BB1/BB2.") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
      legend.position = "none", plot.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 80, 5.5, 5.5))
}

log_event("Rendering 38 neighborhood rank-change figures")
pdf_path <- file.path(figure_dir,
  "all_38_profiles_neighborhood_rank_change.pdf")
grDevices::pdf(pdf_path, width = 15, height = 13, onefile = TRUE)
figure_rows <- vector("list", nrow(inventory))
for (i in seq_len(nrow(inventory))) {
  row <- inventory[i, ]; id <- row$signature_id
  label <- paste(row$cytokine, toupper(row$direction), "-", row$final_size,
    "genes")
  cell_values <- cell_score_matrix[id, ]
  trajectory_data <- signature_trajectory[
    signature_trajectory$feature_set == id, , drop = FALSE]
  driver_data <- driver_change[driver_change$feature_set == id, , drop = FALSE]
  combined <- map_plot(cell_values, "Neighborhood gene-rank change UMAP", label) /
    trajectory_plot(trajectory_data) / driver_plot(driver_data) +
    plot_layout(heights = c(1, 0.9, 0.9)) +
    plot_annotation(title = paste("Coordinated mean change:", label),
      subtitle = paste("k=30 RNA-PCA neighborhoods; within-gene percentile ranks;",
        "each gene compared with five windows earlier before mean aggregation"),
      theme = theme(plot.title = element_text(size = 17, face = "bold")))
  print(combined)
  file <- paste0(sprintf("%02d_", i), slug(paste(row$cytokine,
    row$direction)), "_neighborhood_rank_change.png")
  ggsave(file.path(figure_dir, file), combined, width = 15, height = 13,
    dpi = 220, units = "in", bg = "white")
  local <- trajectory_data[is.finite(trajectory_data$score_mean), ]
  figure_rows[[i]] <- data.frame(signature_id = id,
    cytokine = row$cytokine, direction = row$direction,
    signature_size = as.integer(row$final_size),
    strongest_activation_pt = local$window_midpoint[which.max(local$score_mean)],
    strongest_activation = max(local$score_mean),
    strongest_deactivation_pt = local$window_midpoint[which.min(local$score_mean)],
    strongest_deactivation = min(local$score_mean),
    supported_window_fraction = mean(local$coordinated_support),
    finite_umap_cells = sum(is.finite(cell_values)), file = file,
    stringsAsFactors = FALSE)
}
grDevices::dev.off()
figure_inventory <- do.call(rbind, figure_rows)
write_tsv(figure_inventory, file.path(stage_dir, "figure_inventory.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
first_valid_window <- min(which(is.finite(gene_window_change[1L, ])))
checks <- list(source_has_38_passing_profiles = nrow(inventory) == 38L,
  all_signatures_have_at_least_20_genes = all(inventory$final_size >= 20L),
  saved_native_neighbors_reused = nrow(neighbors) == length(cells) * 50L,
  k30_includes_query = identical(neighbor_index[, 1L], seq_along(cells)),
  all_signature_genes_used = setequal(genes, unique(members$feature_id)),
  gene_percentiles_finite_and_bounded = all(is.finite(gene_percentiles)) &&
    all(gene_percentiles >= 0 & gene_percentiles <= 1),
  five_window_lag = lag_windows == 5L,
  first_five_windows_have_no_change = first_valid_window == lag_windows + 1L,
  cell_scores_finite_and_bounded =
    all(is.finite(cell_score_matrix[, scoreable_cells, drop = FALSE])) &&
    all(abs(cell_score_matrix[, scoreable_cells, drop = FALSE]) <= 1),
  cell_support_finite_and_bounded =
    all(is.finite(cell_support_matrix[, scoreable_cells, drop = FALSE])) &&
    all(cell_support_matrix[, scoreable_cells, drop = FALSE] >= 0 &
      cell_support_matrix[, scoreable_cells, drop = FALSE] <= 1),
  missing_pt_cells_have_missing_score =
    all(is.na(cell_score_matrix[, is.na(metadata$pt), drop = FALSE])),
  exact_source_pt_preserved = identical(as.numeric(metadata$pt),
    as.numeric(cell_scores$pt[cell_scores$key == inventory$signature_id[[1L]]])),
  source_missing_pt_preserved = sum(is.na(metadata$pt)) == 14081L,
  five_drivers_per_profile = all(table(driver_change$feature_set) ==
    5L * n_windows),
  every_figure_exists = all(file.exists(file.path(figure_dir,
    figure_inventory$file))), multipage_pdf_exists = file.exists(pdf_path),
  bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after), result_promotions = 0L)
stopifnot(all(unlist(checks[names(checks) != "result_promotions"])),
  checks$result_promotions == 0L)
write_json(list(status = "passed", checks = checks),
  file.path(stage_dir, "checks.json"))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
provenance_inputs <- c(ucell_paths, driver_paths, neighborhood_paths,
  file.path(bundle_dir, "manifest.json"), script_path,
  "figure_scripts/bdev_reference/run_12_all_cytokine_neighborhood_rank_change_report.sh",
  "config/bdev_reference.yml")
input_manifest <- data.frame(path = normalizePath(provenance_inputs,
  mustWork = TRUE), sha256 = vapply(provenance_inputs, wb_source_sha256,
  character(1)), stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

html_rows <- vapply(seq_len(nrow(figure_inventory)), function(i) {
  x <- figure_inventory[i, ]
  paste0("<h2>", x$cytokine, " ", toupper(x$direction), "</h2>",
    "<p>Strongest activation: ", sprintf("%.3f", x$strongest_activation),
    " at pt ", sprintf("%.2f", x$strongest_activation_pt),
    "; strongest deactivation: ", sprintf("%.3f", x$strongest_deactivation),
    " at pt ", sprintf("%.2f", x$strongest_deactivation_pt), ".</p>",
    '<img src="profile_figures/', x$file, '" alt="', x$cytokine, " ",
    x$direction, ' neighborhood rank change">')
}, character(1))
html <- c('<!doctype html><html><head><meta charset="utf-8">',
  '<title>Neighborhood gene-rank change along Bdev pseudotime</title>',
  '<style>body{font:16px system-ui,sans-serif;max-width:1500px;margin:32px auto;padding:0 22px;color:#172033}img{max-width:100%;height:auto;border:1px solid #ddd}p,li{line-height:1.5}code{overflow-wrap:anywhere}a{color:#1261a0}</style></head><body>',
  '<h1>Neighborhood gene-rank change along Bdev pseudotime</h1>',
  '<p>For every signature gene, SCT expression is averaged over the cell and its 29 nearest cells in the existing 50-dimensional RNA-PCA space, then converted to a within-gene percentile. Each gene is compared with its own median percentile five 1%-step windows earlier. The signature trajectory is the mean of these bounded gene-level changes; positive values indicate net increases and negative values net decreases.</p>',
  '<p>The UMAP score is cell-specific: each cell\'s gene percentiles are compared gene by gene with the corresponding lagged-window medians before taking the signature mean. Directional support is retained separately and requires at least 20% of genes to change by at least 0.05 in the score direction; it is not multiplied into the continuous score. Cells without source <code>pt</code>, and cells without a full five-window lag, are gray. No matched-set calibration is applied.</p>',
  '<p><a href="profile_figures/all_38_profiles_neighborhood_rank_change.pdf">Download 38-page PDF</a> · <a href="signature_neighborhood_rank_change_trajectory.tsv">Signature trajectories</a> · <a href="top5_driver_neighborhood_rank_change.tsv">Driver changes</a> · <a href="neighborhood_rank_change_cell_scores.tsv.gz">Per-cell scores</a></p>',
  html_rows, '</body></html>')
writeLines(html, file.path(stage_dir, "REPORT.html"))
writeLines(c("# Neighborhood gene-rank change report", "",
  "[Open report](REPORT.html)", "",
  "All 38 profiles show a per-cell mean-change UMAP, mean gene-level change along exact Bdev pt, and five driver-gene change curves."),
  file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", promoted = FALSE,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  analysis_scope = paste("user-authorized neighborhood-smoothed gene-level",
    "rank-change analysis; no promotion"),
  source_ucell_report = ucell_report,
  source_driver_report = driver_report,
  source_neighborhood_report = neighborhood_report,
  bundle = bundle_dir, bundle_fingerprint = fingerprint_before,
  profiles = 38L, cytokines = length(unique(inventory$cytokine)),
  neighborhood = list(space = "pca", dimensions = 50L, k = k,
    include_query = TRUE, source = "saved native scWorkbench knn_neighbors"),
  gene_score = list(matrix = "SCT_data",
    smoothing = "mean over k=30 neighborhood",
    transform = "zero-aware within-gene percentile across all Bdev cells"),
  trajectory_change = list(window_width_fraction = 0.05,
    window_step_fraction = 0.01, lag_windows = lag_windows,
    comparison = "per-gene median rank current minus five windows earlier",
    aggregation = "mean across all signature genes",
    median_retained = "diagnostic only",
    band = "gene-wise IQR"),
  cell_change = paste("per-cell current gene percentile minus per-gene",
    "lagged-window median, then mean across signature genes"),
  directional_support = list(gene_change_threshold = 0.05,
    required_fraction = 0.20,
    use = "separate confidence flag; not multiplied into score"),
  matched_set_permutation = FALSE,
  driver_selection_changed = FALSE,
  driver_boundaries = "BB1/BB2 dashed lines at saved median pseudotimes",
  defined_pt_cells = 12642L, missing_pt_cells = 14081L,
  scoreable_umap_cells = length(scoreable_cells), result_promotions = 0L,
  multipage_pdf = file.path("profile_figures", basename(pdf_path)),
  figures = figure_inventory$file,
  report_sha256 = wb_source_sha256(file.path(stage_dir, "REPORT.html"))),
  file.path(stage_dir, "completion.json"))

if (!file.rename(stage_dir, output_dir)) {
  stop("Failed to atomically finalize output directory")
}
log_event("Completed neighborhood gene-rank change report at", output_dir)
