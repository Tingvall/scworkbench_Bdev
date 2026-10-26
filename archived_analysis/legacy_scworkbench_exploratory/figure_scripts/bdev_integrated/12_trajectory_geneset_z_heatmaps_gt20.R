#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L)
  stop("Usage: 12_trajectory_geneset_z_heatmaps_gt20.R <report_root>")
report_root <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- file.path(report_root, "summary_figures", "z_heatmaps_gt20")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

groups <- c("cytokines", "bp", "mf", "cc")
group_labels <- c(cytokines = "Cytokines", bp = "GO BP",
  mf = "GO MF", cc = "GO CC")
read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE, stringsAsFactors = FALSE))
}
cluster_order <- function(x) {
  if (nrow(x) == 1L) return(1L)
  similarity <- cor(t(x), use = "pairwise.complete.obs")
  similarity[!is.finite(similarity)] <- 0
  diag(similarity) <- 1
  distance <- 1 - similarity
  distance[distance < 0] <- 0
  hclust(as.dist(distance), method = "average")$order
}
make_heatmap <- function(long, row_order, title, group, color_limit) {
  row_levels <- rev(row_order)
  long[, plot_row := factor(feature_set, levels = row_levels)]
  tile_width <- min(diff(sort(unique(long$window_center))))
  label_lookup <- unique(long[, .(feature_set, set_name)])
  labels <- setNames(label_lookup$set_name, label_lookup$feature_set)
  p <- ggplot(long, aes(window_center, plot_row, fill = value)) +
    geom_tile(width = tile_width, height = 1) +
    scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
      midpoint = 0, limits = c(-color_limit, color_limit),
      oob = scales::squish, name = "z-score") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = title, x = "Existing Bdev pseudotime (pt)", y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
      axis.ticks.y = element_blank())
  if (group == "cytokines") {
    p <- p + scale_y_discrete(labels = labels, drop = FALSE) +
      theme(axis.text.y = element_text(size = 8))
  } else {
    p <- p + scale_y_discrete(breaks = NULL, drop = FALSE)
  }
  p
}

manifest_rows <- list()
for (g in groups) {
  group_dir <- file.path(report_root, g)
  curves <- read_tsv_gz(file.path(group_dir, "window_curves.tsv.gz"))
  regions <- fread(file.path(group_dir, "significant_positive_regions.tsv"))
  selected <- unique(regions[is.finite(peak_z) & peak_z > 20, feature_set])
  stopifnot(length(selected) > 0L)
  curves <- curves[feature_set %in% selected]
  setorder(curves, feature_set, window_index)
  stopifnot(uniqueN(curves$feature_set) == length(selected),
    !anyDuplicated(curves[, .(feature_set, window_index)]))

  window_map <- unique(curves[, .(window_index, window_center)])
  setorder(window_map, window_index)
  full_wide <- dcast(curves, feature_set ~ window_index,
    value.var = "z_score")
  feature_sets <- full_wide$feature_set
  full_matrix <- as.matrix(full_wide[, -"feature_set"])
  rownames(full_matrix) <- feature_sets
  stopifnot(all(is.finite(full_matrix)))

  mask <- matrix(FALSE, nrow(full_matrix), ncol(full_matrix),
    dimnames = dimnames(full_matrix))
  for (i in seq_along(feature_sets)) {
    fs <- feature_sets[[i]]
    intervals <- regions[feature_set == fs, .(region_left, region_right)]
    for (j in seq_len(nrow(intervals))) {
      mask[i, ] <- mask[i, ] |
        (window_map$window_center >= intervals$region_left[[j]] &
          window_map$window_center <= intervals$region_right[[j]])
    }
  }
  masked_matrix <- full_matrix
  masked_matrix[!mask] <- 0

  order_full <- feature_sets[cluster_order(full_matrix)]
  order_masked <- feature_sets[cluster_order(masked_matrix)]
  set_names <- unique(regions[feature_set %in% selected,
    .(feature_set, set_name)])
  selected_summary <- regions[feature_set %in% selected,
    .(n_final_significant_regions = .N,
      n_regions_peak_z_gt20 = sum(peak_z > 20),
      maximum_peak_z = max(peak_z)), by = feature_set]
  write_order <- function(order, basis) {
    out <- data.table(analysis_group = g, clustering_basis = basis,
      row_position_top_to_bottom = seq_along(order), feature_set = order)
    out <- merge(out, set_names, by = "feature_set", all.x = TRUE, sort = FALSE)
    out <- merge(out, selected_summary, by = "feature_set", all.x = TRUE,
      sort = FALSE)
    setorder(out, row_position_top_to_bottom)
    fwrite(out, file.path(output_dir,
      paste0(g, "_row_order_cluster_", basis, ".tsv")), sep = "\t",
      quote = FALSE, na = "NA")
  }
  write_order(order_full, "complete_z")
  write_order(order_masked, "significant_regions_only")

  base_long <- merge(curves[, .(feature_set, window_index, window_center,
    value = z_score)], set_names, by = "feature_set", all.x = TRUE,
    sort = FALSE)
  setorder(base_long, feature_set, window_index)
  masked_long <- copy(base_long)
  masked_long[, value := as.numeric(t(masked_matrix)[cbind(
    match(window_index, window_map$window_index),
    match(feature_set, feature_sets))])]
  stopifnot(all(is.finite(masked_long$value)))
  color_limit <- max(abs(full_matrix), abs(masked_matrix))

  full_by_full <- make_heatmap(copy(base_long), order_full,
    "Complete z-score profile", g, color_limit)
  masked_by_full <- make_heatmap(copy(masked_long), order_full,
    "Significant regions only — complete-profile order", g, color_limit)
  fig_full <- full_by_full + masked_by_full +
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(title = paste0(group_labels[[g]],
      ": clustering based on complete z-score profiles"),
      subtitle = paste0(length(selected),
        " gene sets with at least one final region having peak z > 20"),
      theme = theme(plot.title = element_text(face = "bold", size = 14)))

  masked_by_masked <- make_heatmap(copy(masked_long), order_masked,
    "Significant regions only", g, color_limit)
  full_by_masked <- make_heatmap(copy(base_long), order_masked,
    "Complete profile — significant-region order", g, color_limit)
  fig_masked <- masked_by_masked + full_by_masked +
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(title = paste0(group_labels[[g]],
      ": clustering based on significant-region-only profiles"),
      subtitle = paste0(length(selected),
        " gene sets with at least one final region having peak z > 20"),
      theme = theme(plot.title = element_text(face = "bold", size = 14)))

  height <- if (g == "cytokines") 6 else 8
  for (format in c("png", "pdf")) {
    device <- if (format == "pdf") cairo_pdf else format
    dpi <- if (format == "png") 220 else 300
    ggsave(file.path(output_dir,
      paste0(g, "_paired_heatmaps_cluster_complete_z.", format)), fig_full,
      width = 14, height = height, dpi = dpi, device = device, bg = "white")
    ggsave(file.path(output_dir,
      paste0(g, "_paired_heatmaps_cluster_significant_regions_only.", format)),
      fig_masked, width = 14, height = height, dpi = dpi, device = device,
      bg = "white")
  }
  manifest_rows[[g]] <- data.table(analysis_group = g,
    group_label = group_labels[[g]], selected_gene_sets = length(selected),
    pseudotime_windows = ncol(full_matrix), color_limit = color_limit,
    row_labels_printed = g == "cytokines",
    clustering_distance = "1 - Pearson correlation",
    clustering_linkage = "average")
}

manifest <- rbindlist(manifest_rows)
fwrite(manifest, file.path(output_dir, "heatmap_manifest.tsv"), sep = "\t",
  quote = FALSE, na = "NA")
