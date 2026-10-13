#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 200)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L)
  stop(paste("Usage: 13_trajectory_geneset_pattern_clusters_gt20.R",
    "<report_root> <bb_up_change_points.tsv>"))
report_root <- normalizePath(args[[1L]], mustWork = TRUE)
bb_path <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- file.path(report_root, "summary_figures",
  "z_heatmap_clusters_gt20")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

groups <- c("bp", "mf", "cc")
group_labels <- c(bp = "GO BP", mf = "GO MF", cc = "GO CC")
n_clusters <- 15L
bb <- fread(bb_path)
stopifnot(nrow(bb) == 2L, identical(as.character(bb$label), c("up1", "up2")),
  all(is.finite(bb$pt_median)))
bb[, boundary_label := c("BB1", "BB2")]
read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE, stringsAsFactors = FALSE))
}
correlation_hclust <- function(x) {
  similarity <- cor(t(x), use = "pairwise.complete.obs")
  similarity[!is.finite(similarity)] <- 0
  diag(similarity) <- 1
  distance <- 1 - similarity
  distance[distance < 0] <- 0
  hclust(as.dist(distance), method = "average")
}
safe_cor <- function(x, y) {
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y)
}
common_keywords <- function(names, maximum = 6L) {
  stop_words <- c("GO", "GOBP", "GOMF", "GOCC", "OF", "IN", "TO", "BY",
    "AND", "OR", "THE", "A", "AN", "PROCESS", "CELLULAR", "REGULATION",
    "POSITIVE", "NEGATIVE", "ACTIVITY", "PROTEIN", "DEPENDENT", "INVOLVED")
  per_set <- lapply(names, function(x) {
    words <- unlist(strsplit(gsub("[^A-Z0-9]+", " ", toupper(x)), " +"))
    unique(words[nchar(words) >= 3L & !words %in% stop_words])
  })
  counts <- sort(table(unlist(per_set)), decreasing = TRUE)
  if (!length(counts)) return("")
  counts <- head(counts, maximum)
  paste0(names(counts), " (", as.integer(counts), ")", collapse = "; ")
}
make_ordered_heatmap <- function(long, row_order, cluster_membership, title,
  color_limit) {
  mapping <- data.table(feature_set = row_order,
    row_position_top_to_bottom = seq_along(row_order))
  d <- merge(long, mapping, by = "feature_set", all.x = TRUE, sort = FALSE)
  n <- length(row_order)
  d[, y := n - row_position_top_to_bottom + 1]
  tile_width <- min(diff(sort(unique(d$window_center))))
  ordered_clusters <- cluster_membership[
    match(row_order, cluster_membership$feature_set), cluster_id]
  runs <- rle(ordered_clusters)
  run_end <- cumsum(runs$lengths)
  run_start <- c(1L, head(run_end, -1L) + 1L)
  mid_top <- (run_start + run_end) / 2
  y_breaks <- n - mid_top + 1
  boundary_y <- n - head(run_end, -1L) + .5
  ggplot(d, aes(window_center, y, fill = value)) +
    geom_tile(width = tile_width, height = 1) +
    geom_hline(yintercept = boundary_y, linewidth = .22, color = "#202020") +
    geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
      linewidth = .55, linetype = "dashed", color = "#374151") +
    geom_text(data = bb, aes(x = pt_median, y = Inf, label = boundary_label),
      inherit.aes = FALSE, vjust = 1.25, size = 3, fontface = "bold",
      color = "#374151") +
    scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
      midpoint = 0, limits = c(-color_limit, color_limit),
      oob = scales::squish, name = "z-score") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = y_breaks, labels = runs$values,
      expand = c(0, 0)) +
    labs(title = title, x = "Existing Bdev pseudotime (pt)",
      y = "Pattern cluster") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold"))
}

membership_list <- list(); summary_list <- list(); centroid_list <- list()
markdown <- c("# GO gene-set pseudotime pattern clusters", "",
  "Gene sets were selected when at least one final region had peak z > 20. Clusters use complete z-score profiles, Pearson-correlation distance, and average-linkage hierarchical clustering. The significant-region-only profile is shown only as a companion in the same order. Cytokines were excluded.", "")

for (g in groups) {
  group_dir <- file.path(report_root, g)
  curves <- read_tsv_gz(file.path(group_dir, "window_curves.tsv.gz"))
  regions <- fread(file.path(group_dir, "significant_positive_regions.tsv"))
  selected <- unique(regions[is.finite(peak_z) & peak_z > 20, feature_set])
  curves <- curves[feature_set %in% selected]
  setorder(curves, feature_set, window_index)
  window_map <- unique(curves[, .(window_index, window_center)])
  setorder(window_map, window_index)
  full_wide <- dcast(curves, feature_set ~ window_index, value.var = "z_score")
  feature_sets <- full_wide$feature_set
  full_matrix <- as.matrix(full_wide[, -"feature_set"])
  rownames(full_matrix) <- feature_sets
  mask <- matrix(FALSE, nrow(full_matrix), ncol(full_matrix),
    dimnames = dimnames(full_matrix))
  for (i in seq_along(feature_sets)) {
    intervals <- regions[feature_set == feature_sets[[i]],
      .(region_left, region_right)]
    for (j in seq_len(nrow(intervals))) mask[i, ] <- mask[i, ] |
      (window_map$window_center >= intervals$region_left[[j]] &
        window_map$window_center <= intervals$region_right[[j]])
  }
  masked_matrix <- full_matrix
  masked_matrix[!mask] <- 0

  hc <- correlation_hclust(full_matrix)
  raw_cluster <- cutree(hc, k = n_clusters)
  ordered_sets <- feature_sets[hc$order]
  raw_in_order <- raw_cluster[ordered_sets]
  raw_levels <- unique(raw_in_order)
  cluster_map <- setNames(sprintf("C%02d", seq_along(raw_levels)), raw_levels)
  cluster_id <- unname(cluster_map[as.character(raw_cluster)])
  membership <- data.table(analysis_group = g, feature_set = feature_sets,
    cluster_id = cluster_id)
  membership[, row_position_top_to_bottom := match(feature_set, ordered_sets)]
  setorder(membership, row_position_top_to_bottom)
  positions <- split(membership$row_position_top_to_bottom,
    membership$cluster_id)
  stopifnot(all(vapply(positions, function(x)
    max(x) - min(x) + 1L == length(x), logical(1))))

  set_names <- unique(regions[feature_set %in% selected,
    .(feature_set, set_name)])
  best_regions <- copy(regions[feature_set %in% selected])
  setorder(best_regions, feature_set, -peak_z)
  region_counts <- best_regions[, .(n_final_significant_regions = .N),
    by = feature_set]
  best_regions <- best_regions[, .SD[1L], by = feature_set][,
    .(feature_set, strongest_region_peak_z = peak_z,
      strongest_region_peak_pt = peak_pseudotime)]
  best_regions <- merge(best_regions, region_counts, by = "feature_set",
    all.x = TRUE, sort = FALSE)
  membership <- merge(membership, set_names, by = "feature_set", all.x = TRUE,
    sort = FALSE)
  membership <- merge(membership, best_regions, by = "feature_set", all.x = TRUE,
    sort = FALSE)
  setorder(membership, row_position_top_to_bottom)

  centroids <- rbindlist(lapply(sprintf("C%02d", seq_len(n_clusters)),
    function(cid) {
      idx <- which(cluster_id == cid)
      data.table(analysis_group = g, cluster_id = cid,
        window_index = window_map$window_index,
        window_center = window_map$window_center,
        mean_masked_z = colMeans(masked_matrix[idx, , drop = FALSE]),
        mean_complete_z = colMeans(full_matrix[idx, , drop = FALSE]))
    }))
  cluster_summaries <- rbindlist(lapply(sprintf("C%02d", seq_len(n_clusters)),
    function(cid) {
      idx <- which(cluster_id == cid)
      centroid <- colMeans(full_matrix[idx, , drop = FALSE])
      correlations <- vapply(idx, function(i)
        safe_cor(full_matrix[i, ], centroid), numeric(1))
      ranking <- order(-correlations,
        -membership$strongest_region_peak_z[match(feature_sets[idx],
          membership$feature_set)], na.last = TRUE)
      representative_ids <- feature_sets[idx][head(ranking, 5L)]
      representative_names <- set_names$set_name[match(representative_ids,
        set_names$feature_set)]
      member_rows <- membership[cluster_id == cid]
      peak_i <- which.max(centroid)
      positive_weight <- pmax(centroid, 0)
      weighted_pt <- if (sum(positive_weight) > 0)
        weighted.mean(window_map$window_center, positive_weight) else NA_real_
      data.table(analysis_group = g, cluster_id = cid,
        n_gene_sets = length(idx), centroid_peak_pt = window_map$window_center[[peak_i]],
        centroid_peak_mean_z = centroid[[peak_i]],
        centroid_positive_weighted_pt = weighted_pt,
        member_peak_pt_q25 = quantile(member_rows$strongest_region_peak_pt, .25),
        member_peak_pt_median = median(member_rows$strongest_region_peak_pt),
        member_peak_pt_q75 = quantile(member_rows$strongest_region_peak_pt, .75),
        maximum_member_peak_z = max(member_rows$strongest_region_peak_z),
        representative_gene_sets = paste(representative_names, collapse = "; "),
        recurring_name_keywords = common_keywords(member_rows$set_name))
    }))
  membership_list[[g]] <- membership
  summary_list[[g]] <- cluster_summaries
  centroid_list[[g]] <- centroids

  base_long <- merge(curves[, .(feature_set, window_index, window_center,
    value = z_score)], set_names, by = "feature_set", all.x = TRUE,
    sort = FALSE)
  setorder(base_long, feature_set, window_index)
  masked_long <- copy(base_long)
  masked_long[, value := as.numeric(t(masked_matrix)[cbind(
    match(window_index, window_map$window_index),
    match(feature_set, feature_sets))])]
  color_limit <- max(abs(full_matrix), abs(masked_matrix))
  left <- make_ordered_heatmap(base_long, ordered_sets, membership,
    "Complete z-score profile", color_limit)
  right <- make_ordered_heatmap(masked_long, ordered_sets, membership,
    "Significant regions only — same order", color_limit)
  paired <- left + right + plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(title = paste0(group_labels[[g]],
      ": 15 complete-z pattern clusters"),
      subtitle = paste0(length(selected),
        " gene sets with at least one final region having peak z > 20"),
      theme = theme(plot.title = element_text(face = "bold", size = 14)))
  ggsave(file.path(output_dir, paste0(g, "_15_complete_z_pattern_clusters_paired.png")),
    paired, width = 14, height = 9, dpi = 220, bg = "white")
  ggsave(file.path(output_dir, paste0(g, "_15_complete_z_pattern_clusters_paired.pdf")),
    paired, width = 14, height = 9, device = cairo_pdf, bg = "white")

  centroid_plot_data <- melt(centroids,
    id.vars = c("analysis_group", "cluster_id", "window_index", "window_center"),
    measure.vars = "mean_complete_z", variable.name = "profile",
    value.name = "value")
  centroid_plot_data[, cluster_id := factor(cluster_id,
    levels = rev(sprintf("C%02d", seq_len(n_clusters))))]
  tile_width <- min(diff(window_map$window_center))
  centroid_limit <- max(abs(centroid_plot_data$value))
  centroid_plot <- ggplot(centroid_plot_data,
    aes(window_center, cluster_id, fill = value)) +
    geom_tile(width = tile_width, height = 1) +
    geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
      linewidth = .55, linetype = "dashed", color = "#374151") +
    geom_text(data = bb, aes(x = pt_median, y = Inf, label = boundary_label),
      inherit.aes = FALSE, vjust = 1.25, size = 3, fontface = "bold",
      color = "#374151") +
    scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
      midpoint = 0, limits = c(-centroid_limit, centroid_limit),
      oob = scales::squish, name = "Mean z") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = paste0(group_labels[[g]],
      ": mean complete z-score profile by cluster"),
      x = "Existing Bdev pseudotime (pt)", y = "Pattern cluster") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))
  ggsave(file.path(output_dir, paste0(g, "_15_complete_z_cluster_centroids.png")),
    centroid_plot, width = 8.5, height = 6, dpi = 220, bg = "white")
  ggsave(file.path(output_dir, paste0(g, "_15_complete_z_cluster_centroids.pdf")),
    centroid_plot, width = 8.5, height = 6, device = cairo_pdf, bg = "white")

  markdown <- c(markdown, paste0("## ", group_labels[[g]]), "",
    "| Cluster | Sets | Centroid peak pt | Member peak pt median [IQR] | Representative gene sets | Recurring keywords |",
    "|---|---:|---:|---|---|---|",
    vapply(seq_len(nrow(cluster_summaries)), function(i) paste0("| ",
      cluster_summaries$cluster_id[[i]], " | ",
      cluster_summaries$n_gene_sets[[i]], " | ",
      signif(cluster_summaries$centroid_peak_pt[[i]], 5), " | ",
      signif(cluster_summaries$member_peak_pt_median[[i]], 5), " [",
      signif(cluster_summaries$member_peak_pt_q25[[i]], 5), ", ",
      signif(cluster_summaries$member_peak_pt_q75[[i]], 5), "] | ",
      cluster_summaries$representative_gene_sets[[i]], " | ",
      cluster_summaries$recurring_name_keywords[[i]], " |"), character(1)), "")
}

all_membership <- rbindlist(membership_list, use.names = TRUE, fill = TRUE)
all_summary <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)
all_centroids <- rbindlist(centroid_list, use.names = TRUE, fill = TRUE)
fwrite(all_membership, file.path(output_dir, "gene_set_cluster_membership.tsv"),
  sep = "\t", quote = FALSE, na = "NA")
fwrite(all_summary, file.path(output_dir, "cluster_summary.tsv"), sep = "\t",
  quote = FALSE, na = "NA")
fwrite(all_centroids, file.path(output_dir, "cluster_centroid_profiles.tsv"),
  sep = "\t", quote = FALSE, na = "NA")
writeLines(markdown, file.path(output_dir, "REPORT.md"))
