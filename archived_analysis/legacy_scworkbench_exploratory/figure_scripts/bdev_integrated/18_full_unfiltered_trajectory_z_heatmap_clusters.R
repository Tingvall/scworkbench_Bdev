#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 200)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 18_full_unfiltered_trajectory_z_heatmap_clusters.R",
  "<report_root> <bb_up_change_points.tsv> <output_dir>"))
report_root <- normalizePath(args[[1L]], mustWork = TRUE)
bb_path <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

groups <- c("cytokines", "bp", "mf", "cc")
group_labels <- c(cytokines = "Cytokines", bp = "GO BP", mf = "GO MF",
  cc = "GO CC")
bb <- fread(bb_path)
stopifnot(nrow(bb) == 2L,
  identical(as.character(bb$label), c("up1", "up2")),
  all(is.finite(bb$pt_median)))
bb[, boundary_label := c("BB1", "BB2")]

read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
safe_cor <- function(x, y) {
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y)
}
correlation_hclust <- function(x) {
  similarity <- cor(t(x), use = "pairwise.complete.obs")
  similarity[!is.finite(similarity)] <- 0
  diag(similarity) <- 1
  distance <- 1 - similarity
  distance[distance < 0] <- 0
  hc <- hclust(as.dist(distance), method = "average")
  rm(similarity, distance)
  hc
}
candidate_k <- function(n) {
  if (n <= 50L) return(3:min(10L, n - 1L))
  if (n <= 600L) return(6:20)
  if (n <= 1200L) return(6:24)
  8:30
}
choose_k <- function(x, hc, seed = 314159L) {
  n <- nrow(x)
  candidates <- candidate_k(n)
  set.seed(seed)
  sampled <- if (n <= 1200L) seq_len(n) else sort(sample.int(n, 1200L))
  sample_similarity <- cor(t(x[sampled, , drop = FALSE]),
    use = "pairwise.complete.obs")
  sample_similarity[!is.finite(sample_similarity)] <- 0
  diag(sample_similarity) <- 1
  sample_distance_matrix <- 1 - sample_similarity
  sample_distance_matrix[sample_distance_matrix < 0] <- 0
  sample_distance <- as.dist(sample_distance_matrix)
  results <- rbindlist(lapply(candidates, function(k) {
    labels <- cutree(hc, k = k)[sampled]
    sizes <- table(labels)
    sil <- cluster::silhouette(labels, sample_distance)
    data.table(k = k, mean_silhouette = mean(sil[, "sil_width"]),
      median_silhouette = median(sil[, "sil_width"]),
      smallest_cluster = min(sizes), largest_cluster = max(sizes))
  }))
  # Avoid solutions supported by tiny clusters. Among effectively tied
  # solutions, prefer the smaller and more interpretable number of clusters.
  eligible <- results[smallest_cluster >= max(2L, floor(length(sampled) * .005))]
  if (!nrow(eligible)) eligible <- results
  best <- max(eligible$mean_silhouette)
  selected <- eligible[mean_silhouette >= best - 0.01][which.min(k), k]
  list(k = selected, diagnostics = results, sample_n = length(sampled))
}
common_keywords <- function(x, maximum = 6L) {
  stop_words <- c("GO", "GOBP", "GOMF", "GOCC", "OF", "IN", "TO", "BY",
    "AND", "OR", "THE", "A", "AN", "PROCESS", "CELLULAR", "REGULATION",
    "POSITIVE", "NEGATIVE", "ACTIVITY", "PROTEIN", "DEPENDENT", "INVOLVED")
  words <- lapply(x, function(z) {
    y <- unlist(strsplit(gsub("[^A-Z0-9]+", " ", toupper(z)), " +"))
    unique(y[nchar(y) >= 3L & !y %in% stop_words])
  })
  counts <- sort(table(unlist(words)), decreasing = TRUE)
  if (!length(counts)) return("")
  counts <- head(counts, maximum)
  paste0(names(counts), " (", as.integer(counts), ")", collapse = "; ")
}

all_membership <- list()
all_centroids <- list()
all_summaries <- list()
all_k_diagnostics <- list()
all_unplottable <- list()
completion_groups <- list()
report <- c("# Complete unfiltered trajectory z-score heatmaps", "",
  paste0("Every gene set with a defined complete z-score profile is included. ",
    "No significance, peak-z, ",
    "GO-hierarchy, or redundancy filter was applied. Rows were clustered from ",
    "their complete 97-window z-score profiles using Pearson-correlation ",
    "distance and average linkage. Cluster number was selected separately per ",
    "group by approximate silhouette width; values were not rescaled or ",
    "truncated analytically. Only the heatmap color display is capped at the ",
    "99.5th percentile of absolute z-score."), "")

for (group_index in seq_along(groups)) {
  g <- groups[[group_index]]
  message("Processing ", g)
  curves <- read_tsv_gz(file.path(report_root, g, "window_curves.tsv.gz"))
  stopifnot(all(c("feature_set", "window_index", "window_center", "z_score") %in%
      names(curves)), all(is.finite(curves$z_score)))
  setorder(curves, feature_set, window_index)
  window_map <- unique(curves[, .(window_index, window_center)])
  setorder(window_map, window_index)
  wide <- dcast(curves, feature_set ~ window_index, value.var = "z_score")
  feature_sets <- wide$feature_set
  permutation_summary <- fread(file.path(report_root, g,
    "permutation_summary.tsv"))
  unplottable <- permutation_summary[!feature_set %chin% feature_sets]
  stopifnot(nrow(permutation_summary) == uniqueN(permutation_summary$feature_set),
    nrow(permutation_summary) == length(feature_sets) + nrow(unplottable))
  all_unplottable[[g]] <- unplottable
  z <- as.matrix(wide[, -"feature_set"])
  rownames(z) <- feature_sets
  stopifnot(ncol(z) == nrow(window_map), !anyNA(z))

  hc <- correlation_hclust(z)
  choice <- choose_k(z, hc, seed = 314159L + group_index)
  k <- choice$k
  raw_cluster <- cutree(hc, k = k)
  ordered_sets <- feature_sets[hc$order]
  raw_levels <- unique(raw_cluster[ordered_sets])
  cluster_map <- setNames(sprintf("C%02d", seq_along(raw_levels)), raw_levels)
  cluster_id <- unname(cluster_map[as.character(raw_cluster)])
  membership <- data.table(analysis_group = g, feature_set = feature_sets,
    cluster_id = cluster_id,
    row_position_top_to_bottom = match(feature_sets, ordered_sets),
    maximum_window_z = apply(z, 1L, max),
    peak_window_pt = window_map$window_center[max.col(z, ties.method = "first")])
  regions <- fread(file.path(report_root, g, "significant_positive_regions.tsv"))
  significant_sets <- unique(regions$feature_set)
  membership[, has_final_significant_region := feature_set %chin% significant_sets]
  setorder(membership, row_position_top_to_bottom)

  centroids <- rbindlist(lapply(sprintf("C%02d", seq_len(k)), function(cid) {
    idx <- which(cluster_id == cid)
    data.table(analysis_group = g, cluster_id = cid,
      window_index = window_map$window_index,
      window_center = window_map$window_center,
      mean_z = colMeans(z[idx, , drop = FALSE]),
      median_z = apply(z[idx, , drop = FALSE], 2L, median))
  }))
  summaries <- rbindlist(lapply(sprintf("C%02d", seq_len(k)), function(cid) {
    idx <- which(cluster_id == cid)
    centroid <- colMeans(z[idx, , drop = FALSE])
    cors <- vapply(idx, function(i) safe_cor(z[i, ], centroid), numeric(1))
    rank_order <- order(-cors, -membership$maximum_window_z[
      match(feature_sets[idx], membership$feature_set)], na.last = TRUE)
    reps <- feature_sets[idx][head(rank_order, 5L)]
    member_sets <- feature_sets[idx]
    peak_i <- which.max(centroid)
    data.table(analysis_group = g, cluster_id = cid,
      n_gene_sets = length(idx),
      n_significant_gene_sets = sum(member_sets %chin% significant_sets),
      centroid_peak_pt = window_map$window_center[[peak_i]],
      centroid_peak_mean_z = centroid[[peak_i]],
      representative_gene_sets = paste(reps, collapse = "; "),
      recurring_name_keywords = common_keywords(member_sets))
  }))

  mapping <- data.table(feature_set = ordered_sets,
    row_position_top_to_bottom = seq_along(ordered_sets))
  plot_data <- merge(curves[, .(feature_set, window_index, window_center,
    value = z_score)], mapping, by = "feature_set", sort = FALSE)
  n_sets <- length(ordered_sets)
  plot_data[, y := n_sets - row_position_top_to_bottom + 1L]
  ordered_cluster <- membership[match(ordered_sets, feature_set), cluster_id]
  runs <- rle(ordered_cluster)
  run_end <- cumsum(runs$lengths)
  run_start <- c(1L, head(run_end, -1L) + 1L)
  y_breaks <- n_sets - (run_start + run_end) / 2 + 1L
  boundary_y <- n_sets - head(run_end, -1L) + 0.5
  color_limit <- as.numeric(quantile(abs(z), .995, na.rm = TRUE))
  tile_width <- min(diff(window_map$window_center))
  y_scale <- if (g == "cytokines") {
    cytokine_labels <- sub(
      "^parse_cytokine_rest_only_refinement_all_top100_v1\\.", "",
      ordered_sets)
    scale_y_continuous(breaks = seq_len(n_sets), labels = rev(cytokine_labels),
      expand = c(0, 0))
  } else {
    scale_y_continuous(breaks = y_breaks, labels = runs$values,
      expand = c(0, 0))
  }
  p <- ggplot(plot_data, aes(window_center, y, fill = value)) +
    geom_raster() +
    geom_hline(yintercept = boundary_y, linewidth = .2, color = "#202020") +
    geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
      linewidth = .5, linetype = "dashed", color = "#374151") +
    geom_text(data = bb, aes(x = pt_median, y = Inf, label = boundary_label),
      inherit.aes = FALSE, vjust = 1.25, size = 3, fontface = "bold",
      color = "#374151") +
    scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
      midpoint = 0, limits = c(-color_limit, color_limit),
      oob = scales::squish, name = "z-score") +
    scale_x_continuous(expand = c(0, 0)) + y_scale +
    labs(title = paste0(group_labels[[g]],
      ": complete unfiltered z-score profiles"),
      subtitle = paste0(n_sets, " defined gene-set profiles; ", k,
        " data-selected pattern clusters"),
      x = "Existing Bdev pseudotime (pt)",
      y = if (g == "cytokines") "Cytokine gene set" else "Pattern cluster") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold"))
  ggsave(file.path(stage_dir, paste0(g, "_full_unfiltered_z_heatmap.png")),
    p, width = 10, height = 10, dpi = 240, bg = "white")
  ggsave(file.path(stage_dir, paste0(g, "_full_unfiltered_z_heatmap.pdf")),
    p, width = 10, height = 10, device = cairo_pdf, bg = "white")

  centroid_data <- copy(centroids)
  centroid_data[, cluster_id := factor(cluster_id,
    levels = rev(sprintf("C%02d", seq_len(k))))]
  centroid_limit <- as.numeric(quantile(abs(centroid_data$mean_z), .995))
  cp <- ggplot(centroid_data, aes(window_center, cluster_id, fill = mean_z)) +
    geom_tile(width = tile_width, height = 1) +
    geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
      linewidth = .5, linetype = "dashed", color = "#374151") +
    geom_text(data = bb, aes(x = pt_median, y = Inf, label = boundary_label),
      inherit.aes = FALSE, vjust = 1.25, size = 3, fontface = "bold",
      color = "#374151") +
    scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
      midpoint = 0, limits = c(-centroid_limit, centroid_limit),
      oob = scales::squish, name = "Mean z") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = paste0(group_labels[[g]], ": cluster mean profiles"),
      subtitle = "Clusters derived from all tested gene sets",
      x = "Existing Bdev pseudotime (pt)", y = "Pattern cluster") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))
  ggsave(file.path(stage_dir, paste0(g, "_full_unfiltered_cluster_centroids.png")),
    cp, width = 9, height = 6.5, dpi = 240, bg = "white")
  ggsave(file.path(stage_dir, paste0(g, "_full_unfiltered_cluster_centroids.pdf")),
    cp, width = 9, height = 6.5, device = cairo_pdf, bg = "white")

  choice$diagnostics[, `:=`(analysis_group = g,
    sample_n = choice$sample_n, selected = k == choice$k)]
  setcolorder(choice$diagnostics, c("analysis_group", "sample_n", "k",
    "selected", "mean_silhouette", "median_silhouette",
    "smallest_cluster", "largest_cluster"))
  all_membership[[g]] <- membership
  all_centroids[[g]] <- centroids
  all_summaries[[g]] <- summaries
  all_k_diagnostics[[g]] <- choice$diagnostics
  completion_groups[[g]] <- list(n_tested_gene_sets = nrow(permutation_summary),
    n_plotted_gene_sets = n_sets, n_undefined_profiles = nrow(unplottable),
    n_windows = ncol(z),
    selected_clusters = k, color_limit_99_5_abs_z = color_limit,
    n_significant_gene_sets = length(significant_sets))
  report <- c(report, paste0("## ", group_labels[[g]]), "",
    paste0("- Tested gene sets: ", nrow(permutation_summary)),
    paste0("- Plotted defined profiles: ", n_sets),
    paste0("- Undefined profiles: ", nrow(unplottable)),
    paste0("- Selected clusters: ", k),
    paste0("- Gene sets with a final significant region (annotation only): ",
      length(significant_sets)),
    paste0("- Display color cap: |z| = ", signif(color_limit, 5)), "",
    "| Cluster | Sets | Significant sets | Centroid peak pt | Representative sets | Recurring keywords |",
    "|---|---:|---:|---:|---|---|",
    vapply(seq_len(nrow(summaries)), function(i) paste0("| ",
      summaries$cluster_id[[i]], " | ", summaries$n_gene_sets[[i]], " | ",
      summaries$n_significant_gene_sets[[i]], " | ",
      signif(summaries$centroid_peak_pt[[i]], 5), " | ",
      summaries$representative_gene_sets[[i]], " | ",
      summaries$recurring_name_keywords[[i]], " |"), character(1)), "")
  rm(curves, wide, z, hc, plot_data)
  invisible(gc())
}

fwrite(rbindlist(all_membership), file.path(stage_dir,
  "gene_set_cluster_membership.tsv.gz"), sep = "\t", quote = FALSE,
  na = "NA", compress = "gzip")
fwrite(rbindlist(all_centroids), file.path(stage_dir,
  "cluster_centroid_profiles.tsv.gz"), sep = "\t", quote = FALSE,
  na = "NA", compress = "gzip")
fwrite(rbindlist(all_summaries), file.path(stage_dir, "cluster_summary.tsv"),
  sep = "\t", quote = FALSE, na = "NA")
fwrite(rbindlist(all_k_diagnostics), file.path(stage_dir,
  "cluster_number_diagnostics.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(rbindlist(all_unplottable, use.names = TRUE, fill = TRUE),
  file.path(stage_dir, "undefined_z_score_profiles.tsv"), sep = "\t",
  quote = FALSE, na = "NA")
writeLines(report, file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", source_report_root = report_root,
  filtering = "none", clustering = list(distance = "1-Pearson correlation",
    linkage = "average", selection = "approximate silhouette"),
  groups = completion_groups), file.path(stage_dir, "completion.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)
file.rename(stage_dir, output_dir)
message("Complete: ", output_dir)
