#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 200)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 19_full_unfiltered_trajectory_z_heatmap_kmeans.R",
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
requested_k <- c(cytokines = 6L, bp = 20L, mf = 12L, cc = 10L)
bb <- fread(bb_path)
stopifnot(nrow(bb) == 2L,
  identical(as.character(bb$label), c("up1", "up2")),
  all(is.finite(bb$pt_median)))
bb[, boundary_label := c("BB1", "BB2")]

read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
row_standardize <- function(x) {
  centers <- rowMeans(x)
  scales <- apply(x, 1L, sd)
  stopifnot(all(is.finite(scales)), all(scales > .Machine$double.eps))
  sweep(sweep(x, 1L, centers, "-"), 1L, scales, "/")
}
correlation_hclust_order <- function(x) {
  if (nrow(x) <= 1L) return(seq_len(nrow(x)))
  similarity <- cor(t(x), use = "pairwise.complete.obs")
  similarity[!is.finite(similarity)] <- 0
  diag(similarity) <- 1
  distance <- 1 - similarity
  distance[distance < 0] <- 0
  hclust(as.dist(distance), method = "average")$order
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

all_membership <- list(); all_centroids <- list(); all_summaries <- list()
all_diagnostics <- list(); all_unplottable <- list(); completion_groups <- list()
report <- c("# Complete unfiltered trajectory z-score heatmaps: k-means", "",
  paste0("Every gene set with a defined complete z-score profile is included. ",
    "No significance, peak-z, GO-hierarchy, or redundancy filter was applied. ",
    "K-means used each gene set's row-standardized 97-window profile, so ",
    "clustering emphasizes trajectory shape rather than absolute z magnitude. ",
    "The heatmaps display the original unscaled z-scores. Cluster counts were ",
    "chosen for interpretable resolution: 20 BP, 12 MF, 10 CC, and 6 cytokine ",
    "clusters. Only the displayed colors are capped at the 99.5th percentile ",
    "of absolute z-score."), "")

for (group_index in seq_along(groups)) {
  g <- groups[[group_index]]; k <- requested_k[[g]]
  message("Processing ", g, " with k=", k)
  curves <- read_tsv_gz(file.path(report_root, g, "window_curves.tsv.gz"))
  stopifnot(all(is.finite(curves$z_score)))
  setorder(curves, feature_set, window_index)
  window_map <- unique(curves[, .(window_index, window_center)])
  setorder(window_map, window_index)
  wide <- dcast(curves, feature_set ~ window_index, value.var = "z_score")
  feature_sets <- wide$feature_set
  z <- as.matrix(wide[, -"feature_set"]); rownames(z) <- feature_sets
  stopifnot(ncol(z) == nrow(window_map), !anyNA(z), nrow(z) > k)
  permutation_summary <- fread(file.path(report_root, g,
    "permutation_summary.tsv"))
  unplottable <- permutation_summary[!feature_set %chin% feature_sets]
  stopifnot(nrow(permutation_summary) == nrow(z) + nrow(unplottable))
  all_unplottable[[g]] <- unplottable

  shape <- row_standardize(z)
  set.seed(314159L + group_index)
  fit <- kmeans(shape, centers = k, iter.max = 200L, nstart = 100L,
    algorithm = "Hartigan-Wong")
  raw_cluster <- fit$cluster
  # Order clusters by their shape centroids, then rows within each cluster by
  # correlation hierarchy. This affects display order, not k-means assignment.
  centroid_order <- correlation_hclust_order(fit$centers)
  raw_to_display <- setNames(sprintf("C%02d", seq_len(k)), centroid_order)
  cluster_id <- unname(raw_to_display[as.character(raw_cluster)])
  ordered_sets <- unlist(lapply(centroid_order, function(raw_id) {
    idx <- which(raw_cluster == raw_id)
    feature_sets[idx][correlation_hclust_order(shape[idx, , drop = FALSE])]
  }), use.names = FALSE)

  regions <- fread(file.path(report_root, g, "significant_positive_regions.tsv"))
  significant_sets <- unique(regions$feature_set)
  membership <- data.table(analysis_group = g, feature_set = feature_sets,
    cluster_id = cluster_id,
    row_position_top_to_bottom = match(feature_sets, ordered_sets),
    maximum_window_z = apply(z, 1L, max),
    peak_window_pt = window_map$window_center[max.col(z, ties.method = "first")],
    has_final_significant_region = feature_sets %chin% significant_sets)
  setorder(membership, row_position_top_to_bottom)

  centroids <- rbindlist(lapply(sprintf("C%02d", seq_len(k)), function(cid) {
    idx <- which(cluster_id == cid)
    data.table(analysis_group = g, cluster_id = cid,
      window_index = window_map$window_index,
      window_center = window_map$window_center,
      mean_original_z = colMeans(z[idx, , drop = FALSE]),
      mean_shape_z = colMeans(shape[idx, , drop = FALSE]))
  }))
  summaries <- rbindlist(lapply(sprintf("C%02d", seq_len(k)), function(cid) {
    idx <- which(cluster_id == cid)
    centroid <- colMeans(shape[idx, , drop = FALSE])
    cors <- as.numeric(shape[idx, , drop = FALSE] %*% centroid /
      ((ncol(shape) - 1L) * sd(centroid)))
    ranking <- order(-cors, -apply(z[idx, , drop = FALSE], 1L, max))
    member_sets <- feature_sets[idx]
    data.table(analysis_group = g, cluster_id = cid,
      n_gene_sets = length(idx),
      fraction_of_group = length(idx) / nrow(z),
      n_significant_gene_sets = sum(member_sets %chin% significant_sets),
      centroid_peak_pt = window_map$window_center[[which.max(centroid)]],
      representative_gene_sets = paste(member_sets[head(ranking, 5L)],
        collapse = "; "), recurring_name_keywords = common_keywords(member_sets))
  }))
  diagnostics <- data.table(analysis_group = g, k = k,
    total_withinss = fit$tot.withinss, totalss = fit$totss,
    between_ss_fraction = fit$betweenss / fit$totss,
    smallest_cluster = min(fit$size), median_cluster = median(fit$size),
    largest_cluster = max(fit$size),
    largest_cluster_fraction = max(fit$size) / nrow(z),
    iterations = fit$iter, ifault = fit$ifault)

  mapping <- data.table(feature_set = ordered_sets,
    row_position_top_to_bottom = seq_along(ordered_sets))
  plot_data <- merge(curves[, .(feature_set, window_center, value = z_score)],
    mapping, by = "feature_set", sort = FALSE)
  n_sets <- length(ordered_sets)
  plot_data[, y := n_sets - row_position_top_to_bottom + 1L]
  ordered_cluster <- membership[match(ordered_sets, feature_set), cluster_id]
  runs <- rle(ordered_cluster); run_end <- cumsum(runs$lengths)
  run_start <- c(1L, head(run_end, -1L) + 1L)
  y_breaks <- n_sets - (run_start + run_end) / 2 + 1L
  boundary_y <- n_sets - head(run_end, -1L) + 0.5
  color_limit <- as.numeric(quantile(abs(z), .995))
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
      ": complete unfiltered z-score profiles — k-means"),
      subtitle = paste0(n_sets, " defined profiles; k = ", k,
        "; clustering uses row-standardized shapes"),
      x = "Existing Bdev pseudotime (pt)",
      y = if (g == "cytokines") "Cytokine gene set" else "K-means cluster") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold"))
  ggsave(file.path(stage_dir, paste0(g, "_full_unfiltered_z_kmeans.png")),
    p, width = 10, height = 10, dpi = 240, bg = "white")
  ggsave(file.path(stage_dir, paste0(g, "_full_unfiltered_z_kmeans.pdf")),
    p, width = 10, height = 10, device = cairo_pdf, bg = "white")

  centroid_plot_data <- copy(centroids)
  centroid_plot_data[, cluster_id := factor(cluster_id,
    levels = rev(sprintf("C%02d", seq_len(k))))]
  shape_limit <- max(abs(centroid_plot_data$mean_shape_z))
  cp <- ggplot(centroid_plot_data,
    aes(window_center, cluster_id, fill = mean_shape_z)) +
    geom_tile(height = 1) +
    geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
      linewidth = .5, linetype = "dashed", color = "#374151") +
    geom_text(data = bb, aes(x = pt_median, y = Inf, label = boundary_label),
      inherit.aes = FALSE, vjust = 1.25, size = 3, fontface = "bold",
      color = "#374151") +
    scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
      midpoint = 0, limits = c(-shape_limit, shape_limit), name = "Mean shape") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = paste0(group_labels[[g]], ": k-means cluster shapes"),
      x = "Existing Bdev pseudotime (pt)", y = "K-means cluster") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))
  ggsave(file.path(stage_dir, paste0(g, "_kmeans_cluster_shapes.png")),
    cp, width = 9, height = 6.5, dpi = 240, bg = "white")
  ggsave(file.path(stage_dir, paste0(g, "_kmeans_cluster_shapes.pdf")),
    cp, width = 9, height = 6.5, device = cairo_pdf, bg = "white")

  all_membership[[g]] <- membership; all_centroids[[g]] <- centroids
  all_summaries[[g]] <- summaries; all_diagnostics[[g]] <- diagnostics
  completion_groups[[g]] <- list(n_tested_gene_sets = nrow(permutation_summary),
    n_plotted_gene_sets = nrow(z), n_undefined_profiles = nrow(unplottable),
    n_windows = ncol(z), k = k, smallest_cluster = min(fit$size),
    median_cluster = median(fit$size), largest_cluster = max(fit$size),
    largest_cluster_fraction = max(fit$size) / nrow(z),
    between_ss_fraction = fit$betweenss / fit$totss)
  report <- c(report, paste0("## ", group_labels[[g]]), "",
    paste0("- Defined profiles: ", nrow(z), " / ", nrow(permutation_summary)),
    paste0("- k: ", k),
    paste0("- Cluster sizes: ", min(fit$size), "–", max(fit$size),
      "; median ", median(fit$size)),
    paste0("- Largest cluster: ", round(100 * max(fit$size) / nrow(z), 1), "%"),
    paste0("- Between-cluster sum-of-squares fraction: ",
      signif(fit$betweenss / fit$totss, 4)), "",
    "| Cluster | Sets | Fraction | Significant sets | Centroid peak pt | Representative sets |",
    "|---|---:|---:|---:|---:|---|",
    vapply(seq_len(nrow(summaries)), function(i) paste0("| ",
      summaries$cluster_id[[i]], " | ", summaries$n_gene_sets[[i]], " | ",
      round(100 * summaries$fraction_of_group[[i]], 1), "% | ",
      summaries$n_significant_gene_sets[[i]], " | ",
      signif(summaries$centroid_peak_pt[[i]], 5), " | ",
      summaries$representative_gene_sets[[i]], " |"), character(1)), "")
  rm(curves, wide, z, shape, fit, plot_data); invisible(gc())
}

fwrite(rbindlist(all_membership), file.path(stage_dir,
  "gene_set_kmeans_membership.tsv.gz"), sep = "\t", quote = FALSE,
  na = "NA", compress = "gzip")
fwrite(rbindlist(all_centroids), file.path(stage_dir,
  "kmeans_centroid_profiles.tsv.gz"), sep = "\t", quote = FALSE,
  na = "NA", compress = "gzip")
fwrite(rbindlist(all_summaries), file.path(stage_dir, "cluster_summary.tsv"),
  sep = "\t", quote = FALSE, na = "NA")
fwrite(rbindlist(all_diagnostics), file.path(stage_dir,
  "kmeans_diagnostics.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(rbindlist(all_unplottable, use.names = TRUE, fill = TRUE),
  file.path(stage_dir, "undefined_z_score_profiles.tsv"), sep = "\t",
  quote = FALSE, na = "NA")
writeLines(report, file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", source_report_root = report_root,
  filtering = "none", clustering = list(method = "kmeans",
    input = "row-standardized complete 97-window z-score profiles",
    algorithm = "Hartigan-Wong", nstart = 100L, iter_max = 200L,
    requested_k = as.list(requested_k)), groups = completion_groups),
  file.path(stage_dir, "completion.json"), pretty = TRUE, auto_unbox = TRUE,
  digits = NA)
stopifnot(file.rename(stage_dir, output_dir))
message("Complete: ", output_dir)
