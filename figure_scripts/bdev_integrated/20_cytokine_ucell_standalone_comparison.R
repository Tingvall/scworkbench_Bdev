#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 200)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 20_cytokine_ucell_standalone_comparison.R",
  "<native_ucell_report> <z_report_root> <bb.tsv> <output_dir>"))
ucell_root <- normalizePath(args[[1L]], mustWork = TRUE)
z_root <- normalizePath(args[[2L]], mustWork = TRUE)
bb_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
row_standardize <- function(x) {
  centers <- rowMeans(x); scales <- apply(x, 1L, sd)
  stopifnot(all(is.finite(scales)))
  variable <- scales > .Machine$double.eps
  out <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
  out[variable, ] <- sweep(sweep(x[variable, , drop = FALSE], 1L,
    centers[variable], "-"), 1L, scales[variable], "/")
  attr(out, "variable_profile") <- variable
  out
}
correlation_order <- function(x) {
  if (nrow(x) <= 1L) return(seq_len(nrow(x)))
  similarity <- cor(t(x)); similarity[!is.finite(similarity)] <- 0
  diag(similarity) <- 1
  distance <- 1 - similarity; distance[distance < 0] <- 0
  hclust(as.dist(distance), method = "average")$order
}
short_label <- function(x) sub(
  "^parse_cytokine_rest_only_refinement_all_top100_v1\\.", "", x)

ucell_completion <- read_json(file.path(ucell_root, "completion.json"),
  simplifyVector = TRUE)
stopifnot(identical(ucell_completion$status, "complete"),
  ucell_completion$passing_signatures == 38L,
  ucell_completion$max_rank == 1500L,
  ucell_completion$defined_pt_cells == 12642L,
  ucell_completion$missing_pt_cells == 14081L)
bb <- fread(bb_path)
stopifnot(nrow(bb) == 2L,
  identical(as.character(bb$label), c("up1", "up2")))
bb[, boundary_label := c("BB1", "BB2")]

ucell <- fread(file.path(ucell_root, "sliding_window_summary.tsv"))
ucell <- ucell[abs(width_fraction - 0.05) < 1e-12 & eligible_window == TRUE,
  .(feature_set, window_index, window_center = pseudotime,
    ucell_median = score_median, ucell_q25 = score_q25,
    ucell_q75 = score_q75, n_cells)]
z_curves <- read_tsv_gz(file.path(z_root, "cytokines", "window_curves.tsv.gz"))[
  , .(feature_set, window_index, z_window_center = window_center,
    trajectory_z = z_score)]
stopifnot(uniqueN(ucell$feature_set) == 38L,
  uniqueN(z_curves$feature_set) == 38L,
  setequal(ucell$feature_set, z_curves$feature_set))
comparison_long <- merge(ucell, z_curves,
  by = c("feature_set", "window_index"), all = TRUE, sort = FALSE)
stopifnot(nrow(comparison_long) == 38L * 97L, !anyNA(comparison_long),
  max(abs(comparison_long$window_center -
    comparison_long$z_window_center)) < 1e-10)
setorder(comparison_long, feature_set, window_index)

ucell_wide <- dcast(comparison_long,
  feature_set ~ window_index, value.var = "ucell_median")
sets <- ucell_wide$feature_set
u <- as.matrix(ucell_wide[, -"feature_set"]); rownames(u) <- sets
shape <- row_standardize(u)
variable_profile <- attr(shape, "variable_profile")
set.seed(314179L)
k <- 6L
fit <- kmeans(shape, centers = k, nstart = 100L, iter.max = 200L,
  algorithm = "Hartigan-Wong")
centroid_order <- correlation_order(fit$centers)
raw_to_display <- setNames(sprintf("C%02d", seq_len(k)), centroid_order)
cluster_id <- unname(raw_to_display[as.character(fit$cluster)])
ordered_sets <- unlist(lapply(centroid_order, function(raw_id) {
  idx <- which(fit$cluster == raw_id)
  sets[idx][correlation_order(shape[idx, , drop = FALSE])]
}), use.names = FALSE)
membership <- data.table(feature_set = sets, cluster_id = cluster_id,
  row_position_top_to_bottom = match(sets, ordered_sets),
  variable_ucell_trajectory = variable_profile)
setorder(membership, row_position_top_to_bottom)

mapping <- membership[, .(feature_set, row_position_top_to_bottom)]
heat <- merge(ucell, mapping, by = "feature_set", sort = FALSE)
heat[, y := length(ordered_sets) - row_position_top_to_bottom + 1L]
ordered_clusters <- membership[match(ordered_sets, feature_set), cluster_id]
runs <- rle(ordered_clusters); boundaries <- length(ordered_sets) -
  head(cumsum(runs$lengths), -1L) + 0.5
upper <- as.numeric(quantile(heat$ucell_median, .995))
p <- ggplot(heat, aes(window_center, y, fill = ucell_median)) +
  geom_raster() +
  geom_hline(yintercept = boundaries, linewidth = .25, color = "#202020") +
  geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
    linewidth = .5, linetype = "dashed", color = "#374151") +
  geom_text(data = bb, aes(x = pt_median, y = Inf, label = boundary_label),
    inherit.aes = FALSE, vjust = 1.25, size = 3, fontface = "bold",
    color = "#374151") +
  scale_fill_gradientn(colors = c("#F7F7F7", "#FEE8C8", "#FDBB84",
    "#E34A33", "#7F0000"), limits = c(0, upper),
    oob = scales::squish, name = "Median\nUCell") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(ordered_sets),
    labels = rev(short_label(ordered_sets)), expand = c(0, 0)) +
  labs(title = "Cytokines: standalone raw-UCell trajectory profiles",
    subtitle = "All 38 sets; six k-means shape clusters; no z-score or combined filter",
    x = "Existing Bdev pseudotime (pt)", y = "Cytokine gene set") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold"))
ggsave(file.path(stage_dir, "cytokine_raw_ucell_heatmap.png"), p,
  width = 10, height = 10, dpi = 240, bg = "white")
ggsave(file.path(stage_dir, "cytokine_raw_ucell_heatmap.pdf"), p,
  width = 10, height = 10, device = cairo_pdf, bg = "white")

set_summary <- comparison_long[, {
  ui <- which.max(ucell_median); zi <- which.max(trajectory_z)
  .(peak_ucell = ucell_median[[ui]], peak_ucell_pt = window_center[[ui]],
    median_ucell_across_pt = median(ucell_median),
    ucell_peak_minus_median = ucell_median[[ui]] - median(ucell_median),
    maximum_trajectory_z = trajectory_z[[zi]],
    maximum_z_pt = window_center[[zi]],
    ucell_at_maximum_z = ucell_median[[zi]],
    peak_pt_difference = window_center[[ui]] - window_center[[zi]])
}, by = feature_set]
set_summary[, `:=`(
  peak_ucell_rank = frank(-peak_ucell, ties.method = "min"),
  maximum_z_rank = frank(-maximum_trajectory_z, ties.method = "min"),
  ucell_effect_rank = frank(-ucell_peak_minus_median, ties.method = "min"))]
set_summary[membership, on = "feature_set", `:=`(
  ucell_cluster_id = i.cluster_id,
  ucell_row_position = i.row_position_top_to_bottom)]
setorder(set_summary, peak_ucell_rank, maximum_z_rank)
fwrite(set_summary, file.path(stage_dir, "cytokine_ucell_vs_z_summary.tsv"),
  sep = "\t", quote = FALSE, na = "NA")
fwrite(comparison_long, file.path(stage_dir,
  "cytokine_ucell_and_z_window_profiles.tsv.gz"), sep = "\t", quote = FALSE,
  na = "NA", compress = "gzip")
fwrite(membership, file.path(stage_dir, "cytokine_ucell_kmeans_membership.tsv"),
  sep = "\t", quote = FALSE, na = "NA")

spearman_peak <- cor(set_summary$peak_ucell,
  set_summary$maximum_trajectory_z, method = "spearman")
spearman_effect <- cor(set_summary$ucell_peak_minus_median,
  set_summary$maximum_trajectory_z, method = "spearman")
scatter <- ggplot(set_summary,
    aes(maximum_trajectory_z, peak_ucell, label = short_label(feature_set))) +
  geom_point(size = 2.2, color = "#7F0000", alpha = .8) +
  ggrepel::geom_text_repel(size = 2.6, max.overlaps = 12, seed = 314179L) +
  labs(title = "Standalone score comparison across cytokine gene sets",
    subtitle = paste0("Spearman rho = ", signif(spearman_peak, 3),
      "; no combined score or joint decision rule"),
    x = "Maximum standardized trajectory z-score",
    y = "Maximum median UCell score") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave(file.path(stage_dir, "cytokine_peak_ucell_vs_maximum_z.png"), scatter,
  width = 9, height = 7, dpi = 220, bg = "white")
ggsave(file.path(stage_dir, "cytokine_peak_ucell_vs_maximum_z.pdf"), scatter,
  width = 9, height = 7, device = cairo_pdf, bg = "white")

top_z <- set_summary[order(-maximum_trajectory_z)][1:10]
top_ucell <- set_summary[order(-peak_ucell)][1:10]
fmt_rows <- function(x) vapply(seq_len(nrow(x)), function(i) paste0("| ",
  short_label(x$feature_set[[i]]), " | ",
  signif(x$maximum_trajectory_z[[i]], 4), " | ",
  signif(x$peak_ucell[[i]], 4), " | ",
  signif(x$ucell_peak_minus_median[[i]], 4), " | "), character(1))
report <- c("# Standalone cytokine UCell trajectory comparison", "",
  paste0("This is a separate raw-UCell analysis and a descriptive comparison. ",
    "No UCell/z combined score, joint cutoff, or combined significance call ",
    "was created."), "",
  "## Inputs and method", "",
  "- Native scWorkbench UCell scores on RNA counts for all 26,723 Bdev cells.",
  "- UCell `max_rank=1500`, no graph smoothing.",
  "- Exact 12,642 cells with existing `pt`; missing `pt` values were preserved.",
  "- Median raw UCell score in the same 5%-wide, 1%-step, minimum-20-cell windows.",
  "- All 38 current cytokine memberships matched the native UCell run exactly (934 gene memberships).",
  "- Heatmap clustering uses k-means on row-standardized UCell trajectory shapes; displayed values remain raw UCell medians.", "",
  paste0("- Constant UCell trajectories retained as all-zero clustering shapes: ",
    sum(!variable_profile), if (any(!variable_profile)) paste0(" (",
      paste(short_label(sets[!variable_profile]), collapse = "; "), ")") else "."), "",
  "## Descriptive comparison", "",
  paste0("- Spearman correlation of maximum raw UCell with maximum trajectory z: ",
    signif(spearman_peak, 4), "."),
  paste0("- Spearman correlation of UCell peak-minus-trajectory-median with maximum trajectory z: ",
    signif(spearman_effect, 4), "."), "",
  "### Top 10 by standardized trajectory z", "",
  "| Gene set | Maximum z | Peak UCell | UCell peak − median |",
  "|---|---:|---:|---:|", fmt_rows(top_z), "",
  "### Top 10 by raw UCell", "",
  "| Gene set | Maximum z | Peak UCell | UCell peak − median |",
  "|---|---:|---:|---:|", fmt_rows(top_ucell), "")
writeLines(report, file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", analysis = "standalone UCell",
  combined_score_or_filter = FALSE, gene_sets = 38L, gene_memberships = 934L,
  defined_pt_cells = 12642L, windows = 97L, kmeans_clusters = 6L,
  constant_ucell_trajectories = sum(!variable_profile),
  spearman_peak_ucell_vs_maximum_z = spearman_peak,
  spearman_ucell_effect_vs_maximum_z = spearman_effect,
  source_ucell_report = ucell_root, source_z_report = z_root),
  file.path(stage_dir, "completion.json"), pretty = TRUE, auto_unbox = TRUE,
  digits = NA)
stopifnot(file.rename(stage_dir, output_dir))
message("Complete: ", output_dir)
