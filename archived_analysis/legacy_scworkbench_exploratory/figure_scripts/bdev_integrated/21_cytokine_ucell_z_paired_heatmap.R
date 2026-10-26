#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(paste(
  "Usage: 21_cytokine_ucell_z_paired_heatmap.R",
  "<comparison_dir> <bb.tsv>"))
comparison_dir <- normalizePath(args[[1L]], mustWork = TRUE)
bb_path <- normalizePath(args[[2L]], mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
short_label <- function(x) sub(
  "^parse_cytokine_rest_only_refinement_all_top100_v1\\.", "", x)

profiles <- read_tsv_gz(file.path(comparison_dir,
  "cytokine_ucell_and_z_window_profiles.tsv.gz"))
membership <- fread(file.path(comparison_dir,
  "cytokine_ucell_kmeans_membership.tsv"))
bb <- fread(bb_path)
stopifnot(nrow(bb) == 2L,
  identical(as.character(bb$label), c("up1", "up2")),
  uniqueN(profiles$feature_set) == 38L,
  uniqueN(profiles$window_index) == 97L,
  nrow(profiles) == 38L * 97L)
bb[, boundary_label := c("BB1", "BB2")]
setorder(membership, row_position_top_to_bottom)
ordered_sets <- membership$feature_set
profiles[membership[, .(feature_set, row_position_top_to_bottom)],
  on = "feature_set", row_position_top_to_bottom :=
    i.row_position_top_to_bottom]
stopifnot(!anyNA(profiles$row_position_top_to_bottom))
profiles[, y := length(ordered_sets) - row_position_top_to_bottom + 1L]
ordered_clusters <- membership$cluster_id
runs <- rle(ordered_clusters)
boundaries <- length(ordered_sets) - head(cumsum(runs$lengths), -1L) + .5
labels <- rev(short_label(ordered_sets))

common_layers <- list(
  geom_hline(yintercept = boundaries, linewidth = .24, color = "#202020"),
  geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
    linewidth = .5, linetype = "dashed", color = "#374151"),
  geom_text(data = bb, aes(x = pt_median, y = Inf, label = boundary_label),
    inherit.aes = FALSE, vjust = 1.25, size = 3, fontface = "bold",
    color = "#374151"),
  scale_x_continuous(expand = c(0, 0)),
  theme_bw(base_size = 10),
  theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold"))
)

ucell_limit <- as.numeric(quantile(profiles$ucell_median, .995))
left <- ggplot(profiles, aes(window_center, y, fill = ucell_median)) +
  geom_raster() + common_layers +
  scale_fill_gradientn(colors = c("#F7F7F7", "#FEE8C8", "#FDBB84",
    "#E34A33", "#7F0000"), limits = c(0, ucell_limit),
    oob = scales::squish, name = "Median\nUCell") +
  scale_y_continuous(breaks = seq_along(ordered_sets), labels = labels,
    expand = c(0, 0)) +
  labs(title = "Raw UCell enrichment",
    subtitle = "Within-cell expression ranks",
    x = "Existing Bdev pseudotime (pt)", y = "Cytokine gene set")

z_limit <- as.numeric(quantile(abs(profiles$trajectory_z), .995))
right <- ggplot(profiles, aes(window_center, y, fill = trajectory_z)) +
  geom_raster() + common_layers +
  scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
    midpoint = 0, limits = c(-z_limit, z_limit),
    oob = scales::squish, name = "Trajectory\nz-score") +
  scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
  labs(title = "Standardized trajectory signal",
    subtitle = "Dataset-relative deviation",
    x = "Existing Bdev pseudotime (pt)", y = NULL)

paired <- left + right + plot_layout(ncol = 2, widths = c(1.18, 1)) +
  plot_annotation(
    title = "Cytokine trajectory heatmap comparison",
    subtitle = paste(
      "Same 38 gene sets, UCell-derived row order, six shape clusters, and 97 windows.",
      "Scores are displayed separately and are not combined."),
    theme = theme(plot.title = element_text(face = "bold", size = 15)))
ggsave(file.path(comparison_dir,
  "cytokine_raw_ucell_vs_trajectory_z_same_order.png"), paired,
  width = 17, height = 10, dpi = 240, bg = "white")
ggsave(file.path(comparison_dir,
  "cytokine_raw_ucell_vs_trajectory_z_same_order.pdf"), paired,
  width = 17, height = 10, device = cairo_pdf, bg = "white")

manifest <- data.table(
  file = c("cytokine_raw_ucell_vs_trajectory_z_same_order.png",
    "cytokine_raw_ucell_vs_trajectory_z_same_order.pdf"),
  format = c("PNG", "PDF"), gene_sets = 38L, windows = 97L,
  row_order = "UCell k-means order",
  ucell_display_cap = ucell_limit, z_display_cap = z_limit,
  combined_score = FALSE)
fwrite(manifest, file.path(comparison_dir, "paired_heatmap_manifest.tsv"),
  sep = "\t", quote = FALSE)
message("Completed paired UCell/z heatmap in: ", comparison_dir)
