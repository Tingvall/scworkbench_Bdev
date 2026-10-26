#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(5L, 6L)) stop(paste(
  "Usage: 23_cytokine_raw_ucell_original_z_passing_heatmap.R",
  "<standalone_comparison_dir> <ucell_permutation_dir>",
  "<row_membership.tsv> <bb.tsv> <output_dir> [raw|per_set_0_1]"))
comparison_dir <- normalizePath(args[[1L]], mustWork = TRUE)
permutation_dir <- normalizePath(args[[2L]], mustWork = TRUE)
membership_path <- normalizePath(args[[3L]], mustWork = TRUE)
bb_path <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
ucell_display <- if (length(args) == 6L) args[[6L]] else "raw"
stopifnot(ucell_display %in% c("raw", "per_set_0_1"))
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
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
significant_sets <- fread(file.path(permutation_dir,
  "significant_gene_sets.tsv"))$feature_set
permutation_completion <- read_json(file.path(permutation_dir,
  "completion.json"), simplifyVector = TRUE)
significant_regions <- fread(file.path(permutation_dir,
  "significant_positive_regions.tsv"))
membership <- fread(membership_path)
setorder(membership, row_position_top_to_bottom)
stopifnot(uniqueN(profiles$feature_set) == 38L,
  uniqueN(profiles$window_index) == 97L,
  nrow(profiles) == 38L * 97L,
  setequal(unique(profiles$feature_set), membership$feature_set),
  all(significant_sets %chin% membership$feature_set))

# Optional min-max scaling is strictly a visualization transform. It is
# calculated independently within each cytokine and never written back to or
# used by the permutation analysis.
profiles[, ucell_visual_0_1 := {
  span <- max(ucell_median) - min(ucell_median)
  if (is.finite(span) && span > .Machine$double.eps)
    (ucell_median - min(ucell_median)) / span else rep(0, .N)
}, by = feature_set]
profiles[, in_significant_positive_region := FALSE]
if (nrow(significant_regions)) for (i in seq_len(nrow(significant_regions))) {
  rr <- significant_regions[i]
  profiles[feature_set == rr$feature_set &
    window_center >= rr$region_left & window_center <= rr$region_right,
    in_significant_positive_region := TRUE]
}
profiles[, significant_region_ucell_0_1 := fifelse(
  in_significant_positive_region, ucell_visual_0_1, NA_real_)]

prior_order <- membership$feature_set
passing_order <- prior_order[prior_order %chin% significant_sets]
nonpassing_order <- prior_order[!prior_order %chin% significant_sets]
display_order <- c(passing_order, nonpassing_order)
is_dynamic <- identical(permutation_completion$analysis,
  "native-UCell trajectory-wide dynamic-range permutation")
status_levels <- c(if (is_dynamic)
  "Passing: dynamic UCell trajectory at BH q < 0.01" else
  "Passing: ≥1 UCell-permutation region at BH q < 0.01", "Not passing")
profiles[, status := factor(ifelse(feature_set %chin% significant_sets,
  status_levels[[1L]], status_levels[[2L]]), levels = status_levels)]
profiles[, feature_label := factor(short_label(feature_set),
  levels = rev(short_label(display_order)))]

bb <- fread(bb_path)
stopifnot(nrow(bb) == 2L,
  identical(as.character(bb$label), c("up1", "up2")))
bb[, boundary_label := c("BB1", "BB2")]

base_layers <- list(
  geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
    linewidth = .45, linetype = "dashed", color = "#374151"),
  geom_text(data = bb, aes(x = pt_median, y = Inf, label = boundary_label),
    inherit.aes = FALSE, vjust = 1.25, size = 3, fontface = "bold",
    color = "#374151"),
  facet_grid(status ~ ., scales = "free_y", space = "free_y"),
  scale_x_continuous(expand = c(0, 0)),
  theme_bw(base_size = 10),
  theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),
    strip.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"))
)

left_value <- if (ucell_display == "per_set_0_1")
  "ucell_visual_0_1" else "ucell_median"
ucell_limit <- if (ucell_display == "per_set_0_1") 1 else
  as.numeric(quantile(profiles$ucell_median, .995))
left <- ggplot(profiles, aes(window_center, feature_label,
    fill = .data[[left_value]])) +
  geom_raster() + base_layers +
  scale_fill_gradientn(colors = c("#FFFFFF", "#FEE8C8", "#FDBB84",
    "#E34A33", "#7F0000"), limits = c(0, ucell_limit),
    oob = scales::squish, name = if (ucell_display == "per_set_0_1")
      "UCell\n0–1" else "Median\nUCell") +
  labs(title = if (ucell_display == "per_set_0_1")
      "UCell enrichment scaled within cytokine" else "Raw UCell enrichment",
    subtitle = if (ucell_display == "per_set_0_1")
      "Each row min–max scaled to 0–1 for visualization only" else
      "Native within-cell rank score",
    x = "Existing Bdev pseudotime (pt)", y = "Cytokine gene set")

# The requested sequential display focuses on positive values from the original
# standardized-gene trajectory analysis. Zero and all negative z-scores map to
# white; no UCell-derived z-score is used anywhere in this figure.
z_limit <- as.numeric(quantile(profiles$trajectory_z[
  is.finite(profiles$trajectory_z) & profiles$trajectory_z > 0], .995))
right <- ggplot(profiles, aes(window_center, feature_label,
    fill = trajectory_z)) +
  geom_raster() + base_layers +
  scale_fill_gradientn(colors = c("#FFFFFF", "#FEE8C8", "#FDBB84",
    "#E34A33", "#7F0000"), limits = c(0, z_limit),
    oob = scales::squish, name = "Original\ntrajectory z") +
  scale_y_discrete(breaks = NULL, drop = TRUE) +
  labs(title = "Original standardized trajectory z-score",
    subtitle = "Positive values shown; z ≤ 0 is white",
    x = "Existing Bdev pseudotime (pt)", y = NULL)

regions_panel <- ggplot(profiles, aes(window_center, feature_label,
    fill = significant_region_ucell_0_1)) +
  geom_raster() + base_layers +
  scale_fill_gradientn(colors = c("#FFFFFF", "#FEE8C8", "#FDBB84",
    "#E34A33", "#7F0000"), limits = c(0, 1),
    na.value = "white", oob = scales::squish,
    name = "UCell\n0–1") +
  scale_y_discrete(breaks = NULL, drop = TRUE) +
  labs(title = "Permutation-supported positive regions",
    subtitle = "Scaled UCell shown only inside supported regions",
    x = "Existing Bdev pseudotime (pt)", y = NULL)

paired <- left + right + regions_panel +
  plot_layout(ncol = 3, widths = c(1.13, 1, 1)) +
  plot_annotation(
    title = paste0("Cytokines: ",
      if (ucell_display == "per_set_0_1") "0–1-scaled UCell" else "raw UCell scores",
      " beside the original trajectory z-score"),
    subtitle = paste(
      "Rows are split by the independent UCell permutation call and kept in the same order.",
      "No UCell z-score is displayed."),
    theme = theme(plot.title = element_text(face = "bold", size = 15)))

file_stem <- if (ucell_display == "per_set_0_1")
  "cytokine_ucell_0_1_vs_original_z_passing_split" else
  "cytokine_raw_ucell_vs_original_z_passing_split"
ggsave(file.path(stage_dir, paste0(file_stem, ".png")), paired,
  width = 25, height = 11, dpi = 240, bg = "white")
ggsave(file.path(stage_dir, paste0(file_stem, ".pdf")), paired,
  width = 25, height = 11, device = cairo_pdf, bg = "white")

row_order <- data.table(feature_set = display_order,
  short_name = short_label(display_order),
  status = ifelse(display_order %chin% significant_sets, "passing",
    "not_passing"),
  row_position_within_status = c(seq_along(passing_order),
    seq_along(nonpassing_order)),
  prior_ucell_row_position = membership$row_position_top_to_bottom[
    match(display_order, membership$feature_set)])
fwrite(row_order, file.path(stage_dir, "heatmap_row_order.tsv"),
  sep = "\t", quote = FALSE)
fwrite(profiles[, .(feature_set, status = as.character(status), window_index,
  window_center, raw_ucell_median = ucell_median,
  ucell_visual_0_1,
  in_significant_positive_region, significant_region_ucell_0_1,
  original_trajectory_z = trajectory_z)],
  file.path(stage_dir, "heatmap_values.tsv.gz"), sep = "\t", quote = FALSE,
  compress = "gzip")

writeLines(c("# UCell versus original trajectory-z heatmap", "",
  paste0("- Passing cytokine sets: ", length(passing_order)),
  paste0("- Non-passing cytokine sets: ", length(nonpassing_order)),
  if (ucell_display == "per_set_0_1")
    "- Left: each cytokine's window-median UCell profile independently min-max scaled to 0–1 for visualization only."
  else "- Left: native median UCell score in each pseudotime window.",
  "- Right: original standardized-gene trajectory z-score from the previous analysis.",
  "- Third panel: 0–1-scaled UCell values shown only in permutation-supported positive regions; all other windows are white.",
  "- No UCell-derived z-score is used.",
  "- All three panels use sequential white-to-red scales, calibrated separately where their units differ.",
  "- Permutation inference uses native, unscaled UCell scores.",
  "- On the original-z panel, z <= 0 is white.",
  "- Rows retain the prior UCell heatmap order within passing status."),
  file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", gene_sets = 38L,
  passing_sets = length(passing_order),
  nonpassing_sets = length(nonpassing_order), windows = 97L,
  left_value = if (ucell_display == "per_set_0_1")
    "per-cytokine min-max-scaled median UCell score (visualization only)" else
    "native median UCell score",
  right_value = "original standardized-gene trajectory z-score",
  ucell_z_score_used = FALSE,
  ucell_rescaled_for_inference = FALSE,
  passing_definition = if (is_dynamic)
    "trajectory-wide UCell range at BH q < 0.01 across 38 cytokines" else
    "at least one UCell permutation region at BH q < 0.01",
  palettes = "separate sequential white-to-red",
  nonpositive_original_z_display = "white",
  ucell_display_cap = ucell_limit, original_z_display_cap = z_limit,
  row_order = "prior UCell heatmap order within passing status",
  source_comparison = comparison_dir,
  source_ucell_permutation_report = permutation_dir),
  file.path(stage_dir, "completion.json"), pretty = TRUE,
  auto_unbox = TRUE, digits = NA)

stopifnot(file.rename(stage_dir, output_dir))
message("Completed paired raw-UCell/original-z heatmap at: ", output_dir)
