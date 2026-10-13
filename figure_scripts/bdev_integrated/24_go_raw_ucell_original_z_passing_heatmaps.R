#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(5L, 6L)) stop(paste(
  "Usage: 24_go_raw_ucell_original_z_passing_heatmaps.R",
  "<go_ucell_report> <original_z_report> <z_kmeans_membership.tsv.gz>",
  "<bb.tsv> <output_dir> [raw|per_set_0_1]"))
ucell_root <- normalizePath(args[[1L]], mustWork = TRUE)
z_root <- normalizePath(args[[2L]], mustWork = TRUE)
order_path <- normalizePath(args[[3L]], mustWork = TRUE)
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
completion <- read_json(file.path(ucell_root, "completion.json"),
  simplifyVector = TRUE)
stopifnot(identical(completion$status, "complete"),
  completion$permutations == 1000L,
  isTRUE(completion$groupwise_multiple_testing))
is_residual <- isTRUE(completion$complexity_residualization_enabled)
if (is_residual && ucell_display != "per_set_0_1")
  stop("Residualized UCell heatmaps require per_set_0_1 display")
z_order_all <- read_tsv_gz(order_path)
bb <- fread(bb_path)
stopifnot(nrow(bb) == 2L,
  identical(as.character(bb$label), c("up1", "up2")))
bb[, boundary_label := c("BB1", "BB2")]

groups <- c("bp", "mf", "cc")
expected <- c(bp = 4087L, mf = 826L, cc = 511L)
labels <- c(bp = "GO Biological Process", mf = "GO Molecular Function",
  cc = "GO Cellular Component")
summaries <- list()

for (g in groups) {
  group_completion <- read_json(file.path(ucell_root, g, "completion.json"),
    simplifyVector = TRUE)
  ucell <- read_tsv_gz(file.path(ucell_root, g, "window_curves.tsv.gz"))[
    , .(feature_set, window_index, window_center,
      ucell_window_median = ucell_median)]
  z <- read_tsv_gz(file.path(z_root, g, "window_curves.tsv.gz"))[
    , .(feature_set, window_index, z_window_center = window_center,
      original_trajectory_z = z_score)]
  d <- merge(ucell, z, by = c("feature_set", "window_index"),
    all.x = TRUE, sort = FALSE)
  stopifnot(uniqueN(d$feature_set) == expected[[g]],
    uniqueN(d$window_index) == 97L,
    nrow(d) == expected[[g]] * 97L,
    all(is.finite(d$ucell_window_median)),
    max(abs(d$window_center - d$z_window_center), na.rm = TRUE) < 1e-10)
  d[, ucell_visual_0_1 := {
    span <- max(ucell_window_median) - min(ucell_window_median)
    if (is.finite(span) && span > .Machine$double.eps)
      (ucell_window_median - min(ucell_window_median)) / span else rep(0, .N)
  }, by = feature_set]
  significant_regions <- fread(file.path(ucell_root, g,
    "significant_positive_regions.tsv"))
  d[, in_significant_positive_region := FALSE]
  if (nrow(significant_regions)) for (i in seq_len(nrow(significant_regions))) {
    rr <- significant_regions[i]
    d[feature_set == rr$feature_set &
      window_center >= rr$region_left & window_center <= rr$region_right,
      in_significant_positive_region := TRUE]
  }
  d[, significant_region_ucell_0_1 := fifelse(
    in_significant_positive_region, ucell_visual_0_1, NA_real_)]

  passing_n <- as.integer(group_completion$passing_gene_sets)
  if (passing_n > 0L) {
    passing_sets <- fread(file.path(ucell_root, g,
      "significant_gene_sets.tsv"))$feature_set
  } else passing_sets <- character()
  stopifnot(length(passing_sets) == passing_n)

  z_order <- z_order_all[analysis_group == g][order(row_position_top_to_bottom)]
  defined_order <- z_order$feature_set
  undefined_sets <- sort(setdiff(unique(d$feature_set), defined_order))
  stopifnot(length(defined_order) + length(undefined_sets) == expected[[g]])
  complete_order <- c(defined_order, undefined_sets)
  passing_order <- complete_order[complete_order %chin% passing_sets]
  nonpassing_order <- complete_order[!complete_order %chin% passing_sets]
  display_order <- c(passing_order, nonpassing_order)
  status_levels <- character()
  if (length(passing_order)) status_levels <- c(status_levels,
    paste0("Passing: ", length(passing_order), " sets at BH q < 0.01"))
  status_levels <- c(status_levels,
    paste0("Not passing: ", length(nonpassing_order), " sets"))
  passing_label <- if (length(passing_order)) status_levels[[1L]] else NA_character_
  nonpassing_label <- tail(status_levels, 1L)
  d[, status := factor(ifelse(feature_set %chin% passing_sets,
    passing_label, nonpassing_label), levels = status_levels)]
  # Gene-set labels are intentionally suppressed for GO; this factor is the
  # exact original-z k-means order within each passing-status section.
  d[, row_id := factor(feature_set, levels = rev(display_order))]

  base_layers <- list(
    geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
      linewidth = .38, linetype = "dashed", color = "#374151"),
    facet_grid(status ~ ., scales = "free_y", space = "free_y"),
    scale_x_continuous(expand = c(0, 0)),
    scale_y_discrete(breaks = NULL, drop = TRUE),
    theme_bw(base_size = 10),
    theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),
      strip.text.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold"))
  )
  left_value <- if (ucell_display == "per_set_0_1")
    "ucell_visual_0_1" else "ucell_window_median"
  ucell_limit <- if (ucell_display == "per_set_0_1") 1 else
    as.numeric(quantile(d$ucell_window_median, .995))
  positive_z <- d$original_trajectory_z[
    is.finite(d$original_trajectory_z) & d$original_trajectory_z > 0]
  z_limit <- as.numeric(quantile(positive_z, .995))
  left <- ggplot(d, aes(window_center, row_id, fill = .data[[left_value]])) +
    geom_raster() + base_layers +
    scale_fill_gradientn(colors = c("#FFFFFF", "#FEE8C8", "#FDBB84",
      "#E34A33", "#7F0000"), limits = c(0, ucell_limit),
      oob = scales::squish, name = if (ucell_display == "per_set_0_1")
        "UCell\n0–1" else "Median\nUCell") +
    labs(title = if (is_residual)
        "RNA-feature-residualized UCell scaled within gene set" else if
        (ucell_display == "per_set_0_1")
        "UCell enrichment scaled within gene set" else "Raw UCell enrichment",
      subtitle = if (ucell_display == "per_set_0_1")
        "Each row min–max scaled to 0–1 for visualization only" else
        "Native-compatible within-cell rank score",
      x = "Existing Bdev pseudotime (pt)", y = "GO gene sets")
  right <- ggplot(d, aes(window_center, row_id,
      fill = original_trajectory_z)) +
    geom_raster() + base_layers +
    scale_fill_gradientn(colors = c("#FFFFFF", "#FEE8C8", "#FDBB84",
      "#E34A33", "#7F0000"), limits = c(0, z_limit),
      oob = scales::squish, na.value = "#D1D5DB",
      name = "Original\ntrajectory z") +
    labs(title = "Original standardized trajectory z-score",
      subtitle = "Positive values shown; z ≤ 0 is white",
      x = "Existing Bdev pseudotime (pt)", y = NULL)
  regions_panel <- ggplot(d, aes(window_center, row_id,
      fill = significant_region_ucell_0_1)) +
    geom_raster() + base_layers +
    scale_fill_gradientn(colors = c("#FFFFFF", "#FEE8C8", "#FDBB84",
      "#E34A33", "#7F0000"), limits = c(0, 1),
      oob = scales::squish, na.value = "white",
      name = if (is_residual) "Residual\n0–1" else "UCell\n0–1") +
    labs(title = "Permutation-supported positive regions",
      subtitle = if (is_residual)
        "Scaled residual-UCell shown only inside supported regions" else
        "Scaled UCell shown only inside supported regions",
      x = "Existing Bdev pseudotime (pt)", y = NULL)
  paired <- left + right + regions_panel +
    plot_layout(ncol = 3, widths = c(1, 1, 1)) +
    plot_annotation(title = paste0(labels[[g]],
      ": UCell dynamics and original trajectory z-score"),
      subtitle = paste(
        "Rows split by UCell-permutation status and ordered by the original z-score k-means heatmap.",
        "Gene-set labels omitted; BB1/BB2 are dashed lines."),
      theme = theme(plot.title = element_text(face = "bold", size = 15)))
  figure_height <- max(10, min(22, expected[[g]] / 220))
  file_tag <- if (ucell_display == "per_set_0_1")
    "ucell_0_1_vs_original_z_dynamic_split" else
    "raw_ucell_vs_original_z_passing_split"
  png_name <- paste0(g, "_", file_tag, ".png")
  pdf_name <- paste0(g, "_", file_tag, ".pdf")
  ggsave(file.path(stage_dir, png_name), paired,
    width = 25, height = figure_height, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, pdf_name), paired,
    width = 25, height = figure_height, device = cairo_pdf, bg = "white")
  row_order <- data.table(analysis_group = g, feature_set = display_order,
    status = ifelse(display_order %chin% passing_sets, "passing",
      "not_passing"),
    row_position_within_status = c(seq_along(passing_order),
      seq_along(nonpassing_order)),
    original_z_row_position = z_order$row_position_top_to_bottom[
      match(display_order, z_order$feature_set)],
    original_z_cluster = z_order$cluster_id[
      match(display_order, z_order$feature_set)],
    original_z_profile_defined = display_order %chin% defined_order)
  fwrite(row_order, file.path(stage_dir, paste0(g, "_heatmap_row_order.tsv")),
    sep = "\t", quote = FALSE, na = "NA")
  summaries[[g]] <- data.table(analysis_group = g,
    tested_sets = expected[[g]], passing_sets = length(passing_order),
    nonpassing_sets = length(nonpassing_order),
    undefined_original_z_profiles = length(undefined_sets),
    ucell_display_cap = ucell_limit, original_z_display_cap = z_limit,
    png = png_name, pdf = pdf_name)
}

summary_table <- rbindlist(summaries)
fwrite(summary_table, file.path(stage_dir, "heatmap_summary.tsv"),
  sep = "\t", quote = FALSE)
writeLines(c(if (is_residual)
    "# GO residualized-UCell versus original-z passing-status heatmaps" else
    "# GO raw-UCell versus original-z passing-status heatmaps", "",
  paste0("One three-panel heatmap was created for BP, MF, and CC: ",
    if (is_residual) "RNA-feature-residualized UCell" else "UCell",
    ", original standardized-gene trajectory z-score, and ",
    if (is_residual) "residualized UCell" else "UCell",
    " restricted to permutation-supported positive regions. No UCell-derived z-score is used."), "",
  "| Group | Passing | Non-passing | Undefined original-z profiles |",
  "|---|---:|---:|---:|",
  vapply(seq_len(nrow(summary_table)), function(i) paste0("| ",
    toupper(summary_table$analysis_group[[i]]), " | ",
    summary_table$passing_sets[[i]], " | ",
    summary_table$nonpassing_sets[[i]], " | ",
    summary_table$undefined_original_z_profiles[[i]], " |"), character(1)), "",
  if (identical(completion$inference_mode, "dynamic_range"))
    "- Passing means a significant trajectory-wide peak-to-trough UCell range at group-wise BH q<0.01."
  else "- Passing means at least one UCell-permutation region at group-wise BH q<0.01.",
  "- Rows are ordered by the full original-z k-means heatmap within passing status.",
  "- GO names are omitted from the figures; exact orders are retained in group-specific TSV files.",
  paste0("- ", if (is_residual) "Residualized UCell" else "UCell",
    " is independently min-max scaled to 0–1 within every gene set for visualization only; inference used unscaled scores."),
  "- All panels use white-to-red scales; original z<=0 and windows outside supported regions are white.",
  "- Undefined original-z profiles are grey and appended after defined profiles within status."),
  file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", groups = summaries,
  ucell_z_score_used = FALSE,
  left_value = if (ucell_display == "per_set_0_1")
    paste("per-gene-set min-max-scaled median",
      if (is_residual) "RNA-feature-residualized UCell" else "UCell",
      "(visualization only)") else
    "raw median UCell",
  right_value = "original standardized-gene trajectory z-score",
  third_value = paste("0-1",
    if (is_residual) "RNA-feature-residualized UCell" else "UCell",
    "inside permutation-supported positive regions"),
  complexity_residualization_enabled = is_residual,
  ucell_rescaled_for_inference = FALSE,
  row_order = "original-z k-means heatmap order within passing status",
  passing_definition = if (identical(completion$inference_mode,
      "dynamic_range"))
    "trajectory-wide UCell range at group-wise BH q<0.01" else
    "at least one UCell-permutation region at group-wise BH q<0.01",
  gene_set_names_displayed = FALSE), file.path(stage_dir, "completion.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)

stopifnot(file.rename(stage_dir, output_dir))
message("Completed GO paired heatmaps at: ", output_dir)
