#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 48_top_adjusted_wt_ko_motif_heatmaps.R",
  "<all_motif_report> <adjusted_score_report> <output_dir>",
  "<top_per_direction> <activity_display_limit>"))
motif_report <- normalizePath(args[[1L]], mustWork = TRUE)
score_report <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
top_n <- as.integer(args[[4L]])
activity_limit <- as.numeric(args[[5L]])
stopifnot(top_n >= 1L, is.finite(activity_limit), activity_limit > 0)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})
write_tsv <- function(x, path) data.table::fwrite(x, path, sep = "\t",
  na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

values <- fread(cmd = paste("zcat", shQuote(file.path(motif_report,
  "all_eligible_motif_pseudotime_values.tsv.gz"))))
scores <- fread(file.path(score_report,
  "trajectory_adjusted_wt_ko_motif_scores.tsv"))
counts <- fread(file.path(motif_report, "pseudotime_bin_cell_counts.tsv"))
stopifnot(setequal(unique(values$region_id), c("Flt3l_peak1", "Il7_peak2")),
  !anyDuplicated(scores[, .(region_id, motif_id)]),
  all(scores$strong_bdev_directional_pattern))

all_selected <- list()
all_heat <- list()
diagnostics <- list()
for (region_key in c("Flt3l_peak1", "Il7_peak2")) {
  region_scores <- scores[region_id == region_key]
  wt_higher <- head(region_scores[adjusted_wt_minus_ko_bdev_sd > 0][
    order(-adjusted_wt_minus_ko_bdev_sd)], top_n)
  ko_higher <- head(region_scores[adjusted_wt_minus_ko_bdev_sd < 0][
    order(adjusted_wt_minus_ko_bdev_sd)], top_n)
  selected <- rbindlist(list(wt_higher, ko_higher), fill = TRUE)
  selected[, `:=`(
    selection_direction = fifelse(adjusted_wt_minus_ko_bdev_sd > 0,
      "WT higher", "KO higher"),
    selection_rank = seq_len(.N)), by = .(
      favored_group)]
  selected[, row_label := sprintf("%s%s  [%+.2f]",
    motif_name, fifelse(strong_consistent_difference, "*", ""),
    adjusted_wt_minus_ko_bdev_sd)]
  selected[, row_order := seq_len(.N)]
  all_selected[[region_key]] <- selected

  region_values <- values[region_id == region_key &
    motif_id %in% selected$motif_id & group %in% c("Bdev", "WT", "KO")]
  activity <- region_values[, .(region_id, motif_id, motif_name, group,
    bin, pt_mid, n_cells, value = bdev_scaled_deviation)]
  wt <- activity[group == "WT", .(motif_id, bin, pt_mid,
    wt = value, wt_n = n_cells)]
  ko <- activity[group == "KO", .(motif_id, bin,
    ko = value, ko_n = n_cells)]
  difference <- merge(wt, ko, by = c("motif_id", "bin"), sort = FALSE)
  difference[, `:=`(region_id = region_key, motif_name = NA_character_,
    group = "WT - KO", n_cells = pmin(wt_n, ko_n), value = wt - ko)]
  difference[, c("wt", "ko", "wt_n", "ko_n") := NULL]
  heat <- rbindlist(list(activity, difference), fill = TRUE)
  heat <- merge(heat, selected[, .(motif_id, row_label, row_order,
    selection_direction, adjusted_wt_minus_ko_bdev_sd,
    strong_consistent_difference)], by = "motif_id", all.x = TRUE,
    sort = FALSE)
  region_counts <- counts[region_id == region_key]
  all_bins <- sort(unique(region_counts$bin))
  bin_mid <- region_counts[, .(pt_mid = unique(pt_mid)), by = bin]
  grid <- CJ(motif_id = selected$motif_id,
    group = c("Bdev", "WT", "KO", "WT - KO"), bin = all_bins,
    unique = TRUE)
  heat <- merge(grid, heat, by = c("motif_id", "group", "bin"),
    all.x = TRUE, sort = FALSE)
  heat <- merge(heat, bin_mid, by = "bin", all.x = TRUE,
    suffixes = c("", ".grid"), sort = FALSE)
  heat[is.na(pt_mid), pt_mid := pt_mid.grid]
  heat[, pt_mid.grid := NULL]
  heat <- merge(heat, selected[, .(motif_id, row_label, row_order,
    selection_direction, adjusted_wt_minus_ko_bdev_sd,
    strong_consistent_difference)], by = "motif_id", all.x = TRUE,
    suffixes = c("", ".selected"), sort = FALSE)
  for (column in c("row_label", "row_order", "selection_direction",
      "adjusted_wt_minus_ko_bdev_sd", "strong_consistent_difference")) {
    selected_column <- paste0(column, ".selected")
    heat[is.na(get(column)), (column) := get(selected_column)]
    heat[, (selected_column) := NULL]
  }
  heat[, `:=`(row_label = factor(row_label,
      levels = rev(selected$row_label)), region_id = region_key)]
  all_heat[[region_key]] <- copy(heat)

  activity_heat <- heat[group %in% c("Bdev", "WT", "KO")]
  activity_labels <- if (region_key == "Flt3l_peak1")
    c(Bdev = "Bdev", WT = "FLT3L WT", KO = "FLT3LKO") else
    c(Bdev = "Bdev", WT = "IL-7 WT", KO = "IL7KO")
  activity_heat[, panel := factor(activity_labels[group],
    levels = unname(activity_labels))]
  difference_heat <- heat[group == "WT - KO"]
  difference_limit <- as.numeric(quantile(abs(difference_heat$value),
    0.98, na.rm = TRUE))
  difference_limit <- max(1, min(6, difference_limit))

  p_activity <- ggplot(activity_heat, aes(bin, row_label, fill = value)) +
    geom_tile() +
    facet_grid(. ~ panel) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-activity_limit, activity_limit),
      oob = scales::squish, na.value = "grey75",
      name = "Activity\n(Bdev SD)") +
    scale_x_continuous(expand = expansion(mult = 0),
      breaks = pretty(all_bins, n = 4)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6.5),
      legend.position = "right") +
    labs(x = "Mapped pseudotime bin", y = NULL)

  p_difference <- ggplot(difference_heat,
      aes(bin, row_label, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-difference_limit, difference_limit),
      oob = scales::squish, na.value = "grey75",
      name = "WT − KO\n(Bdev SD)") +
    scale_x_continuous(expand = expansion(mult = 0),
      breaks = pretty(all_bins, n = 4)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      legend.position = "right") +
    labs(title = "Matched-bin WT − KO", x = "Mapped pseudotime bin",
      y = NULL)

  region_title <- if (region_key == "Flt3l_peak1")
    "Flt3L peak 1: top trajectory-adjusted WT/KO motif differences" else
    "IL-7 peak 2: top trajectory-adjusted WT/KO motif differences"
  figure <- p_activity + p_difference +
    plot_layout(widths = c(3, 1.15)) +
    plot_annotation(title = region_title,
      subtitle = paste0("Top ", nrow(wt_higher), " WT-higher followed by top ",
        nrow(ko_higher),
        " KO-higher motifs; rows selected by the Bdev-distribution-adjusted contrast"),
      caption = paste0("Bracketed row values are adjusted WT−KO effects in Bdev SD; * passes |effect|≥0.5 and ≥70% sign consistency. Gray denotes unsupported bins. Activity clipped at ±",
        activity_limit, "; difference clipped at ±",
        signif(difference_limit, 3), " for display only."),
      theme = theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 9.5),
        plot.caption = element_text(size = 8)))
  stub <- if (region_key == "Flt3l_peak1") "flt3l_peak1" else "il7_peak2"
  ggsave(file.path(stage_dir, "figures", paste0(stub,
    "_top_adjusted_wt_ko_motif_heatmap.png")), figure,
    width = 18, height = 11, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(stub,
    "_top_adjusted_wt_ko_motif_heatmap.pdf")), figure,
    width = 18, height = 11, device = cairo_pdf, bg = "white")
  diagnostics[[region_key]] <- list(
    selected_wt_higher = nrow(wt_higher),
    selected_ko_higher = nrow(ko_higher),
    strong_consistent_wt_higher = sum(wt_higher$strong_consistent_difference),
    strong_consistent_ko_higher = sum(ko_higher$strong_consistent_difference),
    activity_display_limit = activity_limit,
    difference_display_limit = difference_limit)
}

write_tsv(rbindlist(all_selected, fill = TRUE), file.path(stage_dir,
  "selected_top_adjusted_wt_ko_motifs.tsv"))
write_tsv(rbindlist(all_heat, fill = TRUE), file.path(stage_dir,
  "heatmap_values.tsv"))
write_json(list(status = "complete_preliminary_unregistered",
  motif_report = motif_report, adjusted_score_report = score_report,
  top_per_direction = top_n,
  selection = "largest positive and negative trajectory-adjusted WT-minus-KO effects among strong Bdev directional motifs",
  panels = list("Bdev", "WT", "KO", "matched-bin WT-minus-KO"),
  values = "chromVAR deviation centered and scaled by Bdev temporal mean and SD",
  region_diagnostics = diagnostics,
  biological_replicates_per_genotype = 1L,
  inference = "descriptive", bundle_mutated = FALSE,
  capsule_registered = FALSE), file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Top trajectory-adjusted WT versus KO motif heatmaps", "",
  paste0("Each region shows up to ", top_n,
    " WT-higher and ", top_n,
    " KO-higher motifs ranked by the Bdev-distribution-adjusted WT-minus-KO effect."),
  "Rows come from the full strong Bdev directional-motif screen, not the original top-30 display.",
  "The Bdev, WT, and KO panels show Bdev-standardized chromVAR deviations across pseudotime; the final panel shows matched-bin WT-minus-KO.",
  "Stars identify the prespecified strong, sign-consistent descriptive differences.",
  "There is one biological sample per genotype, so the results are descriptive."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed top adjusted WT/KO motif heatmaps at", output_dir, "\n")
