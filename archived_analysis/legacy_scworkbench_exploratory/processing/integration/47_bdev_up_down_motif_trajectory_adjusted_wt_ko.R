#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 47_bdev_up_down_motif_trajectory_adjusted_wt_ko.R",
  "<bdev_motif_report> <output_dir> <minimum_bdev_change_sd>",
  "<minimum_absolute_spearman> <minimum_effect_bdev_sd>",
  "<minimum_sign_consistency>"))
source_report <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
minimum_bdev_change <- as.numeric(args[[3L]])
minimum_absolute_spearman <- as.numeric(args[[4L]])
minimum_effect <- as.numeric(args[[5L]])
minimum_sign_consistency <- as.numeric(args[[6L]])
stopifnot(is.finite(minimum_bdev_change), minimum_bdev_change > 0,
  is.finite(minimum_absolute_spearman), minimum_absolute_spearman > 0,
  minimum_absolute_spearman <= 1, is.finite(minimum_effect),
  minimum_effect > 0, is.finite(minimum_sign_consistency),
  minimum_sign_consistency > 0.5, minimum_sign_consistency <= 1)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
write_tsv <- function(x, path) data.table::fwrite(x, path, sep = "\t",
  na = "NA")
weighted_mean_safe <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}
jackknife_interval <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  n <- length(x)
  if (n < 3L) return(c(se = NA_real_, lower = NA_real_, upper = NA_real_))
  leave_one_out <- vapply(seq_len(n), function(i)
    weighted_mean_safe(x[-i], w[-i]), numeric(1L))
  theta_bar <- mean(leave_one_out)
  se <- sqrt((n - 1) / n * sum((leave_one_out - theta_bar)^2))
  estimate <- weighted_mean_safe(x, w)
  c(se = se, lower = estimate - 1.96 * se, upper = estimate + 1.96 * se)
}

heat <- fread(cmd = paste("zcat", shQuote(file.path(source_report,
  "all_eligible_motif_pseudotime_values.tsv.gz"))))
variability <- fread(file.path(source_report,
  "all_motif_bdev_temporal_variability.tsv"))
counts <- fread(file.path(source_report, "pseudotime_bin_cell_counts.tsv"))
stopifnot(setequal(unique(heat$region_id), c("Flt3l_peak1", "Il7_peak2")),
  all(c("Bdev", "WT", "KO") %in% unique(heat$group)),
  !anyDuplicated(variability[, .(region_id, motif_id)]))

all_scores <- list()
all_bins <- list()
all_plotted <- list()
diagnostics <- list()
for (region_key in c("Flt3l_peak1", "Il7_peak2")) {
  region_variability <- variability[region_id == region_key]
  bdev <- heat[region_id == region_key & group == "Bdev",
    .(motif_id, motif_name, bin, pt_mid,
      bdev_scaled_deviation, bdev_n = n_cells)]
  bdev_bins <- sort(unique(bdev$bin))
  third <- max(1L, floor(length(bdev_bins) / 3L))
  early_bins <- head(bdev_bins, third)
  late_bins <- tail(bdev_bins, third)
  bdev_change <- bdev[, .(
    bdev_early = weighted_mean_safe(
      bdev_scaled_deviation[bin %in% early_bins],
      bdev_n[bin %in% early_bins]),
      bdev_late = weighted_mean_safe(
      bdev_scaled_deviation[bin %in% late_bins],
      bdev_n[bin %in% late_bins]),
    bdev_spearman = cor(pt_mid, bdev_scaled_deviation,
      method = "spearman")), by = .(motif_id, motif_name)]
  bdev_change[, bdev_late_minus_early := bdev_late - bdev_early]
  bdev_change[, bdev_direction := ifelse(bdev_late_minus_early >= 0,
    "Bdev up", "Bdev down")]
  bdev_change[, strong_bdev_directional_pattern :=
    abs(bdev_late_minus_early) >= minimum_bdev_change &
    abs(bdev_spearman) >= minimum_absolute_spearman &
    sign(bdev_late_minus_early) == sign(bdev_spearman)]
  chosen <- bdev_change[strong_bdev_directional_pattern == TRUE]
  chosen[, motif_label := make.unique(sprintf("%s (%s)",
    motif_name, motif_id))]
  stopifnot(nrow(chosen) > 0L, !anyDuplicated(chosen$motif_id))

  wt <- heat[region_id == region_key & group == "WT",
    .(motif_id, bin, pt_mid, wt = bdev_scaled_deviation, wt_n = n_cells)]
  ko <- heat[region_id == region_key & group == "KO",
    .(motif_id, bin, ko = bdev_scaled_deviation, ko_n = n_cells)]
  paired <- merge(wt, ko, by = c("motif_id", "bin"), sort = FALSE)
  bdev_weight <- counts[region_id == region_key & group == "Bdev" &
      included == TRUE, .(bin, bdev_n = n_cells)]
  paired <- merge(paired, bdev_weight, by = "bin", all.x = TRUE,
    sort = FALSE)
  paired[, `:=`(wt_minus_ko = wt - ko, region_id = region_key)]
  paired <- paired[motif_id %in% chosen$motif_id]
  stopifnot(!anyNA(paired$bdev_n))

  scores <- paired[, {
    effect <- weighted_mean_safe(wt_minus_ko, bdev_n)
    jk <- jackknife_interval(wt_minus_ko, bdev_n)
    consistency <- weighted_mean_safe(as.numeric(
      sign(wt_minus_ko) == sign(effect)), bdev_n)
    .(adjusted_wt_minus_ko_bdev_sd = effect,
      unadjusted_wt_minus_ko_bdev_sd =
        weighted_mean_safe(wt, wt_n) - weighted_mean_safe(ko, ko_n),
      jackknife_se = unname(jk[["se"]]),
      jackknife_lower = unname(jk[["lower"]]),
      jackknife_upper = unname(jk[["upper"]]),
      direction_consistency = consistency,
      common_bins = .N,
      common_bdev_cells = sum(bdev_n))
  }, by = motif_id]
  scores <- merge(chosen, scores, by = "motif_id", all.x = TRUE,
    sort = FALSE)
  scores <- merge(scores, region_variability[, .(motif_id,
    bdev_variability_rank, bdev_temporal_sd)], by = "motif_id",
    all.x = TRUE, sort = FALSE)
  scores[, `:=`(
    strong_consistent_difference =
      abs(adjusted_wt_minus_ko_bdev_sd) >= minimum_effect &
      direction_consistency >= minimum_sign_consistency,
    jackknife_interval_excludes_zero =
      jackknife_lower > 0 | jackknife_upper < 0,
    favored_group = fifelse(adjusted_wt_minus_ko_bdev_sd > 0,
      "WT higher", "KO higher"),
    region_id = region_key)]
  scores[, absolute_adjusted_wt_minus_ko_bdev_sd :=
    abs(adjusted_wt_minus_ko_bdev_sd)]
  setorder(scores, bdev_direction,
    -absolute_adjusted_wt_minus_ko_bdev_sd)
  all_scores[[region_key]] <- scores
  paired <- merge(paired, scores[, .(motif_id, motif_name, motif_label,
    bdev_direction)], by = "motif_id", all.x = TRUE, sort = FALSE)
  all_bins[[region_key]] <- paired

  plot_data <- melt(scores,
    id.vars = c("motif_id", "motif_name", "motif_label", "bdev_direction",
      "strong_consistent_difference", "jackknife_lower", "jackknife_upper"),
    measure.vars = c("unadjusted_wt_minus_ko_bdev_sd",
      "adjusted_wt_minus_ko_bdev_sd"),
    variable.name = "estimate_type", value.name = "wt_minus_ko_bdev_sd")
  plot_data[, estimate_type := factor(estimate_type,
    levels = c("unadjusted_wt_minus_ko_bdev_sd",
      "adjusted_wt_minus_ko_bdev_sd"),
    labels = c("Unadjusted", "Trajectory-adjusted"))]
  ordering <- scores[order(bdev_direction,
    adjusted_wt_minus_ko_bdev_sd), motif_label]
  plot_data[, motif_label := factor(motif_label, levels = ordering)]
  scores[order(adjusted_wt_minus_ko_bdev_sd),
    plot_page := ceiling(seq_len(.N) / 50), by = bdev_direction]
  plot_data <- merge(plot_data,
    scores[, .(motif_id, bdev_direction, plot_page)],
    by = c("motif_id", "bdev_direction"), all.x = TRUE, sort = FALSE)
  all_plotted[[region_key]] <- copy(plot_data)

  stub <- if (region_key == "Flt3l_peak1") "flt3l_peak1" else "il7_peak2"
  for (direction in c("Bdev up", "Bdev down")) {
    direction_pages <- sort(unique(
      plot_data[bdev_direction == direction, plot_page]))
    for (page_index in direction_pages) {
      page_data <- plot_data[bdev_direction == direction &
        plot_page == page_index]
      page_levels <- scores[bdev_direction == direction &
        plot_page == page_index][order(adjusted_wt_minus_ko_bdev_sd),
          motif_label]
      page_data[, motif_label := factor(as.character(motif_label),
        levels = page_levels)]
      p <- ggplot(page_data, aes(wt_minus_ko_bdev_sd, motif_label,
          color = estimate_type, shape = estimate_type)) +
        geom_vline(xintercept = 0, color = "grey45", linewidth = 0.45) +
        geom_errorbar(data = page_data[estimate_type ==
            "Trajectory-adjusted"],
          aes(xmin = jackknife_lower, xmax = jackknife_upper),
          orientation = "y", width = 0, color = "#2166AC",
          linewidth = 0.55) +
        geom_point(size = 2.6) +
        geom_point(data = page_data[estimate_type ==
            "Trajectory-adjusted" & strong_consistent_difference == TRUE],
          shape = 21, size = 4.2, stroke = 0.8, color = "black",
          fill = NA, inherit.aes = FALSE,
          aes(x = wt_minus_ko_bdev_sd, y = motif_label)) +
        scale_color_manual(values = c(Unadjusted = "grey65",
          `Trajectory-adjusted` = "#2166AC")) +
        scale_shape_manual(values = c(Unadjusted = 16,
          `Trajectory-adjusted` = 18)) +
        theme_bw(base_size = 10.5) +
        theme(panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 7.5),
          plot.title = element_text(face = "bold"),
          legend.position = "top") +
        labs(title = paste0(if (region_key == "Flt3l_peak1")
            "Flt3L peak 1" else "IL-7 peak 2", ": ", direction,
            " motifs, page ", page_index, "/", length(direction_pages)),
          subtitle = paste0("All strong directional motifs (|change| ≥ ",
            minimum_bdev_change, " SD; |Spearman| ≥ ",
            minimum_absolute_spearman,
            "); matched-bin estimates weighted to Bdev"),
          x = "WT − KO chromVAR deviation (Bdev temporal SD)", y = NULL,
          caption = paste0("Blue intervals: leave-one-bin-out stability, not biological-replicate uncertainty. Black rings: |effect| ≥ ",
            minimum_effect, " and ≥", 100 * minimum_sign_consistency,
            "% sign consistency."), color = NULL, shape = NULL)
      direction_stub <- if (direction == "Bdev up") "bdev_up" else
        "bdev_down"
      file_stub <- sprintf("%s_%s_adjusted_wt_ko_page%02d",
        stub, direction_stub, page_index)
      figure_height <- max(8, 3.5 + 0.24 * length(page_levels))
      ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".png")),
        p, width = 10, height = figure_height, dpi = 300, bg = "white")
      ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".pdf")),
        p, width = 10, height = figure_height, device = cairo_pdf,
        bg = "white")
    }
  }

  diagnostics[[region_key]] <- list(
    selected_bdev_up = sum(scores$bdev_direction == "Bdev up"),
    selected_bdev_down = sum(scores$bdev_direction == "Bdev down"),
    common_bins_range = range(scores$common_bins),
    strong_consistent_differences = sum(scores$strong_consistent_difference),
    jackknife_stable_differences = sum(
      scores$jackknife_interval_excludes_zero))
}

write_tsv(rbindlist(all_scores, fill = TRUE), file.path(stage_dir,
  "trajectory_adjusted_wt_ko_motif_scores.tsv"))
write_tsv(rbindlist(all_bins, fill = TRUE), file.path(stage_dir,
  "matched_bin_wt_ko_motif_differences.tsv"))
write_tsv(rbindlist(all_plotted, fill = TRUE), file.path(stage_dir,
  "figure_values.tsv"))
write_json(list(status = "complete_preliminary_unregistered",
  source_report = source_report,
  motif_selection = "all eligible motifs passing the strong Bdev directional-pattern thresholds",
  bdev_direction = "Bdev-cell-weighted last-third minus first-third mean",
  strong_bdev_directional_thresholds = list(
    minimum_absolute_late_minus_early_bdev_sd = minimum_bdev_change,
    minimum_absolute_spearman = minimum_absolute_spearman,
    change_and_spearman_must_agree = TRUE),
  trajectory_adjustment =
    "matched-bin WT-minus-KO contrasts standardized to the Bdev cell distribution",
  strong_descriptive_thresholds = list(
    minimum_absolute_effect_bdev_sd = minimum_effect,
    minimum_bdev_weighted_sign_consistency = minimum_sign_consistency),
  uncertainty = "leave-one-pseudotime-bin-out stability interval",
  biological_replicates_per_genotype = 1L,
  inference = "descriptive; no biological-replicate significance claim",
  region_diagnostics = diagnostics,
  bundle_mutated = FALSE, capsule_registered = FALSE),
  file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Trajectory-adjusted WT versus KO motif comparison", "",
  "The motif universe and variability selection are inherited unchanged from the Bdev-only motif heatmap analysis.",
  paste0("Within each region, every eligible motif is retained when its absolute Bdev-weighted last-third minus first-third change is at least ",
    minimum_bdev_change, " temporal SD, its absolute Bdev Spearman correlation with pseudotime is at least ",
    minimum_absolute_spearman, ", and the two directions agree."),
  "WT-minus-KO differences are calculated within matched mapped-pseudotime bins and standardized to the Bdev cell distribution, removing differences caused only by unequal WT/KO occupancy along the interval.",
  paste0("Black rings flag descriptive differences with absolute adjusted effect at least ",
    minimum_effect, " Bdev temporal SD and Bdev-weighted sign consistency at least ",
    100 * minimum_sign_consistency, "% across bins."),
  "Leave-one-bin-out intervals quantify sensitivity to pseudotime strata, not uncertainty across biological replicates.",
  "There is one biological sample per genotype, so the analysis is descriptive."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed trajectory-adjusted motif comparison at", output_dir, "\n")
