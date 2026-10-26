#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(paste(
  "Usage: 95_bdev_five_region_ebf1_dynamic_peaks.R BUNDLE",
  "SOURCE_DA_REPORT DENSITY_TABLE OUTPUT_DIR SESSION_ROOT",
  "EXPECTED_FINGERPRINT SCRIPT_VERSION TCF3_COVERAGE"
), call. = FALSE)

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
source_da_dir <- normalizePath(args[[2L]], mustWork = TRUE)
density_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
session_root <- normalizePath(args[[5L]], mustWork = TRUE)
expected_fp <- args[[6L]]
script_version <- args[[7L]]
tcf3_coverage_path <- normalizePath(args[[8L]], mustWork = TRUE)

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
})

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    paste(..., collapse = " "), "\n")
  flush.console()
}

source_peak_path <- file.path(source_da_dir,
  "ebf1_pm2mb_bdev_peak_universe.tsv")
stopifnot(file.exists(source_peak_path))
peaks <- fread(source_peak_path)
stopifnot(nrow(peaks) == 251L, !anyDuplicated(peaks$bdev_feature_id))

# These states partition the union of the current shared WT/KO-supported
# Flt3L (bins 14-25) and IL-7 (bins 22-30) domains. Four-bin flanks match the
# width of the dual state and serve as local Bdev controls.
region_spec <- data.table(
  region_order = 1:5,
  region_id = c("before_flt3l", "flt3l_only", "dual", "il7_only",
    "after_il7"),
  region_label = c("Before Flt3L", "Flt3L only", "Dual",
    "IL-7 only", "After IL-7"),
  first_bin = c(10L, 14L, 22L, 26L, 31L),
  last_bin = c(13L, 21L, 25L, 30L, 34L),
  definition = c("four-bin local flank", "Flt3L shared-support only",
    "overlap of Flt3L and IL-7 shared-support domains",
    "IL-7 shared-support only", "four-bin local flank"))

density <- fread(cmd = paste("zcat", shQuote(density_path)))
bin_grid <- unique(density[, .(bin, bin_left, bin_right, bin_midpoint,
  bdev_n_cells)])
stopifnot(!anyDuplicated(bin_grid$bin), all(1:100 %in% bin_grid$bin))
region_spec[, `:=`(
  pt_left = bin_grid$bin_left[match(first_bin, bin_grid$bin)],
  pt_right = bin_grid$bin_right[match(last_bin, bin_grid$bin)])]
stopifnot(all(is.finite(region_spec$pt_left)),
  all(is.finite(region_spec$pt_right)),
  all(region_spec$pt_right[-nrow(region_spec)] ==
    region_spec$pt_left[-1L]))

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$current_state$fingerprint, expected_fp),
  all(peaks$bdev_feature_id %in%
    wb_feature_ids(bundle, "ATAC_bdev_mm10_peaks")))
bdev_cells <- as.character(wb_cells(bundle, "set:bdev_reference"))
bdev_meta <- as.data.table(wb_meta(bundle, cells = bdev_cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "reference_pt")))
setnames(bdev_meta, names(bdev_meta)[1L], "cell_id")
setnames(bdev_meta, "reference_pt", "pt")
stopifnot(identical(as.character(bdev_meta$cell_id), bdev_cells))

assignments <- rbindlist(lapply(seq_len(nrow(region_spec)), function(i) {
  r <- region_spec[i]
  bdev_meta[is.finite(pt) & pt >= r$pt_left & pt < r$pt_right,
    .(cell_id, pt, region_order = r$region_order,
      region_id = r$region_id, region_label = r$region_label)]
}))
stopifnot(!anyDuplicated(assignments$cell_id), nrow(assignments) > 1000L)
counts <- assignments[, .(bdev_cells = .N, minimum_pt = min(pt),
  maximum_pt = max(pt)), by = .(region_order, region_id, region_label)]
region_spec <- merge(region_spec, counts,
  by = c("region_order", "region_id", "region_label"), all.x = TRUE,
  sort = FALSE)
setorder(region_spec, region_order)
stopifnot(!anyNA(region_spec$bdev_cells), all(region_spec$bdev_cells >= 100L))
write_tsv(region_spec, file.path(stage_dir, "bdev_five_region_definition.tsv"))
fwrite(assignments, file.path(stage_dir,
  "bdev_five_region_cell_assignments.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")

# Reuse native full-matrix log1p-CP10K normalization from analysis 91.
session_id <- "ebf1_pm2mb_shared_support_da_v1"
workspace <- file.path(session_root, session_id)
stopifnot(dir.exists(workspace))
session <- wb_session_resume(bundle, workspace)
stopifnot(is.null(session$durable$promotion_receipt))
find_run <- function(session, marker, operation) {
  ids <- names(session$runs)[vapply(session$runs, function(x) {
    note <- if (is.null(x$notes)) "" else as.character(x$notes)
    identical(x$op_name, operation) && startsWith(note, paste0("[", marker, "]"))
  }, logical(1L))]
  stopifnot(length(ids) <= 1L)
  if (length(ids)) session$runs[[ids[[1L]]]] else NULL
}
normal_run <- find_run(session, "normalize_full_atac", "matrix_log_normalize")
stopifnot(!is.null(normal_run), "matrix" %in% names(normal_run$outputs))
normalized_id <- "ATAC_bdev_mm10_peaks_log1p_cp10k_ebf1_pm2mb_shared_v1"

transition_spec <- data.table(
  transition_order = 1:4,
  transition_id = c("before_to_flt3l_only", "flt3l_only_to_dual",
    "dual_to_il7_only", "il7_only_to_after"),
  group_a_region = region_spec$region_id[2:5],
  group_b_region = region_spec$region_id[1:4],
  transition_label = c("Before → Flt3L only", "Flt3L only → Dual",
    "Dual → IL-7 only", "IL-7 only → After"))

de_list <- list()
for (i in seq_len(nrow(transition_spec))) {
  tr <- transition_spec[i]
  group_a <- assignments[region_id == tr$group_a_region, cell_id]
  group_b <- assignments[region_id == tr$group_b_region, cell_id]
  marker <- paste0("bdev_dynamic_", tr$transition_id)
  de_run <- find_run(session, marker, "de_wilcox")
  if (is.null(de_run)) {
    log_event("Native Bdev DA", tr$transition_label,
      "group A", length(group_a), "group B", length(group_b))
    result <- wb_session_run_op(session, bundle, "de_wilcox",
      inputs = list(matrix = normalized_id,
        cell_set = "set:bdev_reference", min_cells = 10L,
        group_a = as.character(group_a), group_b = as.character(group_b),
        features = as.character(peaks$bdev_feature_id), min_pct = 0,
        p_adjust_method = "BH"),
      notes = paste0("[", marker, "] Ordered Bdev adjacent-state DA; ",
        "group A=", tr$group_a_region, "; group B=", tr$group_b_region,
        "; exact preserved Bdev pt; Ebf1 TSS +/-2 Mb; native full-matrix ",
        "log1p-CP10K; BH within 251 peaks; staged only; no promotion."))
    session <- result$session
    de_run <- result$run
  }
  d <- as.data.table(de_run$outputs$de_table$table)
  stopifnot(nrow(d) == 251L, setequal(d$feature_id,
    peaks$bdev_feature_id))
  d <- merge(d, peaks[, .(feature_id = bdev_feature_id,
    chromosome = mm10_seqname, peak_start_1based = mm10_start_1based,
    peak_end_1based = mm10_end_1based, peak_midpoint_1based,
    midpoint_relative_to_tss_bp)], by = "feature_id", all.x = TRUE)
  d[, `:=`(transition_order = tr$transition_order,
    transition_id = tr$transition_id,
    transition_label = tr$transition_label,
    group_a_region = tr$group_a_region,
    group_b_region = tr$group_b_region,
    detection_difference_later_minus_earlier = pct_a - pct_b,
    direction = fifelse(avg_log1p_fc > 0, "later_higher",
      fifelse(avg_log1p_fc < 0, "earlier_higher", "equal")),
    significant_bh05 = is.finite(p_val_adj) & p_val_adj < 0.05,
    robust_transition = is.finite(p_val_adj) & p_val_adj < 0.05 &
      abs(pct_a - pct_b) >= 0.05 &
      sign(avg_log1p_fc) == sign(pct_a - pct_b))]
  de_list[[tr$transition_id]] <- d
}
transition_da <- rbindlist(de_list, use.names = TRUE)
transition_da <- transition_da[
  order(transition_order, p_val_adj,
    -abs(detection_difference_later_minus_earlier), feature_id)]
write_tsv(transition_da, file.path(stage_dir,
  "all_adjacent_region_ebf1_pm2mb_da.tsv"))
write_tsv(transition_da[significant_bh05 == TRUE], file.path(stage_dir,
  "significant_adjacent_region_ebf1_pm2mb_da_bh05.tsv"))
write_tsv(transition_da[robust_transition == TRUE], file.path(stage_dir,
  "robust_adjacent_region_ebf1_pm2mb_da_bh05_pctdiff05.tsv"))

# Direct per-state accessibility summaries provide the ordered trajectory and
# dynamic range used for ranking, independently of pairwise P-value ordering.
peak_counts <- as.matrix(wb_matrix_read_features(bundle,
  "ATAC_bdev_mm10_peaks", peaks$bdev_feature_id,
  assignments$cell_id, missing = "NA"))
stopifnot(!anyNA(peak_counts), identical(rownames(peak_counts),
  peaks$bdev_feature_id))
region_profiles <- rbindlist(lapply(seq_len(nrow(region_spec)), function(i) {
  r <- region_spec[i]
  idx <- which(assignments$region_id == r$region_id)
  data.table(feature_id = rownames(peak_counts),
    region_order = r$region_order, region_id = r$region_id,
    region_label = r$region_label, bdev_cells = length(idx),
    accessible_fraction = rowMeans(peak_counts[, idx, drop = FALSE] > 0),
    mean_raw_count = rowMeans(peak_counts[, idx, drop = FALSE]))
}))
write_tsv(region_profiles, file.path(stage_dir,
  "bdev_five_region_peak_accessibility_profiles.tsv"))

dynamic <- region_profiles[, .(
  minimum_accessible_fraction = min(accessible_fraction),
  maximum_accessible_fraction = max(accessible_fraction),
  accessible_fraction_range = max(accessible_fraction) -
    min(accessible_fraction),
  maximum_region = region_label[which.max(accessible_fraction)],
  minimum_region = region_label[which.min(accessible_fraction)]),
  by = feature_id]
transition_summary <- transition_da[, .(
  minimum_adjacent_bh_q = if (all(is.na(p_val_adj))) NA_real_ else
    min(p_val_adj, na.rm = TRUE),
  significant_adjacent_transitions = sum(significant_bh05),
  robust_adjacent_transitions = sum(robust_transition),
  largest_adjacent_detection_change = max(
    abs(detection_difference_later_minus_earlier)),
  largest_transition = transition_label[
    which.max(abs(detection_difference_later_minus_earlier))]),
  by = feature_id]
dynamic <- merge(dynamic, transition_summary, by = "feature_id")
dynamic <- merge(dynamic, peaks[, .(feature_id = bdev_feature_id,
  chromosome = mm10_seqname, peak_start_1based = mm10_start_1based,
  peak_end_1based = mm10_end_1based, peak_midpoint_1based,
  midpoint_relative_to_tss_bp)], by = "feature_id")
dynamic[, `:=`(
  statistically_dynamic = significant_adjacent_transitions > 0L,
  robust_dynamic = robust_adjacent_transitions > 0L)]
setorder(dynamic, -robust_dynamic, minimum_adjacent_bh_q,
  -accessible_fraction_range, feature_id)
write_tsv(dynamic, file.path(stage_dir,
  "ebf1_pm2mb_dynamic_peak_summary.tsv"))
write_tsv(dynamic[robust_dynamic == TRUE], file.path(stage_dir,
  "robust_dynamic_ebf1_pm2mb_peaks.tsv"))

transition_counts <- transition_da[, .(
  significant_bh05 = sum(significant_bh05),
  robust_bh05_pctdiff05 = sum(robust_transition),
  later_higher_robust = sum(robust_transition & direction == "later_higher"),
  earlier_higher_robust = sum(robust_transition & direction == "earlier_higher")),
  by = .(transition_order, transition_id, transition_label)]
write_tsv(transition_counts, file.path(stage_dir,
  "dynamic_peak_counts_by_transition.tsv"))

robust_ids <- dynamic[robust_dynamic == TRUE, feature_id]
plot_profiles <- region_profiles[feature_id %in% robust_ids]
plot_profiles <- merge(plot_profiles,
  dynamic[, .(feature_id, minimum_adjacent_bh_q,
    accessible_fraction_range)], by = "feature_id")
top_ids <- dynamic[robust_dynamic == TRUE][order(minimum_adjacent_bh_q,
  -accessible_fraction_range), head(feature_id, 60L)]
plot_profiles <- plot_profiles[feature_id %in% top_ids]
peak_order <- dynamic[feature_id %in% top_ids][order(maximum_region,
  minimum_adjacent_bh_q, feature_id), feature_id]
plot_profiles[, `:=`(
  feature_factor = factor(feature_id, levels = rev(peak_order)),
  region_factor = factor(region_label,
    levels = region_spec$region_label),
  accessible_fraction_display = fifelse(
    accessible_fraction < 0.05, NA_real_, accessible_fraction))]
p_bdev <- ggplot(plot_profiles,
  aes(region_factor, feature_factor, fill = accessible_fraction_display)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradientn(
    colours = c("#deebf7", "#9ecae1", "#4292c6", "#756bb1", "#54278f"),
    limits = c(0.05, max(plot_profiles$accessible_fraction, na.rm = TRUE)),
    breaks = c(0.05, 0.10, 0.20, 0.30, 0.40),
    labels = scales::label_percent(accuracy = 1),
    na.value = "white", name = "Accessible cells\n(white: <5%)") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(
    angle = 35, hjust = 1), axis.title = element_blank()) +
  labs(title = NULL, subtitle = NULL)

tcf3_coverage <- fread(cmd = paste("zcat", shQuote(tcf3_coverage_path)))
stopifnot(all(c("chromosome", "start_1based", "end_1based", "condition",
  "normalized_fragment_coverage") %in% names(tcf3_coverage)),
  setequal(unique(tcf3_coverage$condition), c("control", "TCF3_depleted")))
tcf3_peak_coverage <- rbindlist(lapply(top_ids, function(id) {
  peak <- dynamic[feature_id == id]
  overlap <- tcf3_coverage[
    chromosome == peak$chromosome & end_1based >= peak$peak_start_1based &
      start_1based <= peak$peak_end_1based]
  stopifnot(nrow(overlap) > 0L)
  overlap[, overlap_bp := pmax(0L,
    pmin(end_1based, peak$peak_end_1based) -
      pmax(start_1based, peak$peak_start_1based) + 1L)]
  overlap[, .(mean_fragment_coverage_cpm =
    weighted.mean(normalized_fragment_coverage, overlap_bp)),
    by = condition][, feature_id := id]
}), use.names = TRUE)
tcf3_peak_coverage[, `:=`(
  condition_label = factor(condition,
    levels = c("control", "TCF3_depleted"),
    labels = c("Control/WT", "TCF3 depleted")),
  feature_factor = factor(feature_id, levels = rev(peak_order)))]
write_tsv(tcf3_peak_coverage[, .(feature_id, condition,
  mean_fragment_coverage_cpm)], file.path(stage_dir,
  "top_dynamic_peaks_tcf3_condition_mean_coverage.tsv"))

p_tcf3 <- ggplot(tcf3_peak_coverage,
  aes(condition_label, feature_factor, fill = mean_fragment_coverage_cpm)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradientn(
    colours = c("#fff7ec", "#fdbb84", "#e34a33", "#7f0000"),
    trans = "sqrt", name = "Mean fragment\ncoverage CPM") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(
    angle = 35, hjust = 1), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), axis.title = element_blank())

p_heat <- p_bdev + p_tcf3 +
  plot_layout(widths = c(5, 2), guides = "keep") +
  plot_annotation(
    title = paste0("Top ", length(top_ids),
      " robust dynamic Ebf1-window peaks across five Bdev states"),
    subtitle = paste("Bdev accessible fraction at left; existing bulk",
      "E2A/TCF3 control and depleted coverage at right"))
ggsave(file.path(stage_dir, "figures",
  "top_robust_dynamic_peak_accessibility_heatmap.png"), p_heat,
  width = 11, height = max(6, 0.18 * length(top_ids) + 2.5),
  dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "top_robust_dynamic_peak_accessibility_heatmap.pdf"), p_heat,
  width = 11, height = max(6, 0.18 * length(top_ids) + 2.5),
  device = cairo_pdf, bg = "white")

p_locus <- ggplot(dynamic,
  aes(midpoint_relative_to_tss_bp / 1e6, accessible_fraction_range,
    color = robust_dynamic)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35") +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(values = c(`TRUE` = "#7A3E9D", `FALSE` = "grey72"),
    name = "Robust dynamic") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom") +
  labs(x = "Peak midpoint relative to Ebf1 TSS (Mb)",
    y = "Five-state accessible-fraction range",
    title = "Dynamic Bdev accessibility across Ebf1 TSS ±2 Mb")
ggsave(file.path(stage_dir, "figures",
  "dynamic_peak_range_by_locus_position.png"), p_locus,
  width = 10, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "dynamic_peak_range_by_locus_position.pdf"), p_locus,
  width = 10, height = 5.5, device = cairo_pdf, bg = "white")

summary <- list(
  total_peaks = nrow(dynamic),
  statistically_dynamic = dynamic[, sum(statistically_dynamic)],
  robust_dynamic = dynamic[, sum(robust_dynamic)],
  transition_counts = split(transition_counts, transition_counts$transition_id))
write_json(list(status = "complete_unpromoted_report_only",
  script_version = script_version,
  bundle_fingerprint = expected_fp,
  pseudotime_source = "exact preserved Bdev metadata pt",
  region_basis = "current shared WT/KO-supported Flt3L bins 14-25 and IL-7 bins 22-30",
  local_flanks = "four bins, matching the dual-state width",
  test = "native scWorkbench de_wilcox on four ordered adjacent-state contrasts",
  dynamic_rule = "BH q < 0.05 in at least one native adjacent contrast",
  robust_rule = "native BH q < 0.05, absolute accessible-fraction difference >=0.05, and concordant normalized-effect direction",
  heatmap_scale = "accessible fraction below 0.05 is white; sequential blue-to-purple scale begins at 0.05",
  tcf3_heatmap_columns = paste("weighted mean of existing 50-bp normalized",
    "fragment coverage across each peak; control and TCF3-depleted",
    "condition means over three replicates"),
  summary = summary, bundle_mutated = FALSE, promotion = FALSE),
  file.path(stage_dir, "diagnostics.json"))

region_lines <- apply(region_spec, 1L, function(x) paste0("- `",
  x[["region_label"]], "`: bins ", x[["first_bin"]], "-",
  x[["last_bin"]], ", pt [",
  sprintf("%.3f", as.numeric(x[["pt_left"]])), ", ",
  sprintf("%.3f", as.numeric(x[["pt_right"]])), "); ",
  x[["bdev_cells"]], " Bdev cells."))
transition_lines <- apply(transition_counts, 1L, function(x) paste0("- `",
  x[["transition_label"]], "`: ", x[["significant_bh05"]],
  " BH-significant; ", x[["robust_bh05_pctdiff05"]],
  " robust dynamic peaks."))
writeLines(c(
  "# Five-state Bdev dynamic peaks within Ebf1 TSS ±2 Mb", "",
  "## Suggested ordered Bdev regions", "", region_lines, "",
  "These states are based on the current shared WT/KO-supported analysis domains, not the broader raw WT-enrichment calls. Four-bin flanks match the dual-state width and provide local controls.", "",
  "## Dynamic peak screen", "", transition_lines, "",
  paste0("Across 251 peaks, ", summary$statistically_dynamic,
    " are BH-significant in at least one adjacent-state comparison and ",
    summary$robust_dynamic,
    " also have an absolute accessible-fraction change of at least 0.05 with concordant normalized direction."), "",
  "All four adjacent contrasts use native scWorkbench full-matrix log1p-CP10K and `de_wilcox`, with BH adjustment across the 251 peaks in each contrast.", "",
  "The heatmap is thresholded for readability: accessible fractions below 5% are white, and values from 5% upward use a sequential blue-to-purple scale.", "",
  "Two separately scaled columns at right show the existing bulk E2A/TCF3 control and depleted mean fragment-coverage CPM over the same peaks. These are bulk coverage measurements, not accessible-cell fractions.", "",
  "- [Dynamic peak summary](ebf1_pm2mb_dynamic_peak_summary.tsv)",
  "- [Robust dynamic peaks](robust_dynamic_ebf1_pm2mb_peaks.tsv)",
  "- [All adjacent-region DA](all_adjacent_region_ebf1_pm2mb_da.tsv)",
  "- [Five-region accessibility profiles](bdev_five_region_peak_accessibility_profiles.tsv)",
  "- [TCF3 condition mean coverage over displayed peaks](top_dynamic_peaks_tcf3_condition_mean_coverage.tsv)",
  "- [Top dynamic heatmap](figures/top_robust_dynamic_peak_accessibility_heatmap.png)",
  "- [Locus dynamic-range plot](figures/dynamic_peak_range_by_locus_position.png)"
), file.path(stage_dir, "REPORT.md"))

capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed five-state Bdev Ebf1 dynamic peak report at", output_dir, "\n")
