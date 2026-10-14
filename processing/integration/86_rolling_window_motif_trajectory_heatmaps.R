#!/usr/bin/env Rscript

# Bdev-defined motif analysis over the reviewed centered three-bin windows.
# The Bdev and mapping bundles are opened read-only. Query cells do not affect
# peak/motif eligibility, chromVAR calibration, motif selection, scaling, or
# row order.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(11L, 12L)) stop(paste(
  "Usage: 86_rolling_window_motif_trajectory_heatmaps.R",
  "<bdev_bundle> <mapping_bundle> <joint_peak_matrix>",
  "<window_support_report> <query_support_tsv_gz> <bdev_rdata>",
  "<output_dir> <bdev_fingerprint> <mapping_fingerprint> <top_n> <seed>",
  "[rolling_3bin|single_bin|fine_rolling_half_step]"),
  call. = FALSE)
bdev_dir <- normalizePath(args[[1L]], mustWork = TRUE)
mapping_dir <- normalizePath(args[[2L]], mustWork = TRUE)
joint_dir <- normalizePath(args[[3L]], mustWork = TRUE)
support_dir <- normalizePath(args[[4L]], mustWork = TRUE)
query_path <- normalizePath(args[[5L]], mustWork = TRUE)
bdev_rdata <- normalizePath(args[[6L]], mustWork = TRUE)
output_dir <- normalizePath(args[[7L]], mustWork = FALSE)
bdev_fp <- args[[8L]]
mapping_fp <- args[[9L]]
top_n <- as.integer(args[[10L]])
seed <- as.integer(args[[11L]])
analysis_mode <- if (length(args) == 12L) args[[12L]] else "rolling_3bin"
stopifnot(top_n >= 5L, is.finite(seed))
stopifnot(analysis_mode %in% c("rolling_3bin", "single_bin",
  "fine_rolling_half_step"))
figure_tag <- switch(analysis_mode, rolling_3bin = "rolling_3bin",
  single_bin = "single_bin", fine_rolling_half_step =
    "fine_rolling_half_step")
window_description <- switch(analysis_mode,
  rolling_3bin = "three source bins centered on displayed middle bin",
  single_bin = "one source bin per displayed trajectory position",
  fine_rolling_half_step = paste("three half-width bins (1.5 original-bin total width)",
    "with focal centers spaced at half the original bin width"))
plot_method_label <- switch(analysis_mode,
  rolling_3bin = "rolling-window", single_bin = "single-bin",
  fine_rolling_half_step = "fine-resolution rolling-window")
plot_caption_window <- switch(analysis_mode,
  rolling_3bin = "Each column aggregates three centered source bins.",
  single_bin = "Each column aggregates exactly one source bin.",
  fine_rolling_half_step = paste("Each column aggregates three half-width bins",
    "and centers advance by one half-width bin."))
x_axis_label <- if (analysis_mode == "fine_rolling_half_step")
  "Focal Bdev pseudotime" else "Middle-bin Bdev pseudotime"
if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(BPCells)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(chromVAR)
  library(SummarizedExperiment)
  library(Signac)
  library(SeuratObject)
  library(jsonlite)
})
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json_file <- function(x, path) jsonlite::write_json(x, path,
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(..., collapse = " "), "\n")

bdev <- wb_open_bundle(bdev_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
mapping <- wb_open_bundle(mapping_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(mapping), silent = TRUE), add = TRUE)
bdev_manifest <- wb_manifest(bdev)
mapping_manifest <- wb_manifest(mapping)
stopifnot(identical(bdev_manifest$bundle_id, "bdev_reference"),
  identical(bdev_manifest$current_state$fingerprint, bdev_fp),
  identical(mapping_manifest$bundle_id, "bdev_integrated_ko_mapping_v1"),
  identical(mapping_manifest$current_state$fingerprint, mapping_fp))

regions <- fread(file.path(support_dir, "included_regions.tsv"))
if (analysis_mode == "rolling_3bin") {
  windows <- fread(file.path(support_dir,
    "rolling_3bin_window_cell_counts.tsv"))
  windows <- windows[region_id %in% regions$region_id]
} else if (analysis_mode == "single_bin") {
  source_bins <- fread(file.path(support_dir, "source_bin_cell_counts.tsv"))
  windows <- source_bins[retained_run == TRUE &
      region_id %in% regions$region_id,
    .(target, wt_peak, region_id, center_bin = bin,
      window_first_bin = bin, window_last_bin = bin,
      window_pt_left = bin_left, window_pt_right = bin_right,
      center_pt = bin_midpoint, wt_cells, ko_cells,
      center_bin_wt_cells = wt_cells,
      center_bin_ko_cells = ko_cells,
      ko_failed_source_bins_in_window = as.integer(!ko_pass),
      center_bin_ko_pass = ko_pass)]
} else {
  source_bins <- fread(file.path(support_dir, "source_bin_cell_counts.tsv"))
  original_width <- median(source_bins$bin_right - source_bins$bin_left)
  stopifnot(is.finite(original_width), original_width > 0)
  half_step <- original_width / 2
  # Three fine bins, each half the width of an original source bin.
  full_width <- original_width * 1.5
  windows <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
    reg <- regions[i]
    n_positions <- as.integer(round(
      (reg$center_pt_right - reg$center_pt_left) / half_step)) + 1L
    centers <- reg$center_pt_left + seq.int(0L, n_positions - 1L) *
      half_step
    stopifnot(abs(tail(centers, 1L) - reg$center_pt_right) < 1e-8)
    data.table(target = reg$target, wt_peak = reg$wt_peak,
      region_id = reg$region_id, center_bin = seq_len(n_positions),
      window_first_bin = NA_integer_, window_last_bin = NA_integer_,
      window_pt_left = centers - full_width / 2,
      window_pt_right = centers + full_width / 2,
      center_pt = centers, wt_cells = NA_integer_, ko_cells = NA_integer_,
      center_bin_wt_cells = NA_integer_,
      center_bin_ko_cells = NA_integer_,
      ko_failed_source_bins_in_window = NA_integer_,
      center_bin_ko_pass = NA)
  }))
}
query <- fread(cmd = paste("zcat", shQuote(query_path)))[
  current_accepted == TRUE & trajectory_support >= 0.25 &
    is.finite(mapped_pt)]
stopifnot(nrow(regions) == 2L,
  setequal(regions$region_id,
    c("flt3l_wt_peak1_run1", "il7_wt_peak1_run1")),
  !anyDuplicated(windows[, .(region_id, center_bin)]),
  !anyDuplicated(query$cell_id))

bdev_cells <- wb_cells(bdev, "set:bdev_reference")
bdev_meta <- as.data.table(wb_meta(bdev, cells = bdev_cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt")))
setnames(bdev_meta, names(bdev_meta)[[1L]], "source_cell_id")
stopifnot(identical(as.character(bdev_meta$source_cell_id), bdev_cells))
bdev_meta[, cell_id := paste0("bdev::", source_cell_id)]

log_event("Loading exact source Bdev JASPAR motif matches and GC bias")
source_env <- new.env(parent = emptyenv())
con <- gzfile(bdev_rdata, "rb")
loaded <- load(con, envir = source_env)
close(con)
stopifnot("Bdev" %in% loaded)
source_bdev <- source_env$Bdev
motif_all <- Signac::GetMotifData(source_bdev[["peaks"]], slot = "data")
peak_meta <- as.data.table(source_bdev[["peaks"]][[]],
  keep.rownames = "feature_id")
motif_ids <- colnames(motif_all)
motif_names <- tryCatch(Signac::ConvertMotifID(source_bdev[["peaks"]],
  id = motif_ids), error = function(e) setNames(motif_ids, motif_ids))
if (is.null(names(motif_names))) names(motif_names) <- motif_ids

joint <- BPCells::open_matrix_dir(joint_dir)
peak_ids <- rownames(joint)
joint_cells <- colnames(joint)
stopifnot(length(peak_ids) == 97788L, !anyDuplicated(peak_ids),
  !anyDuplicated(joint_cells), all(peak_ids %in% rownames(motif_all)),
  all(peak_ids %in% peak_meta$feature_id))
motif_all <- motif_all[peak_ids, , drop = FALSE]
peak_bias <- peak_meta$GC.percent[match(peak_ids, peak_meta$feature_id)]
stopifnot(all(is.finite(peak_bias)))
rm(source_bdev, source_env, peak_meta)
invisible(gc())

all_counts <- list()
all_variability <- list()
all_selected <- list()
all_values <- list()
all_eligible <- list()
all_density <- list()
region_diagnostics <- list()
group_levels <- c("Bdev", "WT", "KO")

for (ri in seq_len(nrow(regions))) {
  reg <- regions[ri]
  region_key <- reg$region_id
  target <- reg$target
  win <- windows[region_id == region_key]
  setorder(win, center_bin)
  boundary <- reg$ko_support_boundary_pt
  boundary_definition <- "original middle source-bin KO count"
  if (analysis_mode == "fine_rolling_half_step") {
    fine_bin_width <- median(diff(win$center_pt))
    stopifnot(is.finite(fine_bin_width), fine_bin_width > 0)
    fine_ko_counts <- vapply(win$center_pt, function(focal_pt) {
      query[knockout_target == target & genotype == "KO" &
        mapped_pt >= focal_pt - fine_bin_width / 2 &
        mapped_pt < focal_pt + fine_bin_width / 2, .N]
    }, integer(1L))
    win[, `:=`(center_bin_ko_cells = fine_ko_counts,
      center_bin_ko_pass = fine_ko_counts >= 10L)]
    last_supported <- max(which(win$center_bin_ko_pass))
    stopifnot(last_supported < nrow(win),
      all(!win$center_bin_ko_pass[seq.int(last_supported + 1L, nrow(win))]))
    boundary <- mean(win$center_pt[c(last_supported, last_supported + 1L)])
    boundary_definition <- paste(
      "half-width middle-bin KO count; boundary between last >=10 bin",
      "and first subsequent <10 tail bin")
  }
  log_event("Preparing", region_key, "with", nrow(win), analysis_mode,
    "trajectory positions")

  membership <- rbindlist(lapply(seq_len(nrow(win)), function(wi) {
    w <- win[wi]
    bd <- bdev_meta[is.finite(pt) & pt >= w$window_pt_left &
      pt < w$window_pt_right, .(cell_id, group = "Bdev")]
    qu <- query[knockout_target == target & mapped_pt >= w$window_pt_left &
      mapped_pt < w$window_pt_right,
      .(cell_id, group = as.character(genotype))]
    ans <- rbindlist(list(bd, qu))
    ans[, `:=`(center_bin = w$center_bin,
      center_pt = w$center_pt,
      column_id = paste(group, sprintf("%03d", w$center_bin), sep = "::"))]
    ans
  }))
  stopifnot(all(membership$cell_id %in% joint_cells),
    all(membership$group %in% group_levels),
    !anyDuplicated(membership[, .(cell_id, group, center_bin)]))

  columns <- CJ(group = group_levels, center_bin = win$center_bin,
    unique = TRUE)
  columns <- merge(columns,
    win[, .(center_bin, window_first_bin, window_last_bin,
      window_pt_left, window_pt_right, center_pt,
      center_bin_ko_cells, center_bin_ko_pass)],
    by = "center_bin", all.x = TRUE, sort = FALSE)
  columns[, column_id := paste(group, sprintf("%03d", center_bin),
    sep = "::")]
  observed <- membership[, .(n_cells = .N), by = .(group, center_bin)]
  columns <- merge(columns, observed, by = c("group", "center_bin"),
    all.x = TRUE, sort = FALSE)
  columns[is.na(n_cells), n_cells := 0L]
  columns[, group_order := match(group, group_levels)]
  setorder(columns, group_order, center_bin)
  columns[, `:=`(group_order = NULL, region_id = region_key)]
  if (analysis_mode != "fine_rolling_half_step") {
    stopifnot(identical(columns[group == "WT", n_cells], win$wt_cells),
      identical(columns[group == "KO", n_cells], win$ko_cells))
  }
  stopifnot(all(columns[group == "Bdev", n_cells] > 0L),
    all(columns[group == "WT", n_cells] > 0L))
  all_counts[[region_key]] <- copy(columns)

  analysis_columns <- columns[n_cells > 0L]
  membership <- membership[column_id %in% analysis_columns$column_id]
  design <- sparseMatrix(i = match(membership$cell_id, joint_cells),
    j = match(membership$column_id, analysis_columns$column_id), x = 1,
    dims = c(length(joint_cells), nrow(analysis_columns)),
    dimnames = list(joint_cells, analysis_columns$column_id))
  agg <- as(joint %*% design, "dgCMatrix")
  stopifnot(identical(colnames(agg), analysis_columns$column_id))
  bdev_cols <- which(analysis_columns$group == "Bdev")
  bdev_agg <- agg[, bdev_cols, drop = FALSE]
  peak_keep <- Matrix::rowSums(bdev_agg) >= 10
  motif <- motif_all[peak_keep, , drop = FALSE]
  motif_peak_n <- Matrix::colSums(motif != 0)
  motif_keep <- motif_peak_n >= 10L
  motif <- motif[, motif_keep, drop = FALSE]
  motif_peak_n <- motif_peak_n[motif_keep]
  counts_kept <- agg[peak_keep, , drop = FALSE]
  bdev_counts <- counts_kept[, bdev_cols, drop = FALSE]
  stopifnot(nrow(counts_kept) > 10000L, ncol(motif) > top_n)

  log_event("Computing Bdev-defined chromVAR deviations for", region_key)
  set.seed(seed + ri)
  bg <- chromVAR::getBackgroundPeaks(bdev_counts,
    bias = peak_bias[peak_keep], niterations = 50L, w = 0.1, bs = 50L)
  expectation <- chromVAR::computeExpectations(bdev_counts)
  dev_obj <- chromVAR::computeDeviations(counts_kept, motif,
    background_peaks = bg, expectation = expectation)
  raw <- SummarizedExperiment::assay(dev_obj, "deviations")
  stopifnot(identical(colnames(raw), analysis_columns$column_id),
    all(is.finite(raw)))

  bdev_raw <- raw[, bdev_cols, drop = FALSE]
  center <- rowMeans(bdev_raw)
  scale <- apply(bdev_raw, 1L, sd)
  eligible <- is.finite(scale) & scale > 0
  scaled <- sweep(sweep(raw[eligible, , drop = FALSE], 1L,
    center[eligible], "-"), 1L, scale[eligible], "/")
  variability <- data.table(region_id = region_key,
    motif_id = rownames(raw),
    motif_name = as.character(motif_names[rownames(raw)]),
    bdev_temporal_sd = apply(bdev_raw, 1L, sd),
    bdev_temporal_range = apply(bdev_raw, 1L, function(x) diff(range(x))),
    bdev_accessible_peak_count = as.integer(motif_peak_n[rownames(raw)]))
  variability[is.na(motif_name) | motif_name == "", motif_name := motif_id]
  setorder(variability, -bdev_temporal_sd, -bdev_temporal_range, motif_id)
  variability[, `:=`(bdev_variability_rank = seq_len(.N),
    selected = seq_len(.N) <= top_n)]
  selected <- variability[selected == TRUE & motif_id %in% rownames(scaled)]
  stopifnot(nrow(selected) == top_n)
  all_variability[[region_key]] <- variability
  all_selected[[region_key]] <- selected

  long <- as.data.table(as.table(scaled))
  setnames(long, c("motif_id", "column_id", "bdev_scaled_deviation"))
  long[, `:=`(motif_id = as.character(motif_id),
    column_id = as.character(column_id))]
  raw_long <- as.data.table(as.table(raw[eligible, , drop = FALSE]))
  setnames(raw_long, c("motif_id", "column_id", "raw_chromvar_deviation"))
  raw_long[, `:=`(motif_id = as.character(motif_id),
    column_id = as.character(column_id))]
  long <- merge(long, raw_long, by = c("motif_id", "column_id"),
    sort = FALSE)
  full_grid <- CJ(motif_id = rownames(scaled),
    column_id = columns$column_id, unique = TRUE)
  long <- merge(full_grid, long, by = c("motif_id", "column_id"),
    all.x = TRUE, sort = FALSE)
  long <- merge(long, columns[, .(column_id, group, center_bin, center_pt,
    n_cells, center_bin_ko_cells, center_bin_ko_pass)],
    by = "column_id", sort = FALSE)
  long <- merge(long, variability[, .(motif_id, motif_name,
    bdev_temporal_sd, bdev_variability_rank)], by = "motif_id",
    sort = FALSE)
  long[, region_id := region_key]
  all_eligible[[region_key]] <- copy(long)

  heat <- long[motif_id %in% selected$motif_id]
  labels <- selected[, .(motif_id,
    motif_label = make.unique(sprintf("%s (%s)", motif_name, motif_id)))]
  heat <- merge(heat, labels, by = "motif_id", sort = FALSE)
  bdev_heat <- dcast(heat[group == "Bdev"], motif_id ~ center_bin,
    value.var = "bdev_scaled_deviation")
  row_order <- bdev_heat$motif_id[hclust(dist(as.matrix(
    bdev_heat[, -1L])))$order]
  wt <- heat[group == "WT", .(motif_id, motif_label, center_bin,
    center_pt, wt = bdev_scaled_deviation)]
  ko <- heat[group == "KO", .(motif_id, center_bin,
    ko = bdev_scaled_deviation)]
  diff <- merge(wt, ko, by = c("motif_id", "center_bin"), sort = FALSE)
  diff[, `:=`(group = "WT - KO", n_cells = NA_integer_,
    bdev_scaled_deviation = wt - ko,
    raw_chromvar_deviation = NA_real_)]
  diff[, c("wt", "ko") := NULL]
  heat <- rbindlist(list(heat, diff), fill = TRUE)
  heat[, region_id := region_key]
  heat[, `:=`(display_bdev_scaled_deviation = bdev_scaled_deviation,
    ko_support_displayed =
      !(group %in% c("KO", "WT - KO") & center_pt > boundary))]
  heat[ko_support_displayed == FALSE,
    display_bdev_scaled_deviation := NA_real_]
  group_labels <- if (target == "Flt3l")
    c(Bdev = "Bdev", WT = "Flt3L WT", KO = "FLT3L KO",
      `WT - KO` = "WT - KO") else
    c(Bdev = "Bdev", WT = "IL-7 WT", KO = "IL-7 KO",
      `WT - KO` = "WT - KO")
  heat[, group_label := factor(group_labels[group],
    levels = unname(group_labels))]
  ordered_labels <- labels$motif_label[match(row_order, labels$motif_id)]
  heat[, motif_label := factor(motif_label, levels = rev(ordered_labels))]
  all_values[[region_key]] <- copy(heat)

  density_wide <- dcast(columns,
    center_bin + center_pt ~ group, value.var = "n_cells", fill = 0L)
  bdev_total <- sum(is.finite(bdev_meta$pt))
  wt_total <- query[knockout_target == target & genotype == "WT", .N]
  ko_total <- query[knockout_target == target & genotype == "KO", .N]
  stopifnot(bdev_total > 0L, wt_total > 0L, ko_total > 0L,
    all(density_wide$Bdev > 0L))
  density_wide[, `:=`(
    bdev_fraction = Bdev / bdev_total,
    wt_fraction = WT / wt_total,
    ko_fraction = KO / ko_total)]
  density_wide[, `:=`(
    wt_relative_to_bdev = wt_fraction / bdev_fraction,
    ko_relative_to_bdev = ko_fraction / bdev_fraction)]
  density_plot <- rbindlist(list(
    density_wide[, .(center_bin, center_pt, group = "Bdev",
      series = "Bdev", relative_density = 1)],
    density_wide[, .(center_bin, center_pt, group = "WT",
      series = "WT", relative_density = wt_relative_to_bdev)],
    density_wide[, .(center_bin, center_pt, group = "KO",
      series = "KO", relative_density = ko_relative_to_bdev)],
    density_wide[, .(center_bin, center_pt, group = "WT - KO",
      series = "WT", relative_density = wt_relative_to_bdev)],
    density_wide[, .(center_bin, center_pt, group = "WT - KO",
      series = "KO", relative_density = ko_relative_to_bdev)]))
  density_plot[, `:=`(region_id = region_key,
    group_label = factor(group_labels[group], levels = unname(group_labels)),
    series = factor(series, levels = c("Bdev", "WT", "KO")))]
  all_density[[region_key]] <- merge(density_plot,
    density_wide[, .(center_bin, bdev_cells = Bdev, wt_cells = WT,
      ko_cells = KO, bdev_fraction, wt_fraction, ko_fraction)],
    by = "center_bin", all.x = TRUE, sort = FALSE)

  boundary_data <- data.table(group_label = factor(unname(group_labels),
    levels = unname(group_labels)), boundary = boundary)
  p_heat <- ggplot(heat, aes(center_pt, motif_label,
      fill = display_bdev_scaled_deviation)) +
    geom_tile() +
    geom_vline(data = boundary_data, aes(xintercept = boundary),
      inherit.aes = FALSE, linetype = "dashed", colour = "#D55E00",
      linewidth = 0.55) +
    facet_grid(. ~ group_label) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish,
      na.value = "grey85", name = "Deviation\n(Bdev SD)") +
    scale_x_continuous(expand = expansion(mult = 0)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6.5),
      plot.title = element_text(face = "bold")) +
      labs(x = x_axis_label, y = NULL,
      caption = paste0(plot_caption_window,
        " Orange line: KO middle-bin support boundary (10 cells); KO and WT-KO motifs after it are masked. WT-KO is descriptive; one sample per genotype."))
  p_density <- ggplot(density_plot,
      aes(center_pt, relative_density, fill = series, colour = series,
        group = series)) +
    geom_hline(yintercept = 1, colour = "grey35", linetype = "dotted",
      linewidth = 0.4) +
    geom_area(position = "identity", alpha = 0.32, linewidth = 0) +
    geom_line(linewidth = 0.65) +
    geom_vline(data = boundary_data, aes(xintercept = boundary),
      inherit.aes = FALSE, linetype = "dashed", colour = "#D55E00",
      linewidth = 0.55) +
    facet_grid(. ~ group_label) +
    scale_fill_manual(values = c(Bdev = "grey55", WT = "#0072B2",
      KO = "#D55E00"), drop = FALSE) +
    scale_colour_manual(values = c(Bdev = "grey35", WT = "#0072B2",
      KO = "#D55E00"), drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal(base_size = 8.5) +
    theme(panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.title.x = element_blank(), legend.position = "top",
      plot.margin = margin(5.5, 5.5, 0, 5.5)) +
    labs(y = "Density / Bdev", fill = NULL, colour = NULL)
  p <- p_density / p_heat + patchwork::plot_layout(heights = c(1.8, 8.2)) +
    patchwork::plot_annotation(
      title = paste0(if (target == "Flt3l") "Flt3L" else "IL-7",
        " WT-enriched region: ", plot_method_label, " motif trajectory"),
      subtitle = paste0("Top ", top_n,
        " temporally variable motifs selected and scaled using Bdev only"))
  stub <- if (target == "Flt3l") "flt3l_wt_peak1" else "il7_wt_peak1"
  ggsave(file.path(stage_dir, "figures", paste0(stub,
    "_", figure_tag, "_motif_heatmap.png")), p, width = 18, height = 12,
    dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(stub,
    "_", figure_tag, "_motif_heatmap.pdf")), p, width = 18, height = 12,
    device = cairo_pdf, bg = "white")

  region_diagnostics[[region_key]] <- list(target = target,
    source_bins = c(reg$source_first_bin, reg$source_last_bin),
    trajectory_positions = nrow(win), ko_support_boundary_pt = boundary,
    ko_support_boundary_definition = boundary_definition,
    ko_and_comparison_motifs_after_boundary = "masked_from_display",
    cells_by_group = as.list(columns[, sum(n_cells), by = group][,
      setNames(V1, group)]), bdev_accessible_peaks = sum(peak_keep),
    eligible_motifs = nrow(variability), selected_motifs = top_n)
  rm(design, agg, bdev_agg, counts_kept, bdev_counts, dev_obj, raw,
    bdev_raw, scaled, long, raw_long, heat)
  invisible(gc())
}

write_tsv(rbindlist(all_counts), file.path(stage_dir,
  paste0(figure_tag, "_cell_counts.tsv")))
write_tsv(rbindlist(all_variability), file.path(stage_dir,
  "all_motif_bdev_temporal_variability.tsv"))
write_tsv(rbindlist(all_selected), file.path(stage_dir,
  "selected_bdev_variable_motifs.tsv"))
write_tsv(rbindlist(all_density), file.path(stage_dir,
  paste0(figure_tag, "_cell_density_relative_to_bdev.tsv")))
fwrite(rbindlist(all_values), file.path(stage_dir,
  paste0("selected_motif_", figure_tag, "_trajectory_values.tsv.gz")), sep = "\t",
  na = "NA", compress = "gzip")
fwrite(rbindlist(all_eligible), file.path(stage_dir,
  paste0("all_eligible_motif_", figure_tag, "_trajectory_values.tsv.gz")), sep = "\t",
  na = "NA", compress = "gzip")

diagnostics <- list(status = "complete_unpromoted_report_only",
  bdev_bundle_id = bdev_manifest$bundle_id,
  bdev_bundle_fingerprint = bdev_fp,
  mapping_bundle_id = mapping_manifest$bundle_id,
  mapping_bundle_fingerprint = mapping_fp,
  source_window_report = support_dir, trajectory_support_cutoff = 0.25,
  regions = region_diagnostics, motif_source =
    "exact read-only Bdev JASPAR match matrix",
  selection_calibration_scaling_order = "Bdev only",
  query_role = "display and descriptive WT-KO comparison only",
  analysis_mode = analysis_mode, trajectory_window = window_description,
  density_track = "within-group accepted-cell fraction divided by Bdev defined-pt cell fraction at the same trajectory position",
  biological_replicates_per_genotype = 1L,
  inference = "descriptive", bundle_mutated = FALSE, promotion = FALSE)
write_json_file(diagnostics, file.path(stage_dir, "diagnostics.json"))
report_title <- switch(analysis_mode,
  rolling_3bin = "# Rolling-window motif trajectories",
  single_bin = "# Single-bin motif trajectories",
  fine_rolling_half_step = "# Fine-resolution rolling-window motif trajectories")
display_explanation <- switch(analysis_mode,
  rolling_3bin = "Every displayed trajectory column is a pseudobulk over three centered source bins and is positioned at the middle-bin Bdev pseudotime.",
  single_bin = "Every displayed trajectory column is a pseudobulk from exactly one source bin and is positioned at that bin's Bdev pseudotime midpoint.",
  fine_rolling_half_step = "The original source-bin width is halved. Every displayed column is a centered pseudobulk over three of these fine bins (1.5 original-bin widths total), and centers advance by one fine bin for two-fold finer trajectory resolution.")
dependence_note <- if (analysis_mode == "single_bin")
  "Because there is one biological sample per genotype, the comparison is descriptive rather than replicate-level inference." else
  "Because there is one biological sample per genotype and adjacent rolling windows overlap, the comparison is descriptive rather than replicate-level inference."
writeLines(c(report_title, "",
  "The agreed Flt3L WT peak 1 and IL-7 WT peak 1 regions are shown. IL-7 peak 2 is excluded because it has no KO-supported middle bin.", "",
  paste(display_explanation, if (analysis_mode == "fine_rolling_half_step")
    "The orange vertical line is halfway between the last half-width middle bin with at least 10 KO cells and the first subsequent unsupported tail bin. KO and WT-minus-KO motif tiles after this line are masked; their computed values remain in the exported table with an explicit display-support flag."
  else "The orange vertical line is the reviewed boundary after the last bin with at least 10 KO cells. KO and WT-minus-KO motif tiles after this line are masked; their computed values remain in the exported table with an explicit display-support flag."), "",
  paste0("For each region, Bdev alone determines accessible peaks, chromVAR background and expectation, the top ", top_n, " motif selection, centering, scaling, and row order. WT and KO cells are evaluated only afterward. The fourth panel is the aligned descriptive WT-minus-KO difference in Bdev temporal-SD units."), "",
  "The filled track above each heatmap is cell-density enrichment relative to Bdev: the within-group fraction of accepted mapped cells at that trajectory position divided by the fraction of Bdev defined-pt cells at the same position. One is Bdev-like occupancy; the comparison panel overlays WT and KO.", "",
  dependence_note, "",
  paste0("- [Flt3L figure](figures/flt3l_wt_peak1_", figure_tag,
    "_motif_heatmap.png)"),
  paste0("- [IL-7 figure](figures/il7_wt_peak1_", figure_tag,
    "_motif_heatmap.png)")),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
wb_close_bundle(bdev); bdev <- NULL
wb_close_bundle(mapping); mapping <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed", analysis_mode, "motif analysis at", output_dir, "\n")
