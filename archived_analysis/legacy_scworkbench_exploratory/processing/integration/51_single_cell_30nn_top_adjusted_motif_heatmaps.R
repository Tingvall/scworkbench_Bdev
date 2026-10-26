#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) stop(paste(
  "Usage: 51_single_cell_30nn_top_adjusted_motif_heatmaps.R",
  "<bdev_bundle> <mapping_bundle> <joint_peak_matrix_dir>",
  "<region_report> <bdev_rdata> <motif_selection_report> <output_dir>",
  "<expected_bdev_fingerprint> <expected_mapping_fingerprint>",
  "<neighbor_count> <maximum_neighbor_radius_pt> <seed>"))
bdev_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
mapping_bundle_dir <- normalizePath(args[[2L]], mustWork = TRUE)
joint_matrix_dir <- normalizePath(args[[3L]], mustWork = TRUE)
region_report <- normalizePath(args[[4L]], mustWork = TRUE)
bdev_rdata_path <- normalizePath(args[[5L]], mustWork = TRUE)
selection_report <- normalizePath(args[[6L]], mustWork = TRUE)
output_dir <- normalizePath(args[[7L]], mustWork = FALSE)
expected_bdev_fingerprint <- args[[8L]]
expected_mapping_fingerprint <- args[[9L]]
neighbor_count <- as.integer(args[[10L]])
maximum_neighbor_radius_pt <- as.numeric(args[[11L]])
seed <- as.integer(args[[12L]])
stopifnot(neighbor_count >= 3L, is.finite(maximum_neighbor_radius_pt),
  maximum_neighbor_radius_pt > 0, is.finite(seed))
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
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
write_tsv <- function(x, path) data.table::fwrite(x, path, sep = "\t",
  na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

smooth_30nn <- function(deviation, pt, target_pt, k,
    restrict_to_observed_support = TRUE) {
  stopifnot(ncol(deviation) == length(pt), length(pt) >= k,
    all(is.finite(pt)), all(is.finite(target_pt)))
  ordering <- order(pt)
  sorted_pt <- pt[ordering]
  sorted_deviation <- deviation[, ordering, drop = FALSE]
  result <- matrix(NA_real_, nrow = nrow(deviation),
    ncol = length(target_pt),
    dimnames = list(rownames(deviation), as.character(seq_along(target_pt))))
  radius <- rep(NA_real_, length(target_pt))
  supported <- if (restrict_to_observed_support)
    target_pt >= min(sorted_pt) & target_pt <= max(sorted_pt) else
    rep(TRUE, length(target_pt))
  for (j in which(supported)) {
    insertion <- findInterval(target_pt[[j]], sorted_pt)
    lo <- max(1L, insertion - k)
    hi <- min(length(sorted_pt), insertion + k + 1L)
    candidates <- seq.int(lo, hi)
    candidate_distance <- abs(sorted_pt[candidates] - target_pt[[j]])
    nearest <- candidates[order(candidate_distance, candidates)][seq_len(k)]
    result[, j] <- rowMeans(sorted_deviation[, nearest, drop = FALSE])
    radius[[j]] <- max(abs(sorted_pt[nearest] - target_pt[[j]]))
  }
  list(value = result, radius = radius, supported = supported,
    source_min = min(sorted_pt), source_max = max(sorted_pt))
}

largest_contiguous_run <- function(x) {
  run <- rle(x)
  true_runs <- which(run$values)
  if (!length(true_runs)) stop("No common well-covered pseudotime interval")
  selected_run <- true_runs[which.max(run$lengths[true_runs])]
  end <- cumsum(run$lengths)[selected_run]
  start <- end - run$lengths[selected_run] + 1L
  seq.int(start, end)
}

selected_all <- fread(file.path(selection_report,
  "selected_top_adjusted_wt_ko_motifs.tsv"))
stopifnot(setequal(unique(selected_all$region_id),
    c("Flt3l_peak1", "Il7_peak2")),
  !anyDuplicated(selected_all[, .(region_id, motif_id)]))

bdev <- wb_open_bundle(bdev_bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
mapping <- wb_open_bundle(mapping_bundle_dir, read_only = TRUE,
  validate = FALSE)
on.exit(try(wb_close_bundle(mapping), silent = TRUE), add = TRUE)
bdev_manifest <- wb_manifest(bdev)
mapping_manifest <- wb_manifest(mapping)
stopifnot(identical(bdev_manifest$bundle_id, "bdev_reference"),
  identical(bdev_manifest$current_state$fingerprint,
    expected_bdev_fingerprint),
  identical(mapping_manifest$bundle_id, "bdev_integrated_ko_mapping_v1"),
  identical(mapping_manifest$current_state$fingerprint,
    expected_mapping_fingerprint))

regions <- fread(file.path(region_report, "wt_enriched_pt_regions.tsv"))[
  region_id %in% c("Flt3l_peak1", "Il7_peak2")]
assignments <- fread(cmd = paste("zcat", shQuote(file.path(region_report,
  "cell_region_assignments.tsv.gz"))))[
    region_id %in% c("Flt3l_peak1", "Il7_peak2")]
bdev_cells <- wb_cells(bdev, "set:bdev_reference")
bdev_meta <- as.data.table(wb_meta(bdev, cells = bdev_cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt")))
setnames(bdev_meta, names(bdev_meta)[[1L]], "source_cell_id")
stopifnot(nrow(regions) == 2L,
  !anyDuplicated(assignments[, .(region_id, cell_id)]),
  identical(as.character(bdev_meta$source_cell_id), bdev_cells))

source_env <- new.env(parent = emptyenv())
source_connection <- gzfile(bdev_rdata_path, "rb")
loaded <- load(source_connection, envir = source_env)
close(source_connection)
stopifnot("Bdev" %in% loaded)
source_bdev <- source_env$Bdev
motif_matrix_all <- Signac::GetMotifData(source_bdev[["peaks"]],
  slot = "data")
peak_covariates <- as.data.table(source_bdev[["peaks"]][[]],
  keep.rownames = "feature_id")

joint_counts <- BPCells::open_matrix_dir(joint_matrix_dir)
joint_peak_ids <- rownames(joint_counts)
joint_cell_ids <- colnames(joint_counts)
stopifnot(length(joint_peak_ids) == 97788L,
  all(joint_peak_ids %in% rownames(motif_matrix_all)),
  all(joint_peak_ids %in% peak_covariates$feature_id))
motif_matrix_all <- motif_matrix_all[joint_peak_ids, , drop = FALSE]
peak_bias <- peak_covariates$GC.percent[
  match(joint_peak_ids, peak_covariates$feature_id)]
stopifnot(all(is.finite(peak_bias)))
rm(source_bdev, source_env, peak_covariates)
invisible(gc())

all_values <- list()
all_support <- list()
all_selected <- list()
region_diagnostics <- list()
for (region_index in seq_len(nrow(regions))) {
  region <- regions[region_index]
  region_key <- region$region_id
  selected <- copy(selected_all[region_id == region_key])
  setorder(selected, row_order)
  stopifnot(nrow(selected) == 40L,
    all(selected$motif_id %in% colnames(motif_matrix_all)))
  selected[, display_row := sprintf("%s%s  [%+.2f]", motif_name,
    fifelse(strong_consistent_difference, "*", ""),
    adjusted_wt_minus_ko_bdev_sd)]
  all_selected[[region_key]] <- selected

  bdev_region <- bdev_meta[is.finite(pt) & pt >= region$pt_left &
    pt < region$pt_right]
  setorder(bdev_region, pt, source_cell_id)
  bdev_region[, cell_id := paste0("bdev::", source_cell_id)]
  query_region <- assignments[region_id == region_key]
  query_region[, group := ifelse(genotype == "WT", "WT", "KO")]
  stopifnot(nrow(bdev_region) >= neighbor_count,
    query_region[, min(.N), by = group][, min(V1)] >= neighbor_count,
    all(bdev_region$cell_id %in% joint_cell_ids),
    all(query_region$cell_id %in% joint_cell_ids))

  bdev_counts <- as(joint_counts[, bdev_region$cell_id, drop = FALSE],
    "dgCMatrix")
  peak_keep <- Matrix::rowSums(bdev_counts) >= 10
  counts_kept_bdev <- bdev_counts[peak_keep, , drop = FALSE]
  counts_kept_query <- as(joint_counts[peak_keep, query_region$cell_id,
    drop = FALSE], "dgCMatrix")
  combined_counts <- cbind(counts_kept_bdev, counts_kept_query)
  motif_matrix <- motif_matrix_all[peak_keep, selected$motif_id,
    drop = FALSE]
  stopifnot(all(Matrix::colSums(motif_matrix != 0) >= 10L))

  set.seed(seed + region_index)
  background_peaks <- chromVAR::getBackgroundPeaks(counts_kept_bdev,
    bias = peak_bias[peak_keep], niterations = 50L, w = 0.1, bs = 50L)
  bdev_expectation <- chromVAR::computeExpectations(counts_kept_bdev)
  deviation_object <- chromVAR::computeDeviations(combined_counts,
    motif_matrix, background_peaks = background_peaks,
    expectation = bdev_expectation)
  deviation <- SummarizedExperiment::assay(deviation_object, "deviations")
  stopifnot(all(is.finite(deviation)),
    identical(rownames(deviation), selected$motif_id))

  focal_pt <- bdev_region$pt
  focal_cell_id <- bdev_region$source_cell_id
  bdev_deviation <- deviation[, bdev_region$cell_id, drop = FALSE]
  wt_cells <- query_region[group == "WT"]
  ko_cells <- query_region[group == "KO"]
  wt_deviation <- deviation[, wt_cells$cell_id, drop = FALSE]
  ko_deviation <- deviation[, ko_cells$cell_id, drop = FALSE]
  smooth_bdev <- smooth_30nn(bdev_deviation, bdev_region$pt, focal_pt,
    neighbor_count, restrict_to_observed_support = FALSE)
  smooth_wt <- smooth_30nn(wt_deviation, wt_cells$mapped_pt, focal_pt,
    neighbor_count, restrict_to_observed_support = TRUE)
  smooth_ko <- smooth_30nn(ko_deviation, ko_cells$mapped_pt, focal_pt,
    neighbor_count, restrict_to_observed_support = TRUE)
  common_coverage <- smooth_wt$supported & smooth_ko$supported &
    smooth_wt$radius <= maximum_neighbor_radius_pt &
    smooth_ko$radius <= maximum_neighbor_radius_pt
  displayed_indices <- largest_contiguous_run(common_coverage)
  displayed_common_coverage <- seq_along(focal_pt) %in% displayed_indices

  centers <- rowMeans(smooth_bdev$value)
  scales <- apply(smooth_bdev$value, 1L, sd)
  stopifnot(all(is.finite(scales) & scales > 0))
  standardize <- function(x) sweep(sweep(x, 1L, centers, "-"),
    1L, scales, "/")
  activity <- list(Bdev = standardize(smooth_bdev$value),
    WT = standardize(smooth_wt$value),
    KO = standardize(smooth_ko$value))
  difference <- activity$WT - activity$KO

  make_long <- function(group_name, matrix_value) {
    out <- as.data.table(as.table(matrix_value))
    setnames(out, c("motif_id", "focal_index", "value_bdev_sd"))
    out[, `:=`(motif_id = as.character(motif_id),
      focal_index = as.integer(as.character(focal_index)))]
    out[, `:=`(group = group_name,
      focal_cell_id = focal_cell_id[focal_index],
      focal_pt = focal_pt[focal_index], region_id = region_key)]
    out
  }
  heat <- rbindlist(list(make_long("Bdev", activity$Bdev),
    make_long("WT", activity$WT), make_long("KO", activity$KO),
    make_long("WT - KO", difference)))
  heat <- merge(heat, selected[, .(motif_id, display_row, row_order,
    motif_name, adjusted_wt_minus_ko_bdev_sd,
    strong_consistent_difference)], by = "motif_id", all.x = TRUE,
    sort = FALSE)
  heat[, displayed_common_coverage :=
    displayed_common_coverage[focal_index]]
  all_values[[region_key]] <- copy(heat)

  support <- rbindlist(list(
    data.table(group = "Bdev", focal_index = seq_along(focal_pt),
      neighbor_radius_pt = smooth_bdev$radius,
      supported = smooth_bdev$supported),
    data.table(group = "WT", focal_index = seq_along(focal_pt),
      neighbor_radius_pt = smooth_wt$radius,
      supported = smooth_wt$supported),
    data.table(group = "KO", focal_index = seq_along(focal_pt),
      neighbor_radius_pt = smooth_ko$radius,
      supported = smooth_ko$supported)))
  support[, `:=`(focal_cell_id = focal_cell_id[focal_index],
    focal_pt = focal_pt[focal_index], neighbor_count = neighbor_count,
    region_id = region_key,
    displayed_common_coverage =
      displayed_common_coverage[focal_index])]
  all_support[[region_key]] <- support

  heat[, display_row := factor(display_row,
    levels = rev(selected$display_row))]
  activity_heat <- heat[group %in% c("Bdev", "WT", "KO") &
    displayed_common_coverage]
  activity_labels <- if (region$target == "Flt3l")
    c(Bdev = "Bdev", WT = "FLT3L WT", KO = "FLT3LKO") else
    c(Bdev = "Bdev", WT = "IL-7 WT", KO = "IL7KO")
  activity_heat[, panel := factor(activity_labels[group],
    levels = unname(activity_labels))]
  difference_heat <- heat[group == "WT - KO" &
    displayed_common_coverage]
  difference_heat[, panel := factor("Matched-position WT − KO")]
  difference_limit <- as.numeric(quantile(abs(difference_heat$value_bdev_sd),
    0.98, na.rm = TRUE))
  difference_limit <- max(1, min(6, difference_limit))
  activity_limit <- 2.5
  axis_breaks <- unique(as.integer(round(seq(min(displayed_indices),
    max(displayed_indices),
    length.out = 5))))
  axis_labels <- sprintf("%.2f", focal_pt[axis_breaks])

  p_activity <- ggplot(activity_heat,
      aes(focal_index, display_row, fill = value_bdev_sd)) +
    geom_raster() +
    facet_grid(. ~ panel) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-activity_limit, activity_limit),
      oob = scales::squish, na.value = "grey75",
      name = "30-NN activity\n(Bdev SD)") +
    scale_x_continuous(expand = expansion(mult = 0),
      breaks = axis_breaks, labels = axis_labels) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6.4),
      legend.position = "right") +
    labs(x = "Bdev focal cells ordered by pseudotime (tick labels = pt)",
      y = NULL)
  p_difference <- ggplot(difference_heat,
      aes(focal_index, display_row, fill = value_bdev_sd)) +
    geom_raster() +
    facet_grid(. ~ panel) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-difference_limit, difference_limit),
      oob = scales::squish, na.value = "grey75",
      name = "WT − KO\n(Bdev SD)") +
    scale_x_continuous(expand = expansion(mult = 0),
      breaks = axis_breaks, labels = axis_labels) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      legend.position = "right") +
    labs(x = "Bdev focal cells ordered by pseudotime (tick labels = pt)",
      y = NULL)
  title <- if (region$target == "Flt3l")
    "Flt3L peak 1: single-cell-axis 30-nearest-neighbor motif heatmap" else
    "IL-7 peak 2: single-cell-axis 30-nearest-neighbor motif heatmap"
  figure <- p_activity + p_difference +
    plot_layout(widths = c(3, 1.15)) +
    plot_annotation(title = title,
      caption = paste0(length(displayed_indices), " of ", nrow(bdev_region),
        " Bdev focal cells shown: both query groups have 30-NN radius ≤ ",
        maximum_neighbor_radius_pt,
        " pt. Rows retain the prior top-20 WT-higher/top-20 KO-higher selection and order. Brackets show the prior adjusted effect; * marks its strong-consistency flag. Activity clipped at ±",
        activity_limit, "; difference clipped at ±",
        signif(difference_limit, 3), " for display only."),
      theme = theme(plot.title = element_text(face = "bold", size = 15),
        plot.caption = element_text(size = 8)))
  stub <- if (region$target == "Flt3l") "flt3l_peak1" else "il7_peak2"
  ggsave(file.path(stage_dir, "figures", paste0(stub,
    "_single_cell_axis_30nn_top_adjusted_motif_heatmap.png")),
    figure, width = 18, height = 11, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(stub,
    "_single_cell_axis_30nn_top_adjusted_motif_heatmap.pdf")),
    figure, width = 18, height = 11, device = cairo_pdf, bg = "white")

  region_diagnostics[[region_key]] <- list(
    selected_motifs = nrow(selected),
    bdev_focal_cells = nrow(bdev_region),
    wt_source_cells = nrow(wt_cells), ko_source_cells = nrow(ko_cells),
    neighbor_count = neighbor_count,
    maximum_neighbor_radius_pt = maximum_neighbor_radius_pt,
    displayed_focal_cells = length(displayed_indices),
    displayed_pt_range = focal_pt[range(displayed_indices)],
    wt_supported_focal_cells = sum(smooth_wt$supported),
    ko_supported_focal_cells = sum(smooth_ko$supported),
    wt_source_pt_range = c(smooth_wt$source_min, smooth_wt$source_max),
    ko_source_pt_range = c(smooth_ko$source_min, smooth_ko$source_max),
    bdev_accessible_peaks = sum(peak_keep),
    activity_display_limit = activity_limit,
    difference_display_limit = difference_limit)
  rm(bdev_counts, counts_kept_bdev, counts_kept_query, combined_counts,
    deviation_object, deviation, bdev_deviation, wt_deviation, ko_deviation)
  invisible(gc())
}

write_tsv(rbindlist(all_selected, fill = TRUE), file.path(stage_dir,
  "selected_motifs.tsv"))
fwrite(rbindlist(all_values, fill = TRUE), file.path(stage_dir,
  "single_cell_axis_30nn_heatmap_values.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")
fwrite(rbindlist(all_support, fill = TRUE), file.path(stage_dir,
  "focal_cell_neighbor_support.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")
write_json(list(status = "complete_preliminary_unregistered",
  source_selection_report = selection_report,
  bdev_bundle_fingerprint = expected_bdev_fingerprint,
  mapping_bundle_fingerprint = expected_mapping_fingerprint,
  focal_axis = "every Bdev cell in the respective pseudotime region",
  smoothing = paste0(neighbor_count,
    " nearest pseudotime cells independently within Bdev, WT, and KO"),
  alignment = "all groups evaluated at the same Bdev focal-cell pseudotimes",
  display_coverage = paste0(
    "largest contiguous interval where both WT and KO 30-neighbor radii are <= ",
    maximum_neighbor_radius_pt, " pseudotime units"),
  query_support = "full support audit and values retained outside the displayed interval",
  scaling = "center and SD of the Bdev 30-neighbor smoothed trajectory per motif",
  region_diagnostics = region_diagnostics,
  biological_replicates_per_genotype = 1L,
  inference = "descriptive", bundle_mutated = FALSE,
  capsule_registered = FALSE), file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Single-cell-axis 30-nearest-neighbor motif heatmaps", "",
  "This report retains the prior top-20 adjusted WT-higher and top-20 adjusted KO-higher motif selections separately for Flt3L peak 1 and IL-7 peak 2.",
  paste0("Every Bdev cell in the region supplies one ordered focal pseudotime column. Bdev, WT, and KO are each smoothed independently from their ",
    neighbor_count, " closest pseudotime cells and evaluated at those common focal positions."),
  paste0("Figures are restricted to the largest contiguous interval where both WT and KO 30-neighbor radii are at most ",
    maximum_neighbor_radius_pt, " pseudotime units."),
  "WT-minus-KO is therefore calculated at matched focal pseudotime positions rather than from unequal genotype cell distributions.",
  "The full-axis values and neighbor-support audit remain exported with a displayed-common-coverage flag.",
  "There is one biological sample per genotype, so results are descriptive."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
stopifnot(identical(wb_manifest(bdev)$current_state$fingerprint,
    expected_bdev_fingerprint),
  identical(wb_manifest(mapping)$current_state$fingerprint,
    expected_mapping_fingerprint))
wb_close_bundle(bdev)
bdev <- NULL
wb_close_bundle(mapping)
mapping <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed single-cell-axis 30-NN motif heatmaps at", output_dir, "\n")
