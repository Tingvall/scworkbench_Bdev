#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L) stop(paste(
  "Usage: 89_single_cell_axis_15nn_motif_heatmaps.R",
  "<bdev_bundle> <mapping_bundle> <joint_peak_matrix> <support_report>",
  "<query_support_tsv_gz> <bdev_rdata> <selection_report> <output_dir>",
  "<bdev_fingerprint> <mapping_fingerprint> <neighbor_count>",
  "<maximum_neighbor_radius_pt> <seed>"), call. = FALSE)
bdev_dir <- normalizePath(args[[1L]], mustWork = TRUE)
mapping_dir <- normalizePath(args[[2L]], mustWork = TRUE)
joint_dir <- normalizePath(args[[3L]], mustWork = TRUE)
support_dir <- normalizePath(args[[4L]], mustWork = TRUE)
query_path <- normalizePath(args[[5L]], mustWork = TRUE)
bdev_rdata <- normalizePath(args[[6L]], mustWork = TRUE)
selection_dir <- normalizePath(args[[7L]], mustWork = TRUE)
output_dir <- normalizePath(args[[8L]], mustWork = FALSE)
bdev_fp <- args[[9L]]
mapping_fp <- args[[10L]]
neighbor_count <- as.integer(args[[11L]])
maximum_radius <- as.numeric(args[[12L]])
seed <- as.integer(args[[13L]])
stopifnot(neighbor_count >= 3L, !is.na(maximum_radius),
  maximum_radius > 0, is.finite(seed))
nn_tag <- paste0(neighbor_count, "nn")
apply_radius_mask <- is.finite(maximum_radius)
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
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(..., collapse = " "), "\n")

knn_smooth <- function(deviation, source_pt, focal_pt, k) {
  stopifnot(ncol(deviation) == length(source_pt), length(source_pt) >= k,
    all(is.finite(source_pt)), all(is.finite(focal_pt)))
  ord <- order(source_pt)
  source_pt <- source_pt[ord]
  deviation <- deviation[, ord, drop = FALSE]
  value <- matrix(NA_real_, nrow(deviation), length(focal_pt),
    dimnames = list(rownames(deviation), as.character(seq_along(focal_pt))))
  radius <- rep(NA_real_, length(focal_pt))
  observed <- focal_pt >= min(source_pt) & focal_pt <= max(source_pt)
  for (j in seq_along(focal_pt)) {
    insertion <- findInterval(focal_pt[[j]], source_pt)
    lo <- max(1L, insertion - k)
    hi <- min(length(source_pt), insertion + k + 1L)
    candidates <- seq.int(lo, hi)
    nearest <- candidates[order(abs(source_pt[candidates] - focal_pt[[j]]),
      candidates)][seq_len(k)]
    value[, j] <- rowMeans(deviation[, nearest, drop = FALSE])
    radius[[j]] <- max(abs(source_pt[nearest] - focal_pt[[j]]))
  }
  list(value = value, radius = radius, observed = observed,
    source_range = range(source_pt))
}
knn_radius <- function(source_pt, focal_pt, k) {
  stopifnot(length(source_pt) >= k, all(is.finite(source_pt)),
    all(is.finite(focal_pt)))
  source_pt <- sort(source_pt)
  vapply(focal_pt, function(x) {
    insertion <- findInterval(x, source_pt)
    lo <- max(1L, insertion - k)
    hi <- min(length(source_pt), insertion + k + 1L)
    candidates <- seq.int(lo, hi)
    nearest <- candidates[order(abs(source_pt[candidates] - x),
      candidates)][seq_len(k)]
    max(abs(source_pt[nearest] - x))
  }, numeric(1L))
}

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
source_bins <- fread(file.path(support_dir, "source_bin_cell_counts.tsv"))
original_width <- median(source_bins$bin_right - source_bins$bin_left)
fine_width <- original_width / 2
selected_all <- fread(file.path(selection_dir,
  "selected_bdev_variable_motifs.tsv"))[selected == TRUE]
query <- fread(cmd = paste("zcat", shQuote(query_path)))[
  current_accepted == TRUE & trajectory_support >= 0.25 & is.finite(mapped_pt)]
stopifnot(nrow(regions) == 2L,
  setequal(regions$region_id, c("flt3l_wt_peak1_run1", "il7_wt_peak1_run1")),
  selected_all[, uniqueN(motif_id), by = region_id][, all(V1 == 30L)],
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
joint <- BPCells::open_matrix_dir(joint_dir)
peak_ids <- rownames(joint)
joint_cells <- colnames(joint)
stopifnot(all(peak_ids %in% rownames(motif_all)),
  all(peak_ids %in% peak_meta$feature_id))
motif_all <- motif_all[peak_ids, , drop = FALSE]
peak_bias <- peak_meta$GC.percent[match(peak_ids, peak_meta$feature_id)]
stopifnot(all(is.finite(peak_bias)))
rm(source_bdev, source_env, peak_meta)
invisible(gc())

all_values <- list()
all_support <- list()
all_selected <- list()
all_density <- list()
diagnostics <- list()
for (ri in seq_len(nrow(regions))) {
  reg <- regions[ri]
  region_id <- reg$region_id
  target <- reg$target
  selected <- copy(selected_all[region_id == reg$region_id])
  setorder(selected, bdev_variability_rank)
  stopifnot(nrow(selected) == 30L,
    all(selected$motif_id %in% colnames(motif_all)))
  all_selected[[region_id]] <- selected

  bdev_region <- bdev_meta[is.finite(pt) & pt >= reg$center_pt_left &
    pt <= reg$center_pt_right]
  setorder(bdev_region, pt, source_cell_id)
  focal_pt <- bdev_region$pt
  focal_id <- bdev_region$source_cell_id
  query_target <- query[knockout_target == target]
  wt <- query_target[genotype == "WT"]
  ko <- query_target[genotype == "KO"]
  stopifnot(nrow(bdev_region) >= neighbor_count,
    nrow(wt) >= neighbor_count, nrow(ko) >= neighbor_count,
    all(c(bdev_region$cell_id, query_target$cell_id) %in% joint_cells))

  fine_centers <- seq(reg$center_pt_left, reg$center_pt_right,
    by = fine_width)
  fine_ko_n <- vapply(fine_centers, function(x)
    ko[mapped_pt >= x - fine_width / 2 & mapped_pt < x + fine_width / 2, .N],
    integer(1L))
  last_supported <- max(which(fine_ko_n >= 10L))
  stopifnot(last_supported < length(fine_centers),
    all(fine_ko_n[seq.int(last_supported + 1L, length(fine_centers))] < 10L))
  ko_boundary <- mean(fine_centers[c(last_supported, last_supported + 1L)])

  log_event("Computing per-cell chromVAR for", region_id,
    "with", length(focal_pt), "Bdev focal cells")
  bdev_counts <- as(joint[, bdev_region$cell_id, drop = FALSE], "dgCMatrix")
  peak_keep <- Matrix::rowSums(bdev_counts) >= 10
  combined_ids <- c(bdev_region$cell_id, query_target$cell_id)
  combined <- as(joint[peak_keep, combined_ids, drop = FALSE], "dgCMatrix")
  motif <- motif_all[peak_keep, selected$motif_id, drop = FALSE]
  stopifnot(all(Matrix::colSums(motif != 0) >= 10L))
  set.seed(seed + ri)
  bg <- chromVAR::getBackgroundPeaks(bdev_counts[peak_keep, , drop = FALSE],
    bias = peak_bias[peak_keep], niterations = 50L, w = 0.1, bs = 50L)
  expectation <- chromVAR::computeExpectations(
    bdev_counts[peak_keep, , drop = FALSE])
  dev_obj <- chromVAR::computeDeviations(combined, motif,
    background_peaks = bg, expectation = expectation)
  dev <- SummarizedExperiment::assay(dev_obj, "deviations")
  stopifnot(identical(colnames(dev), combined_ids), all(is.finite(dev)))

  smooth_bdev <- knn_smooth(dev[, bdev_region$cell_id, drop = FALSE],
    bdev_region$pt, focal_pt, neighbor_count)
  smooth_wt <- knn_smooth(dev[, wt$cell_id, drop = FALSE], wt$mapped_pt,
    focal_pt, neighbor_count)
  smooth_ko <- knn_smooth(dev[, ko$cell_id, drop = FALSE], ko$mapped_pt,
    focal_pt, neighbor_count)
  center <- rowMeans(smooth_bdev$value)
  scale <- apply(smooth_bdev$value, 1L, sd)
  stopifnot(all(is.finite(scale) & scale > 0))
  standardize <- function(x) sweep(sweep(x, 1L, center, "-"), 1L, scale, "/")
  activity <- list(Bdev = standardize(smooth_bdev$value),
    WT = standardize(smooth_wt$value), KO = standardize(smooth_ko$value))
  activity[["WT - KO"]] <- activity$WT - activity$KO

  support <- data.table(focal_index = seq_along(focal_pt),
    focal_cell_id = focal_id, focal_pt = focal_pt,
    bdev_radius = smooth_bdev$radius, wt_radius = smooth_wt$radius,
    ko_radius = smooth_ko$radius,
    bdev_supported = if (apply_radius_mask)
      smooth_bdev$radius <= maximum_radius else TRUE,
    wt_supported = if (apply_radius_mask)
      smooth_wt$observed & smooth_wt$radius <= maximum_radius else TRUE,
    ko_supported = if (apply_radius_mask)
      smooth_ko$observed & smooth_ko$radius <= maximum_radius &
        focal_pt <= ko_boundary else focal_pt <= ko_boundary)
  support[, `:=`(comparison_supported = wt_supported & ko_supported,
    region_id = region_id, neighbor_count = neighbor_count,
    maximum_neighbor_radius_pt = if (apply_radius_mask) maximum_radius else Inf,
    ko_support_boundary_pt = ko_boundary)]
  all_support[[region_id]] <- support

  density_bdev_radius <- knn_radius(bdev_meta[is.finite(pt), pt], focal_pt,
    neighbor_count)
  density_wt_radius <- knn_radius(wt$mapped_pt, focal_pt, neighbor_count)
  density_ko_radius <- knn_radius(ko$mapped_pt, focal_pt, neighbor_count)
  # Mapping can assign identical pseudotimes to multiple query cells. Use a
  # small, explicit bandwidth floor so tied neighbors do not imply infinite
  # local density.
  radius_floor <- fine_width / 20
  bdev_total <- bdev_meta[is.finite(pt), .N]
  wt_total <- nrow(wt)
  ko_total <- nrow(ko)
  density <- data.table(focal_index = seq_along(focal_pt), focal_pt,
    bdev_radius = density_bdev_radius, wt_radius = density_wt_radius,
    ko_radius = density_ko_radius)
  density[, `:=`(
    wt_relative_to_bdev =
      pmax(bdev_radius, radius_floor) * bdev_total /
        (pmax(wt_radius, radius_floor) * wt_total),
    ko_relative_to_bdev =
      pmax(bdev_radius, radius_floor) * bdev_total /
        (pmax(ko_radius, radius_floor) * ko_total),
    region_id = region_id)]
  density_plot <- rbindlist(list(
    density[, .(focal_index, focal_pt, group = "Bdev", series = "Bdev",
      relative_density = 1)],
    density[, .(focal_index, focal_pt, group = "WT", series = "WT",
      relative_density = wt_relative_to_bdev)],
    density[, .(focal_index, focal_pt, group = "KO", series = "KO",
      relative_density = ko_relative_to_bdev)],
    density[, .(focal_index, focal_pt, group = "WT - KO", series = "WT",
      relative_density = wt_relative_to_bdev)],
    density[, .(focal_index, focal_pt, group = "WT - KO", series = "KO",
      relative_density = ko_relative_to_bdev)]))
  all_density[[region_id]] <- copy(density)

  make_long <- function(group, value) {
    ans <- as.data.table(as.table(value))
    setnames(ans, c("motif_id", "focal_index", "value_bdev_sd"))
    ans[, `:=`(motif_id = as.character(motif_id),
      focal_index = as.integer(as.character(focal_index)), group = group)]
    ans
  }
  heat <- rbindlist(Map(make_long, names(activity), activity))
  heat <- merge(heat, support[, .(focal_index, focal_cell_id, focal_pt,
    bdev_supported, wt_supported, ko_supported, comparison_supported)],
    by = "focal_index", sort = FALSE)
  heat <- merge(heat, selected[, .(motif_id, motif_name,
    bdev_variability_rank)], by = "motif_id", sort = FALSE)
  heat[, display_supported := fcase(group == "Bdev", bdev_supported,
    group == "WT", wt_supported, group == "KO", ko_supported,
    group == "WT - KO", comparison_supported, default = FALSE)]
  heat[, display_value_bdev_sd := fifelse(display_supported,
    value_bdev_sd, NA_real_)]
  heat[, region_id := region_id]
  row_matrix <- activity$Bdev
  row_order <- rownames(row_matrix)[hclust(dist(row_matrix))$order]
  labels <- selected[, .(motif_id,
    motif_label = make.unique(sprintf("%s (%s)", motif_name, motif_id)))]
  heat <- merge(heat, labels, by = "motif_id", sort = FALSE)
  heat[, motif_label := factor(motif_label,
    levels = rev(labels$motif_label[match(row_order, labels$motif_id)]))]
  group_labels <- if (target == "Flt3l")
    c(Bdev = "Bdev", WT = "Flt3L WT", KO = "FLT3L KO",
      `WT - KO` = "WT - KO") else
    c(Bdev = "Bdev", WT = "IL-7 WT", KO = "IL-7 KO",
      `WT - KO` = "WT - KO")
  heat[, panel := factor(group_labels[group], levels = unname(group_labels))]
  density_plot[, `:=`(panel = factor(group_labels[group],
    levels = unname(group_labels)), series = factor(series,
    levels = c("Bdev", "WT", "KO")))]
  all_values[[region_id]] <- copy(heat)

  boundary_x <- mean(c(max(which(focal_pt <= ko_boundary)),
    min(which(focal_pt > ko_boundary))))
  axis_breaks <- unique(as.integer(round(seq(1, length(focal_pt),
    length.out = 6))))
  axis_labels <- sprintf("%.2f", focal_pt[axis_breaks])
  figure_title <- paste0(if (target == "Flt3l") "Flt3L" else "IL-7",
    " WT-enriched region: individual-cell-axis ", neighbor_count,
    "-NN motif trajectory")
  figure_subtitle <- paste("Every Bdev cell is a focal position; each group uses",
    neighbor_count, "nearest pseudotime cells")
  radius_caption <- if (apply_radius_mask)
    paste0("; any group with ", neighbor_count, "-NN radius > ",
      maximum_radius, " pt is also masked") else
    "; neighbor radius is audited but is not used to mask values"
  p_heat <- ggplot(heat, aes(focal_index, motif_label,
      fill = display_value_bdev_sd)) +
    geom_raster() +
    geom_vline(xintercept = boundary_x, linetype = "dashed",
      colour = "#D55E00", linewidth = 0.55) +
    facet_grid(. ~ panel) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish,
      na.value = "grey85", name = paste0(neighbor_count,
        "-NN activity\n(Bdev SD)")) +
    scale_x_continuous(expand = expansion(mult = 0),
      breaks = axis_breaks, labels = axis_labels) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6.5),
      plot.title = element_text(face = "bold")) +
    labs(x = "Bdev focal cells ordered by pseudotime (labels = pt)", y = NULL,
      caption = paste0("Orange line: fine-bin KO support boundary. KO and WT-KO are masked after it", radius_caption,
        ". One sample per genotype."))
  p_density <- ggplot(density_plot,
      aes(focal_index, relative_density, fill = series, colour = series,
        group = series)) +
    geom_hline(yintercept = 1, colour = "grey35", linetype = "dotted",
      linewidth = 0.4) +
    geom_area(position = "identity", alpha = 0.30, linewidth = 0) +
    geom_line(linewidth = 0.45) +
    geom_vline(xintercept = boundary_x, linetype = "dashed",
      colour = "#D55E00", linewidth = 0.55) +
    facet_grid(. ~ panel) +
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
    labs(y = paste0(neighbor_count, "-NN density / Bdev"),
      fill = NULL, colour = NULL)
  p <- p_density / p_heat + plot_layout(heights = c(1.8, 8.2)) +
    plot_annotation(title = figure_title, subtitle = figure_subtitle)
  stub <- if (target == "Flt3l") "flt3l_wt_peak1" else "il7_wt_peak1"
  ggsave(file.path(stage_dir, "figures", paste0(stub,
    "_single_cell_axis_", nn_tag, "_motif_heatmap.png")), p,
    width = 18, height = 10.5, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(stub,
    "_single_cell_axis_", nn_tag, "_motif_heatmap.pdf")), p,
    width = 18, height = 10.5, device = cairo_pdf, bg = "white")

  diagnostics[[region_id]] <- list(target = target,
    bdev_focal_cells = length(focal_pt), wt_source_cells = nrow(wt),
    ko_source_cells = nrow(ko), neighbor_count = neighbor_count,
    maximum_neighbor_radius_pt = if (apply_radius_mask) maximum_radius else NULL,
    neighbor_radius_used_for_masking = apply_radius_mask,
    ko_support_boundary_pt = ko_boundary,
    supported_focal_cells = list(Bdev = sum(support$bdev_supported),
      WT = sum(support$wt_supported), KO = sum(support$ko_supported),
      comparison = sum(support$comparison_supported)),
    median_neighbor_radius = list(Bdev = median(support$bdev_radius),
      WT = median(support$wt_radius), KO = median(support$ko_radius)),
    selected_motifs = nrow(selected), bdev_accessible_peaks = sum(peak_keep))
  rm(bdev_counts, combined, dev_obj, dev, heat, activity)
  invisible(gc())
}

write_tsv(rbindlist(all_selected), file.path(stage_dir, "selected_motifs.tsv"))
write_tsv(rbindlist(all_density), file.path(stage_dir,
  paste0("single_cell_axis_", nn_tag, "_density_relative_to_bdev.tsv")))
fwrite(rbindlist(all_support), file.path(stage_dir,
  "focal_cell_neighbor_support.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")
fwrite(rbindlist(all_values), file.path(stage_dir,
  paste0("single_cell_axis_", nn_tag, "_motif_values.tsv.gz")),
  sep = "\t", na = "NA",
  compress = "gzip")
write_json(list(status = "complete_unpromoted_report_only",
  source_selection_report = selection_dir,
  selection = "fixed top 30 motifs from the Bdev-only fine rolling analysis",
  focal_axis = "every Bdev cell in each retained region, ordered by preserved pt",
  smoothing = paste(neighbor_count,
    "nearest pseudotime cells independently within Bdev, WT, and KO"),
  mapped_cell_selection = "current_accepted and trajectory_support >= 0.25",
  neighbor_radius_used_for_masking = apply_radius_mask,
  density_track = paste(neighbor_count,
    "neighbor local pseudotime density divided by Bdev local density, normalized by group total; radius floor", fine_width / 20),
  maximum_neighbor_radius_pt = if (apply_radius_mask) maximum_radius else NULL,
  region_diagnostics = diagnostics, biological_replicates_per_genotype = 1L,
  inference = "descriptive", bundle_mutated = FALSE, promotion = FALSE),
  file.path(stage_dir, "diagnostics.json"))
writeLines(c(paste0("# Individual-cell-axis ", neighbor_count,
  "-nearest-neighbor motif trajectories"), "",
  paste0("Every Bdev cell in the retained Flt3L WT peak 1 or IL-7 WT peak 1 region supplies one ordered focal pseudotime position. At each position, Bdev, WT, and KO motif activity is the mean per-cell chromVAR deviation of that group's ", neighbor_count, " nearest pseudotime cells."), "",
  "The motif identities are fixed to the top 30 selected using Bdev in the fine rolling analysis, allowing the smoothing approaches to be compared directly. Per-cell chromVAR calibration and trajectory centering/scaling remain Bdev-only.", "",
  paste0("Filled tracks above the heatmaps show ", neighbor_count, "-neighbor local pseudotime density relative to Bdev, normalized by the total number of cells in each group. The WT-minus-KO panel overlays WT and KO density."), "",
  if (apply_radius_mask) paste0("A group is masked when its ", neighbor_count,
    "th neighbor is more than ", maximum_radius,
    " pseudotime units away. KO and WT-minus-KO are additionally masked after the fine-bin 10-KO-cell support boundary. The complete radius and display-support audit is exported.") else
    "Neighbor radius is audited but does not mask motif values. The 0.25 trajectory-support cutoff is used only to select accepted mapped query cells. KO and WT-minus-KO remain masked after the agreed fine-bin 10-KO-cell support boundary.", "",
  "Because there is one biological sample per genotype and neighboring focal positions reuse cells, comparisons are descriptive.", "",
  paste0("- [Flt3L figure](figures/flt3l_wt_peak1_single_cell_axis_", nn_tag,
    "_motif_heatmap.png)"),
  paste0("- [IL-7 figure](figures/il7_wt_peak1_single_cell_axis_", nn_tag,
    "_motif_heatmap.png)")),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
wb_close_bundle(bdev); bdev <- NULL
wb_close_bundle(mapping); mapping <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed individual-cell-axis", paste0(neighbor_count, "-NN"),
  "motif analysis at", output_dir, "\n")
