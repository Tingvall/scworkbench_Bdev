#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) stop(paste(
  "Usage: 49_selected_motif_bdev_umap_query_violins.R",
  "<bdev_bundle> <mapping_bundle> <joint_peak_matrix_dir>",
  "<region_report> <bdev_rdata> <output_dir>",
  "<expected_bdev_fingerprint> <expected_mapping_fingerprint> <seed>",
  "<motif_preset>"))
bdev_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
mapping_bundle_dir <- normalizePath(args[[2L]], mustWork = TRUE)
joint_matrix_dir <- normalizePath(args[[3L]], mustWork = TRUE)
region_report <- normalizePath(args[[4L]], mustWork = TRUE)
bdev_rdata_path <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
expected_bdev_fingerprint <- args[[7L]]
expected_mapping_fingerprint <- args[[8L]]
seed <- as.integer(args[[9L]])
preset <- args[[10L]]
stopifnot(is.finite(seed))
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
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
write_tsv <- function(x, path) data.table::fwrite(x, path, sep = "\t",
  na = "NA")

if (preset == "snai_tcf3_runx1_srebf2_stat5") {
  motifs <- data.table(
    motif_label = c("SNAI2", "TCF3", "RUNX1", "SREBF2",
      "STAT5A::STAT5B"),
    motif_id = c("MA0745.2", "MA0522.3", "MA0002.2", "MA0596.1",
      "MA0519.1"),
    selection_note = c(
      "most Bdev-variable SNAI-family motif across both regions",
      "exact TCF3 motif", "exact RUNX1 motif", "exact SREBF2 motif",
      "STAT5A::STAT5B family motif"))
} else if (preset == "tcf3_ebf1_zeb_stat5_runx1") {
  motifs <- data.table(
    motif_label = c("TCF3", "EBF1", "ZEB1 (ZEB2 proxy)",
      "STAT5A::STAT5B", "RUNX1"),
    motif_id = c("MA0522.3", "MA0154.4", "MA0103.3", "MA0519.1",
      "MA0002.2"),
    selection_note = c("exact TCF3 motif", "exact EBF1 motif",
      "ZEB1 motif used as ZEB-family proxy because no ZEB2 motif is present",
      "STAT5A::STAT5B family motif", "exact RUNX1 motif"))
} else {
  stop("Unknown motif preset: ", preset)
}
motifs[, motif_order := seq_len(.N)]

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
stopifnot(nrow(regions) == 2L,
  setequal(regions$region_id, c("Flt3l_peak1", "Il7_peak2")),
  !anyDuplicated(assignments[, .(region_id, cell_id)]))

bdev_cells <- wb_cells(bdev, "set:bdev_reference")
bdev_meta <- as.data.table(wb_meta(bdev, cells = bdev_cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt")))
setnames(bdev_meta, names(bdev_meta)[[1L]], "source_cell_id")
stopifnot(identical(as.character(bdev_meta$source_cell_id), bdev_cells),
  all(motifs$motif_id %in% wb_feature_ids(bdev, "chromvar_data")))

# Native source Bdev chromVAR activity on the unchanged ATAC UMAP.
bdev_activity <- as.matrix(wb_matrix_read_features(bdev, "chromvar_data",
  motifs$motif_id, bdev_cells, missing = "NA"))
stopifnot(identical(rownames(bdev_activity), motifs$motif_id),
  identical(colnames(bdev_activity), bdev_cells),
  all(is.finite(bdev_activity)))
space <- wb_space(bdev, "umap.atac", cells = bdev_cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), bdev_cells))
  umap <- data.table(cell_id = bdev_cells, x = as.numeric(space$x),
    y = as.numeric(space$y))
} else {
  stopifnot(identical(rownames(space), bdev_cells), ncol(space) == 2L)
  umap <- data.table(cell_id = bdev_cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]))
}
stopifnot(all(is.finite(umap$x)), all(is.finite(umap$y)))

umap_values <- list()
umap_scales <- list()
umap_plots <- list()
for (i in seq_len(nrow(motifs))) {
  motif <- motifs[i]
  value <- as.numeric(bdev_activity[motif$motif_id, ])
  limit <- as.numeric(quantile(abs(value), 0.99, na.rm = TRUE))
  stopifnot(is.finite(limit), limit > 0)
  d <- copy(umap)
  d[, `:=`(motif_id = motif$motif_id,
    motif_label = motif$motif_label, chromvar_deviation = value,
    display_value = pmax(pmin(value, limit), -limit))]
  setorder(d, display_value)
  umap_values[[i]] <- copy(d)
  umap_scales[[i]] <- data.table(motif_id = motif$motif_id,
    motif_label = motif$motif_label, symmetric_absolute_p99 = limit,
    observed_minimum = min(value), observed_maximum = max(value))
  umap_plots[[i]] <- ggplot(d, aes(x, y, color = display_value)) +
    geom_point(size = 0.28, alpha = 1, stroke = 0) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-limit, limit), oob = scales::squish,
      name = NULL) +
    coord_equal() +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8.5),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.key.width = grid::unit(12, "mm")) +
    labs(title = motif$motif_label,
      subtitle = paste0(motif$motif_id, "; symmetric |p99| = ",
        signif(limit, 3)))
}
umap_figure <- wrap_plots(umap_plots, ncol = 3, guides = "keep") +
  plot_annotation(title = "Selected motif activity on Bdev umap.atac",
    subtitle = paste("Native Bdev chromVAR_data; gene-specific symmetric",
      "blue-white-red scales clipped at absolute p99 for display"),
    caption = "The Bdev bundle and embedding were opened read-only.",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 8.5)))
ggsave(file.path(stage_dir, "figures",
  "selected_motif_activity_bdev_umap.png"), umap_figure,
  width = 14, height = 9.5, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "selected_motif_activity_bdev_umap.pdf"), umap_figure,
  width = 14, height = 9.5, device = cairo_pdf, bg = "white")

# Load exact source motif matches and bias used by the preceding motif report.
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
stopifnot(all(motifs$motif_id %in% colnames(motif_matrix_all)))

joint_counts <- BPCells::open_matrix_dir(joint_matrix_dir)
joint_peak_ids <- rownames(joint_counts)
joint_cell_ids <- colnames(joint_counts)
stopifnot(length(joint_peak_ids) == 97788L,
  all(joint_peak_ids %in% rownames(motif_matrix_all)),
  all(joint_peak_ids %in% peak_covariates$feature_id))
motif_matrix_all <- motif_matrix_all[joint_peak_ids, motifs$motif_id,
  drop = FALSE]
peak_bias <- peak_covariates$GC.percent[
  match(joint_peak_ids, peak_covariates$feature_id)]
stopifnot(all(is.finite(peak_bias)))
rm(source_bdev, source_env, peak_covariates)
invisible(gc())

all_query_values <- list()
all_adjusted <- list()
violin_plots <- list()
region_diagnostics <- list()
for (region_index in seq_len(nrow(regions))) {
  region <- regions[region_index]
  region_key <- region$region_id
  edges <- seq(region$pt_left, region$pt_right,
    length.out = region$n_bins + 1L)
  assign_bin <- function(pt) as.integer(cut(pt, breaks = edges,
    include.lowest = TRUE, right = FALSE, labels = FALSE))
  bdev_region <- bdev_meta[is.finite(pt) & pt >= region$pt_left &
    pt < region$pt_right]
  bdev_region[, `:=`(cell_id = paste0("bdev::", source_cell_id),
    bin = assign_bin(pt))]
  query_region <- assignments[region_id == region_key]
  query_region[, `:=`(mapped_pt = as.numeric(mapped_pt),
    group = ifelse(genotype == "WT", "WT", "KO"),
    bin = assign_bin(mapped_pt))]
  stopifnot(!anyNA(bdev_region$bin), !anyNA(query_region$bin),
    all(bdev_region$cell_id %in% joint_cell_ids),
    all(query_region$cell_id %in% joint_cell_ids))

  bdev_counts <- as(joint_counts[, bdev_region$cell_id, drop = FALSE],
    "dgCMatrix")
  peak_keep <- Matrix::rowSums(bdev_counts) >= 10
  counts_kept_bdev <- bdev_counts[peak_keep, , drop = FALSE]
  counts_kept_query <- as(joint_counts[peak_keep, query_region$cell_id,
    drop = FALSE], "dgCMatrix")
  combined_counts <- cbind(counts_kept_bdev, counts_kept_query)
  motif_matrix <- motif_matrix_all[peak_keep, , drop = FALSE]
  motif_peak_counts <- Matrix::colSums(motif_matrix != 0)
  stopifnot(all(motif_peak_counts >= 10L))
  set.seed(seed + region_index)
  background_peaks <- chromVAR::getBackgroundPeaks(counts_kept_bdev,
    bias = peak_bias[peak_keep], niterations = 50L, w = 0.1, bs = 50L)
  bdev_expectation <- chromVAR::computeExpectations(counts_kept_bdev)
  query_deviation_object <- chromVAR::computeDeviations(combined_counts,
    motif_matrix, background_peaks = background_peaks,
    expectation = bdev_expectation)
  query_deviation <- SummarizedExperiment::assay(query_deviation_object,
    "deviations")[, query_region$cell_id, drop = FALSE]
  stopifnot(identical(rownames(query_deviation), motifs$motif_id),
    identical(colnames(query_deviation), query_region$cell_id),
    all(is.finite(query_deviation)))

  long <- as.data.table(as.table(query_deviation))
  setnames(long, c("motif_id", "cell_id", "chromvar_deviation"))
  long[, `:=`(motif_id = as.character(motif_id),
    cell_id = as.character(cell_id))]
  long <- merge(long, motifs[, .(motif_id, motif_label, motif_order)],
    by = "motif_id", all.x = TRUE, sort = FALSE)
  long <- merge(long, query_region[, .(cell_id, mapped_pt, group, bin)],
    by = "cell_id", all.x = TRUE, sort = FALSE)
  long[, region_id := region_key]

  bin_means <- long[, .(mean_deviation = mean(chromvar_deviation),
    n_cells = .N), by = .(motif_id, group, bin)]
  wt <- bin_means[group == "WT", .(motif_id, bin,
    wt_mean = mean_deviation, wt_n = n_cells)]
  ko <- bin_means[group == "KO", .(motif_id, bin,
    ko_mean = mean_deviation, ko_n = n_cells)]
  paired <- merge(wt, ko, by = c("motif_id", "bin"), sort = FALSE)
  bdev_bin_n <- bdev_region[, .(bdev_n = .N), by = bin]
  paired <- merge(paired, bdev_bin_n, by = "bin", all.x = TRUE,
    sort = FALSE)
  paired[, supported := wt_n >= 20L & ko_n >= 20L]
  adjusted <- paired[supported == TRUE, .(
    adjusted_wt_minus_ko = sum((wt_mean - ko_mean) * bdev_n) /
      sum(bdev_n), common_bins = .N,
    common_bdev_cells = sum(bdev_n)), by = motif_id]
  adjusted <- merge(adjusted,
    motifs[, .(motif_id, motif_label, motif_order)],
    by = "motif_id", all.x = TRUE, sort = FALSE)
  adjusted[, region_id := region_key]
  all_adjusted[[region_key]] <- adjusted
  long <- merge(long, adjusted[, .(motif_id, adjusted_wt_minus_ko)],
    by = "motif_id", all.x = TRUE, sort = FALSE)
  all_query_values[[region_key]] <- copy(long)

  facet_labels <- adjusted[order(motif_order), setNames(
    sprintf("%s\nadjusted WT−KO = %+.3f", motif_label,
      adjusted_wt_minus_ko), motif_label)]
  long[, motif_label_factor := factor(motif_label,
    levels = motifs$motif_label)]
  long[, group_label := if (region$target == "Flt3l")
    factor(group, levels = c("WT", "KO"),
      labels = c("FLT3L WT", "FLT3LKO")) else
    factor(group, levels = c("WT", "KO"),
      labels = c("IL-7 WT", "IL7KO"))]
  p <- ggplot(long, aes(group_label, chromvar_deviation,
      fill = group_label)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.35,
      alpha = 0.75) +
    geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.35,
      fill = "white") +
    facet_wrap(~ motif_label_factor, nrow = 1, scales = "free_y",
      labeller = as_labeller(facet_labels)) +
    scale_fill_manual(values = c("FLT3L WT" = "#D73027",
      "FLT3LKO" = "#4575B4", "IL-7 WT" = "#D73027",
      "IL7KO" = "#4575B4")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none", plot.title = element_text(face = "bold")) +
    labs(title = if (region$target == "Flt3l")
        "Flt3L peak 1: per-cell motif activity" else
        "IL-7 peak 2: per-cell motif activity",
      subtitle = paste0("Raw violins; labels give matched-bin WT−KO weighted to Bdev (",
        nrow(query_region[group == "WT"]), " WT, ",
        nrow(query_region[group == "KO"]), " KO cells)"),
      x = NULL, y = "chromVAR deviation")
  violin_plots[[region_key]] <- p
  region_diagnostics[[region_key]] <- list(
    bdev_cells = nrow(bdev_region),
    wt_cells = nrow(query_region[group == "WT"]),
    ko_cells = nrow(query_region[group == "KO"]),
    bdev_accessible_peaks = sum(peak_keep),
    motif_peak_counts = as.list(setNames(as.integer(motif_peak_counts),
      motifs$motif_label)),
    common_supported_bins = as.list(setNames(adjusted$common_bins,
      adjusted$motif_label)))
  rm(bdev_counts, counts_kept_bdev, counts_kept_query, combined_counts,
    background_peaks, bdev_expectation, query_deviation_object,
    query_deviation)
  invisible(gc())
}

violin_figure <- violin_plots[["Flt3l_peak1"]] /
  violin_plots[["Il7_peak2"]] +
  plot_annotation(title = "Selected motif activity in WT and knockout cells",
    subtitle = paste("Each comparison uses cells in its respective enriched",
      "pseudotime region and a region-specific Bdev-defined chromVAR model"),
    caption = paste("Violins show raw per-cell distributions. Adjusted labels",
      "compare matched pseudotime bins standardized to the Bdev distribution;",
      "one biological sample per genotype."),
    theme = theme(plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 8.5)))
ggsave(file.path(stage_dir, "figures",
  "selected_motif_activity_wt_ko_violins.png"), violin_figure,
  width = 16, height = 10, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "selected_motif_activity_wt_ko_violins.pdf"), violin_figure,
  width = 16, height = 10, device = cairo_pdf, bg = "white")

if (preset == "tcf3_ebf1_zeb_stat5_runx1") {
  query_all <- rbindlist(all_query_values, fill = TRUE)
  adjusted_all <- rbindlist(all_adjusted, fill = TRUE)
  composite_violin_plots <- vector("list", nrow(motifs))
  for (i in seq_len(nrow(motifs))) {
    motif <- motifs[i]
    d <- query_all[motif_id == motif$motif_id]
    d[, region_label := factor(region_id,
      levels = c("Flt3l_peak1", "Il7_peak2"),
      labels = c("Flt3L peak 1", "IL-7 peak 2"))]
    d[, comparison_group := factor(group, levels = c("WT", "KO"))]
    effects <- adjusted_all[motif_id == motif$motif_id]
    flt_effect <- effects[region_id == "Flt3l_peak1",
      adjusted_wt_minus_ko]
    il7_effect <- effects[region_id == "Il7_peak2",
      adjusted_wt_minus_ko]
    composite_violin_plots[[i]] <- ggplot(d,
        aes(comparison_group, chromvar_deviation,
          fill = comparison_group)) +
      geom_violin(scale = "width", trim = TRUE, linewidth = 0.3,
        alpha = 0.75) +
      geom_boxplot(width = 0.11, outlier.shape = NA, linewidth = 0.3,
        fill = "white") +
      facet_wrap(~ region_label, nrow = 1, scales = "free_y") +
      scale_fill_manual(values = c(WT = "#D73027", KO = "#4575B4")) +
      theme_bw(base_size = 9) +
      theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 8.5),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8),
        legend.position = "none") +
      labs(title = motif$motif_label,
        subtitle = sprintf("Adjusted WT−KO: Flt3L %+.3f; IL-7 %+.3f",
          flt_effect, il7_effect),
        x = NULL, y = "chromVAR deviation")
  }
  left_block <- wrap_plots(umap_plots, ncol = 1)
  right_block <- wrap_plots(composite_violin_plots, ncol = 1)
  composite_figure <- (left_block | right_block) +
    plot_layout(widths = c(1, 2.15)) +
    plot_annotation(
      title = "Bdev motif UMAPs and corresponding WT/KO distributions",
      subtitle = paste("UMAPs use native Bdev chromVAR; violins use",
        "region-specific Bdev-calibrated query-cell deviations"),
      caption = paste("ZEB1 is shown as a ZEB-family proxy because this",
        "JASPAR set contains no ZEB2 motif. Violin subtitles give matched-bin",
        "WT−KO effects weighted to Bdev; one sample per genotype."),
      theme = theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 10),
        plot.caption = element_text(size = 8.5)))
  ggsave(file.path(stage_dir, "figures",
    "tcf3_ebf1_zeb_stat5_runx1_umap_left_violins_right.png"),
    composite_figure, width = 16, height = 22, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures",
    "tcf3_ebf1_zeb_stat5_runx1_umap_left_violins_right.pdf"),
    composite_figure, width = 16, height = 22, device = cairo_pdf,
    bg = "white")
}

write_tsv(motifs, file.path(stage_dir, "selected_motifs.tsv"))
fwrite(rbindlist(umap_values), file.path(stage_dir,
  "bdev_umap_motif_activity.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")
write_tsv(rbindlist(umap_scales), file.path(stage_dir,
  "bdev_umap_color_scales.tsv"))
fwrite(rbindlist(all_query_values), file.path(stage_dir,
  "query_cell_motif_activity.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")
write_tsv(rbindlist(all_adjusted), file.path(stage_dir,
  "trajectory_adjusted_wt_ko_motif_effects.tsv"))
write_json(list(status = "complete_preliminary_unregistered",
  bdev_bundle_fingerprint = expected_bdev_fingerprint,
  mapping_bundle_fingerprint = expected_mapping_fingerprint,
  motif_preset = preset,
  motifs = as.list(setNames(motifs$motif_id, motifs$motif_label)),
  motif_selection_notes = as.list(setNames(motifs$selection_note,
    motifs$motif_label)),
  bdev_umap = list(matrix = "native chromvar_data",
    space = "umap.atac", scale = "gene-specific symmetric absolute p99"),
  query_violin = list(
    values = "per-cell chromVAR deviations from region-specific Bdev-defined accessible peaks, background, and expectation",
    raw_distribution = TRUE,
    adjusted_annotation = "matched-bin WT-minus-KO weighted to Bdev cell distribution"),
  region_diagnostics = region_diagnostics,
  biological_replicates_per_genotype = 1L,
  inference = "descriptive", bundle_mutated = FALSE,
  capsule_registered = FALSE), file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Selected motif Bdev UMAP and WT/KO violin plots", "",
  paste0("The selected motifs are ", paste(sprintf("%s (%s)",
    motifs$motif_label, motifs$motif_id), collapse = ", "), "."),
  if (preset == "tcf3_ebf1_zeb_stat5_runx1")
    "ZEB1 is explicitly used as a ZEB-family proxy because the motif set contains no ZEB2 entry." else NULL,
  "Bdev UMAP panels use the native chromvar_data matrix on the unchanged umap.atac embedding.",
  "WT/KO violins use per-cell deviations recomputed independently in Flt3L peak 1 and IL-7 peak 2 with Bdev-defined accessible peaks, background peaks, and expectation.",
  "Violin distributions are raw; facet labels also report the matched-bin WT-minus-KO effect weighted to the Bdev pseudotime distribution.",
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
cat("Completed selected motif UMAP and violins at", output_dir, "\n")
