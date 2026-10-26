#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) stop(paste(
  "Usage: 42_bdev_selected_motif_pt_heatmaps.R",
  "<bdev_bundle> <mapping_bundle> <joint_peak_matrix_dir>",
  "<wt_enriched_report> <bdev_rdata> <output_dir>",
  "<expected_bdev_fingerprint> <expected_mapping_fingerprint>",
  "<top_n_motifs> <seed>"))
bdev_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
mapping_bundle_dir <- normalizePath(args[[2L]], mustWork = TRUE)
joint_matrix_dir <- normalizePath(args[[3L]], mustWork = TRUE)
region_report <- normalizePath(args[[4L]], mustWork = TRUE)
bdev_rdata_path <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
expected_bdev_fingerprint <- args[[7L]]
expected_mapping_fingerprint <- args[[8L]]
top_n <- as.integer(args[[9L]])
seed <- as.integer(args[[10L]])
stopifnot(top_n >= 5L, is.finite(seed))
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
log_event <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(..., collapse = " "), "\n")

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
  all(assignments$support_status %in% c("supported", "borderline")),
  !anyDuplicated(assignments[, .(region_id, cell_id)]))

bdev_cells <- wb_cells(bdev, "set:bdev_reference")
bdev_meta <- as.data.table(wb_meta(bdev, cells = bdev_cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt")))
setnames(bdev_meta, names(bdev_meta)[[1L]], "source_cell_id")
stopifnot(identical(as.character(bdev_meta$source_cell_id), bdev_cells))

log_event("Loading exact source Bdev JASPAR motif matches and peak bias")
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
motif_ids_all <- colnames(motif_matrix_all)
motif_names <- tryCatch(Signac::ConvertMotifID(source_bdev[["peaks"]],
  id = motif_ids_all), error = function(e) setNames(motif_ids_all,
    motif_ids_all))
if (is.null(names(motif_names))) names(motif_names) <- motif_ids_all

joint_counts <- BPCells::open_matrix_dir(joint_matrix_dir)
joint_peak_ids <- rownames(joint_counts)
joint_cell_ids <- colnames(joint_counts)
stopifnot(length(joint_peak_ids) == 97788L,
  all(joint_peak_ids %in% rownames(motif_matrix_all)),
  all(joint_peak_ids %in% peak_covariates$feature_id),
  !anyDuplicated(joint_peak_ids), !anyDuplicated(joint_cell_ids))
motif_matrix_all <- motif_matrix_all[joint_peak_ids, , drop = FALSE]
peak_bias <- peak_covariates$GC.percent[
  match(joint_peak_ids, peak_covariates$feature_id)]
stopifnot(all(is.finite(peak_bias)))
rm(source_bdev, source_env, peak_covariates)
invisible(gc())

all_variability <- list()
all_selected <- list()
all_values <- list()
all_eligible_values <- list()
all_counts <- list()
region_diagnostics <- list()

for (region_index in seq_len(nrow(regions))) {
  region <- regions[region_index]
  region_key <- region$region_id
  log_event("Preparing", region_key)
  edges <- seq(region$pt_left, region$pt_right,
    length.out = region$n_bins + 1L)
  bin_midpoints <- (head(edges, -1L) + tail(edges, -1L)) / 2
  assign_bin <- function(pt) as.integer(cut(pt, breaks = edges,
    include.lowest = TRUE, right = FALSE, labels = FALSE))

  bdev_region <- bdev_meta[is.finite(pt) & pt >= region$pt_left &
    pt < region$pt_right]
  bdev_region[, `:=`(cell_id = paste0("bdev::", source_cell_id),
    group = "Bdev", bin = assign_bin(pt))]
  query_region <- assignments[region_id == region_key]
  query_region[, `:=`(pt = mapped_pt,
    group = ifelse(genotype == "WT", "WT", "KO"),
    bin = assign_bin(mapped_pt))]
  cell_table <- rbindlist(list(
    bdev_region[, .(cell_id, group, bin, pt)],
    query_region[, .(cell_id, group, bin, pt)]))
  stopifnot(!anyNA(cell_table$bin), !anyDuplicated(cell_table$cell_id),
    all(cell_table$cell_id %in% joint_cell_ids))
  group_levels <- c("Bdev", "WT", "KO")
  column_table <- CJ(group = group_levels, bin = seq_len(region$n_bins),
    unique = TRUE)
  column_table[, column_id := paste(group, sprintf("%02d", bin), sep = "::")]
  cell_table <- merge(cell_table, column_table[, .(group, bin, column_id)],
    by = c("group", "bin"), all.x = TRUE, sort = FALSE)
  column_table[, `:=`(pt_left = edges[bin], pt_right = edges[bin + 1L],
    pt_mid = bin_midpoints[bin])]
  count_table <- cell_table[, .(n_cells = .N), by = .(group, bin)]
  column_table <- merge(column_table, count_table,
    by = c("group", "bin"), all.x = TRUE, sort = FALSE)
  column_table[is.na(n_cells), n_cells := 0L]
  column_table[, group_order := match(group, group_levels)]
  setorder(column_table, group_order, bin)
  column_table[, group_order := NULL]
  # Sparse tails are displayed as missing rather than estimating unstable
  # deviations from only a handful of cells. Retain broad Bdev region coverage.
  column_table[, included := n_cells >= 20L]
  stopifnot(column_table[group == "Bdev", sum(included)] >=
    ceiling(2 * region$n_bins / 3))
  column_table[, region_id := region_key]
  all_counts[[region_key]] <- copy(column_table)
  analysis_columns <- column_table[included == TRUE]
  cell_table <- cell_table[column_id %in% analysis_columns$column_id]

  design <- sparseMatrix(
    i = match(cell_table$cell_id, joint_cell_ids),
    j = match(cell_table$column_id, analysis_columns$column_id),
    x = 1, dims = c(length(joint_cell_ids), nrow(analysis_columns)),
    dimnames = list(joint_cell_ids, analysis_columns$column_id))
  log_event("Aggregating common-peak counts into", ncol(design),
    "pseudotime pseudobulks")
  aggregate_counts <- as(joint_counts %*% design, "dgCMatrix")
  stopifnot(identical(rownames(aggregate_counts), joint_peak_ids),
    identical(colnames(aggregate_counts), analysis_columns$column_id))

  bdev_columns <- which(analysis_columns$group == "Bdev")
  bdev_aggregate <- aggregate_counts[, bdev_columns, drop = FALSE]
  # Peak and motif eligibility are determined from Bdev only.
  peak_keep <- Matrix::rowSums(bdev_aggregate) >= 10
  motif_matrix <- motif_matrix_all[peak_keep, , drop = FALSE]
  motif_peak_counts <- Matrix::colSums(motif_matrix != 0)
  motif_keep <- motif_peak_counts >= 10L
  motif_matrix <- motif_matrix[, motif_keep, drop = FALSE]
  motif_peak_counts <- motif_peak_counts[motif_keep]
  counts_kept <- aggregate_counts[peak_keep, , drop = FALSE]
  bdev_counts_kept <- counts_kept[, bdev_columns, drop = FALSE]
  bias_kept <- peak_bias[peak_keep]
  stopifnot(nrow(counts_kept) > 10000L, ncol(motif_matrix) > top_n)

  log_event("Computing Bdev-defined chromVAR background and expectations for",
    region_key, "using", nrow(counts_kept), "Bdev-accessible peaks and",
    ncol(motif_matrix), "eligible motifs")
  set.seed(seed + region_index)
  background_peaks <- chromVAR::getBackgroundPeaks(bdev_counts_kept,
    bias = bias_kept, niterations = 50L, w = 0.1, bs = 50L)
  bdev_expectation <- chromVAR::computeExpectations(bdev_counts_kept)
  deviation_object <- chromVAR::computeDeviations(counts_kept,
    motif_matrix, background_peaks = background_peaks,
    expectation = bdev_expectation)
  raw_deviation <- SummarizedExperiment::assay(deviation_object,
    "deviations")
  stopifnot(identical(colnames(raw_deviation), analysis_columns$column_id),
    all(is.finite(raw_deviation)))

  bdev_deviation <- raw_deviation[, bdev_columns, drop = FALSE]
  bdev_sd <- apply(bdev_deviation, 1L, sd)
  bdev_range <- apply(bdev_deviation, 1L, function(x) diff(range(x)))
  bdev_iqr <- apply(bdev_deviation, 1L, IQR)
  variability <- data.table(region_id = region_key,
    motif_id = rownames(raw_deviation),
    motif_name = as.character(motif_names[rownames(raw_deviation)]),
    bdev_temporal_sd = bdev_sd, bdev_temporal_range = bdev_range,
    bdev_temporal_iqr = bdev_iqr,
    bdev_accessible_peak_count = as.integer(
      motif_peak_counts[rownames(raw_deviation)]))
  variability[is.na(motif_name) | motif_name == "", motif_name := motif_id]
  setorder(variability, -bdev_temporal_sd, -bdev_temporal_range, motif_id)
  variability[, bdev_variability_rank := seq_len(.N)]
  variability[, selected := bdev_variability_rank <= top_n]
  selected <- variability[selected == TRUE]
  all_variability[[region_key]] <- variability
  all_selected[[region_key]] <- selected

  all_centers <- rowMeans(raw_deviation[, bdev_columns, drop = FALSE])
  all_scales <- apply(raw_deviation[, bdev_columns, drop = FALSE], 1L, sd)
  trajectory_eligible <- is.finite(all_scales) & all_scales > 0
  all_scaled <- sweep(sweep(raw_deviation[trajectory_eligible, , drop = FALSE],
    1L, all_centers[trajectory_eligible], "-"), 1L,
    all_scales[trajectory_eligible], "/")
  eligible_long <- as.data.table(as.table(all_scaled))
  setnames(eligible_long,
    c("motif_id", "column_id", "bdev_scaled_deviation"))
  eligible_long[, `:=`(motif_id = as.character(motif_id),
    column_id = as.character(column_id))]
  eligible_long <- merge(eligible_long,
    analysis_columns[, .(column_id, group, bin, pt_mid, n_cells)],
    by = "column_id", all.x = TRUE, sort = FALSE)
  eligible_long <- merge(eligible_long,
    variability[, .(motif_id, motif_name, bdev_temporal_sd,
      bdev_variability_rank)], by = "motif_id", all.x = TRUE, sort = FALSE)
  eligible_raw <- as.data.table(as.table(
    raw_deviation[trajectory_eligible, , drop = FALSE]))
  setnames(eligible_raw,
    c("motif_id", "column_id", "raw_chromvar_deviation"))
  eligible_raw[, `:=`(motif_id = as.character(motif_id),
    column_id = as.character(column_id))]
  eligible_long <- merge(eligible_long, eligible_raw,
    by = c("motif_id", "column_id"), all.x = TRUE, sort = FALSE)
  eligible_long[, region_id := region_key]
  all_eligible_values[[region_key]] <- eligible_long

  selected_raw <- raw_deviation[selected$motif_id, , drop = FALSE]
  centers <- rowMeans(selected_raw[, bdev_columns, drop = FALSE])
  scales <- apply(selected_raw[, bdev_columns, drop = FALSE], 1L, sd)
  stopifnot(all(is.finite(scales) & scales > 0))
  scaled <- sweep(sweep(selected_raw, 1L, centers, "-"), 1L, scales, "/")
  bdev_scaled <- scaled[, bdev_columns, drop = FALSE]
  row_order <- rownames(bdev_scaled)[hclust(dist(bdev_scaled))$order]
  label_table <- selected[, .(motif_id,
    motif_label = sprintf("%s (%s)", motif_name, motif_id))]
  label_table[, motif_label := make.unique(motif_label)]

  heat <- as.data.table(as.table(scaled))
  setnames(heat, c("motif_id", "column_id", "bdev_scaled_deviation"))
  heat[, `:=`(motif_id = as.character(motif_id),
    column_id = as.character(column_id))]
  heat <- merge(heat, analysis_columns[, .(column_id, group, bin, pt_mid,
    n_cells)], by = "column_id", all.x = TRUE, sort = FALSE)
  heat <- merge(heat, label_table, by = "motif_id", all.x = TRUE,
    sort = FALSE)
  raw_long <- as.data.table(as.table(selected_raw))
  setnames(raw_long, c("motif_id", "column_id", "raw_chromvar_deviation"))
  raw_long[, `:=`(motif_id = as.character(motif_id),
    column_id = as.character(column_id))]
  heat <- merge(heat, raw_long, by = c("motif_id", "column_id"),
    all.x = TRUE, sort = FALSE)

  wt <- heat[group == "WT", .(motif_id, motif_label, bin, pt_mid,
    wt = bdev_scaled_deviation)]
  ko <- heat[group == "KO", .(motif_id, bin, ko = bdev_scaled_deviation)]
  difference <- merge(wt, ko, by = c("motif_id", "bin"), sort = FALSE)
  difference[, `:=`(group = "WT - KO", n_cells = NA_integer_,
    bdev_scaled_deviation = wt - ko,
    raw_chromvar_deviation = NA_real_)]
  difference[, c("wt", "ko") := NULL]
  heat <- rbindlist(list(heat, difference), fill = TRUE)
  heat[, region_id := region_key]
  group_labels <- if (region$target == "Flt3l")
    c(Bdev = "Bdev reference", WT = "FLT3L WT", KO = "FLT3LKO",
      `WT - KO` = "WT - FLT3LKO") else
    c(Bdev = "Bdev reference", WT = "IL-7 WT", KO = "IL7KO",
      `WT - KO` = "WT - IL7KO")
  heat[, group_label := group_labels[group]]
  heat[, group_label := factor(group_label, levels = unname(group_labels))]
  ordered_labels <- label_table$motif_label[
    match(row_order, label_table$motif_id)]
  heat[, motif_label := factor(motif_label, levels = rev(ordered_labels))]
  all_values[[region_key]] <- copy(heat)

  plot_limit <- 2.5
  p <- ggplot(heat, aes(pt_mid, motif_label,
      fill = bdev_scaled_deviation)) +
    geom_tile() +
    facet_grid(. ~ group_label) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-plot_limit, plot_limit),
      oob = scales::squish, name = "Deviation\n(Bdev SD)") +
    scale_x_continuous(expand = expansion(mult = 0),
      breaks = pretty(c(region$pt_left, region$pt_right), n = 4)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6.5),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 9.5),
      legend.position = "right") +
    labs(title = paste0(region$cytokine,
        " enriched-region motif accessibility across pseudotime"),
      subtitle = paste0("Top ", top_n,
        " motifs selected, ordered, centered, and scaled using Bdev only; ",
        "query data do not affect selection"),
      x = "Bdev pseudotime", y = NULL,
      caption = paste0(region_key, ": pt [",
        sprintf("%.3f", region$pt_left), ", ",
        sprintf("%.3f", region$pt_right), "). Raw chromVAR deviations from ",
        "Bdev-defined backgrounds; WT-KO is descriptive."))
  file_stub <- if (region$target == "Flt3l")
    "flt3l_peak1_bdev_selected_motif_heatmap" else
    "il7_peak2_bdev_selected_motif_heatmap"
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".png")), p,
    width = 18, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".pdf")), p,
    width = 18, height = 10, device = cairo_pdf, bg = "white")

  group_cell_totals <- column_table[, .(n_cells = sum(n_cells)), by = group]
  region_diagnostics[[region_key]] <- list(
    region_id = region_key, pt_left = region$pt_left,
    pt_right = region$pt_right, bins = region$n_bins,
    cells = as.list(setNames(group_cell_totals$n_cells,
      group_cell_totals$group)),
    bdev_accessible_peaks = sum(peak_keep),
    eligible_motifs = ncol(motif_matrix), selected_motifs = top_n,
    minimum_cells_per_displayed_bin = 20L,
    sparse_query_bins_are_missing = TRUE,
    selection = "top temporal SD of raw chromVAR deviations across Bdev pseudotime bins",
    background_and_expectation = "Bdev-only",
    centering_scaling_and_row_order = "Bdev-only")
  rm(design, aggregate_counts, counts_kept, bdev_counts_kept,
    deviation_object, raw_deviation, selected_raw, scaled, bdev_scaled,
    all_scaled, eligible_long, eligible_raw)
  invisible(gc())
}

write_tsv(rbindlist(all_variability), file.path(stage_dir,
  "all_motif_bdev_temporal_variability.tsv"))
write_tsv(rbindlist(all_selected), file.path(stage_dir,
  "selected_bdev_variable_motifs.tsv"))
fwrite(rbindlist(all_values), file.path(stage_dir,
  "selected_motif_pseudotime_heatmap_values.tsv.gz"), sep = "\t",
  na = "NA", compress = "gzip")
fwrite(rbindlist(all_eligible_values), file.path(stage_dir,
  "all_eligible_motif_pseudotime_values.tsv.gz"), sep = "\t",
  na = "NA", compress = "gzip")
write_tsv(rbindlist(all_counts), file.path(stage_dir,
  "pseudotime_bin_cell_counts.tsv"))

diagnostics <- list(status = "complete_preliminary_unregistered",
  bdev_bundle_id = bdev_manifest$bundle_id,
  bdev_bundle_fingerprint = expected_bdev_fingerprint,
  mapping_bundle_id = mapping_manifest$bundle_id,
  mapping_bundle_fingerprint = expected_mapping_fingerprint,
  peak_matrix = "ATAC_bdev_mm10_peaks",
  common_peak_universe = length(joint_peak_ids),
  motif_source = "exact read-only Bdev JASPAR match matrix",
  top_n_motifs = top_n, selection_uses_query = FALSE,
  region_diagnostics = region_diagnostics,
  biological_replicates_per_genotype = 1L,
  inference = "descriptive",
  bundle_mutated = FALSE, capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Bdev-selected motif activity across cytokine-enriched pseudotime regions", "",
  "Separate heatmaps cover Flt3L peak 1 and IL-7 peak 2.",
  "For each region, Bdev alone determines accessible peaks, chromVAR background and expectation, motif temporal-variability ranking, the top-30 selection, row order, centering, and scaling.",
  "WT and KO cells contribute only after these quantities are fixed, using their mapped pseudotime.",
  "The standardized pseudotime values for every Bdev-eligible motif with nonzero temporal SD are exported for downstream directional screens; the displayed heatmaps remain restricted to the top 30.",
  "The fourth block is the descriptive WT-minus-KO difference on the Bdev-scaled axis.",
  "Because each genotype has one biological sample, these displays are descriptive rather than replicate-level inference."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
wb_close_bundle(bdev)
bdev <- NULL
wb_close_bundle(mapping)
mapping <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed Bdev-selected motif pseudotime heatmaps at", output_dir, "\n")
