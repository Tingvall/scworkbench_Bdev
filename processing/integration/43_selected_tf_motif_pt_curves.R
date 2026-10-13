#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) stop(paste(
  "Usage: 43_selected_tf_motif_pt_curves.R",
  "<bdev_bundle> <mapping_bundle> <joint_peak_matrix_dir>",
  "<wt_enriched_report> <bdev_rdata> <motif_heatmap_report>",
  "<output_dir> <expected_bdev_fingerprint>",
  "<expected_mapping_fingerprint> <minimum_cells> <seed>"))
bdev_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
mapping_bundle_dir <- normalizePath(args[[2L]], mustWork = TRUE)
joint_matrix_dir <- normalizePath(args[[3L]], mustWork = TRUE)
region_report <- normalizePath(args[[4L]], mustWork = TRUE)
bdev_rdata_path <- normalizePath(args[[5L]], mustWork = TRUE)
motif_report <- normalizePath(args[[6L]], mustWork = TRUE)
output_dir <- normalizePath(args[[7L]], mustWork = FALSE)
expected_bdev_fingerprint <- args[[8L]]
expected_mapping_fingerprint <- args[[9L]]
minimum_cells <- as.integer(args[[10L]])
seed <- as.integer(args[[11L]])
stopifnot(minimum_cells >= 5L, is.finite(seed))
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
variability <- fread(file.path(motif_report,
  "all_motif_bdev_temporal_variability.tsv"))
stopifnot(nrow(regions) == 2L,
  setequal(regions$region_id, c("Flt3l_peak1", "Il7_peak2")),
  !anyDuplicated(assignments[, .(region_id, cell_id)]),
  setequal(unique(variability$region_id), regions$region_id))

# The same motif ID is used in both cytokine regions. For each requested TF or
# family, rank matching JASPAR motifs by mean Bdev temporal SD across regions.
requested <- data.table(
  tf = c("EBF1", "PAX5", "TCF3", "CEBP", "HLF", "BACH2", "LEF1"),
  display_order = seq_len(7L))
motif_candidates <- unique(variability[, .(motif_id, motif_name)])
match_tf <- function(tf, motif_name) {
  if (tf == "CEBP") return(grepl("^CEBP[A-Z]?(\\(var\\.[0-9]+\\))?$",
    motif_name, ignore.case = TRUE))
  grepl(paste0("^", tf, "(\\(var\\.[0-9]+\\))?$"), motif_name,
    ignore.case = TRUE)
}
candidate_map <- rbindlist(lapply(requested$tf, function(tf_name) {
  motif_candidates[vapply(motif_name, function(x) match_tf(tf_name, x),
    logical(1L))][, tf := tf_name]
}))
stopifnot(setequal(unique(candidate_map$tf), requested$tf))
candidate_scores <- merge(candidate_map, variability[, .(
  mean_bdev_temporal_sd = mean(bdev_temporal_sd),
  minimum_bdev_temporal_sd = min(bdev_temporal_sd)),
  by = .(motif_id, motif_name)], by = c("motif_id", "motif_name"))
setorder(candidate_scores, tf, -mean_bdev_temporal_sd,
  -minimum_bdev_temporal_sd, motif_id)
selection <- candidate_scores[, .SD[1L], by = tf]
selection <- merge(requested, selection, by = "tf", sort = FALSE)
setorder(selection, display_order)
stopifnot(nrow(selection) == 7L, !anyDuplicated(selection$motif_id))
selection[, selection_rule :=
  "highest mean Bdev temporal SD across Flt3l_peak1 and Il7_peak2 among matching JASPAR motifs"]
write_tsv(selection, file.path(stage_dir, "selected_tf_motifs.tsv"))
write_tsv(candidate_scores, file.path(stage_dir,
  "candidate_tf_motif_ranking.tsv"))

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
stopifnot(all(selection$motif_id %in% colnames(motif_matrix_all)))

joint_counts <- BPCells::open_matrix_dir(joint_matrix_dir)
joint_peak_ids <- rownames(joint_counts)
joint_cell_ids <- colnames(joint_counts)
stopifnot(length(joint_peak_ids) == 97788L,
  all(joint_peak_ids %in% rownames(motif_matrix_all)),
  all(joint_peak_ids %in% peak_covariates$feature_id),
  !anyDuplicated(joint_peak_ids), !anyDuplicated(joint_cell_ids))
motif_matrix_all <- motif_matrix_all[joint_peak_ids,
  selection$motif_id, drop = FALSE]
peak_bias <- peak_covariates$GC.percent[
  match(joint_peak_ids, peak_covariates$feature_id)]
stopifnot(all(is.finite(peak_bias)))
rm(source_bdev, source_env, peak_covariates)
invisible(gc())

all_values <- list()
all_counts <- list()
region_diagnostics <- list()
group_levels <- c("Bdev", "WT", "KO")
group_colors <- c(Bdev = "black", WT = "#D73027", KO = "#2166AC")

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

  column_table <- CJ(group = group_levels, bin = seq_len(region$n_bins),
    unique = TRUE)
  column_table[, column_id := paste(group, sprintf("%02d", bin), sep = "::")]
  column_table[, `:=`(pt_left = edges[bin], pt_right = edges[bin + 1L],
    pt_mid = bin_midpoints[bin])]
  cell_table <- merge(cell_table, column_table[, .(group, bin, column_id)],
    by = c("group", "bin"), all.x = TRUE, sort = FALSE)
  count_table <- cell_table[, .(n_cells = .N), by = .(group, bin)]
  column_table <- merge(column_table, count_table,
    by = c("group", "bin"), all.x = TRUE, sort = FALSE)
  column_table[is.na(n_cells), n_cells := 0L]
  column_table[, `:=`(included = n_cells >= minimum_cells,
    region_id = region_key)]
  column_table[, group_order := match(group, group_levels)]
  setorder(column_table, group_order, bin)
  column_table[, group_order := NULL]
  stopifnot(column_table[group == "Bdev", sum(included)] >=
    ceiling(2 * region$n_bins / 3))
  all_counts[[region_key]] <- copy(column_table)
  analysis_columns <- column_table[included == TRUE]
  analysis_cells <- cell_table[column_id %in% analysis_columns$column_id]

  design <- sparseMatrix(
    i = match(analysis_cells$cell_id, joint_cell_ids),
    j = match(analysis_cells$column_id, analysis_columns$column_id),
    x = 1, dims = c(length(joint_cell_ids), nrow(analysis_columns)),
    dimnames = list(joint_cell_ids, analysis_columns$column_id))
  aggregate_counts <- as(joint_counts %*% design, "dgCMatrix")
  stopifnot(identical(rownames(aggregate_counts), joint_peak_ids),
    identical(colnames(aggregate_counts), analysis_columns$column_id))

  bdev_columns <- which(analysis_columns$group == "Bdev")
  bdev_aggregate <- aggregate_counts[, bdev_columns, drop = FALSE]
  peak_keep <- Matrix::rowSums(bdev_aggregate) >= 10
  counts_kept <- aggregate_counts[peak_keep, , drop = FALSE]
  bdev_counts_kept <- counts_kept[, bdev_columns, drop = FALSE]
  motif_matrix <- motif_matrix_all[peak_keep, , drop = FALSE]
  bias_kept <- peak_bias[peak_keep]
  stopifnot(all(Matrix::colSums(motif_matrix != 0) >= 10L))

  log_event("Computing Bdev-defined chromVAR values for", region_key)
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

  centers <- rowMeans(raw_deviation[, bdev_columns, drop = FALSE])
  scales <- apply(raw_deviation[, bdev_columns, drop = FALSE], 1L, sd)
  stopifnot(all(is.finite(scales) & scales > 0))
  scaled_deviation <- sweep(sweep(raw_deviation, 1L, centers, "-"),
    1L, scales, "/")
  observed <- as.data.table(as.table(scaled_deviation))
  setnames(observed, c("motif_id", "column_id", "bdev_scaled_deviation"))
  observed[, `:=`(motif_id = as.character(motif_id),
    column_id = as.character(column_id))]
  raw_long <- as.data.table(as.table(raw_deviation))
  setnames(raw_long, c("motif_id", "column_id", "raw_chromvar_deviation"))
  raw_long[, `:=`(motif_id = as.character(motif_id),
    column_id = as.character(column_id))]
  observed <- merge(observed, raw_long,
    by = c("motif_id", "column_id"), sort = FALSE)

  plot_grid <- CJ(motif_id = selection$motif_id, group = group_levels,
    bin = seq_len(region$n_bins), unique = TRUE)
  plot_grid <- merge(plot_grid, column_table[, .(group, bin, column_id,
    pt_left, pt_right, pt_mid, n_cells, included)],
    by = c("group", "bin"), all.x = TRUE, sort = FALSE)
  values <- merge(plot_grid, observed,
    by = c("motif_id", "column_id"), all.x = TRUE, sort = FALSE)
  values <- merge(values, selection[, .(motif_id, tf, motif_name,
    display_order)], by = "motif_id", all.x = TRUE, sort = FALSE)
  values[, `:=`(region_id = region_key,
    group = factor(group, levels = group_levels),
    motif_label = factor(sprintf("%s — %s (%s)", tf, motif_name, motif_id),
      levels = sprintf("%s — %s (%s)", selection$tf,
        selection$motif_name, selection$motif_id)))]
  all_values[[region_key]] <- copy(values)

  cytokine_label <- if (region$target == "Flt3l") "Flt3L" else "IL-7"
  knockout_label <- if (region$target == "Flt3l") "FLT3LKO" else "IL7KO"
  p <- ggplot(values, aes(pt_mid, bdev_scaled_deviation,
      color = group, group = group)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.35) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 1.7, na.rm = TRUE) +
    facet_wrap(~ motif_label, ncol = 2, scales = "free_y") +
    scale_color_manual(values = group_colors,
      labels = c(Bdev = "Bdev", WT = "WT", KO = knockout_label),
      drop = FALSE) +
    scale_x_continuous(limits = c(region$pt_left, region$pt_right),
      expand = expansion(mult = c(0.01, 0.01))) +
    theme_minimal(base_size = 10.5) +
    theme(panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "top", legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10)) +
    labs(title = paste0(cytokine_label,
        " enriched-region selected motif trajectories"),
      subtitle = paste0("Same seven JASPAR motifs in both cytokine regions; ",
        "motifs and chromVAR calibration fixed from Bdev only"),
      x = "Bdev pseudotime",
      y = "Motif accessibility (Bdev SD)",
      caption = paste0("Black: Bdev; red: WT; blue: ", knockout_label,
        ". Gaps indicate bins with <", minimum_cells,
        " cells. Facets use independent y-axis ranges."))
  file_stub <- if (region$target == "Flt3l")
    "flt3l_peak1_selected_tf_motif_pt_curves" else
    "il7_peak2_selected_tf_motif_pt_curves"
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".png")), p,
    width = 12, height = 13, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".pdf")), p,
    width = 12, height = 13, device = cairo_pdf, bg = "white")

  region_diagnostics[[region_key]] <- list(
    pt_left = region$pt_left, pt_right = region$pt_right,
    bins = region$n_bins,
    bdev_accessible_peaks = sum(peak_keep),
    displayed_pseudobulks = nrow(analysis_columns),
    omitted_low_cell_pseudobulks = nrow(column_table) - nrow(analysis_columns))
  rm(design, aggregate_counts, counts_kept, bdev_counts_kept,
    deviation_object, raw_deviation, scaled_deviation)
  invisible(gc())
}

write_tsv(rbindlist(all_counts), file.path(stage_dir,
  "pseudotime_bin_cell_counts.tsv"))
fwrite(rbindlist(all_values), file.path(stage_dir,
  "selected_tf_motif_pt_values.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")
diagnostics <- list(status = "complete_preliminary_unregistered",
  bdev_bundle_fingerprint = expected_bdev_fingerprint,
  mapping_bundle_fingerprint = expected_mapping_fingerprint,
  source_motif_report = motif_report,
  motif_selection_uses_query = FALSE,
  same_motifs_in_both_regions = TRUE,
  minimum_cells_per_displayed_pseudobulk = minimum_cells,
  region_diagnostics = region_diagnostics,
  biological_replicates_per_genotype = 1L,
  inference = "descriptive", bundle_mutated = FALSE,
  capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Selected TF motif trajectories", "",
  "The same seven JASPAR motifs are plotted for Flt3L peak 1 and IL-7 peak 2.",
  "For each requested TF/family, the displayed motif has the highest mean Bdev temporal SD across the two regions among matching motifs.",
  "Bdev alone fixes selection, accessible peaks, chromVAR background and expectation, centering, and scaling; WT and KO are display-only afterward.",
  "Black denotes Bdev, red WT, and blue the cytokine-specific knockout. Bins with fewer than 20 cells are gaps.",
  "The curves are descriptive because each genotype has one biological sample."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
wb_close_bundle(bdev)
bdev <- NULL
wb_close_bundle(mapping)
mapping <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed selected TF motif trajectory figures at", output_dir, "\n")
