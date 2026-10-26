#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) stop(paste(
  "Usage: 44_bdev_variable_gene_expression_pt_heatmaps.R",
  "<bdev_bundle> <query_bundle> <wt_enriched_report>",
  "<query_session_workspace> <bdev_session_root> <output_dir>",
  "<expected_bdev_fingerprint> <expected_query_fingerprint>",
  "<top_n_genes> <minimum_cells> <seed>"))
bdev_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
query_bundle_dir <- normalizePath(args[[2L]], mustWork = TRUE)
region_report <- normalizePath(args[[3L]], mustWork = TRUE)
query_session_workspace <- normalizePath(args[[4L]], mustWork = TRUE)
bdev_session_root <- normalizePath(args[[5L]], mustWork = FALSE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
expected_bdev_fingerprint <- args[[7L]]
expected_query_fingerprint <- args[[8L]]
top_n <- as.integer(args[[9L]])
minimum_cells <- as.integer(args[[10L]])
seed <- as.integer(args[[11L]])
stopifnot(top_n >= 5L, minimum_cells >= 5L, is.finite(seed))
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE)
dir.create(bdev_session_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(BPCells)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
write_tsv <- function(x, path) data.table::fwrite(x, path, sep = "\t",
  na = "NA")
log_event <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(..., collapse = " "), "\n")
find_run <- function(session, label, op) {
  marker <- paste0("[", label, "] ")
  ids <- names(session$runs)[vapply(session$runs, function(x) {
    note <- if (is.null(x$notes)) "" else as.character(x$notes)
    identical(x$op_name, op) && startsWith(note, marker)
  }, logical(1L))]
  stopifnot(length(ids) <= 1L)
  if (length(ids)) session$runs[[ids[[1L]]]] else NULL
}
make_design <- function(cell_table, column_table, weighted) {
  cell_table <- cell_table[column_id %in% column_table$column_id]
  column_index <- match(cell_table$column_id, column_table$column_id)
  weights <- if (weighted) 1 / column_table$n_cells[column_index] else
    rep(1, nrow(cell_table))
  sparseMatrix(i = seq_len(nrow(cell_table)), j = column_index, x = weights,
    dims = c(nrow(cell_table), nrow(column_table)),
    dimnames = list(cell_table$source_cell_id, column_table$column_id))
}
aggregate_by_columns <- function(matrix, features, cell_table, column_table,
    weighted = TRUE) {
  kept_cells <- cell_table[column_id %in% column_table$column_id]
  design <- make_design(kept_cells, column_table, weighted)
  matrix_subset <- matrix[features, kept_cells$source_cell_id, drop = FALSE]
  stopifnot(identical(colnames(matrix_subset), rownames(design)))
  result <- as(matrix_subset %*% design, "dgCMatrix")
  stopifnot(identical(rownames(result), features),
    identical(colnames(result), column_table$column_id))
  result
}

bdev <- wb_open_bundle(bdev_bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
query <- wb_open_bundle(query_bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(query), silent = TRUE), add = TRUE)
bdev_manifest <- wb_manifest(bdev)
query_manifest <- wb_manifest(query)
stopifnot(identical(bdev_manifest$bundle_id, "bdev_reference"),
  identical(bdev_manifest$current_state$fingerprint,
    expected_bdev_fingerprint),
  identical(query_manifest$bundle_id, "scmulti_flt3ko_il7ko"),
  identical(query_manifest$current_state$fingerprint,
    expected_query_fingerprint))

regions <- fread(file.path(region_report, "wt_enriched_pt_regions.tsv"))[
  region_id %in% c("Flt3l_peak1", "Il7_peak2")]
assignments <- fread(cmd = paste("zcat", shQuote(file.path(region_report,
  "cell_region_assignments.tsv.gz"))))[
    region_id %in% c("Flt3l_peak1", "Il7_peak2")]
stopifnot(nrow(regions) == 2L,
  setequal(regions$region_id, c("Flt3l_peak1", "Il7_peak2")),
  !anyDuplicated(assignments[, .(region_id, source_cell_id)]))

bdev_cells <- wb_cells(bdev, "set:bdev_reference")
bdev_meta <- as.data.table(wb_meta(bdev, cells = bdev_cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt")))
setnames(bdev_meta, names(bdev_meta)[[1L]], "source_cell_id")
stopifnot(identical(as.character(bdev_meta$source_cell_id), bdev_cells))
query_cells <- wb_cells(query, "set:qc_use_downstream")
stopifnot(all(assignments$source_cell_id %in% query_cells))

# Reuse the native scWorkbench query normalization already staged for the
# accepted WT-enriched-region RNA analysis.
query_session <- wb_session_resume(query, query_session_workspace)
if (!is.null(query_session$durable$promotion_receipt))
  stop("Query report-only session was unexpectedly promoted")
query_run <- find_run(query_session, "normalize_downstream_rna",
  "matrix_log_normalize")
stopifnot(!is.null(query_run), "matrix" %in% names(query_run$outputs))
query_matrix_id <- "RNA_log1p_cp10k_wt_enriched_peak_de_v1"
query_view <- wb_session_view(query_session)
on.exit(try(wb_close_bundle(query_view), silent = TRUE), add = TRUE)

# Stage the identical native operation for Bdev without modifying the bundle.
bdev_session_id <- "bdev_variable_gene_expression_norm_v1"
bdev_workspace <- file.path(bdev_session_root, bdev_session_id)
bdev_session <- if (dir.exists(bdev_workspace))
  wb_session_resume(bdev, bdev_workspace) else
  wb_session_create(bdev, workspace_root = bdev_session_root,
    session_id = bdev_session_id)
if (!is.null(bdev_session$durable$promotion_receipt))
  stop("Bdev report-only session was unexpectedly promoted")
bdev_matrix_id <- "RNA_log1p_cp10k_bdev_variable_gene_expression_v1"
bdev_run <- find_run(bdev_session, "normalize_bdev_rna",
  "matrix_log_normalize")
if (is.null(bdev_run)) {
  log_event("Running native scWorkbench log1p-CP10K normalization for Bdev")
  result <- wb_session_run_op(bdev_session, bdev, "matrix_log_normalize",
    inputs = list(matrix = "RNA_counts", cell_set = "all",
      min_cells = 10L, output_matrix_id = bdev_matrix_id,
      scale_factor = 10000, overwrite = FALSE), notes = paste0(
        "[normalize_bdev_rna] Full-bundle native log1p-CP10K matrix; ",
        "Bdev-only variable-gene selection for report-only heatmaps; ",
        "no promotion."))
  bdev_session <- result$session
  bdev_run <- result$run
} else log_event("Reusing staged native Bdev log1p-CP10K matrix")
stopifnot("matrix" %in% names(bdev_run$outputs))
bdev_view <- wb_session_view(bdev_session)
on.exit(try(wb_close_bundle(bdev_view), silent = TRUE), add = TRUE)

bdev_norm <- wb_matrix(bdev_view, bdev_matrix_id)
query_norm <- wb_matrix(query_view, query_matrix_id)
bdev_counts <- wb_matrix(bdev, "RNA_counts")
common_features <- intersect(rownames(bdev_norm), rownames(query_norm))
stopifnot(length(common_features) > 25000L, !anyDuplicated(common_features),
  all(common_features %in% rownames(bdev_counts)))

all_variability <- list()
all_selected <- list()
all_values <- list()
all_counts <- list()
region_diagnostics <- list()
group_levels <- c("Bdev", "WT", "KO")

for (region_index in seq_len(nrow(regions))) {
  region <- regions[region_index]
  region_key <- region$region_id
  log_event("Selecting Bdev-variable genes for", region_key)
  edges <- seq(region$pt_left, region$pt_right,
    length.out = region$n_bins + 1L)
  bin_midpoints <- (head(edges, -1L) + tail(edges, -1L)) / 2
  assign_bin <- function(pt) as.integer(cut(pt, breaks = edges,
    include.lowest = TRUE, right = FALSE, labels = FALSE))

  bdev_region <- bdev_meta[is.finite(pt) & pt >= region$pt_left &
    pt < region$pt_right]
  bdev_region[, `:=`(group = "Bdev", bin = assign_bin(pt))]
  query_region <- assignments[region_id == region_key]
  query_region[, `:=`(pt = mapped_pt,
    group = ifelse(genotype == "WT", "WT", "KO"),
    bin = assign_bin(mapped_pt))]
  cell_table <- rbindlist(list(
    bdev_region[, .(source_cell_id, group, bin, pt)],
    query_region[, .(source_cell_id, group, bin, pt)]))
  stopifnot(!anyNA(cell_table$bin),
    !anyDuplicated(cell_table[, .(group, source_cell_id)]))

  column_table <- CJ(group = group_levels, bin = seq_len(region$n_bins),
    unique = TRUE)
  column_table[, `:=`(
    column_id = paste(group, sprintf("%02d", bin), sep = "::"),
    pt_left = edges[bin], pt_right = edges[bin + 1L],
    pt_mid = bin_midpoints[bin])]
  cell_table <- merge(cell_table, column_table[, .(group, bin, column_id)],
    by = c("group", "bin"), all.x = TRUE, sort = FALSE)
  count_table <- cell_table[, .(n_cells = .N), by = .(group, bin)]
  column_table <- merge(column_table, count_table,
    by = c("group", "bin"), all.x = TRUE, sort = FALSE)
  column_table[is.na(n_cells), n_cells := 0L]
  column_table[, `:=`(included = n_cells > 0L,
    low_cell_bin = n_cells < minimum_cells, region_id = region_key)]
  column_table[, group_order := match(group, group_levels)]
  setorder(column_table, group_order, bin)
  column_table[, group_order := NULL]
  stopifnot(all(column_table[group == "Bdev", included]))
  all_counts[[region_key]] <- copy(column_table)

  bdev_columns <- column_table[group == "Bdev" & included == TRUE]
  bdev_cell_table <- cell_table[group == "Bdev"]
  bdev_mean <- aggregate_by_columns(bdev_norm, common_features,
    bdev_cell_table, bdev_columns, weighted = TRUE)
  bdev_count_sum <- aggregate_by_columns(bdev_counts, common_features,
    bdev_cell_table, bdev_columns, weighted = FALSE)
  bdev_total_counts <- Matrix::rowSums(bdev_count_sum)
  temporal_sd <- apply(bdev_mean, 1L, sd)
  temporal_range <- apply(bdev_mean, 1L, function(x) diff(range(x)))
  temporal_iqr <- apply(bdev_mean, 1L, IQR)
  nonzero_bins <- Matrix::rowSums(bdev_mean > 0)
  eligible <- bdev_total_counts >= 10 & nonzero_bins >= 3 &
    is.finite(temporal_sd) & temporal_sd > 0
  variability <- data.table(region_id = region_key,
    feature_id = common_features, bdev_temporal_sd = temporal_sd,
    bdev_temporal_range = temporal_range,
    bdev_temporal_iqr = temporal_iqr,
    bdev_total_counts_in_region = as.numeric(bdev_total_counts),
    bdev_nonzero_mean_bins = as.integer(nonzero_bins), eligible = eligible)
  variability[, bdev_variability_rank := NA_integer_]
  eligible_order <- order(-variability$bdev_temporal_sd[eligible],
    -variability$bdev_temporal_range[eligible],
    variability$feature_id[eligible])
  eligible_rows <- which(eligible)[eligible_order]
  variability$bdev_variability_rank[eligible_rows] <-
    seq_along(eligible_rows)
  variability[, selected := eligible & bdev_variability_rank <= top_n]
  setorder(variability, -selected, bdev_variability_rank, feature_id,
    na.last = TRUE)
  selected <- variability[selected == TRUE]
  stopifnot(nrow(selected) == top_n)
  all_variability[[region_key]] <- variability
  all_selected[[region_key]] <- selected

  selected_features <- selected$feature_id
  analysis_columns <- column_table[included == TRUE]
  query_columns <- analysis_columns[group != "Bdev"]
  query_cell_table <- cell_table[group != "Bdev"]
  selected_bdev_mean <- bdev_mean[selected_features, , drop = FALSE]
  query_mean <- aggregate_by_columns(query_norm, selected_features,
    query_cell_table, query_columns, weighted = TRUE)
  mean_expression <- as.matrix(cbind(selected_bdev_mean, query_mean))
  mean_expression <- mean_expression[, analysis_columns$column_id,
    drop = FALSE]
  centers <- rowMeans(selected_bdev_mean)
  scales <- apply(selected_bdev_mean, 1L, sd)
  stopifnot(all(is.finite(scales) & scales > 0))
  scaled_expression <- sweep(sweep(mean_expression, 1L, centers, "-"),
    1L, scales, "/")
  row_order <- rownames(selected_bdev_mean)[
    hclust(dist(sweep(sweep(selected_bdev_mean, 1L, centers, "-"),
      1L, scales, "/")))$order]

  observed <- as.data.table(as.table(scaled_expression))
  setnames(observed, c("feature_id", "column_id",
    "bdev_scaled_mean_expression"))
  observed[, `:=`(feature_id = as.character(feature_id),
    column_id = as.character(column_id))]
  raw_long <- as.data.table(as.table(mean_expression))
  setnames(raw_long, c("feature_id", "column_id", "mean_log1p_cp10k"))
  raw_long[, `:=`(feature_id = as.character(feature_id),
    column_id = as.character(column_id))]
  observed <- merge(observed, raw_long, by = c("feature_id", "column_id"),
    sort = FALSE)
  plot_grid <- CJ(feature_id = selected_features, group = group_levels,
    bin = seq_len(region$n_bins), unique = TRUE)
  plot_grid <- merge(plot_grid, column_table[, .(group, bin, column_id,
    pt_mid, n_cells, included, low_cell_bin)], by = c("group", "bin"),
    all.x = TRUE, sort = FALSE)
  heat <- merge(plot_grid, observed, by = c("feature_id", "column_id"),
    all.x = TRUE, sort = FALSE)

  wt <- heat[group == "WT", .(feature_id, bin, pt_mid,
    wt_n = n_cells, wt = bdev_scaled_mean_expression)]
  ko <- heat[group == "KO", .(feature_id, bin, ko_n = n_cells,
    ko = bdev_scaled_mean_expression)]
  difference <- merge(wt, ko, by = c("feature_id", "bin"),
    all = TRUE, sort = FALSE)
  difference[, `:=`(group = "WT - KO", n_cells = NA_integer_,
    column_id = NA_character_, included = is.finite(wt) & is.finite(ko),
    low_cell_bin = wt_n < minimum_cells | ko_n < minimum_cells,
    mean_log1p_cp10k = NA_real_,
    bdev_scaled_mean_expression = wt - ko)]
  difference[, c("wt_n", "ko_n", "wt", "ko") := NULL]
  heat <- rbindlist(list(heat, difference), fill = TRUE)
  heat[, region_id := region_key]

  group_labels <- if (region$target == "Flt3l")
    c(Bdev = "Bdev reference", WT = "FLT3L WT", KO = "FLT3LKO",
      `WT - KO` = "WT - FLT3LKO") else
    c(Bdev = "Bdev reference", WT = "IL-7 WT", KO = "IL7KO",
      `WT - KO` = "WT - IL7KO")
  heat[, group_label := factor(group_labels[group],
    levels = unname(group_labels))]
  heat[, feature_id := factor(feature_id, levels = rev(row_order))]
  all_values[[region_key]] <- copy(heat)

  display_limit <- 3
  p <- ggplot(heat, aes(pt_mid, feature_id,
      fill = bdev_scaled_mean_expression)) +
    geom_tile(na.rm = FALSE) +
    geom_tile(data = heat[low_cell_bin == TRUE], fill = NA,
      color = "black", linewidth = 0.28) +
    geom_point(data = heat[low_cell_bin == TRUE], shape = 4,
      size = 0.65, stroke = 0.3, color = "black") +
    facet_grid(. ~ group_label) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-display_limit, display_limit),
      oob = scales::squish, name = "Mean expression\n(Bdev SD)") +
    scale_x_continuous(expand = expansion(mult = 0),
      breaks = pretty(c(region$pt_left, region$pt_right), n = 4)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6.5),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 9.5), legend.position = "right") +
    labs(title = paste0(region$cytokine,
        " enriched-region variable-gene expression across pseudotime"),
      subtitle = paste0("Top ", top_n,
        " genes selected, ordered, centered, and scaled using Bdev only; ",
        "WT and KO do not affect selection"),
      x = "Bdev pseudotime", y = NULL,
      caption = paste0("Native scWorkbench log1p-CP10K bin means; all ",
        "nonempty bins are calculated. Black borders/crosses mark <",
        minimum_cells, " cells; gray indicates zero cells. Display ",
        "clipped at ±", display_limit,
        " Bdev SD. Exported values are not clipped; WT-KO is descriptive."))
  file_stub <- if (region$target == "Flt3l")
    "flt3l_peak1_bdev_variable_gene_expression_heatmap" else
    "il7_peak2_bdev_variable_gene_expression_heatmap"
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".png")), p,
    width = 18, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".pdf")), p,
    width = 18, height = 10, device = cairo_pdf, bg = "white")

  region_diagnostics[[region_key]] <- list(
    pt_left = region$pt_left, pt_right = region$pt_right,
    bins = region$n_bins, common_gene_universe = length(common_features),
    eligible_genes = sum(eligible), selected_genes = top_n,
    bdev_cells = nrow(bdev_region), query_cells = nrow(query_region),
    low_cell_pseudobulks = sum(column_table$low_cell_bin),
    zero_cell_pseudobulks = sum(column_table$n_cells == 0L))
  rm(bdev_mean, bdev_count_sum, selected_bdev_mean, query_mean,
    mean_expression, scaled_expression)
  invisible(gc())
}

write_tsv(rbindlist(all_variability), file.path(stage_dir,
  "all_gene_bdev_temporal_variability.tsv"))
write_tsv(rbindlist(all_selected), file.path(stage_dir,
  "selected_bdev_variable_genes.tsv"))
fwrite(rbindlist(all_values), file.path(stage_dir,
  "selected_gene_expression_heatmap_values.tsv.gz"), sep = "\t",
  na = "NA", compress = "gzip")
write_tsv(rbindlist(all_counts), file.path(stage_dir,
  "pseudotime_bin_cell_counts.tsv"))
diagnostics <- list(status = "complete_preliminary_unregistered",
  bdev_bundle_fingerprint = expected_bdev_fingerprint,
  query_bundle_fingerprint = expected_query_fingerprint,
  normalization = "native scWorkbench matrix_log_normalize log1p-CP10K",
  bdev_matrix_id = bdev_matrix_id, query_matrix_id = query_matrix_id,
  top_n_genes = top_n,
  selection = "top temporal SD of mean expression across Bdev pseudotime bins",
  eligibility = "common genes with at least 10 Bdev counts in-region and positive mean expression in at least 3 Bdev bins",
  selection_uses_query = FALSE,
  background_centering_scaling_and_row_order = "Bdev-only",
  low_cell_highlight_threshold = minimum_cells,
  nonempty_bins_are_all_calculated = TRUE,
  display_clip_bdev_sd = 3L, exported_values_clipped = FALSE,
  region_diagnostics = region_diagnostics,
  biological_replicates_per_genotype = 1L,
  inference = "descriptive", bundle_mutated = FALSE,
  capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Bdev-variable gene-expression pseudotime heatmaps", "",
  "Separate heatmaps cover Flt3L peak 1 and IL-7 peak 2.",
  "For each region, the top 30 genes are selected independently by temporal SD of Bdev pseudotime-bin mean expression.",
  "Eligibility requires a common RNA feature, at least 10 Bdev counts in the region, and positive mean expression in at least three retained Bdev bins.",
  "Both datasets use native scWorkbench full-matrix log1p-CP10K normalization.",
  "Bdev alone determines selection, row order, centering, and scaling; WT and KO are evaluated afterward.",
  "Every nonempty bin is calculated. Bins with fewer than 20 cells are marked with a black border and cross; zero-cell bins are gray.",
  "The fourth panel is descriptive WT-minus-KO. There is one biological sample per genotype."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
stopifnot(identical(wb_manifest(bdev)$current_state$fingerprint,
  expected_bdev_fingerprint),
  identical(wb_manifest(query)$current_state$fingerprint,
    expected_query_fingerprint))
wb_close_bundle(bdev_view)
bdev_view <- NULL
wb_close_bundle(query_view)
query_view <- NULL
wb_close_bundle(bdev)
bdev <- NULL
wb_close_bundle(query)
query <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed Bdev-variable gene-expression heatmaps at", output_dir, "\n")
