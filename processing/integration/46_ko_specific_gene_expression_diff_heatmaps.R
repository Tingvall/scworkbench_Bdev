#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) stop(paste(
  "Usage: 46_ko_specific_gene_expression_diff_heatmaps.R",
  "<bdev_bundle> <query_bundle> <wt_enriched_report>",
  "<query_session_workspace> <bdev_session_workspace> <output_dir>",
  "<expected_bdev_fingerprint> <expected_query_fingerprint>",
  "<top_n_genes> <low_cell_threshold> <display_limit>",
  "<bdev_variability_quantile>"))
bdev_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
query_bundle_dir <- normalizePath(args[[2L]], mustWork = TRUE)
region_report <- normalizePath(args[[3L]], mustWork = TRUE)
query_session_workspace <- normalizePath(args[[4L]], mustWork = TRUE)
bdev_session_workspace <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
expected_bdev_fingerprint <- args[[7L]]
expected_query_fingerprint <- args[[8L]]
top_n <- as.integer(args[[9L]])
low_cell_threshold <- as.integer(args[[10L]])
display_limit <- as.numeric(args[[11L]])
bdev_variability_quantile <- as.numeric(args[[12L]])
stopifnot(top_n >= 5L, low_cell_threshold >= 1L,
  is.finite(display_limit), display_limit > 0,
  is.finite(bdev_variability_quantile), bdev_variability_quantile > 0,
  bdev_variability_quantile < 1)
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
  stopifnot(length(ids) == 1L)
  session$runs[[ids[[1L]]]]
}
aggregate_mean <- function(matrix, features, cells, bins, column_ids) {
  stopifnot(length(cells) == length(bins), length(cells) > 0L)
  bin_levels <- as.integer(sub(".*::", "", column_ids))
  kept <- bins %in% bin_levels
  cells <- cells[kept]
  bins <- bins[kept]
  j <- match(bins, bin_levels)
  n <- tabulate(j, nbins = length(column_ids))
  stopifnot(all(n > 0L))
  design <- sparseMatrix(i = seq_along(cells), j = j,
    x = 1 / n[j], dims = c(length(cells), length(column_ids)),
    dimnames = list(cells, column_ids))
  matrix_subset <- matrix[features, cells, drop = FALSE]
  stopifnot(identical(colnames(matrix_subset), rownames(design)))
  as.matrix(matrix_subset %*% design)
}
weighted_mean_rows <- function(x, weights) {
  keep <- is.finite(weights) & weights > 0
  stopifnot(any(keep))
  rowSums(x[, keep, drop = FALSE] *
    matrix(weights[keep], nrow = nrow(x), ncol = sum(keep), byrow = TRUE),
    na.rm = TRUE) / sum(weights[keep])
}

bdev <- wb_open_bundle(bdev_bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
query <- wb_open_bundle(query_bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(query), silent = TRUE), add = TRUE)
bdev_manifest <- wb_manifest(bdev)
query_manifest <- wb_manifest(query)
stopifnot(identical(bdev_manifest$current_state$fingerprint,
    expected_bdev_fingerprint),
  identical(query_manifest$current_state$fingerprint,
    expected_query_fingerprint))

regions <- fread(file.path(region_report, "wt_enriched_pt_regions.tsv"))[
  region_id %in% c("Flt3l_peak1", "Il7_peak2")]
assignments <- fread(cmd = paste("zcat", shQuote(file.path(region_report,
  "cell_region_assignments.tsv.gz"))))[
    region_id %in% c("Flt3l_peak1", "Il7_peak2")]
stopifnot(nrow(regions) == 2L,
  setequal(regions$region_id, c("Flt3l_peak1", "Il7_peak2")))

bdev_cells <- wb_cells(bdev, "set:bdev_reference")
bdev_meta <- as.data.table(wb_meta(bdev, cells = bdev_cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt")))
setnames(bdev_meta, names(bdev_meta)[[1L]], "source_cell_id")
stopifnot(identical(as.character(bdev_meta$source_cell_id), bdev_cells))

query_session <- wb_session_resume(query, query_session_workspace)
query_run <- find_run(query_session, "normalize_downstream_rna",
  "matrix_log_normalize")
stopifnot("matrix" %in% names(query_run$outputs),
  is.null(query_session$durable$promotion_receipt))
query_matrix_id <- "RNA_log1p_cp10k_wt_enriched_peak_de_v1"
query_view <- wb_session_view(query_session)
on.exit(try(wb_close_bundle(query_view), silent = TRUE), add = TRUE)

bdev_session <- wb_session_resume(bdev, bdev_session_workspace)
bdev_run <- find_run(bdev_session, "normalize_bdev_rna",
  "matrix_log_normalize")
stopifnot("matrix" %in% names(bdev_run$outputs),
  is.null(bdev_session$durable$promotion_receipt))
bdev_matrix_id <- "RNA_log1p_cp10k_bdev_variable_gene_expression_v1"
bdev_view <- wb_session_view(bdev_session)
on.exit(try(wb_close_bundle(bdev_view), silent = TRUE), add = TRUE)

bdev_norm <- wb_matrix(bdev_view, bdev_matrix_id)
query_norm <- wb_matrix(query_view, query_matrix_id)
common_features <- intersect(rownames(bdev_norm), rownames(query_norm))
stopifnot(length(common_features) > 25000L, !anyDuplicated(common_features))

all_scores <- list()
all_selected <- list()
all_values <- list()
all_bdev_values <- list()
all_counts <- list()
region_diagnostics <- list()

for (region_index in seq_len(nrow(regions))) {
  region <- regions[region_index]
  region_key <- region$region_id
  log_event("Scoring KO-specific expression divergence for", region_key)
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
    bdev_region[, .(source_cell_id, group, bin)],
    query_region[, .(source_cell_id, group, bin)]))

  column_table <- CJ(group = c("Bdev", "WT", "KO"),
    bin = seq_len(region$n_bins), unique = TRUE)
  column_table[, `:=`(
    column_id = paste(group, sprintf("%02d", bin), sep = "::"),
    pt_mid = bin_midpoints[bin])]
  count_table <- cell_table[, .(n_cells = .N), by = .(group, bin)]
  column_table <- merge(column_table, count_table,
    by = c("group", "bin"), all.x = TRUE, sort = FALSE)
  column_table[is.na(n_cells), n_cells := 0L]
  column_table[, `:=`(calculable = n_cells > 0L,
    low_cell_bin = n_cells < low_cell_threshold,
    region_id = region_key)]
  column_table[, group_order := match(group, c("Bdev", "WT", "KO"))]
  setorder(column_table, group_order, bin)
  column_table[, group_order := NULL]
  stopifnot(all(column_table[group == "Bdev", calculable]))
  all_counts[[region_key]] <- copy(column_table)

  group_means <- list()
  for (group_name in c("Bdev", "WT", "KO")) {
    group_cells <- cell_table[group == group_name]
    group_columns <- column_table[group == group_name & calculable == TRUE]
    source_matrix <- if (group_name == "Bdev") bdev_norm else query_norm
    group_means[[group_name]] <- aggregate_mean(source_matrix,
      common_features, group_cells$source_cell_id, group_cells$bin,
      group_columns$column_id)
  }

  full_matrix <- function(group_name) {
    result <- matrix(NA_real_, nrow = length(common_features),
      ncol = region$n_bins,
      dimnames = list(common_features, as.character(seq_len(region$n_bins))))
    columns <- column_table[group == group_name & calculable == TRUE]
    result[, as.character(columns$bin)] <- group_means[[group_name]]
    result
  }
  bdev_mean <- full_matrix("Bdev")
  wt_mean <- full_matrix("WT")
  ko_mean <- full_matrix("KO")
  wt_difference <- wt_mean - bdev_mean
  ko_difference <- ko_mean - bdev_mean

  bdev_n <- column_table[group == "Bdev", n_cells]
  wt_n <- column_table[group == "WT", n_cells]
  ko_n <- column_table[group == "KO", n_cells]
  reliability <- function(query_n) sqrt(
    pmin(bdev_n, low_cell_threshold) / low_cell_threshold *
      pmin(query_n, low_cell_threshold) / low_cell_threshold)
  wt_weights <- reliability(wt_n)
  ko_weights <- reliability(ko_n)
  wt_mad <- weighted_mean_rows(abs(wt_difference), wt_weights)
  ko_mad <- weighted_mean_rows(abs(ko_difference), ko_weights)
  wt_signed <- weighted_mean_rows(wt_difference, wt_weights)
  ko_signed <- weighted_mean_rows(ko_difference, ko_weights)
  bdev_temporal_sd <- apply(bdev_mean, 1L, sd)
  bdev_temporal_range <- apply(bdev_mean, 1L, function(x) diff(range(x)))
  bdev_temporal_iqr <- apply(bdev_mean, 1L, IQR)
  bdev_nonzero_bins <- rowSums(bdev_mean > 0)
  bdev_maximum_mean_expression <- apply(bdev_mean, 1L, max,
    na.rm = TRUE)
  bdev_variability_base <- is.finite(bdev_temporal_sd) &
    bdev_maximum_mean_expression >= 0.1 & bdev_nonzero_bins >= 3L
  bdev_variability_threshold <- as.numeric(quantile(
    bdev_temporal_sd[bdev_variability_base],
    probs = bdev_variability_quantile, names = FALSE, na.rm = TRUE))
  passes_bdev_variability <- bdev_variability_base &
    bdev_temporal_sd >= bdev_variability_threshold
  maximum_mean_expression <- apply(cbind(bdev_mean, wt_mean, ko_mean),
    1L, max, na.rm = TRUE)
  specificity_score <- ko_mad - wt_mad
  specificity_ratio <- (ko_mad + 0.05) / (wt_mad + 0.05)
  eligible <- is.finite(specificity_score) & passes_bdev_variability &
    ko_mad > wt_mad
  primary <- eligible & ko_mad >= 0.2 & specificity_ratio >= 1.5
  scores <- data.table(region_id = region_key,
    feature_id = common_features,
    wt_bdev_weighted_mean_absolute_difference = wt_mad,
    ko_bdev_weighted_mean_absolute_difference = ko_mad,
    wt_bdev_weighted_signed_difference = wt_signed,
    ko_bdev_weighted_signed_difference = ko_signed,
    ko_specificity_score = specificity_score,
    ko_to_wt_difference_ratio = specificity_ratio,
    maximum_bin_mean_log1p_cp10k = maximum_mean_expression,
    bdev_maximum_bin_mean_log1p_cp10k = bdev_maximum_mean_expression,
    bdev_nonzero_mean_bins = as.integer(bdev_nonzero_bins),
    bdev_temporal_sd = bdev_temporal_sd,
    bdev_temporal_range = bdev_temporal_range,
    bdev_temporal_iqr = bdev_temporal_iqr,
    bdev_variability_threshold = bdev_variability_threshold,
    passes_bdev_variability = passes_bdev_variability,
    eligible = eligible, passes_primary_thresholds = primary)
  setorder(scores, -passes_primary_thresholds, -ko_specificity_score,
    -ko_bdev_weighted_mean_absolute_difference,
    -ko_to_wt_difference_ratio, feature_id)
  scores[, ko_specificity_rank := NA_integer_]
  scores[eligible == TRUE, ko_specificity_rank := seq_len(.N)]
  primary_count <- sum(primary)
  selection_pool <- if (primary_count >= top_n)
    scores[passes_primary_thresholds == TRUE] else scores[eligible == TRUE]
  selected <- head(selection_pool, top_n)
  stopifnot(nrow(selected) == top_n)
  scores[, selected := feature_id %in% selected$feature_id]
  selected[, selection_rank := seq_len(.N)]
  all_scores[[region_key]] <- scores
  all_selected[[region_key]] <- selected

  selected_features <- selected$feature_id
  bdev_counts_by_bin <- bdev_n
  make_long <- function(group_name, difference_matrix,
      query_counts_by_bin) {
    m <- difference_matrix[selected_features, , drop = FALSE]
    out <- as.data.table(as.table(m))
    setnames(out, c("feature_id", "bin", "query_minus_bdev_log1p_cp10k"))
    out[, `:=`(feature_id = as.character(feature_id),
      bin = as.integer(as.character(bin)))]
    out[, `:=`(query_group = group_name,
      pt_mid = bin_midpoints[bin],
      bdev_n = bdev_counts_by_bin[bin],
      query_n = query_counts_by_bin[bin],
      low_cell_bin = bdev_counts_by_bin[bin] < low_cell_threshold |
        query_counts_by_bin[bin] < low_cell_threshold,
      calculable = is.finite(query_minus_bdev_log1p_cp10k),
      region_id = region_key)]
    out
  }
  heat <- rbindlist(list(make_long("WT", wt_difference, wt_n),
    make_long("KO", ko_difference, ko_n)))
  heat <- merge(heat, selected[, .(feature_id, selection_rank,
    wt_bdev_weighted_mean_absolute_difference,
    ko_bdev_weighted_mean_absolute_difference, ko_specificity_score,
    ko_to_wt_difference_ratio)], by = "feature_id", all.x = TRUE,
    sort = FALSE)

  if (region$target == "Flt3l") {
    panel_labels <- c(WT = "FLT3L WT − Bdev",
      KO = "FLT3LKO − Bdev")
    title <- "Flt3L region: KO-specific expression divergence from Bdev"
    file_stub <- "flt3l_peak1_ko_specific_gene_expression_diff_heatmap"
  } else {
    panel_labels <- c(WT = "IL-7 WT − Bdev", KO = "IL7KO − Bdev")
    title <- "IL-7 region: KO-specific expression divergence from Bdev"
    file_stub <- "il7_peak2_ko_specific_gene_expression_diff_heatmap"
  }
  heat[, panel_label := factor(panel_labels[query_group],
    levels = unname(panel_labels))]
  heat[, feature_id := factor(feature_id,
    levels = rev(selected$feature_id))]
  all_values[[region_key]] <- copy(heat)

  bdev_heat <- as.data.table(as.table(
    bdev_mean[selected_features, , drop = FALSE]))
  setnames(bdev_heat, c("feature_id", "bin", "mean_log1p_cp10k"))
  bdev_heat[, `:=`(feature_id = as.character(feature_id),
    bin = as.integer(as.character(bin)))]
  bdev_heat[, `:=`(pt_mid = bin_midpoints[bin],
    n_cells = bdev_counts_by_bin[bin],
    low_cell_bin = bdev_counts_by_bin[bin] < low_cell_threshold,
    region_id = region_key)]
  bdev_heat[, feature_id := factor(feature_id,
    levels = rev(selected$feature_id))]
  all_bdev_values[[region_key]] <- copy(bdev_heat)

  caption_text <- paste0("Candidates first pass a Bdev-only temporal-SD filter ",
    "(top ", round(100 * (1 - bdev_variability_quantile)),
    "%; SD ≥ ", signif(bdev_variability_threshold, 3),
    "), then rank by reliability-weighted mean |KO−Bdev| minus mean ",
    "|WT−Bdev|. Black borders/crosses ",
    "mark either constituent bin with <", low_cell_threshold,
    " cells; gray is not calculable. Values are raw log1p-CP10K differences; ",
    "display clipped at ±", display_limit, ".")
  bdev_display_limit <- as.numeric(quantile(bdev_heat$mean_log1p_cp10k,
    0.99, na.rm = TRUE))
  p_bdev <- ggplot(bdev_heat, aes(pt_mid, feature_id,
      fill = mean_log1p_cp10k)) +
    geom_tile() +
    geom_tile(data = bdev_heat[low_cell_bin == TRUE], fill = NA,
      color = "black", linewidth = 0.28) +
    geom_point(data = bdev_heat[low_cell_bin == TRUE], shape = 4,
      size = 0.65, stroke = 0.3, color = "black") +
    scale_fill_gradient(low = "white", high = "#54278F",
      limits = c(0, bdev_display_limit), oob = scales::squish,
      name = "Bdev mean\nlog1p-CP10K") +
    scale_x_continuous(limits = c(region$pt_left, region$pt_right),
      expand = expansion(mult = 0),
      breaks = pretty(c(region$pt_left, region$pt_right), n = 4)) +
    theme_minimal(base_size = 9.5) +
    theme(panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "right") +
    labs(title = "Bdev expression", x = "Bdev pseudotime", y = NULL)

  p_difference <- ggplot(heat, aes(pt_mid, feature_id,
      fill = query_minus_bdev_log1p_cp10k)) +
    geom_tile(na.rm = FALSE) +
    geom_tile(data = heat[low_cell_bin == TRUE], fill = NA,
      color = "black", linewidth = 0.28) +
    geom_point(data = heat[low_cell_bin == TRUE], shape = 4,
      size = 0.65, stroke = 0.3, color = "black") +
    facet_grid(. ~ panel_label) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-display_limit, display_limit),
      oob = scales::squish, na.value = "grey65",
      name = "Query − Bdev\nlog1p-CP10K") +
    scale_x_continuous(limits = c(region$pt_left, region$pt_right),
      expand = expansion(mult = 0),
      breaks = pretty(c(region$pt_left, region$pt_right), n = 5)) +
    theme_minimal(base_size = 9.5) +
    theme(panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      legend.position = "right") +
    labs(x = "Bdev pseudotime", y = NULL)
  p <- p_bdev + p_difference +
    plot_layout(widths = c(1, 2.15)) +
    plot_annotation(title = title,
      subtitle = paste0("Top ", top_n,
        " Bdev-variable genes with large KO−Bdev but small WT−Bdev deviation"),
      caption = paste(strwrap(caption_text, width = 165), collapse = "\n"),
      theme = theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 9.5),
        plot.caption = element_text(size = 8)))
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".png")), p,
    width = 18, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".pdf")), p,
    width = 18, height = 10, device = cairo_pdf, bg = "white")

  region_diagnostics[[region_key]] <- list(
    common_gene_universe = length(common_features),
    bdev_variability_base_genes = sum(bdev_variability_base),
    bdev_variability_quantile = bdev_variability_quantile,
    bdev_variability_threshold = bdev_variability_threshold,
    bdev_variable_genes = sum(passes_bdev_variability),
    eligible_genes = sum(eligible), primary_threshold_genes = primary_count,
    selected_genes = nrow(selected),
    selection_pool = if (primary_count >= top_n)
      "primary_threshold_passers" else "eligible_ko_greater_than_wt",
    wt_calculable_bins = sum(wt_n > 0L),
    ko_calculable_bins = sum(ko_n > 0L))
  rm(bdev_mean, wt_mean, ko_mean, wt_difference, ko_difference)
  invisible(gc())
}

write_tsv(rbindlist(all_scores), file.path(stage_dir,
  "all_gene_ko_specificity_scores.tsv"))
write_tsv(rbindlist(all_selected), file.path(stage_dir,
  "selected_ko_specific_genes.tsv"))
fwrite(rbindlist(all_values), file.path(stage_dir,
  "ko_specific_gene_expression_diff_values.tsv.gz"), sep = "\t",
  na = "NA", compress = "gzip")
fwrite(rbindlist(all_bdev_values), file.path(stage_dir,
  "selected_gene_bdev_expression_values.tsv.gz"), sep = "\t",
  na = "NA", compress = "gzip")
write_tsv(rbindlist(all_counts), file.path(stage_dir,
  "pseudotime_bin_cell_counts.tsv"))
diagnostics <- list(status = "complete_preliminary_unregistered",
  bdev_bundle_fingerprint = expected_bdev_fingerprint,
  query_bundle_fingerprint = expected_query_fingerprint,
  normalization = "native scWorkbench matrix_log_normalize log1p-CP10K",
  selection = "Bdev-variable genes ranked by reliability-weighted mean absolute KO-Bdev difference minus WT-Bdev difference",
  bdev_variability_filter = list(
    metric = "temporal SD across Bdev pseudotime-bin means",
    base_eligibility = "Bdev maximum bin mean >= 0.1 and at least 3 nonzero Bdev bins",
    retained_quantile = bdev_variability_quantile,
    retained_tail = paste0("top ", 100 * (1 - bdev_variability_quantile), "%")),
  reliability_weight = "geometric mean of Bdev and query cell-count fractions capped at 20",
  primary_thresholds = list(ko_mad_minimum = 0.2,
    ko_to_wt_ratio_minimum = 1.5, maximum_mean_expression_minimum = 0.1),
  all_nonempty_bins_used = TRUE,
  low_cell_bins_downweighted_not_excluded = TRUE,
  display_unit = "raw mean log1p-CP10K difference",
  display_clip = display_limit, exported_values_clipped = FALSE,
  bdev_expression_panel = "raw mean log1p-CP10K with region-specific p99 display cap",
  region_diagnostics = region_diagnostics,
  biological_replicates_per_genotype = 1L,
  inference = "descriptive", bundle_mutated = FALSE,
  capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# KO-specific expression divergence from Bdev", "",
  "Genes are screened independently in Flt3L peak 1 and IL-7 peak 2.",
  paste0("The candidate universe is restricted using Bdev alone to the top ",
    100 * (1 - bdev_variability_quantile),
    "% by temporal SD across Bdev pseudotime-bin means, after requiring a Bdev maximum bin mean of at least 0.1 and expression in at least three Bdev bins."),
  "Selection favors large reliability-weighted mean absolute KO-minus-Bdev differences and small WT-minus-Bdev differences across matched pseudotime bins.",
  "Every nonempty bin contributes; bins below 20 cells are down-weighted rather than excluded and are marked in the heatmaps.",
  "Selection uses raw log1p-CP10K differences to prevent genes with nearly constant Bdev expression from ranking through division by a tiny Bdev SD.",
  "Each figure includes the corresponding Bdev mean-expression heatmap beside WT-minus-Bdev and KO-minus-Bdev.",
  "There is one biological sample per genotype, so results are descriptive."),
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
cat("Completed KO-specific expression divergence analysis at", output_dir,
  "\n")
