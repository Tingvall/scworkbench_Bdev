#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 45_query_vs_bdev_gene_expression_diff_heatmaps.R",
  "<source_gene_heatmap_report> <output_dir> <low_cell_threshold> <display_limit>"))
source_report <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
low_cell_threshold <- as.integer(args[[3L]])
display_limit <- as.numeric(args[[4L]])
stopifnot(low_cell_threshold >= 1L, is.finite(display_limit),
  display_limit > 0)
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

source_values <- fread(cmd = paste("zcat", shQuote(file.path(source_report,
  "selected_gene_expression_heatmap_values.tsv.gz"))))
source_selection <- fread(file.path(source_report,
  "selected_bdev_variable_genes.tsv"))
source_counts <- fread(file.path(source_report,
  "pseudotime_bin_cell_counts.tsv"))
source_diagnostics <- read_json(file.path(source_report, "diagnostics.json"),
  simplifyVector = FALSE)
stopifnot(setequal(unique(source_values$region_id),
    c("Flt3l_peak1", "Il7_peak2")),
  all(c("Bdev", "WT", "KO", "WT - KO") %in%
    unique(source_values$group)),
  all(source_selection[, .N, by = region_id]$N == 30L),
  isTRUE(source_diagnostics$selection_uses_query == FALSE))

all_differences <- list()
region_diagnostics <- list()
region_order <- c("Flt3l_peak1", "Il7_peak2")
for (region_key in region_order) {
  genes <- source_selection[region_id == region_key,
    feature_id]
  values <- source_values[region_id == region_key &
    group %in% c("Bdev", "WT", "KO")]
  counts <- source_counts[region_id == region_key]
  stopifnot(length(genes) == 30L,
    setequal(unique(values$feature_id), genes))

  bdev <- values[group == "Bdev", .(feature_id, bin, pt_mid,
    bdev_n = n_cells,
    bdev_scaled = bdev_scaled_mean_expression)]
  # Preserve the Bdev-only row ordering used by the source analysis.
  bdev_wide <- dcast(bdev, feature_id ~ bin, value.var = "bdev_scaled")
  bdev_matrix <- as.matrix(bdev_wide[, -"feature_id"])
  rownames(bdev_matrix) <- bdev_wide$feature_id
  row_order <- rownames(bdev_matrix)[hclust(dist(bdev_matrix))$order]

  difference_list <- lapply(c("WT", "KO"), function(query_group) {
    query <- values[group == query_group, .(feature_id, bin, pt_mid,
      query_n = n_cells,
      query_scaled = bdev_scaled_mean_expression)]
    out <- merge(bdev, query, by = c("feature_id", "bin", "pt_mid"),
      all = TRUE, sort = FALSE)
    out[, `:=`(query_group = query_group,
      difference_bdev_sd = query_scaled - bdev_scaled,
      low_cell_bin = bdev_n < low_cell_threshold |
        query_n < low_cell_threshold,
      calculable = is.finite(bdev_scaled) & is.finite(query_scaled))]
    out
  })
  difference <- rbindlist(difference_list)
  difference[, region_id := region_key]
  if (region_key == "Flt3l_peak1") {
    panel_labels <- c(WT = "FLT3L WT − Bdev",
      KO = "FLT3LKO − Bdev")
    title <- "Flt3L enriched-region expression difference from Bdev"
    file_stub <- "flt3l_peak1_query_minus_bdev_gene_expression_heatmap"
  } else {
    panel_labels <- c(WT = "IL-7 WT − Bdev",
      KO = "IL7KO − Bdev")
    title <- "IL-7 enriched-region expression difference from Bdev"
    file_stub <- "il7_peak2_query_minus_bdev_gene_expression_heatmap"
  }
  difference[, panel_label := factor(panel_labels[query_group],
    levels = unname(panel_labels))]
  difference[, feature_id := factor(feature_id, levels = rev(row_order))]
  all_differences[[region_key]] <- copy(difference)

  pt_limits <- range(counts$pt_left, counts$pt_right)
  caption_text <- paste0("Positive values indicate higher query expression; ",
    "negative values indicate higher Bdev expression. Black borders/crosses ",
    "mark either constituent bin with <", low_cell_threshold,
    " cells; gray is not calculable. Display clipped at ±", display_limit,
    " Bdev SD; exported differences are not clipped.")
  p <- ggplot(difference, aes(pt_mid, feature_id,
      fill = difference_bdev_sd)) +
    geom_tile(na.rm = FALSE) +
    geom_tile(data = difference[low_cell_bin == TRUE], fill = NA,
      color = "black", linewidth = 0.28) +
    geom_point(data = difference[low_cell_bin == TRUE], shape = 4,
      size = 0.65, stroke = 0.3, color = "black") +
    facet_grid(. ~ panel_label) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-display_limit, display_limit),
      oob = scales::squish, na.value = "grey65",
      name = "Query − Bdev\n(Bdev SD)") +
    scale_x_continuous(limits = pt_limits, expand = expansion(mult = 0),
      breaks = pretty(pt_limits, n = 5)) +
    theme_minimal(base_size = 9.5) +
    theme(panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10), legend.position = "right") +
    labs(title = title,
      subtitle = paste0("The source top-30 genes and Bdev scaling are fixed; ",
        "differences compare matched pseudotime-bin means"),
      x = "Bdev pseudotime", y = NULL,
      caption = paste(strwrap(caption_text, width = 145), collapse = "\n"))
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".png")), p,
    width = 12, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".pdf")), p,
    width = 12, height = 10, device = cairo_pdf, bg = "white")

  region_diagnostics[[region_key]] <- list(
    genes = length(genes), bins = uniqueN(difference$bin),
    calculable_wt_comparisons = difference[query_group == "WT",
      sum(calculable)],
    calculable_ko_comparisons = difference[query_group == "KO",
      sum(calculable)],
    low_cell_comparisons = difference[low_cell_bin == TRUE, .N])
}

differences <- rbindlist(all_differences)
differences[, feature_id := as.character(feature_id)]
differences[, panel_label := as.character(panel_label)]
fwrite(differences, file.path(stage_dir,
  "query_minus_bdev_gene_expression_differences.tsv.gz"), sep = "\t",
  na = "NA", compress = "gzip")
write_tsv(source_selection, file.path(stage_dir,
  "source_selected_bdev_variable_genes.tsv"))
diagnostics <- list(status = "complete_preliminary_unregistered",
  source_report = source_report,
  source_selection_uses_query = FALSE,
  comparisons = c("WT_minus_Bdev", "KO_minus_Bdev"),
  difference_unit = "Bdev temporal standard deviations",
  all_nonempty_source_bins_retained = TRUE,
  low_cell_highlight_threshold = low_cell_threshold,
  display_clip_bdev_sd = display_limit,
  exported_values_clipped = FALSE,
  region_diagnostics = region_diagnostics,
  biological_replicates_per_genotype = 1L,
  inference = "descriptive", bundle_mutated = FALSE,
  capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Query-minus-Bdev gene-expression difference heatmaps", "",
  "The source analysis' independently selected top-30 genes are retained for each cytokine region.",
  "Panels show WT minus matched-bin Bdev and KO minus matched-bin Bdev in Bdev temporal-SD units.",
  "Positive values indicate higher query expression; negative values indicate higher Bdev expression.",
  "Every calculable bin is included. A black border and cross marks a comparison where either constituent bin has fewer than 20 cells; gray is not calculable.",
  "There is one biological sample per genotype, so all query comparisons are descriptive."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed query-minus-Bdev expression heatmaps at", output_dir, "\n")
