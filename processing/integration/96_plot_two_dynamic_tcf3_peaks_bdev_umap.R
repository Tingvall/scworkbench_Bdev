#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 96_plot_two_dynamic_tcf3_peaks_bdev_umap.R BDEV_BUNDLE",
  "SOURCE_VALUES THREE_WAY_TABLE OUTPUT_DIR EXPECTED_BDEV_FP",
  "SCRIPT_VERSION"), call. = FALSE)

bdev_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
source_values_path <- normalizePath(args[[2L]], mustWork = TRUE)
three_way_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_bdev_fp <- args[[5L]]
script_version <- args[[6L]]

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
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
})

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

peak_ids <- c("chr11-45089034-45090047", "chr11-45847272-45848407")
three_way <- fread(three_way_path)
selected <- three_way[match(peak_ids, bdev_feature_id)]
stopifnot(identical(selected$bdev_feature_id, peak_ids),
  all(selected$tcf3_route_pass == TRUE), all(selected$tcf3_best_log2fc < -1),
  all(selected$tcf3_best_padj < 0.05))

bundle <- wb_open_bundle(bdev_bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_bdev_fp))
cells <- as.character(wb_cells(bundle, "all"))
space <- wb_space(bundle, "umap.atac", cells = cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells))
  coords <- data.table(cell_id = cells, umap_1 = as.numeric(space$x),
    umap_2 = as.numeric(space$y))
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  coords <- data.table(cell_id = cells, umap_1 = as.numeric(space[, 1L]),
    umap_2 = as.numeric(space[, 2L]))
}
stopifnot(all(is.finite(coords$umap_1)), all(is.finite(coords$umap_2)))

matrix_id <- "peaks_log1p_cp10k_ebf1_three_way_v1"
source_values <- fread(cmd = paste("zcat", shQuote(source_values_path)))[
  bdev_feature_id %in% peak_ids]
stopifnot(nrow(source_values) == length(cells) * length(peak_ids),
  setequal(source_values$bdev_feature_id, peak_ids),
  setequal(source_values$cell_id, cells),
  !anyDuplicated(source_values[, .(cell_id, bdev_feature_id)]),
  all(is.finite(source_values$log1p_cp10k_accessibility)),
  all(source_values$log1p_cp10k_accessibility >= 0))
accessibility <- matrix(NA_real_, nrow = length(peak_ids), ncol = length(cells),
  dimnames = list(peak_ids, cells))
for (id in peak_ids) {
  x <- source_values[bdev_feature_id == id]
  accessibility[id, x$cell_id] <- x$log1p_cp10k_accessibility
}
stopifnot(!anyNA(accessibility))

shared_cap <- as.numeric(quantile(accessibility[accessibility > 0], 0.99,
  names = FALSE, type = 7))
stopifnot(is.finite(shared_cap), shared_cap > 0)
summary_rows <- list()
plot_rows <- list()
value_rows <- list()
for (i in seq_along(peak_ids)) {
  id <- peak_ids[[i]]
  values <- as.numeric(accessibility[i, ])
  detected <- values > 0
  display <- pmin(values, shared_cap)
  display[!detected] <- NA_real_
  summary_rows[[i]] <- data.table(
    feature_id = id, cells = length(cells), detected_cells = sum(detected),
    detected_percent = 100 * mean(detected), mean_log1p_cp10k = mean(values),
    shared_positive_p99_display_cap = shared_cap,
    tcf3_log2fc_depleted_vs_control = selected$tcf3_best_log2fc[[i]],
    tcf3_bh_q = selected$tcf3_best_padj[[i]],
    placseq_support = selected$placseq_route_pass[[i]])
  label <- sprintf("%s\n%.1f%% accessible; TCF3 depletion log2FC %.2f",
    id, 100 * mean(detected), selected$tcf3_best_log2fc[[i]])
  plot_rows[[i]] <- data.table(umap_1 = coords$umap_1,
    umap_2 = coords$umap_2, accessibility = display,
    panel = label, panel_order = i)
  value_rows[[i]] <- data.table(cell_id = cells, feature_id = id,
    log1p_cp10k_accessibility = values)
}
summary_dt <- rbindlist(summary_rows)
plot_dt <- rbindlist(plot_rows)
value_dt <- rbindlist(value_rows)
plot_dt[, panel := factor(panel, levels = unique(panel[order(panel_order)]))]
setorder(plot_dt, panel_order, accessibility, na.last = FALSE)

plot <- ggplot() +
  geom_point(data = plot_dt[is.na(accessibility)], aes(umap_1, umap_2),
    color = "#D9DDE3", size = 0.12, alpha = 0.5, stroke = 0) +
  geom_point(data = plot_dt[!is.na(accessibility)],
    aes(umap_1, umap_2, color = accessibility),
    size = 0.65, alpha = 0.95, stroke = 0) +
  facet_wrap(~panel, nrow = 1) +
  scale_color_viridis_c(option = "magma", limits = c(0, shared_cap),
    oob = scales::squish, name = "Peak accessibility\nlog1p(CP10K)",
    breaks = c(0, shared_cap / 2, shared_cap),
    labels = c("0+", sprintf("%.2f", shared_cap / 2),
      sprintf("%.2f+", shared_cap))) +
  coord_equal() +
  labs(title = "Bdev accessibility of the two dynamic TCF3-dependent elements",
    subtitle = paste("All 26,723 Bdev cells on the unchanged native ATAC UMAP;",
      "one shared color scale across both peaks"),
    caption = paste("Gray cells have zero accessibility. Colored values are",
      "native log1p(CP10K), capped only for display at the pooled positive p99.")) +
  theme_void(base_size = 10) +
  theme(strip.text = element_text(size = 10, face = "bold", lineheight = 1),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 8, color = "grey30", hjust = 0),
    legend.position = "bottom", legend.key.width = grid::unit(25, "mm"),
    plot.margin = margin(8, 8, 8, 8))

ggsave(file.path(stage_dir, "figures",
  "two_dynamic_tcf3_peaks_bdev_umap.png"), plot,
  width = 11, height = 5.8, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "two_dynamic_tcf3_peaks_bdev_umap.pdf"), plot,
  width = 11, height = 5.8, device = cairo_pdf, bg = "white")
write_tsv(summary_dt, file.path(stage_dir, "peak_accessibility_summary.tsv"))
fwrite(value_dt, file.path(stage_dir,
  "all_cells_two_peak_log1p_cp10k_accessibility.tsv.gz"), sep = "\t",
  na = "NA", compress = "gzip")

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
write_json(list(schema_version = "two_dynamic_tcf3_peaks_bdev_umap.v1",
  status = "complete_unpromoted_report_only", script_version = script_version,
  bdev_bundle_fingerprint = fingerprint_before, embedding = "umap.atac",
  accessibility_matrix = matrix_id, cells = length(cells), peaks = peak_ids,
  source_values = source_values_path,
  shared_positive_p99_display_cap = shared_cap,
  bundle_mutated = FALSE, promotion = FALSE),
  file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Two dynamic TCF3-dependent elements on Bdev UMAP", "",
  "The two Bdev-dynamic elements that overlap significant E2A/TCF3-depletion peaks are displayed across all 26,723 cells on the unchanged native Bdev `umap.atac` embedding.", "",
  "Gray cells have zero peak accessibility. Colored cells show native scWorkbench log1p(CP10K) accessibility. Both panels share one color scale, capped for display at the pooled nonzero p99; uncapped values are exported.", "",
  "![Bdev UMAP](figures/two_dynamic_tcf3_peaks_bdev_umap.png)", "",
  "- [Vector PDF](figures/two_dynamic_tcf3_peaks_bdev_umap.pdf)",
  "- [Peak summary](peak_accessibility_summary.tsv)",
  "- [All-cell uncapped values](all_cells_two_peak_log1p_cp10k_accessibility.tsv.gz)"),
  file.path(stage_dir, "REPORT.md"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed two-peak Bdev UMAP report at", output_dir, "\n")
