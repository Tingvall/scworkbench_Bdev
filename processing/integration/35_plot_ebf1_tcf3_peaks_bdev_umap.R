#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 35_plot_ebf1_tcf3_peaks_bdev_umap.R BDEV_BUNDLE SESSION_WORKSPACE",
  "THREE_WAY_TABLE OUTPUT_DIR EXPECTED_BDEV_FP SCRIPT_VERSION SEED"
), call. = FALSE)

bdev_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
session_workspace <- normalizePath(args[[2L]], mustWork = TRUE)
three_way_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_bdev_fp <- args[[5L]]
script_version <- args[[6L]]
seed <- as.integer(args[[7L]])
stopifnot(is.finite(seed))

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)
set.seed(seed)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

log_event <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(..., collapse = " "), "\n")
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

results <- fread(three_way_path)
required <- c("bdev_feature_id", "tcf3_route_pass", "tcf3_best_log2fc",
  "tcf3_best_padj", "placseq_route_pass", "priority_rank")
stopifnot(all(required %in% names(results)))
selected <- results[tcf3_route_pass == TRUE]
setorder(selected, priority_rank, bdev_feature_id)
stopifnot(nrow(selected) == 8L, !anyDuplicated(selected$bdev_feature_id),
  all(selected$tcf3_best_log2fc < -1),
  all(selected$tcf3_best_padj < 0.05))

bundle <- wb_open_bundle(bdev_bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
bdev_fp <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(bdev_fp, expected_bdev_fp))
cells <- as.character(wb_cells(bundle, "all"))
stopifnot(length(cells) == 26723L, !anyDuplicated(cells))

space_info <- wb_space_info(bundle, "umap.atac")
stopifnot(nrow(space_info) == 1L, as.integer(space_info$dims[[1L]]) == 2L)
space <- wb_space(bundle, "umap.atac", cells = cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells),
    all(c("x", "y") %in% names(space)))
  coords <- data.table(cell_id = cells, umap_1 = as.numeric(space$x),
    umap_2 = as.numeric(space$y))
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  coords <- data.table(cell_id = cells, umap_1 = as.numeric(space[, 1L]),
    umap_2 = as.numeric(space[, 2L]))
}
stopifnot(all(is.finite(coords$umap_1)), all(is.finite(coords$umap_2)))

md <- as.data.table(wb_meta(bundle, cells = cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt")))
setnames(md, names(md)[1L], "cell_id")
stopifnot(identical(as.character(md$cell_id), cells))

session <- wb_session_resume(bundle, session_workspace)
stopifnot(is.null(session$durable$promotion_receipt))
view <- wb_session_view(session)
on.exit(try(wb_close_bundle(view), silent = TRUE), add = TRUE)
matrix_id <- "peaks_log1p_cp10k_ebf1_three_way_v1"
peaks <- as.character(selected$bdev_feature_id)
stopifnot(all(peaks %in% wb_feature_ids(view, matrix_id)))
log_event("Reading", length(peaks), "TCF3-dependent peaks for", length(cells),
  "Bdev cells")
accessibility <- as.matrix(wb_matrix_read_features(view, matrix_id, peaks,
  cells, missing = "NA"))
stopifnot(identical(rownames(accessibility), peaks),
  identical(colnames(accessibility), cells), !anyNA(accessibility),
  all(is.finite(accessibility)), all(accessibility >= 0))

summary_rows <- vector("list", length(peaks))
plot_rows <- vector("list", length(peaks))
value_rows <- vector("list", length(peaks))
for (i in seq_along(peaks)) {
  peak <- peaks[[i]]
  values <- as.numeric(accessibility[i, ])
  detected <- values > 0
  positive <- values[detected]
  cap <- if (length(positive)) as.numeric(quantile(positive, 0.99,
    names = FALSE, type = 7)) else 1
  if (!is.finite(cap) || cap <= 0) cap <- max(c(positive, 1))
  display_value <- pmin(values, cap) / cap
  display_value[!detected] <- NA_real_
  summary_rows[[i]] <- data.table(
    bdev_feature_id = peak,
    panel_order = i,
    cells = length(cells),
    detected_cells = sum(detected),
    detected_percent = 100 * mean(detected),
    mean_log1p_cp10k = mean(values),
    nonzero_median_log1p_cp10k = if (length(positive)) median(positive) else NA_real_,
    nonzero_p99_display_cap = cap,
    maximum_log1p_cp10k = max(values),
    tcf3_log2fc = selected$tcf3_best_log2fc[[i]],
    tcf3_padj = selected$tcf3_best_padj[[i]],
    placseq_support = selected$placseq_route_pass[[i]])
  label <- sprintf("%s\nTCF3 log2FC %.2f; padj %.2g\n%.2f%% detected%s",
    peak, selected$tcf3_best_log2fc[[i]], selected$tcf3_best_padj[[i]],
    100 * mean(detected),
    if (isTRUE(selected$placseq_route_pass[[i]])) "; PLAC-seq+" else "")
  plot_rows[[i]] <- data.table(umap_1 = coords$umap_1,
    umap_2 = coords$umap_2, display_accessibility = display_value,
    panel = label, panel_order = i)
  value_rows[[i]] <- data.table(cell_id = cells, pt = md$pt,
    bdev_feature_id = peak, log1p_cp10k_accessibility = values)
}
summary_dt <- rbindlist(summary_rows)
plot_dt <- rbindlist(plot_rows)
value_dt <- rbindlist(value_rows)
plot_dt[, panel := factor(panel, levels = unique(panel[order(panel_order)]))]
plot_dt[, detected_sort := !is.na(display_accessibility)]
setorder(plot_dt, panel_order, detected_sort, display_accessibility,
  na.last = FALSE)
plot_dt[, detected_sort := NULL]
write_tsv(summary_dt, file.path(stage_dir, "peak_accessibility_summary.tsv"))
fwrite(value_dt, file.path(stage_dir,
  "all_cells_eight_peak_log1p_cp10k_accessibility.tsv.gz"), sep = "\t",
  na = "NA", compress = "gzip")

umap_plot <- ggplot() +
  geom_point(data = plot_dt[is.na(display_accessibility)],
    aes(umap_1, umap_2), color = "#D9DDE3", size = 0.08,
    alpha = 0.45, stroke = 0) +
  geom_point(data = plot_dt[!is.na(display_accessibility)],
    aes(umap_1, umap_2, color = display_accessibility), size = 0.60,
    alpha = 0.95, stroke = 0) +
  facet_wrap(~panel, ncol = 4) +
  scale_color_viridis_c(option = "magma", limits = c(0, 1),
    oob = scales::squish, na.value = "#D9DDE3",
    breaks = c(0, 0.5, 1), labels = c("low", "mid", "p99+"),
    name = "Within-peak\naccessibility") +
  coord_equal() +
  labs(
    title = "Accessibility of eight TCF3-dependent Ebf1-window peaks",
    subtitle = paste0("All 26,723 Bdev cells on the unchanged ATAC UMAP; ",
      "native scWorkbench log1p(CP10K) peak accessibility"),
    caption = paste(
      "Gray indicates zero counts; larger colored points are nonzero cells.",
      "Positive values are capped at each peak's nonzero p99 and scaled to 0–1",
      "for spatial display; raw normalized values are retained in the table.")) +
  theme_void(base_size = 10) +
  theme(
    strip.text = element_text(size = 8, face = "bold", lineheight = 0.95,
      margin = margin(3, 2, 3, 2)),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 8, color = "grey30", hjust = 0),
    plot.margin = margin(8, 8, 8, 8),
    legend.position = "bottom",
    legend.key.width = grid::unit(24, "mm"))

ggsave(file.path(stage_dir, "figures",
  "eight_tcf3_dependent_peaks_bdev_umap.png"), umap_plot,
  width = 14, height = 7.8, dpi = 240, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "eight_tcf3_dependent_peaks_bdev_umap.pdf"), umap_plot,
  width = 14, height = 7.8, device = cairo_pdf, bg = "white")

write_json(list(
  schema_version = "ebf1_tcf3_peaks_bdev_umap.v3",
  status = "complete",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script_version = script_version,
  seed = seed,
  bdev_bundle_fingerprint = bdev_fp,
  embedding = "umap.atac",
  accessibility_matrix = matrix_id,
  cells = length(cells),
  peaks = length(peaks),
  three_way_table_sha256 = digest::digest(three_way_path, algo = "sha256",
    file = TRUE)
), file.path(stage_dir, "completion.json"))

writeLines(c(
  "# Bdev UMAP accessibility of TCF3-dependent Ebf1-window peaks", "",
  paste0("All **", length(peaks), "** Bdev peaks passing the TCF3 route are ",
    "shown across all **26,723** Bdev cells on the unchanged native ",
    "`umap.atac` embedding."), "",
  "Accessibility is the native scWorkbench log1p(CP10K) peak matrix created ",
    "for the three-way Ebf1 screen. Small gray cells have zero accessibility; ",
    "larger colored points are nonzero cells. For ",
    "visualization only, each peak's positive values are capped at its nonzero ",
    "99th percentile and scaled to 0–1; the downloadable table retains the ",
    "uncapped values.", "",
  "![Eight TCF3-dependent peaks](figures/eight_tcf3_dependent_peaks_bdev_umap.png)", "",
  "## Downloads", "",
  "- [PDF figure](figures/eight_tcf3_dependent_peaks_bdev_umap.pdf)",
  "- [Peak accessibility summary](peak_accessibility_summary.tsv)",
  "- [All-cell normalized accessibility](all_cells_eight_peak_log1p_cp10k_accessibility.tsv.gz)"
), file.path(stage_dir, "REPORT.md"))

html <- c(
  "<!doctype html><html><head><meta charset='utf-8'>",
  "<title>TCF3-dependent Ebf1 peaks on Bdev UMAP</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1500px;margin:2rem auto;padding:0 1rem;line-height:1.45}img{max-width:100%;height:auto}code{background:#f2f2f2;padding:.1rem .25rem}</style></head><body>",
  "<h1>Bdev UMAP accessibility of TCF3-dependent Ebf1-window peaks</h1>",
  "<p>All eight peaks are shown across all 26,723 Bdev cells on the unchanged <code>umap.atac</code> embedding. Small gray points denote zero accessibility; larger colored points are nonzero cells. Positive native log1p(CP10K) values are capped at each peak's nonzero p99 and scaled to 0–1 for display.</p>",
  "<img src='figures/eight_tcf3_dependent_peaks_bdev_umap.png'>",
  "<h2>Downloads</h2><ul>",
  "<li><a href='figures/eight_tcf3_dependent_peaks_bdev_umap.pdf'>PDF figure</a></li>",
  "<li><a href='peak_accessibility_summary.tsv'>Peak accessibility summary</a></li>",
  "<li><a href='all_cells_eight_peak_log1p_cp10k_accessibility.tsv.gz'>All-cell normalized accessibility</a></li>",
  "</ul></body></html>")
writeLines(html, file.path(stage_dir, "REPORT.html"))

stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed TCF3-dependent Ebf1 peak UMAP report at", output_dir)
