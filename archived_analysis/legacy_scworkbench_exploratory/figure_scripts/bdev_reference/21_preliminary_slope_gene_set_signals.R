#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "Usage: 21_preliminary_slope_gene_set_signals.R",
  "<bdev_bundle> <w30_s5_report> <smoothed_report> <turning_report>",
  "<lineage_membership> <parse_membership> <output_dir>",
  "<expected_bundle_fingerprint> <loess_span>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
smoothed_report <- normalizePath(args[[3L]], mustWork = TRUE)
turning_report <- normalizePath(args[[4L]], mustWork = TRUE)
lineage_path <- normalizePath(args[[5L]], mustWork = TRUE)
parse_path <- normalizePath(args[[6L]], mustWork = TRUE)
output_dir <- normalizePath(args[[7L]], mustWork = FALSE)
expected_fingerprint <- args[[8L]]
loess_span <- as.numeric(args[[9L]])
stopifnot(is.finite(loess_span), loess_span > 0, loess_span <= 1)
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

read_tsv <- function(path) read.delim(path, check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = "NA")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

slope <- read_tsv(file.path(smoothed_report, "smoothed_slope_windows.tsv"))
trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
turning <- read_tsv(file.path(turning_report, "selected_turning_regions.tsv"))
lineage_source <- read_tsv(lineage_path)
parse_source <- read_tsv(parse_path)
stopifnot(nrow(slope) == 2524L, nrow(trajectory) == 12642L,
  !anyDuplicated(trajectory$cell_id), nrow(turning) == 4L,
  all(table(membership$window_id) == 30L),
  all(c("lineage", "feature_id", "mapped_to_matrix") %in%
    names(lineage_source)),
  all(c("cytokine", "direction", "feature_id", "mapped_to_matrix") %in%
    names(parse_source)))

lineage_source <- lineage_source[lineage_source$mapped_to_matrix, ]
parse_source <- parse_source[parse_source$mapped_to_matrix &
  parse_source$direction == "up" &
  parse_source$cytokine %in% c("Flt3L", "IL-7"), ]
source_gene_sets <- list(
  `B lineage` = unique(lineage_source$feature_id[
    lineage_source$lineage == "B lineage"]),
  `T lineage` = unique(lineage_source$feature_id[
    lineage_source$lineage == "T lineage"]),
  `Flt3L UP` = unique(parse_source$feature_id[
    parse_source$cytokine == "Flt3L"]),
  `IL-7 UP` = unique(parse_source$feature_id[
    parse_source$cytokine == "IL-7"]))
expected_sizes <- c(`B lineage` = 246L, `T lineage` = 233L,
  `Flt3L UP` = 45L, `IL-7 UP` = 20L)
stopifnot(identical(lengths(source_gene_sets), expected_sizes))

region_colors <- c(T1 = "#0072B2", T2 = "#D55E00", T3 = "#009E73",
  T4 = "#CC79A7")
turning$region_id <- factor(turning$region_id, levels = names(region_colors))
stopifnot(!anyNA(turning$region_id))
turning_bands <- turning[, c("region_id", "pt_min", "pt_max")]

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))
available <- wb_feature_ids(bundle, "SCT_data")
gene_sets <- lapply(source_gene_sets, function(x) x[x %in% available])
missing_sets <- lapply(source_gene_sets, function(x) x[!x %in% available])
stopifnot(all(lengths(gene_sets) > 0L))

cell_index <- match(membership$cell_id, trajectory$cell_id)
stopifnot(!anyNA(cell_index))
membership_split <- split(cell_index, factor(membership$window_id,
  levels = slope$window_id))
stopifnot(length(membership_split) == nrow(slope),
  all(lengths(membership_split) == 30L))

smoothing_fallback_sets <- character()
smooth_window_values <- function(values, set_name) {
  fit_data <- data.frame(x = slope$window_index, response = values)
  fit <- tryCatch(loess(response ~ x, data = fit_data, span = loess_span,
      degree = 2L, family = "symmetric", normalize = TRUE,
      control = loess.control(surface = "direct", statistics = "approximate")),
    error = function(e) {
      smoothing_fallback_sets <<- c(smoothing_fallback_sets, set_name)
      loess(response ~ x, data = fit_data, span = loess_span,
        degree = 2L, family = "gaussian", normalize = TRUE,
        control = loess.control(surface = "direct", statistics = "approximate"))
    })
  predicted <- as.numeric(predict(fit,
    newdata = data.frame(x = slope$window_index)))
  stopifnot(length(predicted) == length(values), all(is.finite(predicted)))
  predicted
}

score_rows <- lapply(names(gene_sets), function(set_name) {
  genes <- gene_sets[[set_name]]
  expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", genes,
    trajectory$cell_id, missing = "NA"))
  stopifnot(identical(rownames(expression), genes),
    identical(colnames(expression), trajectory$cell_id),
    all(is.finite(expression)))
  cell_score <- colMeans(expression)
  raw <- vapply(membership_split, function(ix) mean(cell_score[ix]), numeric(1))
  data.frame(window_id = slope$window_id,
    window_index = slope$window_index, pt_median = slope$pt_median,
    gene_set = set_name, mapped_genes = length(genes),
    mean_sct_data_score = raw,
    smoothed_mean_sct_data_score = smooth_window_values(raw, set_name),
    stringsAsFactors = FALSE)
})
scores <- do.call(rbind, score_rows)
set_order <- names(gene_sets)
scores$gene_set <- factor(scores$gene_set, levels = set_order)
stopifnot(nrow(scores) == length(set_order) * nrow(slope),
  all(is.finite(scores$smoothed_mean_sct_data_score)))

coverage <- data.frame(gene_set = set_order,
  source_genes = lengths(source_gene_sets),
  bdev_sct_mapped_genes = lengths(gene_sets),
  bdev_sct_missing_genes = lengths(missing_sets),
  missing_feature_ids = vapply(missing_sets, paste, collapse = ";",
    FUN.VALUE = character(1)),
  source = c(rep(lineage_path, 2L), rep(parse_path, 2L)),
  stringsAsFactors = FALSE)
membership_rows <- do.call(rbind, lapply(set_order, function(set_name)
  data.frame(gene_set = set_name, feature_id = gene_sets[[set_name]],
    stringsAsFactors = FALSE)))
write_tsv(scores, file.path(stage_dir, "gene_set_window_scores.tsv"))
write_tsv(coverage, file.path(stage_dir, "gene_set_coverage.tsv"))
write_tsv(membership_rows, file.path(stage_dir, "gene_set_membership.tsv"))

shared_x <- range(slope$pt_median)
base_theme <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
    plot.margin = margin(3, 8, 3, 4))
slope_plot <- ggplot(slope, aes(pt_median, smoothed_local_slope)) +
  geom_rect(data = turning_bands,
    aes(xmin = pt_min, xmax = pt_max, fill = region_id),
    ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.20) +
  geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.65,
    color = "#4A4A4A") +
  geom_line(color = "#4B1D82", linewidth = 1.05) +
  geom_point(data = turning, aes(pt, smoothed_slope, color = region_id),
    inherit.aes = FALSE, size = 2.3) +
  geom_text(data = turning,
    aes(pt, smoothed_slope, label = region_id, color = region_id),
    inherit.aes = FALSE, vjust = -0.9, size = 3.2, fontface = "bold",
    show.legend = FALSE) +
  scale_x_continuous(limits = shared_x) + scale_y_log10() +
  scale_fill_manual(values = region_colors, guide = "none") +
  scale_color_manual(values = region_colors, guide = "none") +
  base_theme +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  labs(title = "A  RNA/ATAC local slope",
    y = "ΔRNA / ΔATAC\n(log10 scale)")

set_colors <- c(`B lineage` = "#0072B2", `T lineage` = "#6A3D9A",
  `Flt3L UP` = "#D55E00", `IL-7 UP` = "#009E73")
panel_letters <- setNames(LETTERS[2:5], set_order)
make_set_plot <- function(set_name, show_x = FALSE) {
  plot_data <- scores[scores$gene_set == set_name, ]
  p <- ggplot(plot_data, aes(pt_median, smoothed_mean_sct_data_score)) +
    geom_rect(data = turning_bands,
      aes(xmin = pt_min, xmax = pt_max, fill = region_id),
      ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.20) +
    geom_line(color = set_colors[[set_name]], linewidth = 1) +
    scale_x_continuous(limits = shared_x) +
    scale_fill_manual(values = region_colors, guide = "none") +
    base_theme +
    labs(title = paste0(panel_letters[[set_name]], "  ", set_name,
        " (", length(gene_sets[[set_name]]), " genes)"),
      x = if (show_x) "Bdev pseudotime (window median pt)" else NULL,
      y = "Mean SCT_data\ngene-set score")
  if (!show_x) p <- p + theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank())
  p
}
set_plots <- lapply(seq_along(set_order), function(i)
  make_set_plot(set_order[[i]], show_x = i == length(set_order)))
combined <- wrap_plots(c(list(slope_plot), set_plots), ncol = 1,
  heights = c(1.12, rep(1, length(set_order)))) +
  plot_annotation(
    title = "Bdev modality progression and lineage/cytokine gene-set signals",
    subtitle = paste(
      "Arithmetic-mean SCT_data scores in 30-cell windows with 5-cell stride;",
      "robust LOESS span 0.08"),
    caption = paste(
      "T1–T4 are the previously selected persistent slope turning windows.",
      "Gene sets reuse the accepted mapped memberships. Preliminary and read-only."),
    theme = theme(plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10.5),
      plot.caption = element_text(size = 8.5)))

ggsave(file.path(stage_dir, "slope_b_t_lineage_flt3l_il7_up_signals.png"),
  combined, width = 10.5, height = 13, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "slope_b_t_lineage_flt3l_il7_up_signals.pdf"),
  combined, width = 10.5, height = 13, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id,
  bundle_fingerprint = expected_fingerprint,
  expression_matrix = "SCT_data",
  score = "arithmetic mean across mapped gene-set members",
  source_gene_set_sizes = as.list(expected_sizes),
  bdev_sct_mapped_gene_set_sizes = as.list(lengths(gene_sets)),
  bdev_sct_missing_gene_set_sizes = as.list(lengths(missing_sets)),
  window_size = 30L, window_stride = 5L,
  smoothing = list(method = "robust LOESS", span = loess_span,
    fallback_sets = unique(smoothing_fallback_sets)),
  lineage_source = lineage_path, parse_signal_source = parse_path,
  bundle_mutated = FALSE, capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Preliminary lineage and cytokine gene-set trajectories", "",
  "The first panel reuses the smoothed RNA/ATAC local slope and T1–T4 bands.",
  "The remaining panels show arithmetic-mean SCT_data scores for accepted",
  paste0("B-lineage (", length(gene_sets[["B lineage"]]), "), T-lineage (",
    length(gene_sets[["T lineage"]]), "), Flt3L-UP (",
    length(gene_sets[["Flt3L UP"]]), "), and IL-7-UP (",
    length(gene_sets[["IL-7 UP"]]), ") Bdev-mapped genes are scored in the ",
    "same 30-cell/5-cell windows."),
  "Bdev mapping coverage and missing feature IDs are exported explicitly.",
  "The Bdev bundle was opened read-only and no capsule was registered."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed preliminary gene-set signal report at", output_dir, "\n")
