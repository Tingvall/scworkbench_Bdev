#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 26_preliminary_state_top20_separate_figures.R",
  "<bdev_bundle> <window_report> <smoothed_report> <gene_screen_report>",
  "<output_dir> <expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
smoothed_report <- normalizePath(args[[3L]], mustWork = TRUE)
gene_screen_report <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
expected_fingerprint <- args[[6L]]
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(Matrix)
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

trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
windows <- read_tsv(file.path(smoothed_report, "smoothed_slope_windows.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
score_connection <- gzfile(file.path(gene_screen_report,
  "modality_gene_pattern_scores.tsv.gz"), open = "rt")
scores <- read_tsv(score_connection)
close(score_connection)
stopifnot(nrow(trajectory) == 12642L, nrow(windows) == 2524L,
  !anyDuplicated(trajectory$cell_id), !anyDuplicated(windows$window_id),
  all(windows$smoothed_local_slope > 0),
  all(table(membership$window_id) == 30L),
  nrow(scores) == 11745L, !anyDuplicated(scores$feature_id),
  all(is.finite(scores$state_correlation)))

rna_top <- scores[order(-scores$state_correlation, scores$feature_id), ][1:20, ]
atac_top <- scores[order(scores$state_correlation, scores$feature_id), ][1:20, ]
change_top <- scores[scores$change_selected, ]
change_top <- change_top[order(-change_top$change_correlation,
  change_top$feature_id), ]
stopifnot(all(rna_top$state_correlation > 0),
  all(atac_top$state_correlation < 0), nrow(change_top) == 4L,
  all(change_top$change_correlation > 0))
rna_top$rank <- seq_len(nrow(rna_top))
atac_top$rank <- seq_len(nrow(atac_top))
write_tsv(rna_top, file.path(stage_dir, "rna_faster_state_top20.tsv"))
write_tsv(atac_top, file.path(stage_dir, "atac_faster_state_top20.tsv"))
change_top$rank <- seq_len(nrow(change_top))
write_tsv(change_top, file.path(stage_dir,
  "rna_faster_change_associated_four.tsv"))

cell_index <- match(membership$cell_id, trajectory$cell_id)
window_index <- match(membership$window_id, windows$window_id)
stopifnot(!anyNA(cell_index), !anyNA(window_index))
window_design <- sparseMatrix(i = cell_index, j = window_index,
  x = rep(1 / 30, nrow(membership)),
  dims = c(nrow(trajectory), nrow(windows)))
membership_split <- split(cell_index, factor(membership$window_id,
  levels = windows$window_id))
stopifnot(length(membership_split) == nrow(windows),
  all(lengths(membership_split) == 30L))

smoothing_fraction <- 0.08
smoothing_half_width <- floor(nrow(windows) * smoothing_fraction / 2)
smoothing_rows <- unlist(lapply(seq_len(nrow(windows)), function(j)
  seq.int(max(1L, j - smoothing_half_width),
    min(nrow(windows), j + smoothing_half_width))), use.names = FALSE)
smoothing_columns <- rep(seq_len(nrow(windows)),
  vapply(seq_len(nrow(windows)), function(j)
    min(nrow(windows), j + smoothing_half_width) -
      max(1L, j - smoothing_half_width) + 1L, numeric(1)))
smoothing_counts <- tabulate(smoothing_columns, nbins = nrow(windows))
smoothing_design <- sparseMatrix(i = smoothing_rows, j = smoothing_columns,
  x = 1 / smoothing_counts[smoothing_columns],
  dims = c(nrow(windows), nrow(windows)))
stopifnot(all(abs(colSums(window_design) - 1) < 1e-12),
  all(abs(colSums(smoothing_design) - 1) < 1e-12))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))
top_features <- unique(c(rna_top$feature_id, atac_top$feature_id,
  change_top$feature_id))
stopifnot(all(top_features %in% wb_feature_ids(bundle, "SCT_data")))
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data",
  top_features, trajectory$cell_id, missing = "NA"))
stopifnot(identical(rownames(expression), top_features),
  identical(colnames(expression), trajectory$cell_id),
  all(is.finite(expression)))

raw_mean <- as.matrix(expression %*% window_design)
raw_q25 <- vapply(membership_split, function(ix)
  apply(expression[, ix, drop = FALSE], 1L, quantile, probs = 0.25,
    names = FALSE, type = 7), numeric(length(top_features)))
raw_q75 <- vapply(membership_split, function(ix)
  apply(expression[, ix, drop = FALSE], 1L, quantile, probs = 0.75,
    names = FALSE, type = 7), numeric(length(top_features)))
rownames(raw_q25) <- rownames(raw_q75) <- top_features
smoothed_mean <- as.matrix(raw_mean %*% smoothing_design)
smoothed_q25 <- as.matrix(raw_q25 %*% smoothing_design)
smoothed_q75 <- as.matrix(raw_q75 %*% smoothing_design)
rownames(smoothed_mean) <- rownames(smoothed_q25) <-
  rownames(smoothed_q75) <- top_features
profile_center <- rowMeans(smoothed_mean)
profile_sd <- sqrt(rowSums((smoothed_mean - profile_center)^2) /
  (ncol(smoothed_mean) - 1L))
stopifnot(all(is.finite(profile_sd)), all(profile_sd > 0))
mean_z <- (smoothed_mean - profile_center) / profile_sd
q25_z <- (smoothed_q25 - profile_center) / profile_sd
q75_z <- (smoothed_q75 - profile_center) / profile_sd
rna_target_z <- as.numeric(scale(log2(windows$smoothed_local_slope)))
atac_target_z <- -rna_target_z

make_plot_data <- function(top, target) {
  rows <- lapply(seq_len(nrow(top)), function(i) {
    gene <- top$feature_id[[i]]
    data.frame(pt = windows$pt_median,
      gene_panel = sprintf("%02d  %s  (r = %.3f)", i, gene,
        top$state_correlation[[i]]),
      mean_z = mean_z[gene, ], q25_z = q25_z[gene, ],
      q75_z = q75_z[gene, ], target_z = target,
      stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  result$gene_panel <- factor(result$gene_panel,
    levels = unique(result$gene_panel))
  result
}
rna_plot_data <- make_plot_data(rna_top, rna_target_z)
atac_plot_data <- make_plot_data(atac_top, atac_target_z)

make_figure <- function(data, title, subtitle, gene_color, target_label) {
  ggplot(data, aes(pt)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "#888888",
      linewidth = 0.3) +
    geom_ribbon(aes(ymin = q25_z, ymax = q75_z, fill = "Gene IQR"),
      alpha = 0.20, color = NA) +
    geom_line(aes(y = mean_z, color = "Gene mean", linewidth = "Gene mean")) +
    geom_line(aes(y = target_z, color = target_label,
      linewidth = target_label)) +
    facet_wrap(~ gene_panel, ncol = 4L) +
    scale_fill_manual(values = c("Gene IQR" = gene_color)) +
    scale_color_manual(values = c("Gene mean" = gene_color,
      setNames("#222222", target_label))) +
    scale_linewidth_manual(values = c("Gene mean" = 0.8,
      setNames(0.55, target_label))) +
    theme_classic(base_size = 9) +
    theme(strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 8.5),
      panel.spacing = grid::unit(5, "mm"), legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10)) +
    guides(fill = guide_legend(order = 1), color = guide_legend(order = 2),
      linewidth = "none") +
    labs(title = title, subtitle = subtitle,
      x = "Bdev pseudotime", y = "Standardized smoothed SCT_data",
      fill = NULL, color = NULL)
}

rna_plot <- make_figure(rna_plot_data,
  "Top 20 genes following the RNA-faster state",
  paste("Ribbon: within-window IQR; colored line: smoothed mean; black line:",
    "standardized log2(delta RNA / delta ATAC)"),
  "#C51B5A", "RNA-faster target")
atac_plot <- make_figure(atac_plot_data,
  "Top 20 genes following the ATAC-faster state",
  paste("Ribbon: within-window IQR; colored line: smoothed mean; black line:",
    "standardized inverse log2(delta ATAC / delta RNA)"),
  "#2166AC", "ATAC-faster target")

change_expression_rows <- lapply(seq_len(nrow(change_top)), function(i) {
  gene <- change_top$feature_id[[i]]
  data.frame(pt = windows$pt_median,
    gene_panel = sprintf("%02d  %s  (change r = %.3f)", i, gene,
      change_top$change_correlation[[i]]),
    mean_z = mean_z[gene, ], q25_z = q25_z[gene, ],
    q75_z = q75_z[gene, ], stringsAsFactors = FALSE)
})
change_expression_data <- do.call(rbind, change_expression_rows)
change_expression_data$gene_panel <- factor(change_expression_data$gene_panel,
  levels = unique(change_expression_data$gene_panel))
change_derivative <- (smoothed_mean[, 3:nrow(windows), drop = FALSE] -
  smoothed_mean[, 1:(nrow(windows) - 2L), drop = FALSE]) / 2
change_derivative_z <- t(apply(change_derivative, 1L, function(x)
  as.numeric(scale(x))))
rownames(change_derivative_z) <- rownames(change_derivative)
change_target_z <- as.numeric(scale(
  log2(windows$smoothed_local_slope)[2:(nrow(windows) - 1L)]))
change_derivative_rows <- lapply(seq_len(nrow(change_top)), function(i) {
  gene <- change_top$feature_id[[i]]
  data.frame(pt = windows$pt_median[2:(nrow(windows) - 1L)],
    gene_panel = sprintf("%02d  %s  (change r = %.3f)", i, gene,
      change_top$change_correlation[[i]]),
    derivative_z = change_derivative_z[gene, ], target_z = change_target_z,
    stringsAsFactors = FALSE)
})
change_derivative_data <- do.call(rbind, change_derivative_rows)
change_derivative_data$gene_panel <- factor(change_derivative_data$gene_panel,
  levels = levels(change_expression_data$gene_panel))

change_expression_plot <- ggplot(change_expression_data, aes(pt)) +
  geom_ribbon(aes(ymin = q25_z, ymax = q75_z), fill = "#E66101",
    alpha = 0.20, color = NA) +
  geom_line(aes(y = mean_z), color = "#E66101", linewidth = 0.85) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "#888888",
    linewidth = 0.3) +
  facet_wrap(~ gene_panel, ncol = 4L) + theme_classic(base_size = 9.5) +
  theme(strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  labs(title = "Expression level and within-window IQR",
    y = "Standardized smoothed SCT_data")
change_derivative_plot <- ggplot(change_derivative_data, aes(pt)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "#888888",
    linewidth = 0.3) +
  geom_line(aes(y = derivative_z, color = "Expression derivative"),
    linewidth = 0.85) +
  geom_line(aes(y = target_z, color = "RNA-faster target"),
    linewidth = 0.6) +
  facet_wrap(~ gene_panel, ncol = 4L) +
  scale_color_manual(values = c("Expression derivative" = "#E66101",
    "RNA-faster target" = "#222222")) +
  theme_classic(base_size = 9.5) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
    legend.position = "bottom") +
  labs(title = "Screened change profile versus modality target",
    x = "Bdev pseudotime", y = "Standardized derivative", color = NULL)
change_plot <- change_expression_plot / change_derivative_plot +
  plot_layout(heights = c(1, 1.08)) +
  plot_annotation(
    title = "Four genes associated with RNA-faster expression change",
    subtitle = paste(
      "IQR is shown for expression levels; the lower tier shows the",
      "window-index derivative used by the change screen"),
    theme = theme(plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10)))

ggsave(file.path(stage_dir, "rna_faster_state_top20_with_iqr.png"), rna_plot,
  width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "rna_faster_state_top20_with_iqr.pdf"), rna_plot,
  width = 14, height = 12, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "atac_faster_state_top20_with_iqr.png"), atac_plot,
  width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "atac_faster_state_top20_with_iqr.pdf"), atac_plot,
  width = 14, height = 12, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "rna_faster_change_four_with_iqr.png"), change_plot,
  width = 14, height = 7.8, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "rna_faster_change_four_with_iqr.pdf"), change_plot,
  width = 14, height = 7.8, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  source_gene_screen = gene_screen_report, expression_matrix = "SCT_data",
  trajectory_cells = nrow(trajectory), windows = nrow(windows),
  window_size = 30L, stride = 5L, smoothing_fraction = smoothing_fraction,
  smoothing_half_width_windows = smoothing_half_width,
  interval = "within-window cell-level IQR, then fixed 8%-window smoothing",
  rna_state_definition = "top 20 positive state correlations",
  atac_state_definition = paste(
    "top 20 negative RNA/ATAC state correlations, displayed against",
    "the sign-reversed ATAC/RNA target"),
  rna_top20 = as.list(rna_top$feature_id),
  atac_top20 = as.list(atac_top$feature_id),
  change_associated_genes = as.list(change_top$feature_id),
  change_figure = paste(
    "upper tier: smoothed expression mean and within-window IQR;",
    "lower tier: standardized window-index derivative and RNA-faster target"),
  bundle_mutated = FALSE, capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Separate top-20 modality-state figures", "",
  "The RNA-faster figure contains the twenty largest positive state-following",
  "correlations. The ATAC-faster figure contains the twenty most negative",
  "RNA/ATAC correlations and uses the exact sign-reversed ATAC/RNA target.",
  "Each panel shows the fixed 8%-window-smoothed mean and the smoothed",
  "within-window cell-level IQR from the exact 30-cell windows.",
  "A companion figure shows the four RNA-faster change-associated genes:",
  "Angpt1, Adgrl4, Gcnt2, and Mecom. Its upper tier shows expression plus IQR;",
  "its lower tier shows the derivative used by the change screen and the target.",
  "The bundle was opened read-only."), file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed separate top-20 state figures at", output_dir, "\n")
