#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 39_preliminary_adjusted_selected_genes_bb.R",
  "<bundle> <window_report> <ratio_report> <direction_screen>",
  "<output_dir> <expected_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
ratio_report <- normalizePath(args[[3L]], mustWork = TRUE)
screen_dir <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
expected_fingerprint <- args[[6L]]
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(Matrix)
  library(matrixStats)
  library(ggplot2)
  library(splines)
  library(jsonlite)
})
read_tsv <- function(path) read.delim(path, check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = "NA")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")

genes <- c("Ebf1", "Tcf3", "Stat5a", "Stat5b", "Ezh2", "Nipbl", "Wapl")
trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
windows <- read_tsv(file.path(window_report, "window_distance_and_slope.tsv"))
bb <- read_tsv(file.path(window_report, "bb1_bb2_regions.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
ratio <- read_tsv(file.path(ratio_report,
  "ratio_total_signal_fitted_and_residual_windows.tsv"))
score_connection <- gzfile(file.path(screen_dir,
  "residual_ratio_weighted_gene_scores.tsv.gz"), "rt")
scores <- read_tsv(score_connection)
close(score_connection)
segment_columns <- paste0("segment_", 1:8,
  "_direction_similarity_direct")
stopifnot(all(genes %in% scores$feature_id),
  all(segment_columns %in% names(scores)),
  identical(bb$boundary, c("BB1", "BB2")), nrow(windows) == nrow(ratio))
audit <- scores[match(genes, scores$feature_id), ]
segment_matrix <- as.matrix(audit[, segment_columns])
audit$direct_segments <- rowSums(segment_matrix > 0)
audit$inverse_segments <- rowSums(segment_matrix < 0)
audit$preferred_direction <- ifelse(
  audit$evidence_weighted_score_direct >=
    audit$evidence_weighted_score_inverse, "direct", "inverse")
write_tsv(audit, file.path(stage_dir,
  "selected_adjusted_gene_direction_scores.tsv"))
write_tsv(bb, file.path(stage_dir, "bb1_bb2_regions.tsv"))

cell_index <- match(membership$cell_id, trajectory$cell_id)
window_index <- match(membership$window_id, windows$window_id)
stopifnot(!anyNA(cell_index), !anyNA(window_index))
window_design <- sparseMatrix(i = cell_index, j = window_index,
  x = rep(1 / 30, nrow(membership)),
  dims = c(nrow(trajectory), nrow(windows)))
membership_split <- split(cell_index, factor(membership$window_id,
  levels = windows$window_id))
smoothing_half_width <- floor(nrow(windows) * 0.03 / 2)
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
target <- as.numeric(t(ratio$full_residual) %*% smoothing_design)
target_z <- as.numeric(scale(target))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint),
  all(genes %in% wb_feature_ids(bundle, "SCT_data")))
metadata <- wb_meta(bundle, cells = trajectory$cell_id,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt", "nCount_RNA"))
stopifnot(identical(as.character(metadata$cell_id), trajectory$cell_id))
depth_z <- as.numeric(scale(log1p(metadata$nCount_RNA)))
depth_design <- cbind(intercept = 1, ns_pt = ns(trajectory$pt, df = 8),
  log_RNA_depth_z = depth_z)
depth_solver <- solve(crossprod(depth_design))
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", genes,
  trajectory$cell_id, missing = "NA"))
coefficients <- depth_solver %*% crossprod(depth_design, t(expression))
depth_beta <- as.numeric(coefficients[nrow(coefficients), ])
adjusted <- expression - tcrossprod(depth_beta, depth_z)
profile <- as.matrix((adjusted %*% window_design) %*% smoothing_design)
q25 <- q75 <- matrix(NA_real_, nrow = length(genes), ncol = nrow(windows))
for (j in seq_len(nrow(windows))) {
  q <- rowQuantiles(adjusted[, membership_split[[j]], drop = FALSE],
    probs = c(0.25, 0.75), drop = FALSE)
  q25[, j] <- q[, 1L]
  q75[, j] <- q[, 2L]
}
q25 <- as.matrix(q25 %*% smoothing_design)
q75 <- as.matrix(q75 %*% smoothing_design)
rownames(profile) <- rownames(q25) <- rownames(q75) <- genes
center <- rowMeans(profile)
spread <- rowSds(profile)
profile_z <- (profile - center) / spread
q25_z <- (q25 - center) / spread
q75_z <- (q75 - center) / spread

rows <- list(data.frame(pt = windows$pt_median,
  panel = "00 Corrected RNA/ATAC ratio", mean = target_z,
  q25 = NA_real_, q75 = NA_real_, direction = "target"))
for (i in seq_along(genes)) {
  gene <- genes[[i]]
  score <- audit[audit$feature_id == gene, ]
  rows[[length(rows) + 1L]] <- data.frame(pt = windows$pt_median,
    panel = sprintf("%02d %s | direct %d/8 | inverse %d/8", i, gene,
      score$direct_segments, score$inverse_segments),
    mean = profile_z[gene, ], q25 = q25_z[gene, ], q75 = q75_z[gene, ],
    direction = score$preferred_direction)
}
plot_data <- do.call(rbind, rows)
plot_data$panel <- factor(plot_data$panel, levels = unique(plot_data$panel))
bb_fill <- c(BB1 = "#2C7FB8", BB2 = "#D95F0E")
p <- ggplot(plot_data, aes(pt, mean)) +
  geom_rect(data = bb, aes(xmin = pt_min, xmax = pt_max,
      ymin = -Inf, ymax = Inf, fill = boundary), inherit.aes = FALSE,
    alpha = 0.14, color = NA) +
  geom_vline(data = bb, aes(xintercept = pt_median, color = boundary),
    inherit.aes = FALSE, linetype = "dashed", linewidth = 0.4,
    show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "#999999", linetype = "dotted",
    linewidth = 0.25) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = direction), alpha = 0.16,
    color = NA, na.rm = TRUE) +
  geom_line(aes(color = direction), linewidth = 0.8) +
  facet_wrap(~panel, ncol = 2, scales = "fixed") +
  scale_color_manual(values = c(target = "#111111", direct = "#C51B5A",
    inverse = "#2166AC", BB1 = bb_fill[["BB1"]], BB2 = bb_fill[["BB2"]]),
    guide = "none") +
  scale_fill_manual(values = c(BB1 = bb_fill[["BB1"]],
    BB2 = bb_fill[["BB2"]], target = "#111111", direct = "#C51B5A",
    inverse = "#2166AC"), breaks = c("BB1", "BB2"),
    name = "Source bb_up region") +
  coord_cartesian(ylim = c(-3.5, 3.5)) +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = 8.5),
    panel.spacing = grid::unit(5, "mm"), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15)) +
  labs(title = "Corrected ratio and selected depth-adjusted genes",
    subtitle = paste("Full pseudotime; exact BB1/BB2 source ranges; identical",
      "3% smoothing; ribbons are adjusted-expression IQR"),
    x = "Bdev pseudotime", y = "Standardized adjusted profile")
ggsave(file.path(stage_dir, "adjusted_selected_genes_bb.png"), p,
  width = 14, height = 8.5, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "adjusted_selected_genes_bb.pdf"), p,
  width = 14, height = 8.5, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  genes = as.list(genes), depth_adjustment =
    "SCT_data minus conditional log1p(nCount_RNA) coefficient after ns(pt,df=8)",
  smoothing_fraction = 0.03,
  highlighted_regions = "source bb_up BB1 and BB2 pt min-to-max with medians",
  figure_width = 14, figure_height = 8.5, bundle_mutated = FALSE,
  capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"), pretty = TRUE,
  auto_unbox = TRUE)
writeLines(c("# Selected adjusted genes with BB regions", "",
  paste("Genes:", paste(genes, collapse = ", ")), "",
  "The compact landscape figure shows the full corrected ratio and standardized depth-adjusted expression with exact source BB1/BB2 shading."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed selected adjusted BB figure at", output_dir, "\n")
