#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 37_preliminary_residual_ratio_full_direction_matches.R",
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

trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
windows <- read_tsv(file.path(window_report, "window_distance_and_slope.tsv"))
con <- gzfile(file.path(window_report, "window_membership.tsv.gz"), "rt")
membership <- read_tsv(con)
close(con)
ratio <- read_tsv(file.path(ratio_report,
  "ratio_total_signal_fitted_and_residual_windows.tsv"))
score_connection <- gzfile(file.path(screen_dir,
  "residual_ratio_weighted_gene_scores.tsv.gz"), "rt")
scores <- read_tsv(score_connection)
close(score_connection)

segment_columns <- paste0("segment_", 1:8,
  "_direction_similarity_direct")
stopifnot(all(segment_columns %in% names(scores)), nrow(scores) == 11745L,
  nrow(trajectory) == 12642L, nrow(windows) == 2524L,
  identical(windows$window_id, ratio$window_id))
segment_matrix <- as.matrix(scores[, segment_columns])
scores$concordant_segments_direct <- rowSums(segment_matrix > 0)
scores$concordant_segments_inverse <- rowSums(segment_matrix < 0)
scores$weakest_segment_direct <- rowMins(segment_matrix)
scores$weakest_segment_inverse <- -rowMaxs(segment_matrix)
direct <- scores[scores$concordant_segments_direct == 8L, ]
inverse <- scores[scores$concordant_segments_inverse == 8L, ]
direct <- direct[order(-direct$evidence_weighted_score_direct,
  direct$feature_id), ]
inverse <- inverse[order(-inverse$evidence_weighted_score_inverse,
  inverse$feature_id), ]
stopifnot(nrow(direct) > 0L, nrow(inverse) > 0L)
write_tsv(direct, file.path(stage_dir,
  "full_direction_concordant_direct.tsv"))
write_tsv(inverse, file.path(stage_dir,
  "full_direction_concordant_inverse.tsv"))

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
  identical(manifest$current_state$fingerprint, expected_fingerprint))
metadata <- wb_meta(bundle, cells = trajectory$cell_id,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt", "nCount_RNA"))
stopifnot(identical(as.character(metadata$cell_id), trajectory$cell_id))
depth_z <- as.numeric(scale(log1p(metadata$nCount_RNA)))
depth_design <- cbind(intercept = 1, ns_pt = ns(trajectory$pt, df = 8),
  log_RNA_depth_z = depth_z)
depth_solver <- solve(crossprod(depth_design))

top_direct <- head(direct, 15L)
top_inverse <- head(inverse, 20L)
plot_features <- unique(c(top_direct$feature_id, top_inverse$feature_id))
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data",
  plot_features, trajectory$cell_id, missing = "NA"))
coefficients <- depth_solver %*% crossprod(depth_design, t(expression))
depth_beta <- as.numeric(coefficients[nrow(coefficients), ])
adjusted <- expression - tcrossprod(depth_beta, depth_z)
profile <- as.matrix((adjusted %*% window_design) %*% smoothing_design)
q25 <- q75 <- matrix(NA_real_, nrow = length(plot_features),
  ncol = nrow(windows))
for (j in seq_len(nrow(windows))) {
  q <- rowQuantiles(adjusted[, membership_split[[j]], drop = FALSE],
    probs = c(0.25, 0.75), drop = FALSE)
  q25[, j] <- q[, 1L]
  q75[, j] <- q[, 2L]
}
q25 <- as.matrix(q25 %*% smoothing_design)
q75 <- as.matrix(q75 %*% smoothing_design)
rownames(profile) <- rownames(q25) <- rownames(q75) <- plot_features
center <- rowMeans(profile)
spread <- rowSds(profile)
profile <- (profile - center) / spread
q25 <- (q25 - center) / spread
q75 <- (q75 - center) / spread

make_data <- function(top, direction) {
  score_column <- paste0("evidence_weighted_score_", direction)
  weakest_column <- paste0("weakest_segment_", direction)
  out <- do.call(rbind, lapply(seq_len(nrow(top)), function(i) {
    gene <- top$feature_id[[i]]
    data.frame(pt = windows$pt_median,
      panel = sprintf("%02d %s | score %.3f | weakest segment %.2f", i,
        gene, top[[score_column]][[i]], top[[weakest_column]][[i]]),
      mean = profile[gene, ], q25 = q25[gene, ], q75 = q75[gene, ],
      target = if (direction == "direct") target_z else -target_z)
  }))
  out$panel <- factor(out$panel, levels = unique(out$panel))
  out
}
make_plot <- function(top, direction) {
  color <- if (direction == "direct") "#C51B5A" else "#2166AC"
  ggplot(make_data(top, direction), aes(pt)) +
    geom_hline(yintercept = 0, color = "#999999", linetype = "dotted",
      linewidth = 0.25) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = color,
      alpha = 0.18, color = NA) +
    geom_line(aes(y = mean), color = color, linewidth = 0.8) +
    geom_line(aes(y = target), color = "#111111", linewidth = 0.55) +
    facet_wrap(~panel, ncol = 4) +
    coord_cartesian(ylim = c(-3.5, 3.5)) +
    theme_classic(base_size = 9) +
    theme(strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 7.2),
      panel.spacing = grid::unit(5, "mm"),
      plot.title = element_text(face = "bold", size = 15)) +
    labs(title = paste("Full-trajectory", direction,
        "up/down-concordant genes"),
      subtitle = paste("All 8 balanced trajectory segments agree in direction;",
        "magnitude is standardized; ribbon is adjusted-expression IQR"),
      x = "Bdev pseudotime", y = "Standardized profile")
}
direct_plot <- make_plot(top_direct, "direct")
inverse_plot <- make_plot(top_inverse, "inverse")
ggsave(file.path(stage_dir, "full_direction_concordant_direct.png"),
  direct_plot, width = 14, height = 10, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "full_direction_concordant_direct.pdf"),
  direct_plot, width = 14, height = 10, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "full_direction_concordant_inverse.png"),
  inverse_plot, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "full_direction_concordant_inverse.pdf"),
  inverse_plot, width = 14, height = 12, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  source_screen = screen_dir, direct_all_eight_count = nrow(direct),
  inverse_all_eight_count = nrow(inverse), magnitude_invariant = TRUE,
  definition = "normalized clipped slope similarity has the requested sign in all eight equal-interest trajectory segments",
  bundle_mutated = FALSE, capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"), pretty = TRUE,
  auto_unbox = TRUE)
writeLines(c("# Full-trajectory direction-concordant gene matches", "",
  paste("Direct genes concordant in all eight segments:", nrow(direct)),
  paste("Inverse genes concordant in all eight segments:", nrow(inverse)), "",
  paste("Top direct:", paste(top_direct$feature_id, collapse = ", ")), "",
  paste("Top inverse:", paste(top_inverse$feature_id, collapse = ", ")), "",
  "Magnitude and offset are standardized. These are direction-concordant candidates, not claims of exact curve reproduction or statistical significance."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed full-direction match report at", output_dir, "\n")
