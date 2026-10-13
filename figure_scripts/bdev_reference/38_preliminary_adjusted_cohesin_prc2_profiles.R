#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 38_preliminary_adjusted_cohesin_prc2_profiles.R",
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

cohesin <- c("Wapl", "Smc1a", "Smc3", "Rad21", "Pds5a", "Ctcf",
  "Stag2", "Pds5b", "Stag1", "Mau2", "Nipbl")
prc2 <- c("Ezh2", "Eed", "Suz12")
genes <- unique(c("Ebf1", cohesin, prc2))

trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
windows <- read_tsv(file.path(window_report, "window_distance_and_slope.tsv"))
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
  all(segment_columns %in% names(scores)), nrow(windows) == nrow(ratio))
audit <- scores[match(genes, scores$feature_id), ]
segment_matrix <- as.matrix(audit[, segment_columns])
audit$direct_segments <- rowSums(segment_matrix > 0)
audit$inverse_segments <- rowSums(segment_matrix < 0)
audit$preferred_direction <- ifelse(
  audit$evidence_weighted_score_direct >=
    audit$evidence_weighted_score_inverse, "direct", "inverse")
write_tsv(audit, file.path(stage_dir,
  "adjusted_cohesin_prc2_ebf1_direction_scores.tsv"))

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

make_data <- function(selected) {
  selected <- unique(c("Ebf1", selected))
  rows <- list(data.frame(pt = windows$pt_median,
    panel = "00 Corrected RNA/ATAC ratio", mean = target_z,
    q25 = NA_real_, q75 = NA_real_, direction = "target"))
  for (i in seq_along(selected)) {
    gene <- selected[[i]]
    score <- audit[audit$feature_id == gene, ]
    rows[[length(rows) + 1L]] <- data.frame(pt = windows$pt_median,
      panel = sprintf("%02d %s | direct %d/8 | inverse %d/8 | depth beta %.3f",
        i, gene, score$direct_segments, score$inverse_segments,
        score$depth_beta), mean = profile_z[gene, ], q25 = q25_z[gene, ],
      q75 = q75_z[gene, ], direction = score$preferred_direction)
  }
  out <- do.call(rbind, rows)
  out$panel <- factor(out$panel, levels = unique(out$panel))
  out
}
make_plot <- function(selected, title) {
  plot_data <- make_data(selected)
  ggplot(plot_data, aes(pt, mean, color = direction, fill = direction)) +
    geom_hline(yintercept = 0, color = "#999999", linetype = "dotted",
      linewidth = 0.25) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.16,
      color = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~panel, ncol = 1, scales = "fixed") +
    scale_color_manual(values = c(target = "#111111", direct = "#C51B5A",
      inverse = "#2166AC"), guide = "none") +
    scale_fill_manual(values = c(target = "#111111", direct = "#C51B5A",
      inverse = "#2166AC"), guide = "none") +
    coord_cartesian(ylim = c(-3.5, 3.5)) +
    theme_classic(base_size = 9) +
    theme(strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0, size = 8),
      panel.spacing = grid::unit(2.5, "mm"),
      plot.title = element_text(face = "bold", size = 15)) +
    labs(title = title,
      subtitle = paste("RNA-depth component removed conditional on ns(pt, df=8);",
        "identical 3% smoothing; ribbon is adjusted-expression IQR"),
      x = "Bdev pseudotime", y = "Standardized adjusted profile")
}
cohesin_plot <- make_plot(cohesin,
  "Adjusted Ebf1 and cohesin-associated expression")
prc2_plot <- make_plot(prc2,
  "Adjusted Ebf1 and PRC2-associated expression")
ggsave(file.path(stage_dir, "adjusted_ebf1_cohesin_profiles.png"),
  cohesin_plot, width = 12, height = 18, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "adjusted_ebf1_cohesin_profiles.pdf"),
  cohesin_plot, width = 12, height = 18, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "adjusted_ebf1_prc2_profiles.png"),
  prc2_plot, width = 12, height = 9, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "adjusted_ebf1_prc2_profiles.pdf"),
  prc2_plot, width = 12, height = 9, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  genes = as.list(genes), depth_adjustment =
    "SCT_data minus conditional log1p(nCount_RNA) coefficient after ns(pt,df=8)",
  smoothing_fraction = 0.03, bundle_mutated = FALSE,
  capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"), pretty = TRUE,
  auto_unbox = TRUE)
writeLines(c("# Adjusted cohesin, PRC2, and Ebf1 profiles", "",
  "The corrected ratio and standardized depth-adjusted expression profiles use identical 3% trajectory smoothing.",
  "Ebf1 is included in both figures. Cohesin and PRC2 membership follows the prior focused audit."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed adjusted cohesin/PRC2 profiles at", output_dir, "\n")
