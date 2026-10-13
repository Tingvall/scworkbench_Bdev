#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 28_preliminary_epigenetic_pattern_candidates.R",
  "<bdev_bundle> <window_report> <smoothed_report> <turning_report>",
  "<turn_pattern_report> <output_dir> <expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
smoothed_report <- normalizePath(args[[3L]], mustWork = TRUE)
turning_report <- normalizePath(args[[4L]], mustWork = TRUE)
turn_pattern_report <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
expected_fingerprint <- args[[7L]]
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(Matrix)
  library(ggplot2)
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
turns <- read_tsv(file.path(turning_report, "selected_turning_regions.tsv"))
turns <- turns[match(paste0("T", 1:4), turns$region_id), ]
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
score_connection <- gzfile(file.path(turn_pattern_report,
  "modality_turn_pattern_scores.tsv.gz"), open = "rt")
scores <- read_tsv(score_connection)
close(score_connection)
stopifnot(nrow(trajectory) == 12642L, nrow(windows) == 2524L,
  !anyDuplicated(trajectory$cell_id), !anyDuplicated(windows$window_id),
  identical(turns$region_id, paste0("T", 1:4)),
  nrow(scores) == 11745L, !anyDuplicated(scores$feature_id))

direct_genes <- c("Baz2b", "Smarca4", "Ncoa2", "Aff3", "Kansl1l",
  "Nipbl", "Arid1b", "Aff1", "Med13", "Bcl7a")
inverse_genes <- c("Set", "Nap1l1", "Cbx5", "Cbx3", "H2afz", "Hmgn2",
  "Smc2")
cohesin_prc2_genes <- c("Nipbl", "Mau2", "Wapl", "Ctcf", "Rad21",
  "Smc1a", "Smc3", "Stag1", "Stag2", "Pds5a", "Pds5b", "Ezh2",
  "Eed", "Suz12")
candidate_definition <- rbind(
  data.frame(figure_group = "Direct chromatin candidates",
    direction = "direct", feature_id = direct_genes),
  data.frame(figure_group = "Inverse chromatin candidates",
    direction = "inverse", feature_id = inverse_genes),
  data.frame(figure_group = "Cohesin and PRC2 audit",
    direction = "direct", feature_id = cohesin_prc2_genes))
stopifnot(all(candidate_definition$feature_id %in% scores$feature_id))

candidate_rows <- lapply(seq_len(nrow(candidate_definition)), function(i) {
  definition <- candidate_definition[i, ]
  score <- scores[match(definition$feature_id, scores$feature_id), ]
  direction <- definition$direction
  data.frame(figure_group = definition$figure_group, direction = direction,
    feature_id = definition$feature_id,
    role_order = i,
    expression_dynamic_range_p95_p05 = score$expression_dynamic_range_p95_p05,
    weighted_sign_concordance = score[[paste0(direction,
      "_weighted_sign_concordance")]],
    weighted_derivative_cosine = score[[paste0(direction,
      "_weighted_derivative_cosine")]],
    matched_flanks_of_8 = score[[paste0(direction, "_matched_flanks_of_8")]],
    complete_turns_of_4 = score[[paste0(direction,
      "_complete_turns_of_4")]],
    composite_score = score[[paste0(direction, "_composite_score")]],
    selected = score[[paste0(direction, "_selected")]],
    stringsAsFactors = FALSE)
})
candidates <- do.call(rbind, candidate_rows)
candidates$status <- ifelse(candidates$selected, "PASS", "NOT PASS")
write_tsv(candidates, file.path(stage_dir,
  "epigenetic_pattern_candidate_status.tsv"))

cell_index <- match(membership$cell_id, trajectory$cell_id)
window_index <- match(membership$window_id, windows$window_id)
stopifnot(!anyNA(cell_index), !anyNA(window_index))
window_design <- sparseMatrix(i = cell_index, j = window_index,
  x = rep(1 / 30, nrow(membership)),
  dims = c(nrow(trajectory), nrow(windows)))
membership_split <- split(cell_index, factor(membership$window_id,
  levels = windows$window_id))

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

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))
plot_features <- unique(candidate_definition$feature_id)
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data",
  plot_features, trajectory$cell_id, missing = "NA"))
stopifnot(identical(rownames(expression), plot_features),
  identical(colnames(expression), trajectory$cell_id), all(is.finite(expression)))

raw_mean <- as.matrix(expression %*% window_design)
raw_q25 <- vapply(membership_split, function(ix)
  apply(expression[, ix, drop = FALSE], 1L, quantile, probs = 0.25,
    names = FALSE, type = 7), numeric(length(plot_features)))
raw_q75 <- vapply(membership_split, function(ix)
  apply(expression[, ix, drop = FALSE], 1L, quantile, probs = 0.75,
    names = FALSE, type = 7), numeric(length(plot_features)))
rownames(raw_q25) <- rownames(raw_q75) <- plot_features
smoothed_mean <- as.matrix(raw_mean %*% smoothing_design)
smoothed_q25 <- as.matrix(raw_q25 %*% smoothing_design)
smoothed_q75 <- as.matrix(raw_q75 %*% smoothing_design)
rownames(smoothed_mean) <- rownames(smoothed_q25) <-
  rownames(smoothed_q75) <- plot_features
profile_center <- rowMeans(smoothed_mean)
profile_sd <- sqrt(rowSums((smoothed_mean - profile_center)^2) /
  (ncol(smoothed_mean) - 1L))
stopifnot(all(is.finite(profile_sd)), all(profile_sd > 0))
mean_z <- (smoothed_mean - profile_center) / profile_sd
q25_z <- (smoothed_q25 - profile_center) / profile_sd
q75_z <- (smoothed_q75 - profile_center) / profile_sd
direct_target_z <- as.numeric(scale(log2(windows$smoothed_local_slope)))
inverse_target_z <- -direct_target_z

make_plot_data <- function(group_name) {
  subset <- candidates[candidates$figure_group == group_name, ]
  rows <- lapply(seq_len(nrow(subset)), function(i) {
    row <- subset[i, ]
    direction_target <- if (row$direction == "direct") direct_target_z else
      inverse_target_z
    data.frame(pt = windows$pt_median,
      gene_panel = sprintf("%s  %s | range %.2f | flanks %d/8 | turns %d/4",
        row$status, row$feature_id, row$expression_dynamic_range_p95_p05,
        row$matched_flanks_of_8, row$complete_turns_of_4),
      status = row$status, mean_z = mean_z[row$feature_id, ],
      q25_z = q25_z[row$feature_id, ], q75_z = q75_z[row$feature_id, ],
      target_z = direction_target, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$gene_panel <- factor(out$gene_panel, levels = unique(out$gene_panel))
  out$status <- factor(out$status, levels = c("PASS", "NOT PASS"))
  out
}

make_figure <- function(group_name, title, target_label, ncol) {
  data <- make_plot_data(group_name)
  ggplot(data, aes(pt)) +
    geom_hline(yintercept = 0, color = "#999999", linetype = "dotted",
      linewidth = 0.25) +
    geom_ribbon(aes(ymin = q25_z, ymax = q75_z, fill = status),
      alpha = 0.20, color = NA) +
    geom_line(aes(y = mean_z, color = status), linewidth = 0.85) +
    geom_line(aes(y = target_z), color = "#111111", linewidth = 0.55) +
    geom_vline(data = turns, aes(xintercept = pt), inherit.aes = FALSE,
      color = "#555555", linetype = "dashed", linewidth = 0.22) +
    facet_wrap(~ gene_panel, ncol = ncol) +
    coord_cartesian(ylim = c(-3.5, 3.5)) +
    scale_color_manual(values = c("PASS" = "#008837",
      "NOT PASS" = "#8C8C8C"), drop = FALSE) +
    scale_fill_manual(values = c("PASS" = "#00A64F",
      "NOT PASS" = "#AFAFAF"), drop = FALSE) +
    theme_classic(base_size = 9) +
    theme(strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 8.1),
      panel.spacing = grid::unit(5, "mm"), legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 15)) +
    labs(title = title,
      subtitle = paste0(
        "Colored line/ribbon: standardized gene mean and within-window IQR; ",
        "black: ", target_label, "; dashed: T1-T4"),
      x = "Bdev pseudotime", y = "Standardized profile", color = NULL,
      fill = NULL)
}

direct_plot <- make_figure("Direct chromatin candidates",
  "Chromatin regulators ranked toward the direct RNA/ATAC pattern",
  "direct RNA/ATAC target", 2L)
inverse_plot <- make_figure("Inverse chromatin candidates",
  "Chromatin regulators ranked toward the inverse ATAC/RNA pattern",
  "inverse ATAC/RNA target", 2L)
cohesin_prc2_plot <- make_figure("Cohesin and PRC2 audit",
  "Focused audit of cohesin-associated and PRC2 genes",
  "direct RNA/ATAC target", 2L)

ggsave(file.path(stage_dir, "direct_chromatin_candidates_pass_fail.png"),
  direct_plot, width = 11, height = 13, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "direct_chromatin_candidates_pass_fail.pdf"),
  direct_plot, width = 11, height = 13, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "inverse_chromatin_candidates_pass_fail.png"),
  inverse_plot, width = 11, height = 10, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "inverse_chromatin_candidates_pass_fail.pdf"),
  inverse_plot, width = 11, height = 10, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "cohesin_prc2_candidates_pass_fail.png"),
  cohesin_prc2_plot, width = 11, height = 17, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "cohesin_prc2_candidates_pass_fail.pdf"),
  cohesin_prc2_plot, width = 11, height = 17, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  source_turn_pattern_report = turn_pattern_report,
  expression_matrix = "SCT_data", trajectory_cells = nrow(trajectory),
  windows = nrow(windows), window_size = 30L, stride = 5L,
  smoothing_fraction = smoothing_fraction,
  interval = "within-window cell-level IQR, then fixed 8%-window smoothing",
  direct_candidates = as.list(direct_genes),
  inverse_candidates = as.list(inverse_genes),
  cohesin_prc2_audit = as.list(cohesin_prc2_genes),
  passing_candidates = as.list(unique(candidates$feature_id[candidates$selected])),
  bundle_mutated = FALSE, capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Epigenetic pattern candidates: pass/fail plots", "",
  paste("Figures show standardized smoothed expression with smoothed",
    "within-window cell-level IQR against the relevant direct RNA/ATAC or",
    "inverse ATAC/RNA target. Every facet states PASS or NOT PASS and records",
    "the expression range, matched event flanks, and complete turns."), "",
  paste("Baz2b is the only strict passing epigenetic candidate. Smarca4 is a",
    "direct near miss with five of eight flanks and two complete turns.",
    "None of the inverse candidates or focused cohesin/PRC2 audit genes pass."),
  "The bundle was opened read-only."), file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed epigenetic pass/fail plots at", output_dir, "\n")
