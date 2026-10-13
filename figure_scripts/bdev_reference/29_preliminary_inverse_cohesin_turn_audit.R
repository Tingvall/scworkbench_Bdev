#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 29_preliminary_inverse_cohesin_turn_audit.R",
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
event_connection <- gzfile(file.path(turn_pattern_report,
  "gene_turn_event_matches.tsv.gz"), open = "rt")
events <- read_tsv(event_connection)
close(event_connection)
event_definition <- read_tsv(file.path(turn_pattern_report,
  "turn_event_definition.tsv"))
stopifnot(nrow(trajectory) == 12642L, nrow(windows) == 2524L,
  identical(turns$region_id, paste0("T", 1:4)),
  nrow(scores) == 11745L, !anyDuplicated(scores$feature_id))

cohesin_genes <- c("Wapl", "Smc1a", "Smc3", "Rad21", "Pds5a", "Ctcf",
  "Stag2", "Pds5b", "Stag1", "Mau2", "Nipbl")
stopifnot(all(cohesin_genes %in% scores$feature_id))
cohesin_scores <- scores[match(cohesin_genes, scores$feature_id), ]
cohesin_scores <- cohesin_scores[order(-cohesin_scores$inverse_composite_score), ]
cohesin_scores$inverse_rank_within_set <- seq_len(nrow(cohesin_scores))
cohesin_scores$inverse_interpretation <- ifelse(
  cohesin_scores$inverse_selected, "STRICT PASS",
  ifelse(cohesin_scores$inverse_weighted_sign_concordance >= 0.65 &
    cohesin_scores$inverse_weighted_derivative_cosine >= 0.20,
    "INVERSE TREND; NOT STRICT PASS", "NOT INVERSE PATTERN"))
write_tsv(cohesin_scores, file.path(stage_dir,
  "inverse_cohesin_candidate_scores.tsv"))

cohesin_events <- events[events$feature_id %in% cohesin_genes, ]
cohesin_events$feature_id <- factor(cohesin_events$feature_id,
  levels = rev(cohesin_scores$feature_id))
cohesin_events$region_id <- factor(cohesin_events$region_id,
  levels = paste0("T", 1:4))
cohesin_events <- cohesin_events[order(cohesin_events$feature_id,
  cohesin_events$region_id), ]
write_tsv(cohesin_events, file.path(stage_dir,
  "inverse_cohesin_turn_flank_calls.tsv"))

# Build an explicit 12-column matrix: two flank calls and one complete-turn
# call per T1-T4. Flank calls already include the 5%-of-gene-range local
# amplitude requirement from the source screen.
heatmap_rows <- do.call(rbind, lapply(seq_len(nrow(cohesin_events)), function(i) {
  row <- cohesin_events[i, ]
  rbind(
    data.frame(feature_id = as.character(row$feature_id),
      region_id = as.character(row$region_id), component = "pre flank",
      pass = row$inverse_pre_match),
    data.frame(feature_id = as.character(row$feature_id),
      region_id = as.character(row$region_id), component = "post flank",
      pass = row$inverse_post_match),
    data.frame(feature_id = as.character(row$feature_id),
      region_id = as.character(row$region_id), component = "complete turn",
      pass = row$inverse_complete_turn))
}))
heatmap_rows$column <- paste(heatmap_rows$region_id,
  ifelse(heatmap_rows$component == "pre flank", "pre",
    ifelse(heatmap_rows$component == "post flank", "post", "turn")))
column_levels <- unlist(lapply(paste0("T", 1:4), function(x)
  paste(x, c("pre", "post", "turn"))))
heatmap_rows$column <- factor(heatmap_rows$column, levels = column_levels)
heatmap_rows$feature_id <- factor(heatmap_rows$feature_id,
  levels = rev(cohesin_scores$feature_id))
heatmap_rows$status <- factor(ifelse(heatmap_rows$pass, "PASS", "FAIL"),
  levels = c("PASS", "FAIL"))

heatmap_plot <- ggplot(heatmap_rows, aes(column, feature_id, fill = status)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(pass, "PASS", "FAIL")), color = "white",
    size = 3, fontface = "bold") +
  geom_vline(xintercept = c(3.5, 6.5, 9.5), color = "#222222",
    linewidth = 0.8) +
  scale_fill_manual(values = c("PASS" = "#1B9E77", "FAIL" = "#D73027"),
    drop = FALSE) +
  scale_x_discrete(labels = gsub(" ", "\n", column_levels),
    position = "top") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
    axis.title = element_blank(), axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15)) +
  labs(title = "Inverse ATAC/RNA flank and complete-turn calls for cohesin genes",
    subtitle = paste(
      "A whole turn passes only when both its pre- and post-flanks pass;",
      "flank calls include the local-amplitude requirement"), fill = NULL)
ggsave(file.path(stage_dir, "inverse_cohesin_turn_flank_pass_fail.png"),
  heatmap_plot, width = 13, height = 7.5, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "inverse_cohesin_turn_flank_pass_fail.pdf"),
  heatmap_plot, width = 13, height = 7.5, device = cairo_pdf, bg = "white")

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
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data",
  cohesin_scores$feature_id, trajectory$cell_id, missing = "NA"))
stopifnot(identical(rownames(expression), cohesin_scores$feature_id),
  all(is.finite(expression)))
raw_mean <- as.matrix(expression %*% window_design)
raw_q25 <- vapply(membership_split, function(ix)
  apply(expression[, ix, drop = FALSE], 1L, quantile, probs = 0.25,
    names = FALSE, type = 7), numeric(nrow(expression)))
raw_q75 <- vapply(membership_split, function(ix)
  apply(expression[, ix, drop = FALSE], 1L, quantile, probs = 0.75,
    names = FALSE, type = 7), numeric(nrow(expression)))
rownames(raw_q25) <- rownames(raw_q75) <- cohesin_scores$feature_id
smoothed_mean <- as.matrix(raw_mean %*% smoothing_design)
smoothed_q25 <- as.matrix(raw_q25 %*% smoothing_design)
smoothed_q75 <- as.matrix(raw_q75 %*% smoothing_design)
rownames(smoothed_mean) <- rownames(smoothed_q25) <-
  rownames(smoothed_q75) <- cohesin_scores$feature_id
profile_center <- rowMeans(smoothed_mean)
profile_sd <- sqrt(rowSums((smoothed_mean - profile_center)^2) /
  (ncol(smoothed_mean) - 1L))
mean_z <- (smoothed_mean - profile_center) / profile_sd
q25_z <- (smoothed_q25 - profile_center) / profile_sd
q75_z <- (smoothed_q75 - profile_center) / profile_sd
inverse_target_z <- -as.numeric(scale(log2(windows$smoothed_local_slope)))

profile_rows <- lapply(seq_len(nrow(cohesin_scores)), function(i) {
  row <- cohesin_scores[i, ]
  gene <- row$feature_id
  data.frame(pt = windows$pt_median,
    gene_panel = sprintf(
      "%02d  %s | range %.2f | inverse flanks %d/8 | turns %d/4",
      i, gene, row$expression_dynamic_range_p95_p05,
      row$inverse_matched_flanks_of_8, row$inverse_complete_turns_of_4),
    feature_id = gene, mean_z = mean_z[gene, ], q25_z = q25_z[gene, ],
    q75_z = q75_z[gene, ], target_z = inverse_target_z,
    emphasis = ifelse(gene == "Wapl", "Wapl inverse near-hit", "Other cohesin"),
    stringsAsFactors = FALSE)
})
profile_data <- do.call(rbind, profile_rows)
profile_data$gene_panel <- factor(profile_data$gene_panel,
  levels = unique(profile_data$gene_panel))

# Position pre/post markers halfway along the corresponding 100-window flank.
event_marker_rows <- lapply(seq_len(nrow(cohesin_events)), function(i) {
  row <- cohesin_events[i, ]
  definition <- event_definition[event_definition$region_id ==
    as.character(row$region_id), ]
  score_row <- cohesin_scores[cohesin_scores$feature_id ==
    as.character(row$feature_id), ]
  panel <- sprintf("%02d  %s | range %.2f | inverse flanks %d/8 | turns %d/4",
    score_row$inverse_rank_within_set, score_row$feature_id,
    score_row$expression_dynamic_range_p95_p05,
    score_row$inverse_matched_flanks_of_8,
    score_row$inverse_complete_turns_of_4)
  rbind(
    data.frame(gene_panel = panel,
      pt = mean(windows$pt_median[c(definition$left_window_index,
        definition$window_index)]), label = paste0(row$region_id, "pre"),
      status = ifelse(row$inverse_pre_match, "PASS", "FAIL")),
    data.frame(gene_panel = panel,
      pt = mean(windows$pt_median[c(definition$window_index,
        definition$right_window_index)]), label = paste0(row$region_id, "post"),
      status = ifelse(row$inverse_post_match, "PASS", "FAIL")))
})
event_markers <- do.call(rbind, event_marker_rows)
event_markers$gene_panel <- factor(event_markers$gene_panel,
  levels = levels(profile_data$gene_panel))
event_markers$status <- factor(event_markers$status, levels = c("PASS", "FAIL"))

trajectory_plot <- ggplot(profile_data, aes(pt)) +
  geom_hline(yintercept = 0, color = "#999999", linetype = "dotted",
    linewidth = 0.25) +
  geom_ribbon(aes(ymin = q25_z, ymax = q75_z), fill = "#9ECAE1",
    alpha = 0.22, color = NA) +
  geom_line(data = profile_data[profile_data$feature_id != "Wapl", ],
    aes(y = mean_z), color = "#737373", linewidth = 0.75) +
  geom_line(data = profile_data[profile_data$feature_id == "Wapl", ],
    aes(y = mean_z), color = "#2166AC", linewidth = 1.0) +
  geom_line(aes(y = target_z), color = "#111111", linewidth = 0.55) +
  geom_vline(data = turns, aes(xintercept = pt), inherit.aes = FALSE,
    color = "#555555", linetype = "dashed", linewidth = 0.22) +
  geom_point(data = event_markers, aes(x = pt, y = 3.0, fill = status),
    inherit.aes = FALSE, shape = 21, size = 3.1, color = "white",
    stroke = 0.5) +
  facet_wrap(~ gene_panel, ncol = 2L) +
  coord_cartesian(ylim = c(-3.5, 3.5)) +
  scale_fill_manual(values = c("PASS" = "#1B9E77", "FAIL" = "#D73027"),
    drop = FALSE) +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8.0),
    panel.spacing = grid::unit(5, "mm"), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15)) +
  labs(title = "Cohesin-associated genes against the inverse ATAC/RNA pattern",
    subtitle = paste(
      "Blue: Wapl; gray: other genes; black: inverse target;",
      "top markers give each T1-T4 pre/post flank call"),
    x = "Bdev pseudotime", y = "Standardized profile", fill = "Flank call")
ggsave(file.path(stage_dir, "inverse_cohesin_trajectories_with_flank_calls.png"),
  trajectory_plot, width = 12, height = 17, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "inverse_cohesin_trajectories_with_flank_calls.pdf"),
  trajectory_plot, width = 12, height = 17, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  source_turn_pattern_report = turn_pattern_report,
  direction = "inverse ATAC/RNA", expression_matrix = "SCT_data",
  trajectory_cells = nrow(trajectory), windows = nrow(windows),
  window_size = 30L, stride = 5L, smoothing_fraction = smoothing_fraction,
  flank_call = paste(
    "source inverse direction match plus local change at least 5%",
    "of the gene p95-p05 range"),
  strict_thresholds = list(dynamic_range = 0.50,
    sign_concordance = 0.65, derivative_cosine = 0.20,
    matched_flanks = 6L, complete_turns = 2L),
  highlighted_near_hit = "Wapl", bundle_mutated = FALSE,
  capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Inverse cohesin turn/flank audit", "",
  paste("The trajectory figure compares eleven cohesin-associated genes with",
    "the inverse ATAC/RNA target and marks every T1-T4 pre/post flank as pass",
    "or fail. The companion heatmap gives all eight flank calls plus the four",
    "complete-turn calls."), "",
  paste("Wapl is the strongest event-level inverse near-hit: T2 and T3 are",
    "complete passes and the T4 pre-flank also passes, giving five of eight",
    "flanks and two of four complete turns. It misses strict selection because",
    "its range is 0.461 rather than 0.50 and it has five rather than six",
    "matching flanks."), "The bundle was opened read-only."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed inverse cohesin turn audit at", output_dir, "\n")
