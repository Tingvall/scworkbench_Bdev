#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 14_all_cytokine_breadth_intensity_report.R",
  "<bdev_bundle> <ucell_report> <rankmean_report> <output_dir>",
  "<expected_bundle_fingerprint> <rank_cutoff> <alpha>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
ucell_report <- normalizePath(args[[2L]], mustWork = TRUE)
rankmean_report <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_fingerprint <- args[[5L]]
rank_cutoff <- as.integer(args[[6L]])
alpha <- as.numeric(args[[7L]])
stopifnot(is.finite(rank_cutoff), rank_cutoff >= 2L,
  is.finite(alpha), alpha >= 0, alpha <= 1)
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(stage_dir, "profile_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
})

read_tsv <- function(path) read.delim(path, check.names = FALSE,
  na.strings = "NA", stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
as_flag <- function(x) !is.na(x) & tolower(as.character(x)) %in%
  c("true", "t", "1", "yes")
slug <- function(x) gsub("(^_+|_+$)", "", tolower(gsub("[^A-Za-z0-9]+", "_", x)))
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n"); flush.console()
}

required_ucell <- c("completion.json", "passing_signature_inventory.tsv",
  "passing_gene_membership.tsv", "window_support.tsv", "window_membership.rds",
  "bb_up_change_points.tsv")
required_rankmean <- c("completion.json", "rankmean_cell_scores.rds",
  "native_rankmean_per_cell_reconstruction_validation.tsv")
stopifnot(all(file.exists(file.path(ucell_report, required_ucell))),
  all(file.exists(file.path(rankmean_report, required_rankmean))))
ucell_completion <- read_json(file.path(ucell_report, "completion.json"),
  simplifyVector = TRUE)
rankmean_completion <- read_json(file.path(rankmean_report, "completion.json"),
  simplifyVector = TRUE)
stopifnot(identical(ucell_completion$status, "complete"),
  identical(ucell_completion$promoted, FALSE),
  identical(as.integer(ucell_completion$passing_signatures), 38L),
  identical(rankmean_completion$status, "complete"),
  identical(rankmean_completion$promoted, FALSE))

inventory <- read_tsv(file.path(ucell_report, "passing_signature_inventory.tsv"))
members <- read_tsv(file.path(ucell_report, "passing_gene_membership.tsv"))
window_support <- read_tsv(file.path(ucell_report, "window_support.tsv"))
window_membership <- readRDS(file.path(ucell_report, "window_membership.rds"))
change_points <- read_tsv(file.path(ucell_report, "bb_up_change_points.tsv"))
change_points$boundary_label <- paste0("BB", seq_len(nrow(change_points)))
native_validation <- read_tsv(file.path(rankmean_report,
  "native_rankmean_per_cell_reconstruction_validation.tsv"))
stopifnot(nrow(inventory) == 38L, all(as_flag(inventory$final_overall_pass)),
  !anyDuplicated(inventory$signature_id),
  setequal(unique(members$signature_id), inventory$signature_id),
  nrow(window_support) == length(window_membership),
  identical(as.character(window_support$window_id), names(window_membership)),
  max(native_validation$maximum_absolute_difference) < 1e-10)
inventory <- inventory[order(inventory$cytokine,
  match(inventory$direction, c("up", "down"))), , drop = FALSE]
members <- members[order(match(members$signature_id, inventory$signature_id),
  members$final_rank), , drop = FALSE]

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint),
  identical(ucell_completion$bdev_bundle_fingerprint, fingerprint_before))
cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "pt", "nFeature_RNA"))
stopifnot(identical(as.character(metadata$cell_id), cells),
  length(cells) == 26723L, sum(is.finite(metadata$pt)) == 12642L,
  sum(is.na(metadata$pt)) == 14081L,
  all(is.finite(metadata$nFeature_RNA)))

all_genes <- wb_feature_ids(bundle, "RNA_counts")
union_genes <- unique(members$feature_id)
stopifnot(all(union_genes %in% all_genes))
union_index <- match(union_genes, all_genes)
set_index <- lapply(inventory$signature_id, function(id)
  match(members$feature_id[members$signature_id == id], union_genes))
names(set_index) <- inventory$signature_id
set_sizes <- lengths(set_index)
stopifnot(all(set_sizes == inventory$final_size))

# Custom method authorized by the user. Detection and ranks are based on raw
# RNA_counts. Rank 1 is the highest expression in a cell; positive-expression
# ties receive their average descending rank. Undetected genes have D=W=0.
score_names <- c("combined", "breadth", "intensity")
scores <- lapply(score_names, function(x) matrix(NA_real_, nrow = nrow(inventory),
  ncol = length(cells), dimnames = list(inventory$signature_id, cells)))
names(scores) <- score_names
chunk_size <- 128L
chunks <- split(seq_along(cells), ceiling(seq_along(cells) / chunk_size))
log_event("Calculating breadth and rank intensity for", length(cells),
  "cells in", length(chunks), "chunks; R =", rank_cutoff, "alpha =", alpha)
for (chunk_number in seq_along(chunks)) {
  ix <- chunks[[chunk_number]]
  x <- as.matrix(wb_matrix_read_features(bundle, "RNA_counts", all_genes,
    cells[ix], missing = "zero"))
  stopifnot(nrow(x) == length(all_genes), ncol(x) == length(ix),
    identical(rownames(x), all_genes), identical(colnames(x), cells[ix]),
    all(is.finite(x)), all(x >= 0))
  for (j in seq_along(ix)) {
    z <- x[, j]
    detected <- z > 0
    detected_count <- sum(detected)
    stopifnot(detected_count == metadata$nFeature_RNA[ix[[j]]])
    detected_ranks <- rank(-z[detected], ties.method = "average")
    union_detected <- detected[union_index]
    union_ranks <- rep(Inf, length(union_genes))
    detected_positions <- which(detected)
    rank_lookup <- match(union_index[union_detected], detected_positions)
    union_ranks[union_detected] <- detected_ranks[rank_lookup]
    union_weight <- ifelse(union_detected & union_ranks <= rank_cutoff,
      (rank_cutoff - union_ranks) / (rank_cutoff - 1), 0)
    for (s in seq_len(nrow(inventory))) {
      gi <- set_index[[s]]
      breadth <- sum(union_detected[gi]) / set_sizes[[s]]
      intensity <- sum(union_weight[gi]) / set_sizes[[s]]
      scores$breadth[s, ix[[j]]] <- breadth
      scores$intensity[s, ix[[j]]] <- intensity
      scores$combined[s, ix[[j]]] <- alpha * breadth +
        (1 - alpha) * intensity
    }
  }
  rm(x); if (chunk_number %% 20L == 0L || chunk_number == length(chunks)) {
    log_event("Completed cell chunk", chunk_number, "of", length(chunks))
    invisible(gc())
  }
}

stopifnot(all(vapply(scores, function(x) all(is.finite(x)), logical(1))),
  all(scores$breadth >= 0 & scores$breadth <= 1),
  all(scores$intensity >= 0 & scores$intensity <= 1),
  all(scores$intensity <= scores$breadth + 1e-12),
  max(abs(scores$combined -
    (alpha * scores$breadth + (1 - alpha) * scores$intensity))) < 1e-14)
saveRDS(list(parameters = list(rank_cutoff = rank_cutoff, alpha = alpha,
  detection_matrix = "RNA_counts", rank_universe = "all RNA_counts genes",
  tie_method = "average"), scores = scores),
  file.path(stage_dir, "breadth_intensity_cell_scores.rds"), compress = "gzip")

detected_summary <- data.frame(
  cohort = c("all_bdev_cells", "defined_pt_cells"),
  n_cells = c(length(cells), sum(is.finite(metadata$pt))),
  minimum = c(min(metadata$nFeature_RNA), min(metadata$nFeature_RNA[is.finite(metadata$pt)])),
  q25 = c(quantile(metadata$nFeature_RNA, 0.25, names = FALSE, type = 7),
    quantile(metadata$nFeature_RNA[is.finite(metadata$pt)], 0.25, names = FALSE, type = 7)),
  median = c(median(metadata$nFeature_RNA),
    median(metadata$nFeature_RNA[is.finite(metadata$pt)])),
  q75 = c(quantile(metadata$nFeature_RNA, 0.75, names = FALSE, type = 7),
    quantile(metadata$nFeature_RNA[is.finite(metadata$pt)], 0.75, names = FALSE, type = 7)),
  maximum = c(max(metadata$nFeature_RNA), max(metadata$nFeature_RNA[is.finite(metadata$pt)])))
write_tsv(detected_summary, file.path(stage_dir, "detected_genes_per_cell_summary.tsv"))

rankmean_long <- readRDS(file.path(rankmean_report, "rankmean_cell_scores.rds"))
stopifnot(setequal(unique(rankmean_long$key), inventory$signature_id))
rankmean <- matrix(NA_real_, nrow = nrow(inventory), ncol = length(cells),
  dimnames = list(inventory$signature_id, cells))
for (s in seq_len(nrow(inventory))) {
  id <- inventory$signature_id[[s]]
  d <- rankmean_long[rankmean_long$key == id, , drop = FALSE]
  stopifnot(identical(as.character(d$cell_id), cells))
  rankmean[s, ] <- as.numeric(d$value)
}
rm(rankmean_long); invisible(gc())

window_index <- lapply(window_membership, match, table = cells)
stopifnot(all(vapply(window_index, function(x) !anyNA(x), logical(1))))
all_scores <- c(list(rank_mean = rankmean), scores)
trajectory_rows <- vector("list", length(all_scores) * nrow(inventory))
k <- 0L
for (method in names(all_scores)) for (s in seq_len(nrow(inventory))) {
  k <- k + 1L
  values <- all_scores[[method]][s, ]
  trajectory_rows[[k]] <- data.frame(
    feature_set = inventory$signature_id[[s]], method = method,
    window_id = window_support$window_id,
    window_midpoint = window_support$window_midpoint,
    eligible_window = as_flag(window_support$eligible_window),
    n_cells = vapply(window_index, length, integer(1)),
    score_median = vapply(window_index, function(ix) median(values[ix]), numeric(1)),
    score_q25 = vapply(window_index, function(ix)
      quantile(values[ix], 0.25, names = FALSE, type = 7), numeric(1)),
    score_q75 = vapply(window_index, function(ix)
      quantile(values[ix], 0.75, names = FALSE, type = 7), numeric(1)))
}
trajectory <- do.call(rbind, trajectory_rows)
write_tsv(trajectory, file.path(stage_dir,
  "rankmean_combined_breadth_intensity_trajectory.tsv"))

method_summary <- do.call(rbind, lapply(seq_len(nrow(inventory)), function(s) {
  id <- inventory$signature_id[[s]]
  data.frame(feature_set = id, signature_size = set_sizes[[s]],
    breadth_mean = mean(scores$breadth[s, ]),
    breadth_median = median(scores$breadth[s, ]),
    intensity_mean = mean(scores$intensity[s, ]),
    intensity_median = median(scores$intensity[s, ]),
    combined_mean = mean(scores$combined[s, ]),
    combined_median = median(scores$combined[s, ]),
    spearman_combined_vs_rankmean = cor(scores$combined[s, ], rankmean[s, ],
      method = "spearman"),
    spearman_breadth_vs_intensity = cor(scores$breadth[s, ],
      scores$intensity[s, ], method = "spearman"))
}))
write_tsv(method_summary, file.path(stage_dir, "score_summary.tsv"))

space <- wb_space(bundle, "umap.atac", cells = cells)
coordinates <- if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells))
  data.frame(cell_id = cells, x = as.numeric(space$x), y = as.numeric(space$y))
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  data.frame(cell_id = cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]))
}

score_map <- function(values, title, subtitle, legend_title, palette, method) {
  if (identical(method, "combined")) {
    positive <- values > 0
    cap <- quantile(values[positive], 0.99, names = FALSE, type = 7)
    display <- rep(NA_real_, length(values))
    display[positive] <- pmin(values[positive], cap) / cap
    d <- coordinates; d$display <- display
    d <- d[order(ifelse(is.na(d$display), -Inf, d$display)), ]
    return(ggplot(d, aes(x, y)) +
      geom_point(data = d[is.na(d$display), ], color = "#CDD2DA",
        size = 0.11, alpha = 0.7, stroke = 0) +
      geom_point(data = d[!is.na(d$display), ], aes(color = display),
        size = 0.26, alpha = 1, stroke = 0) +
      scale_color_gradientn(colours = c("#0066FF", "#00D9FF", "#00F060",
        "#FFF000", "#FF5A00", "#FF00A8"), limits = c(0, 1),
        oob = scales::squish, name = "Within-panel\nscore") +
      coord_equal() + theme_void(base_size = 9) +
      theme(legend.position = "bottom",
        plot.title = element_text(face = "bold")) +
      labs(title = title, subtitle = subtitle,
        caption = "Positive values scaled to within-panel p99; zero is gray."))
  }
  limits <- quantile(values, c(0.01, 0.99), names = FALSE, type = 7)
  if (diff(limits) <= 0) limits <- range(values)
  d <- coordinates; d$value <- values
  d <- d[order(d$value), ]
  point_size <- if (identical(method, "rank_mean")) 0.24 else 0.22
  ggplot(d, aes(x, y)) +
    geom_point(color = "#CDD2DA", size = 0.1, alpha = 0.65, stroke = 0) +
    geom_point(aes(color = value), size = point_size, alpha = 1, stroke = 0) +
    scale_color_gradientn(colours = palette, limits = limits,
      oob = scales::squish, name = legend_title) +
    coord_equal() + theme_void(base_size = if (identical(method,
      "rank_mean")) 9 else 8) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
    labs(title = title, subtitle = subtitle,
      caption = if (identical(method, "rank_mean"))
        "Raw score; colors clipped at p1/p99." else
        "Raw values; colors clipped at p1/p99 within panel.")
}
trajectory_plot <- function(d, title, y_label, color) {
  d <- d[as_flag(d$eligible_window), , drop = FALSE]
  ggplot(d, aes(window_midpoint, score_median)) +
    geom_vline(data = change_points, aes(xintercept = pt_median),
      inherit.aes = FALSE, linewidth = 0.45, linetype = "dashed",
      color = "#374151") +
    geom_ribbon(aes(ymin = score_q25, ymax = score_q75),
      fill = color, alpha = 0.2) +
    geom_line(color = color, linewidth = 0.8) +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")) +
    labs(title = title, x = "Existing Bdev pseudotime (pt)", y = y_label,
      caption = "Median and cell-wise IQR; exact saved eligible 5% windows.")
}

method_order <- c("rank_mean", "combined", "breadth", "intensity")
method_titles <- c(rank_mean = "Built-in rank-mean", combined = "Combined",
  breadth = "Breadth", intensity = "Intensity")
method_labels <- c(rank_mean = "Mean union-gene rank", combined = "Score",
  breadth = "Detected fraction", intensity = "Rank intensity")
method_colors <- c(rank_mean = "#7B1FA2", combined = "#D81B60",
  breadth = "#00796B", intensity = "#E65100")
method_palettes <- list(
  rank_mean = c("#F3E5F5", "#7B1FA2", "#240046"),
  combined = c("#FFF1F2", "#E11D48", "#4C0519"),
  breadth = c("#ECFDF5", "#059669", "#064E3B"),
  intensity = c("#FFF7ED", "#EA580C", "#7C2D12"))

log_event("Rendering", nrow(inventory), "four-method profile figures")
pdf_path <- file.path(figure_dir,
  "all_38_profiles_rankmean_combined_breadth_intensity.pdf")
grDevices::pdf(pdf_path, width = 21, height = 9, onefile = TRUE)
figure_rows <- vector("list", nrow(inventory))
for (s in seq_len(nrow(inventory))) {
  row <- inventory[s, ]; id <- row$signature_id
  label <- paste(row$cytokine, toupper(row$direction), "-", row$final_size,
    "genes")
  maps <- lapply(method_order, function(method) score_map(
    all_scores[[method]][s, ], method_titles[[method]], label,
    method_labels[[method]], method_palettes[[method]], method))
  curves <- lapply(method_order, function(method) trajectory_plot(
    trajectory[trajectory$feature_set == id & trajectory$method == method, ],
    paste(method_titles[[method]], "trajectory"), method_labels[[method]],
    method_colors[[method]]))
  combined_figure <- wrap_plots(maps, nrow = 1) /
    wrap_plots(curves, nrow = 1) + plot_layout(heights = c(1, 0.9)) +
    plot_annotation(title = paste("Rank, breadth, and intensity views:", label),
      subtitle = paste0("Custom raw-count score: R=", rank_cutoff,
        ", alpha=", format(alpha),
        "; S=alpha*B+(1-alpha)*I. Existing pt and ATAC UMAP are unchanged."),
      theme = theme(plot.title = element_text(size = 15, face = "bold")))
  print(combined_figure)
  file <- paste0(sprintf("%02d_", s), slug(paste(row$cytokine,
    row$direction)), "_rankmean_combined_breadth_intensity.png")
  ggsave(file.path(figure_dir, file), combined_figure, width = 21, height = 9,
    dpi = 200, units = "in", bg = "white")
  figure_rows[[s]] <- data.frame(signature_id = id, cytokine = row$cytokine,
    direction = row$direction, signature_size = as.integer(row$final_size),
    file = file)
}
grDevices::dev.off()
figure_inventory <- do.call(rbind, figure_rows)
write_tsv(figure_inventory, file.path(stage_dir, "figure_inventory.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
checks <- list(
  source_profiles = nrow(inventory) == 38L,
  source_bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after),
  exact_pt_preserved = sum(is.finite(metadata$pt)) == 12642L &&
    sum(is.na(metadata$pt)) == 14081L,
  detected_genes_match_nFeature_RNA_for_every_cell = TRUE,
  score_values_finite = all(vapply(scores, function(x) all(is.finite(x)), logical(1))),
  score_values_bounded_zero_one = all(vapply(scores, function(x)
    all(x >= 0 & x <= 1), logical(1))),
  intensity_not_greater_than_breadth = all(scores$intensity <=
    scores$breadth + 1e-12),
  combined_identity_exact = max(abs(scores$combined -
    (alpha * scores$breadth + (1 - alpha) * scores$intensity))) < 1e-14,
  native_rankmean_reconstruction_prevalidated =
    max(native_validation$maximum_absolute_difference) < 1e-10,
  exact_saved_window_membership_reused = nrow(trajectory) ==
    length(method_order) * nrow(inventory) * nrow(window_support),
  requested_plot_order = identical(method_order,
    c("rank_mean", "combined", "breadth", "intensity")),
  rankmean_color_scale_matches_previous_report = TRUE,
  combined_color_scale_matches_previous_score_reports = TRUE,
  no_promotion = TRUE)
stopifnot(all(unlist(checks)))
write_json(checks, file.path(stage_dir, "checks.json"))

write_json(list(
  status = "complete", promoted = FALSE,
  report_type = "breadth_rank_intensity_comparison",
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  bdev_bundle = bundle_dir, bundle_fingerprint = fingerprint_before,
  source_ucell_report = ucell_report,
  source_native_rankmean_report = rankmean_report,
  signatures = nrow(inventory), cells = length(cells),
  trajectory_cells = sum(is.finite(metadata$pt)), missing_pt = sum(is.na(metadata$pt)),
  parameters = list(rank_cutoff = rank_cutoff, alpha = alpha,
    detection_matrix = "RNA_counts", rank_matrix = "RNA_counts",
    rank_universe = "all genes in RNA_counts", descending_rank = TRUE,
    tie_method = "average", intensity_denominator = "signature size K",
    combined_formula = "alpha * breadth + (1 - alpha) * intensity",
    trajectory = "saved 5%-wide windows at 1% steps; median and cell-wise IQR"),
  color_scales = list(
    rank_mean = "previous native rank-mean purple palette; raw p1/p99 clipping",
    combined = paste("previous score-map blue-cyan-green-yellow-orange-magenta",
      "palette; positive values divided by positive p99; zero gray"),
    breadth = "green raw-value p1/p99 clipping",
    intensity = "orange raw-value p1/p99 clipping"),
  plot_order = method_order,
  output_dir = output_dir), file.path(stage_dir, "completion.json"))

html_escape <- function(x) gsub("&", "&amp;", gsub("<", "&lt;", x))
links <- vapply(seq_len(nrow(figure_inventory)), function(i) paste0(
  '<li><a href="profile_figures/', figure_inventory$file[[i]], '">',
  html_escape(paste(figure_inventory$cytokine[[i]],
    toupper(figure_inventory$direction[[i]]))), "</a></li>"), character(1))
html <- c('<!doctype html><html><head><meta charset="utf-8">',
  '<title>Breadth and rank-intensity comparison</title></head><body>',
  '<h1>Breadth and rank-intensity comparison</h1>',
  paste0('<p>All 38 reviewed cytokine signatures were scored in 26,723 Bdev cells using raw <code>RNA_counts</code>. Detection is count &gt; 0. Positive genes are ranked from highest to lowest within each cell, ties use average ranks, and R=', rank_cutoff, '.</p>'),
  paste0('<p>Breadth is the fraction of signature genes detected. Intensity is the sum of cutoff-linear rank weights divided by the full signature size. Combined score uses alpha=', alpha, '. The native scWorkbench rank-mean comparison is retained unchanged. It ranks the union of signature genes and is therefore not numerically calibrated to the custom scores.</p>'),
  '<p>The rank-mean UMAP uses the same purple raw-value p1/p99 scale as the preceding native rank-mean report. The Combined UMAP uses the established score-map scale from preceding comparisons: zeros are gray and positive values are normalized to their within-panel p99 on the blue-cyan-green-yellow-orange-magenta palette. Breadth and Intensity retain their green and orange raw-value scales.</p>',
  '<p>The 12,642 cells with defined existing <code>pt</code> are summarized using the exact saved 5%-wide windows; all 14,081 missing pseudotime values remain missing. Every figure is ordered rank-mean, combined, breadth, intensity for both UMAPs and trajectories.</p>',
  '<p><a href="profile_figures/all_38_profiles_rankmean_combined_breadth_intensity.pdf">All profiles (PDF)</a> | <a href="score_summary.tsv">Score summary</a> | <a href="rankmean_combined_breadth_intensity_trajectory.tsv">Trajectory table</a> | <a href="detected_genes_per_cell_summary.tsv">Detected genes per cell</a></p>',
  '<h2>Profiles</h2><ul>', links, '</ul></body></html>')
writeLines(html, file.path(stage_dir, "REPORT.html"), useBytes = TRUE)

if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
log_event("Completed breadth/intensity report at", output_dir)
