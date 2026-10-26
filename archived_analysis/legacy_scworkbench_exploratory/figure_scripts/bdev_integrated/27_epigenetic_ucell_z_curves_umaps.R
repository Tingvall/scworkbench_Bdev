#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "Usage: 27_epigenetic_ucell_z_curves_umaps.R",
  "<bundle_dir> <memberships.tsv.gz> <selection.tsv>",
  "<ucell_rankings.rds> <ucell_report> <z_report> <bb.tsv>",
  "<output_dir> <expected_sets>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
membership_path <- normalizePath(args[[2L]], mustWork = TRUE)
selection_path <- normalizePath(args[[3L]], mustWork = TRUE)
ranking_path <- normalizePath(args[[4L]], mustWork = TRUE)
ucell_root <- normalizePath(args[[5L]], mustWork = TRUE)
z_root <- normalizePath(args[[6L]], mustWork = TRUE)
bb_path <- normalizePath(args[[7L]], mustWork = TRUE)
output_dir <- normalizePath(args[[8L]], mustWork = FALSE)
expected_sets <- as.integer(args[[9L]])
stopifnot(is.finite(expected_sets), expected_sets > 0L)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
  library(scales)
})

read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
slugify <- function(x) gsub("(^_+|_+$)", "",
  tolower(gsub("[^A-Za-z0-9]+", "_", x)))
interval_mask <- function(pt, regions) {
  out <- rep(FALSE, length(pt))
  if (nrow(regions)) for (i in seq_len(nrow(regions)))
    out <- out | (is.finite(pt) & pt >= regions$region_left[[i]] &
      pt <= regions$region_right[[i]])
  out
}
safe_limits <- function(x) {
  lim <- as.numeric(quantile(x[is.finite(x)], c(.01, .99), names = FALSE))
  if (!all(is.finite(lim)) || diff(lim) <= .Machine$double.eps) {
    mid <- mean(x[is.finite(x)])
    lim <- mid + c(-.5, .5)
  }
  lim
}

selection <- fread(selection_path)
stopifnot(nrow(selection) == expected_sets,
  !anyDuplicated(selection$feature_set),
  setequal(selection$analysis_group, c("bp", "mf", "cc")))
setorder(selection, analysis_group, rank_within_group)
members_all <- read_tsv_gz(membership_path)
members <- members_all[feature_set %chin% selection$feature_set]
stopifnot(setequal(unique(members$feature_set), selection$feature_set))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
cells <- as.character(wb_cells(bundle, "all"))
metadata <- as.data.table(wb_meta(bundle, cells = cells))
umap <- as.data.table(wb_space(bundle, "umap.atac", cells = cells))
stopifnot(length(cells) == 26723L,
  identical(as.character(metadata$cell_id), cells),
  identical(as.character(umap$cell_id), cells),
  all(is.finite(metadata$nCount_RNA)), all(metadata$nCount_RNA > 0))
pt_all <- as.numeric(metadata$pt)
pt_defined <- is.finite(pt_all)
pt <- pt_all[pt_defined]
pt_cells <- cells[pt_defined]
stopifnot(length(pt_cells) == 12642L, sum(!pt_defined) == 14081L)

# Reconstruct the exact in-memory window grid. Reading decimal boundaries back
# from TSV can move cells lying exactly on a boundary; the original analyses
# used these unrounded values before writing their reports.
pt_limits <- range(pt)
pt_span <- diff(pt_limits)
exact_centers <- pt_limits[[1L]] + seq(0, 1, by = 0.01) * pt_span
exact_width <- 0.05 * pt_span
exact_windows <- data.table(grid_window_index = seq_along(exact_centers),
  window_center_exact = exact_centers,
  window_left_exact = pmax(pt_limits[[1L]], exact_centers - exact_width / 2),
  window_right_exact = pmin(pt_limits[[2L]], exact_centers + exact_width / 2))
exact_windows[, n_cells := vapply(seq_len(.N), function(i)
  sum(pt >= window_left_exact[[i]] & pt <= window_right_exact[[i]]),
  integer(1))]
exact_windows <- exact_windows[n_cells >= 20L]
exact_windows[, window_index := seq_len(.N)]
stopifnot(nrow(exact_windows) == 97L)

bb <- fread(bb_path)
stopifnot(nrow(bb) == 2L,
  identical(as.character(bb$label), c("up1", "up2")))
bb[, boundary_label := c("BB1", "BB2")]

# Reconstruct exact native-compatible UCell values for only the nine selected
# signatures from the validated sparse rank matrix used by the GO analysis.
rankings <- readRDS(ranking_path)
stopifnot(identical(colnames(rankings), pt_cells))
rank_presence <- rankings
rank_presence@x[] <- 1
ranking_gene_index <- setNames(seq_len(nrow(rankings)), rownames(rankings))
signatures <- split(members$feature_id, members$feature_set)
signatures <- lapply(signatures[selection$feature_set], unique)
stopifnot(all(unlist(signatures, use.names = FALSE) %in% rownames(rankings)))
sizes <- lengths(signatures)
membership_matrix <- Matrix::sparseMatrix(
  i = rep(seq_along(signatures), sizes),
  j = unname(ranking_gene_index[unlist(signatures, use.names = FALSE)]),
  x = 1, dims = c(length(signatures), nrow(rankings)))
stored_sum <- as.matrix(membership_matrix %*% rankings)
stored_n <- as.matrix(membership_matrix %*% rank_presence)
full_rank_sum <- stored_sum + (sizes - stored_n) * 1500
rank_sum_min <- sizes * (sizes + 1) / 2
denominator <- sizes * 1500 - rank_sum_min
ucell_scores <- 1 - sweep(sweep(full_rank_sum, 1L, rank_sum_min, "-"),
  1L, denominator, "/")
ucell_scores[ucell_scores < 0 & ucell_scores > -1e-12] <- 0
rownames(ucell_scores) <- names(signatures)
colnames(ucell_scores) <- pt_cells
stopifnot(all(is.finite(ucell_scores)), all(ucell_scores >= 0),
  all(ucell_scores <= 1))
rm(rankings, rank_presence, membership_matrix, stored_sum, stored_n,
  full_rank_sum); invisible(gc())

# Reconstruct the earlier depth-adjusted, per-gene-standardized cell scores.
log_depth <- log1p(as.numeric(metadata$nCount_RNA[pt_defined]))
depth_z <- (log_depth - mean(log_depth)) / sd(log_depth)
z_scores <- matrix(NA_real_, nrow = nrow(selection), ncol = length(pt_cells),
  dimnames = list(selection$feature_set, pt_cells))
for (g in c("bp", "mf", "cc")) {
  selected_sets <- selection[analysis_group == g, feature_set]
  mem <- members[feature_set %chin% selected_sets]
  diagnostics <- read_tsv_gz(file.path(z_root, g,
    "gene_standardization_diagnostics.tsv.gz"))
  diagnostics <- diagnostics[standardizable == TRUE]
  genes <- sort(intersect(unique(mem$feature_id), diagnostics$feature_id))
  stopifnot(length(genes) > 0L)
  expr <- as.matrix(wb_matrix_read_features(bundle, "RNA_log1p_cp10k",
    genes, pt_cells, missing = "NA"))
  diagnostics <- diagnostics[match(genes, feature_id)]
  standardized <- (expr -
    tcrossprod(diagnostics$conditional_depth_beta, depth_z) -
    diagnostics$adjusted_center) / diagnostics$adjusted_scale
  rownames(standardized) <- genes
  for (fs in selected_sets) {
    fs_genes <- intersect(unique(mem[feature_set == fs, feature_id]), genes)
    stopifnot(length(fs_genes) > 0L)
    z_scores[fs, ] <- colMeans(standardized[fs_genes, , drop = FALSE])
  }
  rm(expr, standardized); invisible(gc())
}
stopifnot(all(is.finite(z_scores)))

ucell_curves <- rbindlist(lapply(c("bp", "mf", "cc"), function(g)
  read_tsv_gz(file.path(ucell_root, g, "window_curves.tsv.gz"))[
    feature_set %chin% selection$feature_set]))
z_curves <- rbindlist(lapply(c("bp", "mf", "cc"), function(g)
  read_tsv_gz(file.path(z_root, g, "window_curves.tsv.gz"))[
    feature_set %chin% selection$feature_set]))
ucell_regions <- rbindlist(lapply(c("bp", "mf", "cc"), function(g)
  fread(file.path(ucell_root, g, "significant_positive_regions.tsv"))[
    feature_set %chin% selection$feature_set]), fill = TRUE)
z_regions <- rbindlist(lapply(c("bp", "mf", "cc"), function(g)
  fread(file.path(z_root, g, "significant_positive_regions.tsv"))[
    feature_set %chin% selection$feature_set]), fill = TRUE)
stopifnot(setequal(unique(ucell_regions$feature_set), selection$feature_set),
  setequal(unique(z_regions$feature_set), selection$feature_set))

snap_regions_to_exact_grid <- function(regions) {
  regions <- copy(regions)
  regions[, region_left := vapply(region_left, function(x)
    exact_windows$window_left_exact[[which.min(abs(
      exact_windows$window_left_exact - x))]], numeric(1))]
  regions[, region_right := vapply(region_right, function(x)
    exact_windows$window_right_exact[[which.min(abs(
      exact_windows$window_right_exact - x))]], numeric(1))]
  regions
}
ucell_regions <- snap_regions_to_exact_grid(ucell_regions)
z_regions <- snap_regions_to_exact_grid(z_regions)

# Validate reconstructed cell scores against every stored curve value.
validation <- rbindlist(lapply(selection$feature_set, function(fs) {
  uc <- ucell_curves[feature_set == fs][order(window_index)]
  zc <- z_curves[feature_set == fs][order(window_index)]
  uc_reproduced <- vapply(seq_len(nrow(uc)), function(i)
    median(ucell_scores[fs, pt >= exact_windows$window_left_exact[[i]] &
      pt <= exact_windows$window_right_exact[[i]]]), numeric(1))
  z_reproduced <- vapply(seq_len(nrow(zc)), function(i)
    mean(z_scores[fs, pt >= exact_windows$window_left_exact[[i]] &
      pt <= exact_windows$window_right_exact[[i]]]), numeric(1))
  data.table(feature_set = fs,
    maximum_ucell_window_difference = max(abs(uc_reproduced -
      uc$ucell_median)),
    maximum_z_score_window_difference = max(abs(z_reproduced -
      zc$score_mean)))
}))
fwrite(validation, file.path(stage_dir, "score_reconstruction_validation.tsv"),
  sep = "\t", quote = FALSE)
message("Maximum UCell curve difference: ",
  format(max(validation$maximum_ucell_window_difference), scientific = TRUE))
message("Maximum original-z curve difference: ",
  format(max(validation$maximum_z_score_window_difference), scientific = TRUE))
stopifnot(max(validation$maximum_ucell_window_difference) < 1e-12,
  max(validation$maximum_z_score_window_difference) < 1e-10)

base_curve_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"))
base_umap <- data.table(cell_id = cells, x = umap$x, y = umap$y,
  pt = pt_all, pt_defined = pt_defined)
figure_index <- list()

make_umap <- function(score, regions, title, legend_title) {
  mask <- interval_mask(pt_all, regions)
  score_all <- rep(NA_real_, length(cells))
  score_all[pt_defined] <- score
  lim <- safe_limits(score)
  d <- copy(base_umap)
  d[, `:=`(score = score_all, significant = mask)]
  ggplot() +
    geom_point(data = d[!pt_defined], aes(x, y), color = "#D9DDE2",
      size = .16, alpha = .92, stroke = 0) +
    geom_point(data = d[pt_defined & !significant], aes(x, y),
      color = "#AEB5BD", size = .18, alpha = .92, stroke = 0) +
    geom_point(data = d[significant == TRUE][order(score)],
      aes(x, y, color = score),
      size = .46, alpha = .97, stroke = 0) +
    scale_color_gradientn(colors = c("#FECACA", "#EF4444", "#7F0000"),
      limits = lim, oob = scales::squish, name = legend_title) +
    coord_equal() + theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
      legend.position = "right") +
    labs(title = title,
      subtitle = paste0(sum(mask), " significant-region cells; ",
        sum(pt_defined & !mask), " other trajectory cells; ",
        sum(!pt_defined), " cells without pt"),
      caption = "Very light grey: pt missing; grey: not significant; red intensity: cell score")
}

for (i in seq_len(nrow(selection))) {
  fs <- selection$feature_set[[i]]
  group <- selection$analysis_group[[i]]
  uc_curve <- ucell_curves[feature_set == fs][order(window_index)]
  z_curve <- z_curves[feature_set == fs][order(window_index)]
  uc_reg <- ucell_regions[feature_set == fs][order(region_left)]
  z_reg <- z_regions[feature_set == fs][order(region_left)]

  p_uc_curve <- ggplot(uc_curve, aes(window_center, ucell_median)) +
    geom_rect(data = uc_reg, aes(xmin = region_left, xmax = region_right,
      ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
      fill = "#DC2626", alpha = .13) +
    geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
      linewidth = .35, linetype = "dashed", color = "#4B5563") +
    geom_line(color = "#252525", linewidth = .75) +
    geom_point(data = uc_reg, aes(peak_pseudotime, peak_ucell_median),
      inherit.aes = FALSE, color = "#B91C1C", size = 1.7) +
    base_curve_theme +
    labs(title = "UCell trajectory-dynamics analysis",
      subtitle = paste0(nrow(uc_reg),
        " permutation-supported positive region(s)"),
      x = "Existing Bdev pseudotime (pt)", y = "Median UCell")

  p_z_curve <- ggplot(z_curve, aes(window_center, z_score)) +
    geom_rect(data = z_reg, aes(xmin = region_left, xmax = region_right,
      ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
      fill = "#DC2626", alpha = .13) +
    geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
      linewidth = .35, linetype = "dashed", color = "#4B5563") +
    geom_hline(yintercept = 0, linewidth = .3, color = "#6B7280") +
    geom_line(color = "#252525", linewidth = .75) +
    geom_point(data = z_reg, aes(peak_pseudotime, peak_z),
      inherit.aes = FALSE, color = "#B91C1C", size = 1.7) +
    base_curve_theme +
    labs(title = "Earlier standardized-expression z analysis",
      subtitle = paste0(nrow(z_reg), " retained positive region(s)"),
      x = "Existing Bdev pseudotime (pt)", y = "Trajectory z-score")

  p_uc_umap <- make_umap(as.numeric(ucell_scores[fs, ]), uc_reg,
    "UCell significant regions on UMAP", "Cell\nUCell")
  p_z_umap <- make_umap(as.numeric(z_scores[fs, ]), z_reg,
    "Original-z significant regions on UMAP", "Standardized\nset score")
  combined <- (p_uc_curve | p_z_curve) / (p_uc_umap | p_z_umap) +
    plot_layout(heights = c(.85, 1.3)) +
    plot_annotation(title = selection$set_name[[i]],
      subtitle = paste0(toupper(group),
        "; each method uses its own significant positive regions"),
      theme = theme(plot.title = element_text(face = "bold", size = 15)))

  stem <- paste0(group, "_", slugify(fs), "_ucell_vs_original_z")
  ggsave(file.path(stage_dir, paste0(stem, ".png")), combined,
    width = 15, height = 11, dpi = 220, bg = "white")
  ggsave(file.path(stage_dir, paste0(stem, ".pdf")), combined,
    width = 15, height = 11, device = cairo_pdf, bg = "white")
  figure_index[[i]] <- data.table(analysis_group = group,
    rank_within_group = selection$rank_within_group[[i]],
    feature_set = fs, set_name = selection$set_name[[i]],
    ucell_regions = nrow(uc_reg),
    ucell_significant_cells = sum(interval_mask(pt_all, uc_reg)),
    original_z_regions = nrow(z_reg),
    original_z_significant_cells = sum(interval_mask(pt_all, z_reg)),
    png = paste0(stem, ".png"), pdf = paste0(stem, ".pdf"))
}
figure_index <- rbindlist(figure_index)
fwrite(figure_index, file.path(stage_dir, "figure_index.tsv"), sep = "\t",
  quote = FALSE)

writeLines(c("# Epigenetic GO: UCell versus original-z curves and UMAPs", "",
  paste0("- Gene sets: ", nrow(selection)),
  "- Each gene set has four panels: UCell pt curve, original standardized-expression z curve, UCell-region UMAP, and original-z-region UMAP.",
  "- Each method uses its own final significant positive regions.",
  "- UMAP includes all 26,723 cells: cells without pt are very light grey, trajectory cells outside that method's regions are grey, and significant-region cells are colored by the method-specific cell score.",
  "- UCell cell scores were reconstructed from the validated maxRank=1500 rank matrix.",
  "- Original-z cell scores exactly reproduce the earlier per-gene depth adjustment, standardization, and equal-weight gene-set mean.",
  "- Reconstructed scores matched all stored window curves within the recorded tolerances."),
  file.path(stage_dir, "REPORT.md"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
write_json(list(status = "complete", gene_sets = nrow(selection),
  figures = nrow(figure_index), panels_per_figure = 4L,
  all_cells_in_umap = length(cells), defined_pt_cells = sum(pt_defined),
  missing_pt_cells = sum(!pt_defined),
  maximum_ucell_window_difference =
    max(validation$maximum_ucell_window_difference),
  maximum_z_score_window_difference =
    max(validation$maximum_z_score_window_difference),
  bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = fingerprint_after,
  bundle_unchanged = TRUE, source_selection = selection_path,
  source_ucell_report = ucell_root, source_z_report = z_root,
  source_rankings = ranking_path), file.path(stage_dir, "completion.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)

wb_close_bundle(bundle)
stopifnot(file.rename(stage_dir, output_dir))
message("Completed epigenetic UCell/z comparison figures at: ", output_dir)
