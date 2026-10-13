#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 62_trajectory_geneset_standardized_min20.R",
  "<bundle_dir> <inventory.tsv.gz> <memberships.tsv.gz>",
  "<bp|mf|cc|cytokines> <output_dir> <seed> <permutations>"))

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
inventory_path <- normalizePath(args[[2L]], mustWork = TRUE)
membership_path <- normalizePath(args[[3L]], mustWork = TRUE)
group <- args[[4L]]
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
base_seed <- as.integer(args[[6L]])
n_permutations <- as.integer(args[[7L]])
stopifnot(group %in% c("bp", "mf", "cc", "cytokines"),
  is.finite(base_seed), identical(n_permutations, 1000L))
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)
figure_dir <- file.path(stage_dir, "figures")
dir.create(figure_dir)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
  library(splines)
})
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt"); on.exit(close(con), add = TRUE)
  write.table(as.data.frame(x), con, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ..., "\n")
  flush.console()
}
slugify <- function(x) gsub("(^_+|_+$)", "",
  tolower(gsub("[^A-Za-z0-9]+", "_", x)))

inventory <- read_tsv_gz(inventory_path)[analysis_group == group &
  eligible == TRUE]
members <- read_tsv_gz(membership_path)[analysis_group == group]
stopifnot(nrow(inventory) > 0L, !anyDuplicated(inventory$feature_set),
  setequal(inventory$feature_set, unique(members$feature_set)))
setorder(inventory, feature_set); setorder(members, feature_set, feature_id)
inventory[, original_inventory_index := .I]
mapped_counts <- members[, .(mapped_genes_total = uniqueN(feature_id)),
  by = feature_set]
inventory[mapped_counts, on = "feature_set",
  mapped_genes_total := i.mapped_genes_total]
stopifnot(all(is.finite(inventory$mapped_genes_total)),
  all(inventory$mapped_genes_total > 0L))

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
pt_all <- as.numeric(metadata$pt); pt_defined <- is.finite(pt_all)
pt <- pt_all[pt_defined]; pt_cells <- cells[pt_defined]
stopifnot(length(pt) == 12642L, sum(!pt_defined) == 14081L)

depth_spline_df <- 8L
log_depth_all <- log1p(as.numeric(metadata$nCount_RNA))
depth_center <- mean(log_depth_all[pt_defined])
depth_scale <- sd(log_depth_all[pt_defined])
depth_z_all <- (log_depth_all - depth_center) / depth_scale
depth_z <- depth_z_all[pt_defined]
depth_design <- cbind(1, ns(pt, df = depth_spline_df), depth_z)
depth_qr <- qr(depth_design)

width_fraction <- 0.05; step_fraction <- 0.01; min_cells <- 20L
candidate_z <- qnorm(0.99); q_cutoff <- 0.01
prominence_cutoff <- 0.10; merge_gap_fraction <- width_fraction
minimum_supporting_genes <- 3L; minimum_supporting_fraction <- 0.20
single_gene_fraction_cutoff <- 0.50
pt_limits <- range(pt); pt_span <- diff(pt_limits)
centers <- pt_limits[[1L]] + seq(0, 1, by = step_fraction) * pt_span
window_width <- width_fraction * pt_span
windows <- data.table(grid_window_index = seq_along(centers),
  window_center = centers,
  window_left = pmax(pt_limits[[1L]], centers - window_width / 2),
  window_right = pmin(pt_limits[[2L]], centers + window_width / 2))
membership_index <- lapply(seq_len(nrow(windows)), function(i)
  which(pt >= windows$window_left[[i]] & pt <= windows$window_right[[i]]))
windows[, n_cells := lengths(membership_index)]
windows <- windows[n_cells >= min_cells]
membership_index <- lapply(seq_len(nrow(windows)), function(i)
  which(pt >= windows$window_left[[i]] & pt <= windows$window_right[[i]]))
windows[, window_index := seq_len(.N)]
weight_matrix <- matrix(0, nrow(windows), length(pt))
for (i in seq_len(nrow(windows)))
  weight_matrix[i, membership_index[[i]]] <- 1 / length(membership_index[[i]])
segment_id <- cumsum(c(TRUE, diff(windows$grid_window_index) != 1L))
segments <- split(seq_len(nrow(windows)), segment_id)
stopifnot(min(windows$n_cells) >= min_cells)

local_maxima <- function(values, idx) {
  if (length(idx) == 1L) return(idx)
  idx[vapply(idx, function(i) {
    (i == idx[[1L]] || values[[i]] >= values[[i - 1L]]) &&
      (i == idx[[length(idx)]] || values[[i]] >= values[[i + 1L]])
  }, logical(1))]
}
peak_prominence <- function(values, index) {
  peak <- values[[index]]
  lh <- if (index > 1L) which(values[seq.int(index - 1L, 1L)] > peak) else integer()
  rh <- if (index < length(values))
    which(values[seq.int(index + 1L, length(values))] > peak) else integer()
  lb <- if (length(lh)) index - lh[[1L]] else 1L
  rb <- if (length(rh)) index + rh[[1L]] else length(values)
  peak - max(min(values[seq.int(lb, index)]),
    min(values[seq.int(index, rb)]))
}

# Transform each unique mapped gene once: remove only its conditional depth
# component, then center and scale the adjusted profile across pt-defined cells.
genes <- sort(unique(members$feature_id))
log_event(group, "reading", length(genes), "unique mapped genes")
gene_matrix <- as.matrix(wb_matrix_read_features(bundle,
  "RNA_log1p_cp10k", genes, pt_cells, missing = "NA"))
stopifnot(identical(rownames(gene_matrix), genes), all(is.finite(gene_matrix)))
gene_beta <- gene_center <- gene_scale <- rep(NA_real_, length(genes))
block_size <- 256L
for (start in seq.int(1L, length(genes), by = block_size)) {
  idx <- seq.int(start, min(length(genes), start + block_size - 1L))
  block <- gene_matrix[idx, , drop = FALSE]
  coefs <- qr.coef(depth_qr, t(block))
  beta <- as.numeric(coefs[ncol(depth_design), ])
  adjusted <- block - tcrossprod(beta, depth_z)
  center <- rowMeans(adjusted)
  centered <- adjusted - center
  scale <- sqrt(rowSums(centered^2) / (ncol(centered) - 1L))
  valid <- is.finite(scale) & scale > .Machine$double.eps
  if (any(valid))
    gene_matrix[idx[valid], ] <- centered[valid, , drop = FALSE] / scale[valid]
  if (any(!valid)) gene_matrix[idx[!valid], ] <- NA_real_
  gene_beta[idx] <- beta; gene_center[idx] <- center; gene_scale[idx] <- scale
  if (max(idx) %% 2048L < block_size || max(idx) == length(genes))
    log_event(group, "standardized", max(idx), "of", length(genes), "genes")
}
names(gene_beta) <- names(gene_center) <- names(gene_scale) <- genes
standardizable <- is.finite(gene_scale) & gene_scale > .Machine$double.eps
gene_transform <- data.table(analysis_group = group, feature_id = genes,
  conditional_depth_beta = gene_beta, adjusted_center = gene_center,
  adjusted_scale = gene_scale, standardizable = standardizable)
write_tsv_gz(gene_transform,
  file.path(stage_dir, "gene_standardization_diagnostics.tsv.gz"))

valid_genes <- genes[standardizable]
valid_gene_matrix <- gene_matrix[standardizable, , drop = FALSE]
rm(gene_matrix); invisible(gc())
members_valid <- members[feature_id %chin% valid_genes]
standardizable_counts <- members_valid[,
  .(standardizable_genes = uniqueN(feature_id)), by = feature_set]
inventory[standardizable_counts, on = "feature_set",
  standardizable_genes := i.standardizable_genes]
inventory[is.na(standardizable_genes), standardizable_genes := 0L]

# Sparse set-by-gene averaging matrix gives one equal-weight contribution per
# standardized gene, irrespective of its original abundance or variance.
set_index <- setNames(seq_len(nrow(inventory)), inventory$feature_set)
gene_index <- setNames(seq_along(valid_genes), valid_genes)
mv <- members_valid[feature_set %chin% inventory$feature_set]
mv[, `:=`(set_i = unname(set_index[feature_set]),
  gene_i = unname(gene_index[feature_id]))]
mv[, weight := 1 / uniqueN(feature_id), by = feature_set]
averaging_matrix <- sparseMatrix(i = mv$set_i, j = mv$gene_i,
  x = mv$weight, dims = c(nrow(inventory), length(valid_genes)))
score_matrix <- as.matrix(averaging_matrix %*% valid_gene_matrix)
gene_window_matrix <- valid_gene_matrix %*% t(weight_matrix)
rownames(gene_window_matrix) <- valid_genes
rm(valid_gene_matrix, averaging_matrix); invisible(gc())

analyze_one <- function(ii) {
  fs <- inventory$feature_set[[ii]]
  scores <- score_matrix[ii, ]
  if (!all(is.finite(scores)) || sd(scores) <= .Machine$double.eps) {
    return(list(curve = NULL, candidates = NULL,
      permutation = data.table(analysis_group = group, feature_set = fs,
        candidate_peaks = 0L, permutations_used = 0L,
        strongest_peak_z = NA_real_, strongest_peak_exceedances = NA_integer_,
        stop_reason = "no_standardizable_or_variable_genes")))
  }
  curve <- as.numeric(weight_matrix %*% scores)
  null_mean <- mean(scores)
  null_sd <- sd(scores) * sqrt((1 - windows$n_cells / length(scores)) /
    windows$n_cells)
  z <- (curve - null_mean) / null_sd
  curve_row <- data.table(analysis_group = group, feature_set = fs, windows,
    score_mean = curve, null_mean = null_mean, null_sd = null_sd, z_score = z)
  peaks <- unlist(lapply(segments, function(idx) local_maxima(curve, idx)),
    use.names = FALSE)
  peaks <- peaks[z[peaks] >= candidate_z]
  if (!length(peaks)) return(list(curve = curve_row, candidates = NULL,
    permutation = data.table(analysis_group = group, feature_set = fs,
      candidate_peaks = 0L, permutations_used = 0L,
      strongest_peak_z = NA_real_, strongest_peak_exceedances = NA_integer_,
      stop_reason = "no_peak_above_candidate_z")))

  exceed <- integer(length(peaks))
  set.seed(base_seed + inventory$original_inventory_index[[ii]] - 1L)
  for (start in seq.int(1L, n_permutations, by = 100L)) {
    take <- min(100L, n_permutations - start + 1L)
    permuted <- replicate(take, sample(scores, length(scores), replace = FALSE))
    if (is.null(dim(permuted))) permuted <- matrix(permuted, ncol = 1L)
    null_curves <- weight_matrix %*% permuted
    null_z <- sweep(sweep(null_curves, 1L, null_mean, "-"), 1L,
      null_sd, "/")
    null_max <- apply(null_z, 2L, max)
    exceed <- exceed + vapply(z[peaks], function(v) sum(null_max >= v),
      integer(1))
  }
  score_range <- diff(range(curve))
  candidate_row <- rbindlist(lapply(seq_along(peaks), function(jj) {
    pidx <- peaks[[jj]]; prom <- peak_prominence(curve, pidx)
    data.table(analysis_group = group, feature_set = fs,
      set_name = inventory$set_name[[ii]],
      description = inventory$description[[ii]], peak_window = pidx,
      peak_pseudotime = windows$window_center[[pidx]],
      region_left = windows$window_left[[pidx]],
      region_right = windows$window_right[[pidx]],
      n_cells = windows$n_cells[[pidx]], peak_score = curve[[pidx]],
      peak_z = z[[pidx]], peak_prominence = prom,
      relative_peak_prominence = prom / score_range,
      exceedances = exceed[[jj]], permutations_used = n_permutations,
      p_region = (1 + exceed[[jj]]) / (1 + n_permutations))
  }))
  list(curve = curve_row, candidates = candidate_row,
    permutation = data.table(analysis_group = group, feature_set = fs,
      candidate_peaks = length(peaks), permutations_used = n_permutations,
      strongest_peak_z = max(z[peaks]),
      strongest_peak_exceedances = min(exceed), stop_reason = "complete"))
}

workers_requested <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", "8")))
if (!is.finite(workers_requested) || workers_requested < 1L)
  stop("N_WORKERS must be a positive integer")
workers <- min(workers_requested, parallel::detectCores(), nrow(inventory))
log_event(group, "starting", n_permutations, "fixed permutations for",
  nrow(inventory), "sets with", workers, "workers")
analysis_results <- vector("list", nrow(inventory)); batch_size <- 50L
for (batch_start in seq.int(1L, nrow(inventory), by = batch_size)) {
  batch <- seq.int(batch_start,
    min(nrow(inventory), batch_start + batch_size - 1L))
  batch_results <- parallel::mclapply(batch, analyze_one,
    mc.cores = min(workers, length(batch)), mc.preschedule = TRUE)
  if (any(vapply(batch_results, inherits, logical(1), "try-error")))
    stop("Numeric analysis failed in batch beginning ", batch_start)
  analysis_results[batch] <- batch_results
  log_event(group, "completed", max(batch), "of", nrow(inventory), "sets")
}
curve_rows <- Filter(Negate(is.null), lapply(analysis_results, `[[`, "curve"))
candidate_rows <- Filter(Negate(is.null),
  lapply(analysis_results, `[[`, "candidates"))
curves <- rbindlist(curve_rows)
candidates <- rbindlist(candidate_rows)
permutation_summary <- rbindlist(lapply(analysis_results, `[[`, "permutation"))
rm(analysis_results); invisible(gc())
if (!nrow(candidates)) stop("No candidate peaks found for ", group)

candidates[, q_region := p.adjust(p_region, method = "BH")]
setorder(candidates, q_region, p_region, -peak_z, feature_set, peak_pseudotime)
candidates[, peak_id := sprintf(paste0(group, "_peak_%05d"), seq_len(.N))]
candidates[inventory[, .(feature_set, mapped_genes_total,
  standardizable_genes)], on = "feature_set",
  `:=`(mapped_genes_total = i.mapped_genes_total,
    standardizable_genes = i.standardizable_genes)]

# Contribution checks are evaluated at each individual discovered peak before
# any merging. Positive standardized gene-window means count as support.
support_rows <- lapply(seq_len(nrow(candidates)), function(i) {
  fs <- candidates$feature_set[[i]]
  fs_genes <- members_valid[feature_set == fs, feature_id]
  contribution <- gene_window_matrix[fs_genes,
    candidates$peak_window[[i]], drop = TRUE]
  abs_total <- sum(abs(contribution))
  top_i <- if (length(contribution)) which.max(abs(contribution)) else NA_integer_
  data.table(supporting_genes = sum(contribution > 0),
    supporting_fraction_mapped = sum(contribution > 0) /
      candidates$mapped_genes_total[[i]],
    maximum_single_gene_absolute_fraction = if (abs_total > 0)
      max(abs(contribution)) / abs_total else NA_real_,
    top_absolute_contributor = if (length(contribution)) fs_genes[[top_i]] else NA_character_,
    top_absolute_contribution = if (length(contribution)) contribution[[top_i]] else NA_real_)
})
candidates <- cbind(candidates, rbindlist(support_rows))
candidates[, required_supporting_genes := pmax(minimum_supporting_genes,
  ceiling(minimum_supporting_fraction * mapped_genes_total))]
candidates[, `:=`(
  significant = q_region < q_cutoff,
  passes_prominence = is.finite(relative_peak_prominence) &
    relative_peak_prominence >= prominence_cutoff,
  passes_multi_gene_support = supporting_genes >= required_supporting_genes,
  passes_single_gene_dominance =
    is.finite(maximum_single_gene_absolute_fraction) &
    maximum_single_gene_absolute_fraction <= single_gene_fraction_cutoff)]
candidates[, retained_for_reporting := significant & passes_prominence &
  passes_multi_gene_support & passes_single_gene_dominance]

significant_individual <- candidates[significant == TRUE]
prominence_retained <- candidates[significant & passes_prominence]
support_retained <- candidates[significant & passes_prominence &
  passes_multi_gene_support]
retained <- candidates[retained_for_reporting == TRUE]
merge_gap_pt <- merge_gap_fraction * pt_span
merge_peaks <- function(x) {
  x <- copy(x[order(region_left, region_right)])
  gid <- integer(nrow(x)); gid[[1L]] <- 1L
  if (nrow(x) > 1L) for (i in 2:nrow(x)) {
    gap <- x$region_left[[i]] - x$region_right[[i - 1L]]
    gid[[i]] <- gid[[i - 1L]] + as.integer(gap >= merge_gap_pt)
  }
  rbindlist(lapply(split(seq_len(nrow(x)), gid), function(idx) {
    g <- x[idx]; representative <- which.max(g$peak_score)
    out <- copy(g[representative])
    out[, `:=`(region_left = min(g$region_left),
      region_right = max(g$region_right), n_constituent_peaks = nrow(g),
      constituent_peak_pts = paste(signif(g$peak_pseudotime, 7), collapse = ";"),
      constituent_peak_ids = paste(g$peak_id, collapse = ";"))]
    out
  }))
}
if (nrow(retained)) {
  final_regions <- rbindlist(lapply(split(seq_len(nrow(retained)),
    retained$feature_set), function(idx) merge_peaks(retained[idx])))
  setorder(final_regions, q_region, p_region, -peak_z, feature_set, region_left)
  final_regions[, region_id := sprintf(paste0(group, "_region_%05d"),
    seq_len(.N))]
} else final_regions <- data.table()

write_tsv_gz(curves, file.path(stage_dir, "window_curves.tsv.gz"))
write_tsv_gz(candidates, file.path(stage_dir, "all_candidate_peaks.tsv.gz"))
write_tsv(significant_individual,
  file.path(stage_dir, "significant_individual_peaks_unfiltered.tsv"))
write_tsv(prominence_retained,
  file.path(stage_dir, "prominence_retained_individual_peaks.tsv"))
write_tsv(support_retained,
  file.path(stage_dir, "multi_gene_support_retained_individual_peaks.tsv"))
write_tsv(retained,
  file.path(stage_dir, "dominance_retained_constituent_peaks.tsv"))
write_tsv(final_regions, file.path(stage_dir, "significant_positive_regions.tsv"))
write_tsv(permutation_summary,
  file.path(stage_dir, "permutation_summary.tsv"))
write_tsv(inventory[, .(analysis_group, feature_set, set_name,
  mapped_genes_total, standardizable_genes)],
  file.path(stage_dir, "gene_set_mapped_gene_counts.tsv"))

# Preserve the earlier reporting convention: UMAP/curve figures for only the
# 20 most significant final gene sets in each group.
significant_sets <- unique(final_regions$feature_set)
sets_to_plot <- head(significant_sets, 20L)
figure_index <- list()
base_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"))
for (ii in seq_along(sets_to_plot)) {
  fs <- sets_to_plot[[ii]]
  fs_genes <- members_valid[feature_set == fs, feature_id]
  expr <- as.matrix(wb_matrix_read_features(bundle, "RNA_log1p_cp10k",
    fs_genes, cells, missing = "NA"))
  adjusted <- expr - tcrossprod(gene_beta[fs_genes], depth_z_all)
  standardized <- (adjusted - gene_center[fs_genes]) / gene_scale[fs_genes]
  score <- colMeans(standardized)
  pos <- final_regions[feature_set == fs][order(region_left)]
  curve <- curves[feature_set == fs]
  mask <- rep(FALSE, length(cells))
  for (j in seq_len(nrow(pos))) mask <- mask |
    (is.finite(pt_all) & pt_all >= pos$region_left[[j]] &
      pt_all <= pos$region_right[[j]])
  limits <- as.numeric(quantile(score, c(.01, .99), names = FALSE))
  scaled <- pmin(1, pmax(0, (score - limits[[1L]]) / diff(limits)))
  palette <- grDevices::colorRampPalette(c("#FECACA", "#B91C1C"))(256L)
  ci <- pmax(1L, pmin(256L, 1L + floor(scaled * 255L)))
  d <- data.table(x = umap$x, y = umap$y, scaled = scaled,
    status = fifelse(!pt_defined, "missing_pt",
      fifelse(mask, "significant", "not_significant")),
    point_color = fifelse(!pt_defined, "#CDD2D8",
      fifelse(mask, palette[ci], "#AEB5BD")))
  p_umap <- ggplot() +
    geom_point(data = d[status == "missing_pt"], aes(x, y), color = "#CDD2D8",
      size = .17, alpha = .90, stroke = 0) +
    geom_point(data = d[status == "not_significant"], aes(x, y),
      color = "#AEB5BD", size = .19, alpha = .92, stroke = 0) +
    geom_point(data = d[status == "significant"][order(scaled)],
      aes(x, y, color = point_color), size = .48, alpha = .96, stroke = 0) +
    scale_color_identity() + coord_equal() + theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold")) +
    labs(title = paste(inventory[feature_set == fs, set_name],
      "standardized score"),
      subtitle = "Red: retained region; greys: trajectory/non-trajectory",
      caption = "Darker red = higher standardized gene-set score.")
  p_curve <- ggplot(curve, aes(window_center, score_mean)) +
    geom_rect(data = pos, aes(xmin = region_left, xmax = region_right,
      ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
      fill = "#D81B60", alpha = .13) +
    geom_line(color = "#252525", linewidth = .8) +
    geom_point(data = pos, aes(peak_pseudotime, peak_score),
      inherit.aes = FALSE, color = "#D81B60", size = 2) + base_theme +
    labs(title = "Positive trajectory enrichment",
      subtitle = "Per-gene depth adjusted and standardized",
      x = "Existing Bdev pseudotime (pt)", y = "Mean standardized score",
      caption = "5% windows; 1% steps; minimum 20 cells; 1,000 permutations")
  filename <- paste0(sprintf("%04d", ii), "_", slugify(fs), ".png")
  ggsave(file.path(figure_dir, filename), p_umap | p_curve,
    width = 12.5, height = 5.6, dpi = 180, bg = "white")
  figure_index[[ii]] <- data.table(analysis_group = group, feature_set = fs,
    set_name = inventory[feature_set == fs, set_name], n_regions = nrow(pos),
    significant_cells = sum(mask), figure = file.path("figures", filename))
}
figure_index <- if (length(figure_index)) rbindlist(figure_index) else data.table()
write_tsv(figure_index, file.path(stage_dir, "umap_figure_index.tsv"))

filter_counts <- data.table(analysis_group = group,
  tested_gene_sets = nrow(inventory), candidate_peaks = nrow(candidates),
  significant_peaks = nrow(significant_individual),
  prominence_retained_peaks = nrow(prominence_retained),
  multi_gene_support_retained_peaks = nrow(support_retained),
  dominance_retained_constituent_peaks = nrow(retained),
  final_gene_sets = uniqueN(final_regions$feature_set),
  final_regions = nrow(final_regions))
write_tsv(filter_counts, file.path(stage_dir, "filter_counts.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
completion <- list(status = "complete", analysis_group = group,
  tested_sets = nrow(inventory), candidate_peaks = nrow(candidates),
  significant_individual_peaks = nrow(significant_individual),
  prominence_retained_individual_peaks = nrow(prominence_retained),
  multi_gene_support_retained_individual_peaks = nrow(support_retained),
  dominance_retained_constituent_peaks = nrow(retained),
  significant_positive_regions = nrow(final_regions),
  sets_with_significant_positive_regions = uniqueN(final_regions$feature_set),
  figures_generated_for_top_significant_sets = length(sets_to_plot),
  score = "mean of individually depth-adjusted and standardized genes",
  depth_adjustment = paste0("per gene: subtract conditional scaled ",
    "log1p(nCount_RNA) coefficient after ns(pt,df=8)"),
  minimum_window_cells = min_cells, window_width_fraction = width_fraction,
  window_step_fraction = step_fraction, permutations = n_permutations,
  permutation_null = "maximum positive window z-score across trajectory",
  multiple_testing = paste("BH within", group, "candidate peaks only"),
  q_cutoff = q_cutoff, relative_prominence_cutoff = prominence_cutoff,
  minimum_supporting_genes = minimum_supporting_genes,
  minimum_supporting_fraction_of_mapped_genes = minimum_supporting_fraction,
  maximum_single_gene_absolute_contribution_fraction =
    single_gene_fraction_cutoff,
  merge_gap_fraction_strictly_less_than = merge_gap_fraction,
  pseudotime_source = "exact Bdev metadata pt", defined_pt_cells = length(pt),
  missing_pt_cells = sum(!pt_defined), promoted = FALSE,
  random_seed = base_seed, workers = workers,
  bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = fingerprint_after, bundle_unchanged = TRUE)
write_json(completion, file.path(stage_dir, "completion.json"))

writeLines(c(paste0("# ", group,
  " standardized-gene trajectory enrichment"), "",
  paste0("- Tested gene sets: ", nrow(inventory)),
  paste0("- Significant individual peaks: ", nrow(significant_individual)),
  paste0("- After 10% prominence: ", nrow(prominence_retained)),
  paste0("- After multi-gene support: ", nrow(support_retained)),
  paste0("- After <=50% single-gene dominance: ", nrow(retained)),
  paste0("- Final merged regions: ", nrow(final_regions), " across ",
    uniqueN(final_regions$feature_set), " gene sets"),
  "- Scores: per-gene depth adjustment, per-gene standardization, then averaging.",
  "- Windows: 5% width, 1% step, minimum 20 cells.",
  "- Testing: 1,000 fixed permutations; group-wise BH q<0.01.",
  "- Support: >=3 genes and >=20% of all mapped genes positive at the peak.",
  "- Dominance: no gene >50% of total absolute peak contribution.",
  "- Peak discovery precedes downstream prominence/support/dominance/merge filters.",
  "- Merge: uncovered gap strictly less than one 5% window."),
  file.path(stage_dir, "REPORT.md"))

wb_close_bundle(bundle)
stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed", group, "standardized min20 report at", output_dir)
