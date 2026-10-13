#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "Usage: 70_diagnose_ucell_boundary_peak_alignment.R",
  "<bundle_dir> <memberships.tsv.gz> <selection.tsv> <rankings.rds>",
  "<ucell_report> <z_report> <bb.tsv> <output_dir> <expected_sets>"))
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
  library(jsonlite)
})

read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA")

selection <- fread(selection_path)
members_all <- read_tsv_gz(membership_path)
members <- members_all[feature_set %chin% selection$feature_set]
stopifnot(nrow(selection) == expected_sets,
  setequal(selection$feature_set, unique(members$feature_set)))
bb <- fread(bb_path)
stopifnot(nrow(bb) == 2L)

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
cells <- as.character(wb_cells(bundle, "all"))
metadata <- as.data.table(wb_meta(bundle, cells = cells))
pt_defined <- is.finite(metadata$pt)
pt <- as.numeric(metadata$pt[pt_defined])
pt_cells <- cells[pt_defined]
stopifnot(length(pt_cells) == 12642L)

rankings <- readRDS(ranking_path)
stopifnot(identical(colnames(rankings), pt_cells))
rank_presence <- rankings
rank_presence@x[] <- 1
gene_index <- setNames(seq_len(nrow(rankings)), rownames(rankings))
signatures <- split(members$feature_id, members$feature_set)
signatures <- lapply(signatures[selection$feature_set], unique)
sizes <- lengths(signatures)
mm <- Matrix::sparseMatrix(i = rep(seq_along(signatures), sizes),
  j = unname(gene_index[unlist(signatures, use.names = FALSE)]), x = 1,
  dims = c(length(signatures), nrow(rankings)))
stored_sum <- as.matrix(mm %*% rankings)
stored_n <- as.matrix(mm %*% rank_presence)
rank_min <- sizes * (sizes + 1) / 2
scores <- 1 - sweep(sweep(stored_sum + (sizes - stored_n) * 1500,
  1L, rank_min, "-"), 1L, sizes * 1500 - rank_min, "/")
scores[scores < 0 & scores > -1e-12] <- 0
rownames(scores) <- names(signatures)
colnames(scores) <- pt_cells
stopifnot(all(is.finite(scores)), all(scores >= 0), all(scores <= 1))
rm(rankings, rank_presence, mm, stored_sum, stored_n); invisible(gc())

pt_limits <- range(pt)
pt_span <- diff(pt_limits)
window_width <- .05 * pt_span
half_window <- window_width / 2

# Compare final peak localization across all GO sets, not only the selected 9.
peak_proximity <- rbindlist(lapply(c("bp", "mf", "cc"), function(g) {
  rbindlist(lapply(c("ucell", "original_z"), function(method) {
    root <- if (method == "ucell") ucell_root else z_root
    x <- fread(file.path(root, g, "significant_positive_regions.tsv"))
    after_half <- vapply(x$peak_pseudotime, function(p)
      any(p >= bb$pt_median & p <= bb$pt_median + half_window), logical(1))
    after_full <- vapply(x$peak_pseudotime, function(p)
      any(p >= bb$pt_median & p <= bb$pt_median + window_width), logical(1))
    data.table(analysis_group = g, method = method, regions = nrow(x),
      within_half_window_after_boundary = sum(after_half),
      fraction_within_half_window_after_boundary = mean(after_half),
      within_full_window_after_boundary = sum(after_full),
      fraction_within_full_window_after_boundary = mean(after_full))
  }))
}))
write_tsv(peak_proximity, file.path(stage_dir,
  "all_go_peak_boundary_proximity.tsv"))

selected_offsets <- rbindlist(lapply(c("bp", "mf", "cc"), function(g) {
  rbindlist(lapply(c("ucell", "original_z"), function(method) {
    root <- if (method == "ucell") ucell_root else z_root
    x <- fread(file.path(root, g, "significant_positive_regions.tsv"))[
      feature_set %chin% selection$feature_set]
    rbindlist(lapply(seq_len(nrow(x)), function(i) {
      delta <- x$peak_pseudotime[[i]] - bb$pt_median
      j <- which.min(abs(delta))
      data.table(analysis_group = g, feature_set = x$feature_set[[i]],
        method = method, peak_pseudotime = x$peak_pseudotime[[i]],
        nearest_boundary = bb$label[[j]], offset_from_boundary = delta[[j]],
        within_half_window_after = delta[[j]] >= 0 & delta[[j]] <= half_window,
        within_full_window_after = delta[[j]] >= 0 & delta[[j]] <= window_width)
    }))
  }))
}))
write_tsv(selected_offsets, file.path(stage_dir,
  "selected_set_peak_offsets.tsv"))

# Cell-level dependence on library complexity. UCell is rank based, but sparse
# maxRank scoring can still correlate with the number of detected features.
complexity <- data.table(
  log1p_nCount_RNA = log1p(as.numeric(metadata$nCount_RNA[pt_defined])),
  nFeature_RNA = as.numeric(metadata$nFeature_RNA[pt_defined]),
  log1p_nCount_ATAC = log1p(as.numeric(metadata$nCount_ATAC[pt_defined])),
  nFeature_ATAC = as.numeric(metadata$nFeature_ATAC[pt_defined]))
complexity_cor <- rbindlist(lapply(rownames(scores), function(fs)
  rbindlist(lapply(names(complexity), function(v) data.table(
    feature_set = fs, covariate = v,
    spearman_rho = cor(as.numeric(scores[fs, ]), complexity[[v]],
      method = "spearman", use = "complete.obs"))))))
write_tsv(complexity_cor, file.path(stage_dir,
  "selected_set_cell_complexity_correlations.tsv"))

# Pairwise membership overlap identifies whether synchronized curves could be
# caused simply by the same genes appearing in multiple GO categories.
overlap <- rbindlist(lapply(seq_len(length(signatures) - 1L), function(i)
  rbindlist(lapply(seq.int(i + 1L, length(signatures)), function(j) {
    a <- signatures[[i]]; b <- signatures[[j]]
    data.table(feature_set_1 = names(signatures)[[i]],
      feature_set_2 = names(signatures)[[j]], shared_genes = length(intersect(a,b)),
      jaccard = length(intersect(a,b)) / length(union(a,b)),
      overlap_smaller_set = length(intersect(a,b)) / min(length(a), length(b)))
  }))))
setorder(overlap, -jaccard, -overlap_smaller_set)
write_tsv(overlap, file.path(stage_dir,
  "selected_set_pairwise_gene_overlap.tsv"))

# Sensitivity of the local maximum near each boundary to window width and to
# median versus mean aggregation. This is descriptive and does not redefine
# significance.
specs <- data.table(spec = c("median_3pct", "median_5pct", "median_7pct",
  "mean_5pct"), width_fraction = c(.03, .05, .07, .05),
  statistic = c("median", "median", "median", "mean"))
sensitivity <- rbindlist(lapply(seq_len(nrow(specs)), function(si) {
  s <- specs[si]
  centers <- pt_limits[[1L]] + seq(0, 1, by = .01) * pt_span
  ww <- s$width_fraction * pt_span
  windows <- data.table(center = centers,
    left = pmax(pt_limits[[1L]], centers - ww / 2),
    right = pmin(pt_limits[[2L]], centers + ww / 2))
  windows[, idx := lapply(seq_len(.N), function(i)
    which(pt >= left[[i]] & pt <= right[[i]]))]
  windows[, n_cells := lengths(idx)]
  windows <- windows[n_cells >= 20L]
  rbindlist(lapply(rownames(scores), function(fs) {
    values <- vapply(windows$idx, function(idx) if (s$statistic == "median")
      median(scores[fs, idx]) else mean(scores[fs, idx]), numeric(1))
    rbindlist(lapply(seq_len(nrow(bb)), function(bi) {
      eligible <- which(windows$center >= bb$pt_median[[bi]] - 2 &
        windows$center <= bb$pt_median[[bi]] + 2)
      best <- eligible[[which.max(values[eligible])]]
      data.table(feature_set = fs, spec = s$spec,
        width_fraction = s$width_fraction, statistic = s$statistic,
        boundary = bb$label[[bi]], boundary_pt = bb$pt_median[[bi]],
        local_max_center = windows$center[[best]],
        offset_from_boundary = windows$center[[best]] - bb$pt_median[[bi]],
        local_max_value = values[[best]])
    }))
  }))
}))
write_tsv(sensitivity, file.path(stage_dir,
  "selected_set_window_sensitivity.tsv"))

# Reproduce stored selected UCell curves at the original 5% setting.
stored_curves <- rbindlist(lapply(c("bp", "mf", "cc"), function(g)
  read_tsv_gz(file.path(ucell_root, g, "window_curves.tsv.gz"))[
    feature_set %chin% selection$feature_set]))
reproduction <- sensitivity[spec == "median_5pct",
  .(feature_set, boundary, local_max_center, local_max_value)]

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))

fmt_pct <- function(x) paste0(round(100 * x, 1), "%")
prox_rows <- vapply(seq_len(nrow(peak_proximity)), function(i) paste0("| ",
  toupper(peak_proximity$analysis_group[[i]]), " | ",
  peak_proximity$method[[i]], " | ", peak_proximity$regions[[i]], " | ",
  fmt_pct(peak_proximity$fraction_within_half_window_after_boundary[[i]]),
  " | ", fmt_pct(peak_proximity$fraction_within_full_window_after_boundary[[i]]),
  " |"), character(1))
median_abs_cor <- complexity_cor[, .(median_absolute_rho = median(abs(spearman_rho)),
  maximum_absolute_rho = max(abs(spearman_rho))), by = covariate]
cor_rows <- vapply(seq_len(nrow(median_abs_cor)), function(i) paste0("| ",
  median_abs_cor$covariate[[i]], " | ",
  signif(median_abs_cor$median_absolute_rho[[i]], 3), " | ",
  signif(median_abs_cor$maximum_absolute_rho[[i]], 3), " |"), character(1))
writeLines(c("# Diagnosis of UCell peak alignment after BB boundaries", "",
  paste0("The 5% window width is ", signif(window_width, 5),
    " pt and its half-width is ", signif(half_window, 5), " pt."), "",
  "## Peak proximity across all final GO regions", "",
  "| Group | Method | Regions | Within half-window after BB1/BB2 | Within full window after BB1/BB2 |",
  "|---|---|---:|---:|---:|", prox_rows, "",
  "## Cell-level UCell association with complexity", "",
  "| Covariate | Median absolute Spearman rho across 9 sets | Maximum absolute rho |",
  "|---|---:|---:|", cor_rows, "",
  paste0("Maximum pairwise gene-set Jaccard overlap: ",
    signif(max(overlap$jaccard), 4), "."),
  "Window-width/statistic sensitivity, exact selected peak offsets, gene overlap, and per-set complexity correlations are provided as TSV files.",
  "These diagnostics are descriptive and do not change any significance calls."),
  file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", selected_sets = nrow(selection),
  window_width_pt = window_width, half_window_width_pt = half_window,
  maximum_pairwise_jaccard = max(overlap$jaccard),
  bundle_unchanged = TRUE, source_ucell_report = ucell_root,
  source_z_report = z_root), file.path(stage_dir, "completion.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)

wb_close_bundle(bundle)
stopifnot(file.rename(stage_dir, output_dir))
message("Completed UCell boundary-alignment diagnosis at: ", output_dir)
