#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "Usage: 71_diagnose_ucell_rna_complexity_residualization.R",
  "<bundle_dir> <inventory.tsv.gz> <memberships.tsv.gz> <selection.tsv>",
  "<rankings.rds> <ucell_report> <bb.tsv> <output_dir> <permutations>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
inventory_path <- normalizePath(args[[2L]], mustWork = TRUE)
membership_path <- normalizePath(args[[3L]], mustWork = TRUE)
selection_path <- normalizePath(args[[4L]], mustWork = TRUE)
ranking_path <- normalizePath(args[[5L]], mustWork = TRUE)
ucell_root <- normalizePath(args[[6L]], mustWork = TRUE)
bb_path <- normalizePath(args[[7L]], mustWork = TRUE)
output_dir <- normalizePath(args[[8L]], mustWork = FALSE)
n_permutations <- as.integer(args[[9L]])
stopifnot(identical(n_permutations, 1000L))
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(Matrix)
  library(data.table)
  library(matrixStats)
  library(ggplot2)
  library(jsonlite)
})

read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA")
fmt_pct <- function(x) ifelse(is.finite(x), paste0(round(100 * x, 1), "%"), "NA")

inventory <- read_tsv_gz(inventory_path)[eligible == TRUE]
members_all <- read_tsv_gz(membership_path)
selection <- fread(selection_path)
bb <- fread(bb_path)
stopifnot(nrow(selection) == 9L, nrow(bb) == 2L,
  all(c("analysis_group", "feature_set") %in% names(selection)),
  !anyDuplicated(selection$feature_set))
members <- members_all[feature_set %chin% selection$feature_set]
stopifnot(setequal(selection$feature_set, unique(members$feature_set)))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
cells <- as.character(wb_cells(bundle, "all"))
metadata <- as.data.table(wb_meta(bundle, cells = cells))
stopifnot(identical(as.character(metadata$cell_id), cells))
pt_defined <- is.finite(metadata$pt)
pt <- as.numeric(metadata$pt[pt_defined])
pt_cells <- cells[pt_defined]
rna_features <- as.numeric(metadata$nFeature_RNA[pt_defined])
stopifnot(length(pt_cells) == 12642L, all(is.finite(rna_features)))

rankings <- readRDS(ranking_path)
stopifnot(identical(colnames(rankings), pt_cells))
rank_presence <- rankings
rank_presence@x[] <- 1
gene_index <- setNames(seq_len(nrow(rankings)), rownames(rankings))
signatures <- split(members$feature_id, members$feature_set)
signatures <- lapply(signatures[selection$feature_set], unique)
sizes <- lengths(signatures)
membership <- Matrix::sparseMatrix(i = rep(seq_along(signatures), sizes),
  j = unname(gene_index[unlist(signatures, use.names = FALSE)]), x = 1,
  dims = c(length(signatures), nrow(rankings)))
stored_sum <- as.matrix(membership %*% rankings)
stored_n <- as.matrix(membership %*% rank_presence)
rank_min <- sizes * (sizes + 1) / 2
raw_scores <- 1 - sweep(sweep(stored_sum + (sizes - stored_n) * 1500,
  1L, rank_min, "-"), 1L, sizes * 1500 - rank_min, "/")
raw_scores[raw_scores < 0 & raw_scores > -1e-12] <- 0
rownames(raw_scores) <- names(signatures)
colnames(raw_scores) <- pt_cells
stopifnot(all(is.finite(raw_scores)), all(raw_scores >= 0),
  all(raw_scores <= 1))
rm(rankings, rank_presence, membership, stored_sum, stored_n)
invisible(gc())

# Diagnostic adjustment only: remove the smooth marginal association between
# UCell and RNA feature complexity without using pseudotime in the model.
complexity_x <- log1p(rna_features)
residual_scores <- matrix(NA_real_, nrow(raw_scores), ncol(raw_scores),
  dimnames = dimnames(raw_scores))
adjustment_qc <- rbindlist(lapply(rownames(raw_scores), function(fs) {
  y <- as.numeric(raw_scores[fs, ])
  fit <- lm(y ~ splines::ns(complexity_x, df = 4L))
  residual_scores[fs, ] <<- residuals(fit)
  data.table(analysis_group = selection$analysis_group[
      match(fs, selection$feature_set)], feature_set = fs,
    model = "lm(UCell ~ ns(log1p(nFeature_RNA), df=4))",
    model_r_squared = summary(fit)$r.squared,
    raw_spearman_nFeature_RNA = cor(y, rna_features,
      method = "spearman"),
    residual_spearman_nFeature_RNA = cor(residuals(fit), rna_features,
      method = "spearman"),
    raw_sd = sd(y), residual_sd = sd(residuals(fit)))
}))
write_tsv(adjustment_qc, file.path(stage_dir, "complexity_adjustment_qc.tsv"))

width_fraction <- 0.05
step_fraction <- 0.01
min_cells <- 20L
pt_limits <- range(pt)
pt_span <- diff(pt_limits)
window_width <- width_fraction * pt_span
half_window <- window_width / 2
centers <- pt_limits[[1L]] + seq(0, 1, by = step_fraction) * pt_span
windows <- data.table(grid_window_index = seq_along(centers),
  window_center = centers,
  window_left = pmax(pt_limits[[1L]], centers - window_width / 2),
  window_right = pmin(pt_limits[[2L]], centers + window_width / 2))
membership_index <- lapply(seq_len(nrow(windows)), function(i)
  which(pt >= windows$window_left[[i]] & pt <= windows$window_right[[i]]))
windows[, n_cells := lengths(membership_index)]
keep <- windows$n_cells >= min_cells
windows <- windows[keep]
membership_index <- membership_index[keep]
windows[, window_index := seq_len(.N)]
stopifnot(nrow(windows) == 97L)

local_maxima_plateau <- function(values) {
  rr <- rle(values)
  ends <- cumsum(rr$lengths)
  starts <- ends - rr$lengths + 1L
  left <- c(-Inf, head(rr$values, -1L))
  right <- c(tail(rr$values, -1L), -Inf)
  keep <- rr$values >= left & rr$values >= right &
    (rr$values > left | rr$values > right)
  as.integer(floor((starts[keep] + ends[keep]) / 2))
}

analyze_scores <- function(scores, fs, method, seed) {
  baseline <- median(scores)
  observed <- vapply(membership_index,
    function(idx) median(scores[idx]), numeric(1))
  effect <- observed - baseline
  observed_range <- max(observed) - min(observed)
  peak_index <- local_maxima_plateau(observed)
  peak_index <- peak_index[effect[peak_index] > .Machine$double.eps]
  null_max <- numeric(n_permutations)
  null_range <- numeric(n_permutations)
  set.seed(seed)
  for (start in seq.int(1L, n_permutations, by = 100L)) {
    take <- min(100L, n_permutations - start + 1L)
    permuted <- replicate(take, sample(scores, length(scores), replace = FALSE))
    if (is.null(dim(permuted))) permuted <- matrix(permuted, ncol = 1L)
    batch_min <- rep(Inf, take)
    batch_max <- rep(-Inf, take)
    for (idx in membership_index) {
      wm <- matrixStats::colMedians(permuted[idx, , drop = FALSE])
      batch_min <- pmin(batch_min, wm)
      batch_max <- pmax(batch_max, wm - baseline)
    }
    pos <- start:(start + take - 1L)
    null_max[pos] <- batch_max
    null_range[pos] <- batch_max + baseline - batch_min
  }
  peaks <- if (length(peak_index)) data.table(
    analysis_group = selection$analysis_group[match(fs, selection$feature_set)],
    feature_set = fs, method = method, peak_window = peak_index,
    peak_pseudotime = windows$window_center[peak_index],
    region_left = windows$window_left[peak_index],
    region_right = windows$window_right[peak_index],
    peak_value = observed[peak_index], trajectory_median = baseline,
    peak_effect = effect[peak_index],
    p_region = (1 + vapply(effect[peak_index], function(x)
      sum(null_max >= x), integer(1))) / (1 + n_permutations)) else data.table()
  list(curve = data.table(
      analysis_group = selection$analysis_group[match(fs, selection$feature_set)],
      feature_set = fs, method = method, windows,
      window_median = observed, trajectory_median = baseline,
      effect = effect),
    peaks = peaks,
    test = data.table(
      analysis_group = selection$analysis_group[match(fs, selection$feature_set)],
      feature_set = fs, method = method,
      observed_dynamic_range = observed_range,
      null_range_mean = mean(null_range), null_range_sd = sd(null_range),
      dynamic_exceedances = sum(null_range >= observed_range),
      p_dynamic = (1 + sum(null_range >= observed_range)) /
        (1 + n_permutations)))
}

# Use the exact original per-set seeds so the raw-score branch is a strict
# reproduction control for the full GO analysis.
group_seed_offset <- c(bp = 0L, mf = 100000L, cc = 200000L)
seed_for <- function(fs, g) {
  inv <- copy(inventory[analysis_group == g])
  setorder(inv, feature_set)
  314159L + group_seed_offset[[g]] + match(fs, inv$feature_set) - 1L
}
results <- list(); cursor <- 0L
for (fs in selection$feature_set) {
  g <- selection$analysis_group[match(fs, selection$feature_set)]
  seed <- seed_for(fs, g)
  stopifnot(is.finite(seed))
  for (method in c("raw_ucell", "complexity_residual")) {
    cursor <- cursor + 1L
    score <- if (method == "raw_ucell") raw_scores[fs, ] else
      residual_scores[fs, ]
    results[[cursor]] <- analyze_scores(as.numeric(score), fs, method, seed)
  }
}
curves <- rbindlist(lapply(results, `[[`, "curve"))
peaks <- rbindlist(lapply(results, `[[`, "peaks"), fill = TRUE)
tests <- rbindlist(lapply(results, `[[`, "test"))
tests[, q_dynamic_selected3 := p.adjust(p_dynamic, method = "BH"),
  by = .(analysis_group, method)]
peaks[tests[, .(analysis_group, feature_set, method, p_dynamic,
    q_dynamic_selected3)], on = .(analysis_group, feature_set, method),
  `:=`(p_dynamic = i.p_dynamic,
    q_dynamic_selected3 = i.q_dynamic_selected3)]
peaks[, supported_peak := q_dynamic_selected3 < 0.01 & p_region < 0.01]

merge_gap_pt <- width_fraction * pt_span
merge_peaks <- function(x) {
  x <- copy(x[order(region_left, region_right)])
  gid <- integer(nrow(x)); gid[[1L]] <- 1L
  if (nrow(x) > 1L) for (i in 2:nrow(x)) {
    gap <- x$region_left[[i]] - x$region_right[[i - 1L]]
    gid[[i]] <- gid[[i - 1L]] + as.integer(gap >= merge_gap_pt)
  }
  rbindlist(lapply(split(seq_len(nrow(x)), gid), function(ii) {
    z <- x[ii]
    representative <- order(z$p_region, -z$peak_effect)[[1L]]
    out <- copy(z[representative])
    out[, `:=`(region_left = min(z$region_left),
      region_right = max(z$region_right),
      n_constituent_peaks = nrow(z))]
    out
  }))
}
supported <- peaks[supported_peak == TRUE]
regions <- if (nrow(supported)) rbindlist(lapply(
  split(seq_len(nrow(supported)), interaction(supported$feature_set,
    supported$method, drop = TRUE)), function(ii) merge_peaks(supported[ii]))) else
  data.table()

add_boundary <- function(x) {
  if (!nrow(x)) return(x)
  proximity <- rbindlist(lapply(seq_len(nrow(x)), function(i) {
    delta <- x$peak_pseudotime[[i]] - bb$pt_median
    eligible <- which(delta >= 0)
    nearest <- if (length(eligible)) eligible[[which.min(delta[eligible])]] else
      which.min(abs(delta))
    data.table(nearest_boundary = bb$label[[nearest]],
      offset_from_boundary = delta[[nearest]],
      within_half_window_after = delta[[nearest]] >= 0 &&
        delta[[nearest]] <= half_window,
      within_full_window_after = delta[[nearest]] >= 0 &&
        delta[[nearest]] <= window_width)
  }))
  cbind(x, proximity)
}
regions <- add_boundary(regions)

# Strictly validate the raw-score reconstruction and permutation control.
stored_tests <- rbindlist(lapply(c("bp", "mf", "cc"), function(g)
  fread(file.path(ucell_root, g, "gene_set_dynamic_tests.tsv"))[
    feature_set %chin% selection$feature_set]))
raw_test_check <- merge(tests[method == "raw_ucell"],
  stored_tests[, .(feature_set, stored_range = observed_dynamic_range,
    stored_p = p_dynamic)], by = "feature_set")
raw_test_check[, `:=`(
  absolute_range_difference = abs(observed_dynamic_range - stored_range),
  absolute_p_difference = abs(p_dynamic - stored_p))]
stopifnot(max(raw_test_check$absolute_range_difference) < 1e-12,
  max(raw_test_check$absolute_p_difference) < 1e-12)

stored_curves <- rbindlist(lapply(c("bp", "mf", "cc"), function(g)
  read_tsv_gz(file.path(ucell_root, g, "window_curves.tsv.gz"))[
    feature_set %chin% selection$feature_set]))
raw_curve_check <- merge(curves[method == "raw_ucell",
    .(feature_set, grid_window_index, reconstructed = window_median)],
  stored_curves[, .(feature_set, grid_window_index,
    stored = ucell_median)], by = c("feature_set", "grid_window_index"))
max_curve_difference <- max(abs(raw_curve_check$reconstructed -
  raw_curve_check$stored))
stopifnot(max_curve_difference < 1e-12)

write_tsv(tests, file.path(stage_dir, "dynamic_tests.tsv"))
write_tsv(peaks, file.path(stage_dir, "all_candidate_peaks.tsv"))
write_tsv(regions, file.path(stage_dir, "supported_merged_regions.tsv"))
write_tsv(curves, file.path(stage_dir, "window_curves.tsv.gz"))
write_tsv(raw_test_check, file.path(stage_dir,
  "raw_analysis_reproduction_check.tsv"))

region_summary <- rbindlist(lapply(c("raw_ucell", "complexity_residual"),
  function(m) {
    z <- regions[method == m]
    data.table(method = m,
      dynamically_significant_sets_selected3 = tests[method == m &
        q_dynamic_selected3 < .01, uniqueN(feature_set)],
      supported_individual_peaks = peaks[method == m & supported_peak == TRUE,
        .N], merged_regions = nrow(z),
      within_half_window_after_boundary = if (nrow(z))
        sum(z$within_half_window_after) else 0L,
      fraction_within_half_window_after_boundary = if (nrow(z))
        mean(z$within_half_window_after) else NA_real_,
      within_full_window_after_boundary = if (nrow(z))
        sum(z$within_full_window_after) else 0L,
      fraction_within_full_window_after_boundary = if (nrow(z))
        mean(z$within_full_window_after) else NA_real_)
  }))
write_tsv(region_summary, file.path(stage_dir, "boundary_pileup_summary.tsv"))

comparison <- dcast(tests, analysis_group + feature_set ~ method,
  value.var = c("observed_dynamic_range", "p_dynamic",
    "q_dynamic_selected3"))
comparison <- merge(comparison, adjustment_qc, by = c("analysis_group",
  "feature_set"))
raw_region_counts <- regions[method == "raw_ucell", .(raw_regions = .N,
  raw_full_window_after = sum(within_full_window_after)), by = feature_set]
resid_region_counts <- regions[method == "complexity_residual", .(
  residual_regions = .N,
  residual_full_window_after = sum(within_full_window_after)), by = feature_set]
comparison <- merge(comparison, raw_region_counts, by = "feature_set",
  all.x = TRUE)
comparison <- merge(comparison, resid_region_counts, by = "feature_set",
  all.x = TRUE)
for (v in c("raw_regions", "raw_full_window_after", "residual_regions",
    "residual_full_window_after")) set(comparison, which(is.na(comparison[[v]])),
      v, 0L)
write_tsv(comparison, file.path(stage_dir, "per_set_comparison.tsv"))

# Curves are standardized only within this diagnostic figure so raw UCell and
# residual units can be visually compared; all inference used unscaled values.
plot_curves <- copy(curves)
plot_curves[, display_z := {
  s <- sd(window_median)
  if (is.finite(s) && s > 0) (window_median - mean(window_median)) / s else 0
}, by = .(feature_set, method)]
plot_curves[, method_label := fifelse(method == "raw_ucell",
  "Unadjusted UCell", "UCell residualized for RNA features")]
plot_curves[, set_label := selection$set_name[match(feature_set,
  selection$feature_set)]]
p <- ggplot(plot_curves, aes(window_center, display_z,
    color = method_label)) +
  geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
    linetype = 2, linewidth = .3, color = "grey45") +
  geom_line(linewidth = .55) +
  facet_wrap(~set_label, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("Unadjusted UCell" = "#252525",
    "UCell residualized for RNA features" = "#D7301F")) +
  labs(title = "RNA-complexity sensitivity of selected UCell trajectories",
    subtitle = "Display-standardized curves; dashed lines mark BB1 and BB2",
    x = "Existing Bdev pseudotime (pt)",
    y = "Within-curve standardized window median", color = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", strip.text = element_text(size = 7))
ggsave(file.path(stage_dir, "raw_vs_complexity_residualized_curves.png"), p,
  width = 12, height = 9, dpi = 220)
ggsave(file.path(stage_dir, "raw_vs_complexity_residualized_curves.pdf"), p,
  width = 12, height = 9)

summary_rows <- vapply(seq_len(nrow(region_summary)), function(i) paste0("| ",
  region_summary$method[[i]], " | ",
  region_summary$dynamically_significant_sets_selected3[[i]], " | ",
  region_summary$supported_individual_peaks[[i]], " | ",
  region_summary$merged_regions[[i]], " | ",
  fmt_pct(region_summary$fraction_within_half_window_after_boundary[[i]]),
  " | ", fmt_pct(region_summary$fraction_within_full_window_after_boundary[[i]]),
  " |"), character(1))
med_raw_cor <- median(abs(adjustment_qc$raw_spearman_nFeature_RNA))
med_resid_cor <- median(abs(adjustment_qc$residual_spearman_nFeature_RNA))
med_r2 <- median(adjustment_qc$model_r_squared)
writeLines(c("# UCell RNA-complexity residualization diagnostic", "",
  "This nine-set sensitivity analysis asks whether the repeated placement of UCell regions immediately after BB1/BB2 is explained by RNA feature complexity.", "",
  "## Method", "",
  "- Cell-level native-compatible UCell scores were reconstructed from the validated maxRank=1500 ranking matrix.",
  "- For each gene set, a natural-spline model `UCell ~ ns(log1p(nFeature_RNA), df=4)` was fitted without pseudotime. Model residuals were analyzed.",
  "- Unadjusted and residual scores used identical 5%-wide, 1%-step, minimum-20-cell windows and 1,000 permutations.",
  "- Dynamic tests and positive-region localization match the original algorithm. BH values here are calculated among only the three selected sets in each GO group and are diagnostic, not replacements for the original full-GO q-values.", "",
  "## Boundary pile-up", "",
  "| Method | Dynamic sets | Supported peaks | Merged regions | Within half-window after BB | Within full window after BB |",
  "|---|---:|---:|---:|---:|---:|", summary_rows, "",
  paste0("Median absolute Spearman correlation with `nFeature_RNA` changed from ",
    signif(med_raw_cor, 3), " to ", signif(med_resid_cor, 3), "."),
  paste0("The complexity model's median R-squared across the nine sets was ",
    signif(med_r2, 3), "."),
  paste0("The unadjusted analysis reproduced stored window curves with maximum absolute difference ",
    format(max_curve_difference, scientific = TRUE, digits = 3),
    " and reproduced every dynamic empirical p-value exactly."), "",
  "Interpretation should distinguish persistence of trajectory dynamics from exact peak localization: residualization can remove genuine biology when RNA complexity changes with cell state, so this is a sensitivity analysis rather than a corrected primary result."),
  file.path(stage_dir, "REPORT.md"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
write_json(list(status = "complete", selected_sets = 9L,
  adjustment = "lm(UCell ~ ns(log1p(nFeature_RNA), df=4))",
  permutations = n_permutations, windows = nrow(windows),
  raw_curve_maximum_absolute_difference = max_curve_difference,
  raw_dynamic_p_values_reproduced_exactly = TRUE,
  median_absolute_raw_spearman_nFeature_RNA = med_raw_cor,
  median_absolute_residual_spearman_nFeature_RNA = med_resid_cor,
  bundle_unchanged = TRUE), file.path(stage_dir, "completion.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)

wb_close_bundle(bundle)
stopifnot(file.rename(stage_dir, output_dir))
message("Completed RNA-complexity residualization diagnostic at: ", output_dir)
