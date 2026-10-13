#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 72_diagnose_cytokine_ucell_z_pt_shift.R",
  "<bundle_dir> <native_ucell_report> <curve_comparison_dir>",
  "<z_cytokine_dir> <output_dir> <max_lag_windows>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
native_dir <- normalizePath(args[[2L]], mustWork = TRUE)
comparison_dir <- normalizePath(args[[3L]], mustWork = TRUE)
z_dir <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
max_lag <- as.integer(args[[6L]])
stopifnot(is.finite(max_lag), max_lag >= 1L)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(jsonlite)
})

read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
short_label <- function(x) sub(
  "^parse_cytokine_rest_only_refinement_all_top100_v1\\.", "", x)

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
pt_table <- fread(file.path(native_dir, "source_pseudotime.tsv"))
setnames(pt_table, c("cell_id", "pt")); setorder(pt_table, cell_id)
metadata <- as.data.table(wb_meta(bundle, cells = pt_table$cell_id))
stopifnot(identical(as.character(metadata$cell_id), pt_table$cell_id))

scores <- read_tsv_gz(file.path(native_dir, "ucell_cell_scores.tsv.gz"))[
  cell_id %chin% pt_table$cell_id, .(feature_set = key, cell_id,
    raw_ucell = value)]
wide <- dcast(scores, cell_id ~ feature_set, value.var = "raw_ucell")
setorder(wide, cell_id)
stopifnot(identical(wide$cell_id, pt_table$cell_id))
feature_sets <- sort(setdiff(names(wide), "cell_id"))
raw <- t(as.matrix(wide[, ..feature_sets])); rownames(raw) <- feature_sets
complexity_x <- log1p(as.numeric(metadata$nFeature_RNA))
residual <- matrix(NA_real_, nrow(raw), ncol(raw), dimnames = dimnames(raw))
for (fs in feature_sets) residual[fs, ] <- residuals(lm(
  as.numeric(raw[fs, ]) ~ splines::ns(complexity_x, df = 4L)))

pt <- pt_table$pt; pt_limits <- range(pt); pt_span <- diff(pt_limits)
centers <- pt_limits[[1L]] + seq(0, 1, by = .01) * pt_span
window_width <- .05 * pt_span
windows <- data.table(grid_window_index = seq_along(centers),
  window_center = centers,
  window_left = pmax(pt_limits[[1L]], centers - window_width / 2),
  window_right = pmin(pt_limits[[2L]], centers + window_width / 2))
idx <- lapply(seq_len(nrow(windows)), function(i)
  which(pt >= windows$window_left[[i]] & pt <= windows$window_right[[i]]))
windows[, n_cells := lengths(idx)]
keep <- windows$n_cells >= 20L; windows <- windows[keep]; idx <- idx[keep]
windows[, window_index := seq_len(.N)]
stopifnot(nrow(windows) == 97L)

ucell_curves <- rbindlist(lapply(feature_sets, function(fs) data.table(
  feature_set = fs, windows,
  raw_median = vapply(idx, function(ii) median(raw[fs, ii]), numeric(1)),
  raw_mean = vapply(idx, function(ii) mean(raw[fs, ii]), numeric(1)),
  residual_median = vapply(idx, function(ii)
    median(residual[fs, ii]), numeric(1)),
  residual_mean = vapply(idx, function(ii)
    mean(residual[fs, ii]), numeric(1)))))
z <- read_tsv_gz(file.path(z_dir, "window_curves.tsv.gz"))
stopifnot(setequal(unique(z$feature_set), feature_sets),
  nrow(z) == 38L * 97L)

# Positive shift means the UCell curve must be moved toward later pseudotime
# to align with the expression curve. Ties prefer the smallest absolute shift.
best_lag <- function(u, target) {
  candidates <- rbindlist(lapply(-max_lag:max_lag, function(s) {
    ii <- seq_along(u); jj <- ii + s
    keep <- jj >= 1L & jj <= length(target)
    ux <- u[ii[keep]]; tx <- target[jj[keep]]
    r <- if (sd(ux) > .Machine$double.eps &&
      sd(tx) > .Machine$double.eps) cor(ux, tx, method = "pearson") else NA_real_
    data.table(shift_windows = s,
      correlation = r)
  }))
  candidates <- candidates[is.finite(correlation)]
  if (!nrow(candidates)) return(data.table(shift_windows = NA_integer_,
    correlation = NA_real_))
  candidates[, absolute_shift := abs(shift_windows)]
  setorder(candidates, -correlation, absolute_shift)
  candidates[1L]
}

merged <- merge(ucell_curves, z[, .(feature_set, grid_window_index,
    expression_mean = score_mean, expression_z = z_score)],
  by = c("feature_set", "grid_window_index"))
methods <- c("raw_median", "raw_mean", "residual_median", "residual_mean")
targets <- c("expression_z", "expression_mean")
lag_results <- rbindlist(lapply(feature_sets, function(fs) {
  d <- merged[feature_set == fs][order(window_index)]
  rbindlist(lapply(methods, function(m) rbindlist(lapply(targets, function(t) {
    best <- best_lag(d[[m]], d[[t]])
    no_shift <- if (sd(d[[m]]) > .Machine$double.eps &&
      sd(d[[t]]) > .Machine$double.eps) cor(d[[m]], d[[t]]) else NA_real_
    data.table(feature_set = fs, short_name = short_label(fs),
      ucell_summary = m, expression_target = t,
      no_shift_correlation = no_shift,
      best_shift_windows = best$shift_windows,
      best_shift_pt = best$shift_windows * .01 * pt_span,
      best_shift_correlation = best$correlation)
  }))))
}))
lag_summary <- lag_results[, .(
  cytokines = .N,
  cytokines_with_defined_lag = sum(is.finite(best_shift_windows)),
  median_shift_windows = as.numeric(median(best_shift_windows, na.rm = TRUE)),
  q25_shift_windows = unname(quantile(best_shift_windows, .25, na.rm = TRUE)),
  q75_shift_windows = unname(quantile(best_shift_windows, .75, na.rm = TRUE)),
  median_shift_pt = median(best_shift_pt, na.rm = TRUE),
  median_absolute_shift_pt = median(abs(best_shift_pt), na.rm = TRUE),
  fraction_zero_shift = mean(best_shift_windows == 0L, na.rm = TRUE),
  fraction_ucell_shifted_later = mean(best_shift_windows > 0L, na.rm = TRUE),
  median_no_shift_correlation = median(no_shift_correlation, na.rm = TRUE),
  median_best_shift_correlation = median(best_shift_correlation, na.rm = TRUE)),
  by = .(ucell_summary, expression_target)]
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA")
write_tsv(lag_results, file.path(stage_dir,
  "per_cytokine_best_lag.tsv"))
write_tsv(lag_summary, file.path(stage_dir, "lag_summary.tsv"))

# Verify both analyses use the same exact window grid. This excludes plotting
# offsets or different window centers as the cause of the apparent shift.
grid <- merge(unique(ucell_curves[, .(grid_window_index,
    ucell_center = window_center, ucell_left = window_left,
    ucell_right = window_right, ucell_n = n_cells)]),
  unique(z[, .(grid_window_index, z_center = window_center,
    z_left = window_left, z_right = window_right, z_n = n_cells)]),
  by = "grid_window_index")
grid_qc <- data.table(windows = nrow(grid),
  maximum_center_difference = max(abs(grid$ucell_center - grid$z_center)),
  maximum_left_difference = max(abs(grid$ucell_left - grid$z_left)),
  maximum_right_difference = max(abs(grid$ucell_right - grid$z_right)),
  cell_counts_identical = identical(grid$ucell_n, grid$z_n),
  pt_step = .01 * pt_span, window_width = .05 * pt_span)
write_tsv(grid_qc, file.path(stage_dir, "window_grid_qc.tsv"))
stopifnot(grid_qc$maximum_center_difference < 1e-10,
  grid_qc$maximum_left_difference < 1e-10,
  grid_qc$maximum_right_difference < 1e-10,
  grid_qc$cell_counts_identical)

row_for <- function(m, t = "expression_z")
  lag_summary[ucell_summary == m & expression_target == t]
fmt <- function(x) signif(x, 3)
raw_med <- row_for("raw_median")
raw_mean <- row_for("raw_mean")
resid_med <- row_for("residual_median")
resid_mean <- row_for("residual_mean")
writeLines(c("# Cytokine UCell versus expression-z pseudotime-shift diagnosis", "",
  "Positive lag means the UCell curve must be shifted toward later pseudotime to align with the expression curve. Lags were searched over the same 97-window grid.", "",
  "## Results against the displayed expression z-score", "",
  "| UCell curve | Median signed shift (pt) | Median absolute shift (pt) | Shifted later | Median correlation: unshifted -> aligned |",
  "|---|---:|---:|---:|---:|",
  paste0("| Raw window median | ", fmt(raw_med$median_shift_pt), " | ",
    fmt(raw_med$median_absolute_shift_pt), " | ",
    round(100 * raw_med$fraction_ucell_shifted_later, 1), "% | ",
    fmt(raw_med$median_no_shift_correlation), " -> ",
    fmt(raw_med$median_best_shift_correlation), " |"),
  paste0("| Raw window mean | ", fmt(raw_mean$median_shift_pt), " | ",
    fmt(raw_mean$median_absolute_shift_pt), " | ",
    round(100 * raw_mean$fraction_ucell_shifted_later, 1), "% | ",
    fmt(raw_mean$median_no_shift_correlation), " -> ",
    fmt(raw_mean$median_best_shift_correlation), " |"),
  paste0("| Residual window median | ", fmt(resid_med$median_shift_pt), " | ",
    fmt(resid_med$median_absolute_shift_pt), " | ",
    round(100 * resid_med$fraction_ucell_shifted_later, 1), "% | ",
    fmt(resid_med$median_no_shift_correlation), " -> ",
    fmt(resid_med$median_best_shift_correlation), " |"),
  paste0("| Residual window mean | ", fmt(resid_mean$median_shift_pt), " | ",
    fmt(resid_mean$median_absolute_shift_pt), " | ",
    round(100 * resid_mean$fraction_ucell_shifted_later, 1), "% | ",
    fmt(resid_mean$median_no_shift_correlation), " -> ",
    fmt(resid_mean$median_best_shift_correlation), " |"), "",
  "## Technical interpretation", "",
  "- The UCell and expression analyses have identical window centers, bounds, cells, and pseudotime values, so the shift is not a plotting or grid-registration error.",
  "- UCell is a nonlinear within-cell rank statistic truncated at maxRank=1500; the expression analysis averages depth-adjusted, individually standardized log-expression genes. Their maxima need not occur in the same cells or windows.",
  "- UCell inference uses the window median, whereas the expression trajectory uses the window mean. Comparing mean-UCell with median-UCell quantifies this aggregation contribution.",
  "- The displayed expression z-score additionally divides the window mean by a window-specific null SD that depends on cell count. Results against the unstandardized expression window mean are included in the TSV to identify any display-z contribution.",
  "- Complexity residualization is a sensitivity analysis and may remove biological state variation correlated with RNA complexity."),
  file.path(stage_dir, "REPORT.md"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
write_json(list(status = "complete", cytokine_sets = 38L,
  lag_search_windows = c(-max_lag, max_lag),
  lag_step_pt = .01 * pt_span, window_grid_identical = TRUE,
  bundle_unchanged = TRUE), file.path(stage_dir, "completion.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)
wb_close_bundle(bundle)
stopifnot(file.rename(stage_dir, output_dir))
message("Completed cytokine UCell/z shift diagnosis at: ", output_dir)
