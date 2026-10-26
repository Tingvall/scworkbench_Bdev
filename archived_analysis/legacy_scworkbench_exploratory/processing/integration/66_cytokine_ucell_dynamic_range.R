#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 66_cytokine_ucell_dynamic_range.R",
  "<native_ucell_report> <output_dir> <seed> <permutations>"))
source_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
base_seed <- as.integer(args[[3L]])
n_permutations <- as.integer(args[[4L]])
stopifnot(is.finite(base_seed), identical(n_permutations, 1000L))
complexity_residualize <- identical(tolower(Sys.getenv(
  "CYTOKINE_UCELL_COMPLEXITY_RESIDUALIZE", "false")), "true")
bundle_dir <- Sys.getenv("CYTOKINE_UCELL_BUNDLE", "")
if (complexity_residualize)
  bundle_dir <- normalizePath(bundle_dir, mustWork = TRUE)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(jsonlite)
  library(matrixStats)
})

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA")
write_tsv_gz <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA", compress = "gzip")
read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
short_label <- function(x) sub(
  "^parse_cytokine_rest_only_refinement_all_top100_v1\\.", "", x)
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ...,
    "\n")
  flush.console()
}

completion <- read_json(file.path(source_dir, "completion.json"),
  simplifyVector = TRUE)
stopifnot(identical(completion$status, "complete"),
  completion$passing_signatures == 38L,
  completion$defined_pt_cells == 12642L,
  completion$missing_pt_cells == 14081L,
  completion$max_rank == 1500L)

pt_table <- fread(file.path(source_dir, "source_pseudotime.tsv"))
setnames(pt_table, c("cell_id", "pt"))
stopifnot(nrow(pt_table) == 12642L, !anyDuplicated(pt_table$cell_id),
  all(is.finite(pt_table$pt)))
setorder(pt_table, cell_id)
pt <- pt_table$pt

log_event("Reading retained native scWorkbench UCell scores")
scores_long <- read_tsv_gz(file.path(source_dir, "ucell_cell_scores.tsv.gz"))[
  , .(key, cell_id, value)]
scores_long <- scores_long[cell_id %chin% pt_table$cell_id]
stopifnot(uniqueN(scores_long$key) == 38L,
  uniqueN(scores_long$cell_id) == 12642L,
  nrow(scores_long) == 38L * 12642L,
  !anyDuplicated(scores_long, by = c("key", "cell_id")),
  all(is.finite(scores_long$value)), all(scores_long$value >= 0))
score_wide <- dcast(scores_long, cell_id ~ key, value.var = "value")
setorder(score_wide, cell_id)
stopifnot(identical(score_wide$cell_id, pt_table$cell_id))
feature_sets <- sort(setdiff(names(score_wide), "cell_id"))
score_matrix <- t(as.matrix(score_wide[, ..feature_sets]))
rownames(score_matrix) <- feature_sets
rm(scores_long, score_wide); invisible(gc())

validation_score_matrix <- score_matrix
bundle <- NULL
fingerprint_before <- NA_character_
fingerprint_after <- NA_character_
adjustment_qc <- NULL
if (complexity_residualize) {
  log_event("Residualizing each cytokine UCell score for RNA feature complexity")
  bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
  on.exit(if (!is.null(bundle))
    try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
  metadata <- as.data.table(wb_meta(bundle, cells = pt_table$cell_id))
  stopifnot(identical(as.character(metadata$cell_id), pt_table$cell_id),
    all(is.finite(metadata$nFeature_RNA)))
  complexity_x <- log1p(as.numeric(metadata$nFeature_RNA))
  residual_matrix <- matrix(NA_real_, nrow(score_matrix), ncol(score_matrix),
    dimnames = dimnames(score_matrix))
  adjustment_qc <- rbindlist(lapply(feature_sets, function(fs) {
    y <- as.numeric(score_matrix[fs, ])
    fit <- lm(y ~ splines::ns(complexity_x, df = 4L))
    adjusted <- residuals(fit)
    residual_matrix[fs, ] <<- adjusted
    data.table(feature_set = fs, short_name = short_label(fs),
      model = "lm(UCell ~ ns(log1p(nFeature_RNA), df=4))",
      model_r_squared = summary(fit)$r.squared,
      raw_spearman_nFeature_RNA = cor(y, metadata$nFeature_RNA,
        method = "spearman"),
      residual_spearman_nFeature_RNA = cor(adjusted,
        metadata$nFeature_RNA, method = "spearman"))
  }))
  score_matrix <- residual_matrix
  rm(residual_matrix, metadata); invisible(gc())
}

width_fraction <- 0.05
step_fraction <- 0.01
min_cells <- 20L
set_q_cutoff <- 0.01
region_fwer_cutoff <- 0.01
merge_gap_fraction <- width_fraction
pt_limits <- range(pt)
pt_span <- diff(pt_limits)
centers <- pt_limits[[1L]] + seq(0, 1, by = step_fraction) * pt_span
window_width <- width_fraction * pt_span
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
stopifnot(nrow(windows) == 97L, min(windows$n_cells) >= min_cells)

native <- fread(file.path(source_dir, "sliding_window_summary.tsv"))[
  abs(width_fraction - 0.05) < 1e-12 & eligible_window == TRUE,
  .(feature_set, grid_window_index = window_index,
    native_window_center = pseudotime, native_n_cells = n_cells,
    native_median = score_median)]
stopifnot(nrow(native) == 38L * 97L,
  setequal(native$feature_set, feature_sets))

# A tied plateau contributes one midpoint. The first and last plateaus are
# eligible extrema, so monotonic high-to-low and low-to-high profiles count as
# dynamic and can localize enrichment at a trajectory endpoint.
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

analyze_one <- function(ii) {
  fs <- feature_sets[[ii]]
  scores <- as.numeric(score_matrix[ii, ])
  baseline <- median(scores)
  observed <- vapply(membership_index,
    function(idx) median(scores[idx]), numeric(1))
  effect <- observed - baseline
  observed_range <- max(observed) - min(observed)
  peak_index <- local_maxima_plateau(observed)
  peak_index <- peak_index[effect[peak_index] > .Machine$double.eps]

  null_range <- numeric(n_permutations)
  null_max_positive <- numeric(n_permutations)
  batch_size <- 100L
  set.seed(base_seed + ii - 1L)
  for (start in seq.int(1L, n_permutations, by = batch_size)) {
    take <- min(batch_size, n_permutations - start + 1L)
    permuted <- replicate(take, sample(scores, length(scores), replace = FALSE))
    if (is.null(dim(permuted))) permuted <- matrix(permuted, ncol = 1L)
    batch_min <- rep(Inf, take)
    batch_max <- rep(-Inf, take)
    for (idx in membership_index) {
      window_median <- matrixStats::colMedians(permuted[idx, , drop = FALSE])
      batch_min <- pmin(batch_min, window_median)
      batch_max <- pmax(batch_max, window_median)
    }
    pos <- start:(start + take - 1L)
    null_range[pos] <- batch_max - batch_min
    null_max_positive[pos] <- batch_max - baseline
  }

  curve <- data.table(analysis_group = "cytokines", feature_set = fs,
    windows, ucell_median = observed, trajectory_median = baseline,
    ucell_effect = effect)
  test <- data.table(analysis_group = "cytokines", feature_set = fs,
    short_name = short_label(fs), observed_dynamic_range = observed_range,
    null_range_mean = mean(null_range), null_range_sd = sd(null_range),
    null_range_q95 = unname(quantile(null_range, .95)),
    null_range_q99 = unname(quantile(null_range, .99)),
    dynamic_exceedances = sum(null_range >= observed_range),
    permutations_used = n_permutations,
    p_dynamic = (1 + sum(null_range >= observed_range)) /
      (1 + n_permutations), candidate_positive_peaks = length(peak_index),
    positive_null_max_q99 = unname(quantile(null_max_positive, .99)))

  candidates <- if (length(peak_index)) data.table(
    analysis_group = "cytokines", feature_set = fs,
    peak_window = peak_index,
    peak_pseudotime = windows$window_center[peak_index],
    region_left = windows$window_left[peak_index],
    region_right = windows$window_right[peak_index],
    n_cells = windows$n_cells[peak_index],
    peak_ucell_median = observed[peak_index], trajectory_median = baseline,
    peak_ucell_effect = effect[peak_index],
    region_exceedances = vapply(effect[peak_index],
      function(x) sum(null_max_positive >= x), integer(1)),
    permutations_used = n_permutations) else NULL
  if (!is.null(candidates)) candidates[, p_region_fwer :=
    (1 + region_exceedances) / (1 + permutations_used)]
  list(curve = curve, test = test, candidates = candidates)
}

workers_requested <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", "8")))
if (!is.finite(workers_requested) || workers_requested < 1L)
  stop("N_WORKERS must be a positive integer")
workers <- min(workers_requested, parallel::detectCores(), length(feature_sets))
log_event("Starting", n_permutations, "dynamic-range permutations for",
  length(feature_sets), "cytokine sets with", workers, "workers")
results <- parallel::mclapply(seq_along(feature_sets), analyze_one,
  mc.cores = workers, mc.preschedule = TRUE)
if (any(vapply(results, inherits, logical(1), "try-error")))
  stop("At least one cytokine permutation worker failed")
curves <- rbindlist(lapply(results, `[[`, "curve"))
tests <- rbindlist(lapply(results, `[[`, "test"))
candidate_rows <- Filter(Negate(is.null), lapply(results, `[[`, "candidates"))
candidates <- if (length(candidate_rows)) rbindlist(candidate_rows) else
  data.table()

validation_curves <- rbindlist(lapply(seq_along(feature_sets), function(ii)
  data.table(feature_set = feature_sets[[ii]], windows,
    native_reconstructed = vapply(membership_index, function(idx)
      median(validation_score_matrix[ii, idx]), numeric(1)))))
native_check <- merge(validation_curves[, .(feature_set, grid_window_index,
  window_center, n_cells, native_reconstructed)], native,
  by = c("feature_set", "grid_window_index"), all = TRUE)
stopifnot(nrow(native_check) == 38L * 97L, !anyNA(native_check),
  max(abs(native_check$window_center - native_check$native_window_center)) < 1e-10,
  identical(native_check$n_cells, native_check$native_n_cells),
  max(abs(native_check$native_reconstructed - native_check$native_median)) < 1e-12)

tests[, q_dynamic := p.adjust(p_dynamic, method = "BH")]
tests[, dynamic_significant := q_dynamic < set_q_cutoff]
setorder(tests, q_dynamic, p_dynamic, -observed_dynamic_range, feature_set)
candidates[tests[, .(feature_set, p_dynamic, q_dynamic,
  dynamic_significant)], on = "feature_set",
  `:=`(p_dynamic = i.p_dynamic, q_dynamic = i.q_dynamic,
    dynamic_significant = i.dynamic_significant)]
candidates[, region_supported := dynamic_significant &
  p_region_fwer < region_fwer_cutoff]
setorder(candidates, q_dynamic, p_region_fwer, -peak_ucell_effect,
  feature_set, peak_pseudotime)
candidates[, peak_id := sprintf("ucell_dynamic_cytokine_peak_%04d", seq_len(.N))]
supported <- candidates[region_supported == TRUE]

merge_gap_pt <- merge_gap_fraction * pt_span
merge_peaks <- function(x) {
  x <- copy(x[order(region_left, region_right)])
  gid <- integer(nrow(x)); gid[[1L]] <- 1L
  if (nrow(x) > 1L) for (i in 2:nrow(x)) {
    gap <- x$region_left[[i]] - x$region_right[[i - 1L]]
    gid[[i]] <- gid[[i - 1L]] + as.integer(gap >= merge_gap_pt)
  }
  rbindlist(lapply(split(seq_len(nrow(x)), gid), function(idx) {
    g <- x[idx]
    representative <- order(g$p_region_fwer, -g$peak_ucell_effect)[[1L]]
    out <- copy(g[representative])
    out[, `:=`(region_left = min(g$region_left),
      region_right = max(g$region_right),
      n_constituent_peaks = nrow(g),
      constituent_peak_pts = paste(signif(g$peak_pseudotime, 7), collapse = ";"),
      constituent_peak_ids = paste(g$peak_id, collapse = ";"))]
    out
  }))
}
if (nrow(supported)) {
  final_regions <- rbindlist(lapply(split(seq_len(nrow(supported)),
    supported$feature_set), function(idx) merge_peaks(supported[idx])))
  setorder(final_regions, q_dynamic, p_region_fwer, -peak_ucell_effect,
    feature_set, region_left)
  final_regions[, region_id := sprintf("ucell_dynamic_cytokine_region_%04d",
    seq_len(.N))]
} else final_regions <- data.table()

region_counts <- if (nrow(final_regions)) final_regions[, .(
  n_positive_regions = .N,
  best_region_fwer_p = min(p_region_fwer),
  positive_region_peak_pts = paste(signif(peak_pseudotime[order(region_left)],
    6), collapse = ";")), by = feature_set] else data.table(
      feature_set = character(), n_positive_regions = integer(),
      best_region_fwer_p = numeric(), positive_region_peak_pts = character())
significant_sets <- merge(tests[dynamic_significant == TRUE], region_counts,
  by = "feature_set", all.x = TRUE, sort = FALSE)
significant_sets[is.na(n_positive_regions), n_positive_regions := 0L]
setorder(significant_sets, q_dynamic, p_dynamic, -observed_dynamic_range,
  feature_set)

write_tsv_gz(curves, file.path(stage_dir, "window_curves.tsv.gz"))
write_tsv(tests, file.path(stage_dir, "gene_set_dynamic_tests.tsv"))
write_tsv_gz(candidates,
  file.path(stage_dir, "all_candidate_positive_peaks.tsv.gz"))
write_tsv(supported,
  file.path(stage_dir, "permutation_supported_positive_peaks.tsv"))
write_tsv(final_regions,
  file.path(stage_dir, "significant_positive_regions.tsv"))
write_tsv(significant_sets,
  file.path(stage_dir, "significant_gene_sets.tsv"))
if (complexity_residualize)
  write_tsv(adjustment_qc,
    file.path(stage_dir, "complexity_adjustment_qc.tsv"))

fmt_set <- function(i) paste0("| ", significant_sets$short_name[[i]], " | ",
  signif(significant_sets$observed_dynamic_range[[i]], 5), " | ",
  signif(significant_sets$p_dynamic[[i]], 4), " | ",
  signif(significant_sets$q_dynamic[[i]], 4), " | ",
  significant_sets$n_positive_regions[[i]], " | ",
  ifelse(is.na(significant_sets$positive_region_peak_pts[[i]]), "—",
    significant_sets$positive_region_peak_pts[[i]]), " |")
set_rows <- if (nrow(significant_sets))
  vapply(seq_len(nrow(significant_sets)), fmt_set, character(1)) else
  "| None | NA | NA | NA | 0 | — |"
score_description <- if (complexity_residualize)
  "RNA-feature-residualized native scWorkbench UCell on RNA counts (`max_rank=1500`), using `UCell ~ ns(log1p(nFeature_RNA), df=4)` without pseudotime." else
  "Native scWorkbench UCell on RNA counts (`max_rank=1500`, no smoothing)."
report <- c("# Cytokine UCell trajectory-dynamics permutation analysis", "",
  paste0("- Dynamically significant gene sets: ", nrow(significant_sets),
    " / 38"),
  paste0("- Permutation-supported positive peaks: ", nrow(supported)),
  paste0("- Final merged positive regions: ", nrow(final_regions)),
  paste0("- Score: ", score_description),
  "- Dynamic statistic: peak-to-trough range of the 97 window-median UCell scores; endpoints are eligible extrema.",
  "- Dynamic null: shuffle the same cells' UCell scores across fixed pseudotime positions and recalculate the full trajectory range.",
  "- Set inference: one permutation p-value per cytokine; BH across 38 cytokines; q < 0.01.",
  "- Positive-region localization: all positive local maxima, without a score cutoff, are compared with the permutation maximum positive deviation. Peaks require within-set trajectory-wide p < 0.01 and a significant parent set.",
  "- Windows: 5% width, 1% step, minimum 20 cells; 97 retained windows.",
  "- Supported peak windows merge only when the uncovered gap is strictly less than one 5% window.",
  "- UCell scores were not rescaled for inference.", "",
  "## Dynamically significant gene sets", "",
  "| Gene set | UCell range | Dynamic p | Dynamic q | Positive regions | Peak pseudotime(s) |",
  "|---|---:|---:|---:|---:|---|", set_rows)
writeLines(report, file.path(stage_dir, "REPORT.md"))

completion_record <- list(status = "complete",
  analysis = if (complexity_residualize)
    "RNA-feature-residualized native-UCell trajectory-wide dynamic-range permutation" else
    "native-UCell trajectory-wide dynamic-range permutation",
  tested_gene_sets = 38L,
  significant_dynamic_gene_sets = nrow(significant_sets),
  permutation_supported_positive_peaks = nrow(supported),
  significant_positive_regions = nrow(final_regions),
  score_source = if (complexity_residualize)
    "residual of native scWorkbench UCell RNA_counts max_rank=1500 from ns(log1p(nFeature_RNA), df=4)" else
    "native scWorkbench UCell RNA_counts max_rank=1500",
  complexity_residualization_enabled = complexity_residualize,
  complexity_residualization_model = if (complexity_residualize)
    "lm(UCell ~ ns(log1p(nFeature_RNA), df=4))" else NULL,
  score_rescaled_for_inference = FALSE,
  window_statistic = if (complexity_residualize)
    "median RNA-feature-residualized UCell" else "median UCell",
  dynamic_test_statistic = "maximum window median minus minimum window median",
  endpoint_extrema_included = TRUE,
  permutations = n_permutations,
  set_multiple_testing = "BH across 38 cytokines",
  set_q_cutoff = set_q_cutoff,
  positive_region_test = "local positive peak versus permutation maximum positive deviation",
  positive_region_fwer_p_cutoff = region_fwer_cutoff,
  arbitrary_candidate_score_cutoff = FALSE,
  window_width_fraction = width_fraction,
  window_step_fraction = step_fraction,
  minimum_window_cells = min_cells,
  retained_windows = nrow(windows),
  merge_gap_fraction_strictly_less_than = merge_gap_fraction,
  defined_pt_cells = nrow(pt_table),
  missing_pt_cells_preserved = completion$missing_pt_cells,
  native_window_medians_reproduced = TRUE,
  random_seed = base_seed, workers = workers,
  source_ucell_report = source_dir, promoted = FALSE)
if (complexity_residualize) {
  fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
  stopifnot(identical(fingerprint_before, fingerprint_after))
  completion_record$median_absolute_raw_spearman_nFeature_RNA <-
    median(abs(adjustment_qc$raw_spearman_nFeature_RNA))
  completion_record$median_absolute_residual_spearman_nFeature_RNA <-
    median(abs(adjustment_qc$residual_spearman_nFeature_RNA))
  completion_record$bundle_fingerprint_before <- fingerprint_before
  completion_record$bundle_fingerprint_after <- fingerprint_after
  completion_record$bundle_unchanged <- TRUE
  wb_close_bundle(bundle)
  bundle <- NULL
}
write_json(completion_record,
  file.path(stage_dir, "completion.json"), pretty = TRUE, auto_unbox = TRUE,
  digits = NA)

stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed cytokine UCell dynamic-range report at", output_dir)
