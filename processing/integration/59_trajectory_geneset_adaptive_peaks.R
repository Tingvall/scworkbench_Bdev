#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 59_trajectory_geneset_adaptive_peaks.R",
  "<bundle_dir> <inventory.tsv.gz> <memberships.tsv.gz>",
  "<bp|mf|cc|cytokines> <output_dir> <seed> <max_permutations>"))

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
inventory_path <- normalizePath(args[[2L]], mustWork = TRUE)
membership_path <- normalizePath(args[[3L]], mustWork = TRUE)
group <- args[[4L]]
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
base_seed <- as.integer(args[[6L]])
max_permutations <- as.integer(args[[7L]])
resume_stage <- identical(Sys.getenv("RESUME_STAGE", "0"), "1")
stopifnot(group %in% c("bp", "mf", "cc", "cytokines"),
  is.finite(base_seed), max_permutations %in% c(100L, 1000L, 10000L))
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir) && !resume_stage)
  stop("Incomplete output exists: ", stage_dir)
if (!dir.exists(stage_dir) && resume_stage)
  stop("RESUME_STAGE=1 but incomplete output does not exist: ", stage_dir)
if (!dir.exists(stage_dir)) dir.create(stage_dir, recursive = TRUE)
figure_dir <- file.path(stage_dir, "figures")
if (!dir.exists(figure_dir)) dir.create(figure_dir)
if (resume_stage) unlink(list.files(figure_dir, full.names = TRUE), recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
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
  as.data.table(read.delim(con, check.names = FALSE, stringsAsFactors = FALSE))
}
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ..., "\n")
  flush.console()
}
slugify <- function(x) gsub("(^_+|_+$)", "",
  tolower(gsub("[^A-Za-z0-9]+", "_", x)))

inventory <- read_tsv_gz(inventory_path)[analysis_group == group & eligible == TRUE]
members <- read_tsv_gz(membership_path)[analysis_group == group]
stopifnot(nrow(inventory) > 0L, !anyDuplicated(inventory$feature_set),
  setequal(inventory$feature_set, unique(members$feature_set)))
setorder(inventory, feature_set); setorder(members, feature_set, feature_id)

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
pt_only_qr <- qr(depth_design[, -ncol(depth_design), drop = FALSE])
partial_depth <- qr.resid(pt_only_qr, depth_z)

width_fraction <- 0.05; step_fraction <- 0.01; min_cells <- 30L
candidate_alpha <- 0.01; candidate_z <- qnorm(1 - candidate_alpha)
q_cutoff <- 0.01; prominence_cutoff <- 0.10
merge_gap_fraction <- width_fraction
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

permutation_stages <- c(100L, 1000L, 10000L)
permutation_stages <- permutation_stages[permutation_stages <= max_permutations]
stop_confidence_tail <- 0.001
stop_p_threshold <- q_cutoff
if (!resume_stage) {
workers_requested <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", "8")))
if (!is.finite(workers_requested) || workers_requested < 1L)
  stop("N_WORKERS must be a positive integer")
workers <- min(workers_requested, parallel::detectCores(), nrow(inventory))

log_event("Starting", group, "adaptive analysis for", nrow(inventory),
  "sets with", workers, "numeric-analysis workers")
raw_score_matrix <- matrix(NA_real_, nrow(inventory), length(pt),
  dimnames = list(inventory$feature_set, pt_cells))
for (ii in seq_len(nrow(inventory))) {
  feature_set <- inventory$feature_set[[ii]]
  fs_id <- feature_set
  genes <- members[feature_set == fs_id, feature_id]
  expression <- as.matrix(wb_matrix_read_features(bundle,
    "RNA_log1p_cp10k", genes, pt_cells, missing = "NA"))
  raw_score_matrix[ii, ] <- colMeans(expression, na.rm = TRUE)
  if (ii %% 100L == 0L || ii == nrow(inventory))
    log_event(group, "read expression for", ii, "of", nrow(inventory), "sets")
}

analyze_one <- function(ii) {
  feature_set <- inventory$feature_set[[ii]]
  raw_scores <- raw_score_matrix[ii, ]
  stopifnot(length(raw_scores) == length(pt), all(is.finite(raw_scores)))
  beta <- unname(qr.coef(depth_qr, raw_scores)[ncol(depth_design)])
  scores <- raw_scores - beta * depth_z
  raw_curve <- as.numeric(weight_matrix %*% raw_scores)
  curve <- as.numeric(weight_matrix %*% scores)
  population_sd <- sd(scores)
  null_mean <- mean(scores)
  raw_sd <- sd(raw_scores)
  raw_resid <- qr.resid(pt_only_qr, raw_scores)
  adj_resid <- qr.resid(pt_only_qr, scores)
  safe_cor <- function(x, y) if (sd(x) > 0 && sd(y) > 0) cor(x, y) else NA_real_
  depth_row <- data.table(analysis_group = group,
    feature_set = feature_set, depth_beta_per_sd_log1p_nCount_RNA = beta,
    conditional_depth_effect_score_sd = if (raw_sd > 0) abs(beta) / raw_sd else NA_real_,
    partial_depth_correlation_raw = safe_cor(raw_resid, partial_depth),
    partial_depth_correlation_adjusted = safe_cor(adj_resid, partial_depth))
  if (!is.finite(population_sd) || population_sd <= .Machine$double.eps) {
    curve_row <- data.table(analysis_group = group,
      feature_set = feature_set, windows, score_mean = curve,
      raw_score_mean = raw_curve, null_mean = null_mean, null_sd = NA_real_,
      z_score = NA_real_)
    adaptive_row <- data.table(analysis_group = group,
      feature_set = feature_set, candidate_peaks = 0L,
      permutations_used = 0L, strongest_peak_z = NA_real_,
      strongest_peak_exceedances = NA_integer_, stopped_early = TRUE,
      stop_reason = "zero_variance_score")
    return(list(curve = curve_row, candidates = NULL, depth = depth_row,
      adaptive = adaptive_row, beta = beta))
  }
  null_sd <- population_sd * sqrt((1 - windows$n_cells / length(scores)) /
    windows$n_cells)
  stopifnot(all(is.finite(null_sd)), all(null_sd > 0))
  z <- (curve - null_mean) / null_sd
  curve_row <- data.table(analysis_group = group,
    feature_set = feature_set, windows, score_mean = curve,
    raw_score_mean = raw_curve, null_mean = null_mean, null_sd = null_sd,
    z_score = z)
  peaks <- unlist(lapply(segments, function(idx) local_maxima(curve, idx)),
    use.names = FALSE)
  peaks <- peaks[z[peaks] >= candidate_z]
  if (!length(peaks)) {
    candidate_row <- NULL
    adaptive_row <- data.table(analysis_group = group,
      feature_set = feature_set, candidate_peaks = 0L,
      permutations_used = 0L, strongest_peak_z = NA_real_,
      strongest_peak_exceedances = NA_integer_, stopped_early = TRUE,
      stop_reason = "no_peak_above_candidate_z")
  } else {
    exceed <- integer(length(peaks)); used <- 0L; stopped <- FALSE
    set.seed(base_seed + ii - 1L)
    for (target in permutation_stages) {
      remaining <- target - used
      for (start in seq.int(1L, remaining, by = 100L)) {
        take <- min(100L, remaining - start + 1L)
        permuted <- replicate(take, sample(scores, length(scores), replace = FALSE))
        if (is.null(dim(permuted))) permuted <- matrix(permuted, ncol = 1L)
        null_curves <- weight_matrix %*% permuted
        null_z <- sweep(sweep(null_curves, 1L, null_mean, "-"), 1L,
          null_sd, "/")
        null_max <- apply(null_z, 2L, max)
        exceed <- exceed + vapply(z[peaks], function(v) sum(null_max >= v),
          integer(1))
      }
      used <- target
      strongest_k <- min(exceed)
      lower <- if (strongest_k > 0L)
        qbeta(stop_confidence_tail, strongest_k, used - strongest_k + 1L) else 0
      if (target < max(permutation_stages) && lower > stop_p_threshold) {
        stopped <- TRUE; break
      }
    }
    score_range <- diff(range(curve))
    candidate_row <- rbindlist(lapply(seq_along(peaks), function(jj) {
      pidx <- peaks[[jj]]
      prom <- peak_prominence(curve, pidx)
      data.table(
        analysis_group = group, feature_set = feature_set,
        set_name = inventory$set_name[[ii]],
        description = inventory$description[[ii]],
        peak_window = pidx, peak_pseudotime = windows$window_center[[pidx]],
        region_left = windows$window_left[[pidx]],
        region_right = windows$window_right[[pidx]],
        peak_score = curve[[pidx]], raw_peak_score = raw_curve[[pidx]],
        peak_z = z[[pidx]], peak_prominence = prom,
        relative_peak_prominence = prom / score_range,
        exceedances = exceed[[jj]], permutations_used = used,
        p_region = (1 + exceed[[jj]]) / (1 + used))
    }))
    adaptive_row <- data.table(analysis_group = group,
      feature_set = feature_set, candidate_peaks = length(peaks),
      permutations_used = used, strongest_peak_z = max(z[peaks]),
      strongest_peak_exceedances = min(exceed), stopped_early = stopped,
      stop_reason = if (stopped) "99.9pct_lower_bound_above_0.01" else
        "maximum_permutations_reached")
  }
  list(curve = curve_row, candidates = candidate_row, depth = depth_row,
    adaptive = adaptive_row, beta = beta)
}
analyze_one_safe <- function(ii) tryCatch(analyze_one(ii), error = function(e)
  structure(list(index = ii, feature_set = inventory$feature_set[[ii]],
    message = conditionMessage(e)), class = "analysis_error"))

analysis_results <- vector("list", nrow(inventory))
batch_size <- 50L
for (batch_start in seq.int(1L, nrow(inventory), by = batch_size)) {
  batch <- seq.int(batch_start, min(nrow(inventory),
    batch_start + batch_size - 1L))
  batch_results <- if (workers == 1L) lapply(batch, analyze_one_safe) else
    parallel::mclapply(batch, analyze_one_safe,
      mc.cores = min(workers, length(batch)), mc.preschedule = TRUE)
  failed <- which(vapply(batch_results, inherits, logical(1), "analysis_error"))
  if (length(failed)) {
    details <- vapply(batch_results[failed], function(x) paste0("index ",
      x$index, " (", x$feature_set, "): ", x$message), character(1))
    stop("Numeric analysis failed: ", paste(details, collapse = "; "))
  }
  analysis_results[batch] <- batch_results
  log_event(group, "completed numeric analysis for", max(batch), "of",
    nrow(inventory), "sets")
}
curves <- rbindlist(lapply(analysis_results, `[[`, "curve"))
depth_diagnostics <- rbindlist(lapply(analysis_results, `[[`, "depth"))
adaptive_summary <- rbindlist(lapply(analysis_results, `[[`, "adaptive"))
candidate_rows <- Filter(Negate(is.null),
  lapply(analysis_results, `[[`, "candidates"))
candidates <- if (length(candidate_rows)) rbindlist(candidate_rows) else data.table()
depth_beta_by_set <- setNames(vapply(analysis_results, `[[`, numeric(1),
  "beta"), inventory$feature_set)
rm(raw_score_matrix, analysis_results)
invisible(gc())
if (!nrow(candidates)) stop("No candidate peaks found for ", group)
candidates[, q_region := p.adjust(p_region, method = "BH")]
setorder(candidates, q_region, p_region, -peak_z, feature_set, peak_pseudotime)
candidates[, peak_id := sprintf(paste0(group, "_peak_%05d"), seq_len(.N))]
candidates[, significant := q_region < q_cutoff]
candidates[, passes_prominence := is.finite(relative_peak_prominence) &
  relative_peak_prominence >= prominence_cutoff]
candidates[, retained_for_reporting := significant & passes_prominence]
significant_individual <- candidates[significant == TRUE]
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
  final_regions[, region_id := sprintf(paste0(group, "_region_%05d"), seq_len(.N))]
} else final_regions <- data.table()

write_tsv_gz(curves, file.path(stage_dir, "window_curves.tsv.gz"))
write_tsv_gz(candidates, file.path(stage_dir, "all_candidate_peaks.tsv.gz"))
write_tsv(significant_individual,
  file.path(stage_dir, "significant_individual_peaks_unfiltered.tsv"))
write_tsv(retained,
  file.path(stage_dir, "prominence_retained_constituent_peaks.tsv"))
write_tsv(final_regions, file.path(stage_dir, "significant_positive_regions.tsv"))
write_tsv(depth_diagnostics,
  file.path(stage_dir, "rna_depth_adjustment_diagnostics.tsv"))
write_tsv(adaptive_summary,
  file.path(stage_dir, "adaptive_permutation_summary.tsv"))
} else {
  curves <- read_tsv_gz(file.path(stage_dir, "window_curves.tsv.gz"))
  candidates <- read_tsv_gz(file.path(stage_dir, "all_candidate_peaks.tsv.gz"))
  significant_individual <- fread(file.path(stage_dir,
    "significant_individual_peaks_unfiltered.tsv"))
  retained <- fread(file.path(stage_dir,
    "prominence_retained_constituent_peaks.tsv"))
  final_regions <- fread(file.path(stage_dir, "significant_positive_regions.tsv"))
  depth_diagnostics <- fread(file.path(stage_dir,
    "rna_depth_adjustment_diagnostics.tsv"))
  adaptive_summary <- fread(file.path(stage_dir,
    "adaptive_permutation_summary.tsv"))
  depth_beta_by_set <- setNames(
    depth_diagnostics$depth_beta_per_sd_log1p_nCount_RNA,
    depth_diagnostics$feature_set)
  log_event("Resuming", group, "from completed statistical tables")
}

significant_sets <- unique(final_regions$feature_set)
sets_to_plot <- head(significant_sets, 20L)
figure_index <- list()
base_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"))
for (ii in seq_along(sets_to_plot)) {
  feature_set <- sets_to_plot[[ii]]
  fs_id <- feature_set
  genes <- members[feature_set == fs_id, feature_id]
  expr <- as.matrix(wb_matrix_read_features(bundle, "RNA_log1p_cp10k",
    genes, cells, missing = "NA"))
  raw_score <- colMeans(expr, na.rm = TRUE)
  score <- raw_score - depth_beta_by_set[[feature_set]] * depth_z_all
  pos <- final_regions[feature_set == fs_id][order(region_left)]
  curve <- curves[feature_set == fs_id]
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
    labs(title = paste(inventory[feature_set == fs_id, set_name],
      "depth-adjusted score"),
      subtitle = "Red: significant region; greys: trajectory/non-trajectory",
      caption = "Darker red = higher adjusted score.")
  p_curve <- ggplot(curve, aes(window_center, score_mean)) +
    geom_rect(data = pos, aes(xmin = region_left, xmax = region_right,
      ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
      fill = "#D81B60", alpha = .13) +
    geom_line(aes(y = raw_score_mean), color = "#8B929A", linewidth = .65,
      linetype = "22") + geom_line(color = "#252525", linewidth = .8) +
    geom_point(data = pos, aes(peak_pseudotime, peak_score),
      inherit.aes = FALSE, color = "#D81B60", size = 2) + base_theme +
    labs(title = "Positive trajectory enrichment",
      subtitle = "Black: adjusted; dashed grey: raw; pink: retained",
      x = "Existing Bdev pseudotime (pt)", y = "Mean gene-set score",
      caption = "5% windows; 1% steps; group-wise q<0.01; prominence>=10%")
  filename <- paste0(sprintf("%04d", ii), "_", slugify(feature_set), ".png")
  ggsave(file.path(figure_dir, filename), p_umap | p_curve,
    width = 12.5, height = 5.6, dpi = 180, bg = "white")
  figure_index[[ii]] <- data.table(analysis_group = group,
    feature_set = feature_set,
    set_name = inventory[feature_set == fs_id, set_name],
    n_regions = nrow(pos), significant_cells = sum(mask),
    nonsignificant_trajectory_cells = sum(pt_defined & !mask),
    missing_pt_cells = sum(!pt_defined), total_cells = length(cells),
    figure = file.path("figures", filename))
  if (ii %% 25L == 0L || ii == length(sets_to_plot))
    log_event(group, "rendered", ii, "of", length(sets_to_plot), "figures")
}
figure_index <- if (length(figure_index)) rbindlist(figure_index) else data.table()
write_tsv(figure_index, file.path(stage_dir, "umap_figure_index.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
completion <- list(status = "complete", analysis_group = group,
  tested_sets = nrow(inventory),
  zero_variance_untestable_sets = sum(adaptive_summary$stop_reason ==
    "zero_variance_score"), candidate_peaks = nrow(candidates),
  significant_individual_peaks_before_prominence = nrow(significant_individual),
  prominence_retained_constituent_peaks = nrow(retained),
  significant_positive_regions = nrow(final_regions),
  sets_with_significant_positive_regions = length(significant_sets),
  figures_generated_for_top_significant_sets = length(sets_to_plot),
  multiple_testing = paste("BH within", group, "candidate peaks only"),
  q_cutoff = q_cutoff, relative_prominence_cutoff = prominence_cutoff,
  merge_gap_fraction_strictly_less_than = merge_gap_fraction,
  adaptive_permutations = list(stages = as.list(permutation_stages),
    maximum = max_permutations, stop_confidence = 1 - stop_confidence_tail,
    stop_rule = "stop if binomial lower confidence bound for strongest peak exceeds 0.01"),
  depth_adjustment = "subtract conditional scaled log1p(nCount_RNA) coefficient after ns(pt,df=8)",
  pseudotime_source = "exact Bdev metadata pt", defined_pt_cells = length(pt),
  missing_pt_cells = sum(!pt_defined), promoted = FALSE,
  bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = fingerprint_after, bundle_unchanged = TRUE)
write_json(completion, file.path(stage_dir, "completion.json"))

top <- head(final_regions, 25L)
table_lines <- if (nrow(top)) c(
  "| Gene set | Peak pt | Peak z | Prominence | p | q |",
  "|---|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(top)), function(i) paste0("| ", top$set_name[[i]],
    " | ", signif(top$peak_pseudotime[[i]], 5), " | ",
    signif(top$peak_z[[i]], 4), " | ",
    signif(top$relative_peak_prominence[[i]], 4), " | ",
    signif(top$p_region[[i]], 4), " | ", signif(top$q_region[[i]], 4),
    " |"), character(1))) else "No regions passed."
writeLines(c(paste0("# ", group, " adaptive trajectory enrichment"), "",
  paste0("- Tested gene sets: ", nrow(inventory)),
  paste0("- Zero-variance, untestable gene sets: ",
    sum(adaptive_summary$stop_reason == "zero_variance_score")),
  paste0("- Candidate peaks: ", nrow(candidates)),
  paste0("- Significant individual peaks before prominence: ",
    nrow(significant_individual)),
  paste0("- Final regions: ", nrow(final_regions), " across ",
    length(significant_sets), " gene sets"),
  paste0("- Figures: top ", length(sets_to_plot),
    " gene sets ranked by best q, p, then peak z"),
  paste0("- Multiple testing: BH separately within ", group, " (q < 0.01)"),
  "- Positive enrichment only; 10% relative prominence; strict <5% gap merge.",
  "", "## Top final regions", "", table_lines),
  file.path(stage_dir, "REPORT.md"))

wb_close_bundle(bundle)
stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed", group, "adaptive report at", output_dir)
