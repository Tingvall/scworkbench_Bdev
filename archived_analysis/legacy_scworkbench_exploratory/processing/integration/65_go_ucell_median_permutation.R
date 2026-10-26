#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 65_go_ucell_median_permutation.R",
  "<bundle_dir> <inventory.tsv.gz> <memberships.tsv.gz>",
  "<native_cytokine_ucell_report> <output_dir> <seed> <permutations>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
inventory_path <- normalizePath(args[[2L]], mustWork = TRUE)
membership_path <- normalizePath(args[[3L]], mustWork = TRUE)
cytokine_report <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
base_seed <- as.integer(args[[6L]])
n_permutations <- as.integer(args[[7L]])
stopifnot(is.finite(base_seed), identical(n_permutations, 1000L))
inference_mode <- Sys.getenv("GO_UCELL_INFERENCE_MODE", "candidate_peak_bh")
stopifnot(inference_mode %in% c("candidate_peak_bh", "dynamic_range"))
complexity_residualization <- identical(tolower(Sys.getenv(
  "GO_UCELL_COMPLEXITY_RESIDUALIZE", "false")), "true")
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(UCell)
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
write_tsv_gz <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA", compress = "gzip")
write_json_file <- function(x, path) write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ...,
    "\n")
  flush.console()
}

inventory_all <- read_tsv_gz(inventory_path)[eligible == TRUE]
members_all <- read_tsv_gz(membership_path)
groups <- c("bp", "mf", "cc")
expected_counts <- c(bp = 4087L, mf = 826L, cc = 511L)
stopifnot(all(groups %chin% inventory_all$analysis_group))
for (g in groups) {
  inv <- inventory_all[analysis_group == g]
  mem <- members_all[analysis_group == g]
  stopifnot(nrow(inv) == expected_counts[[g]], !anyDuplicated(inv$feature_set),
    setequal(inv$feature_set, unique(mem$feature_set)))
}

workers_requested <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", "8")))
if (!is.finite(workers_requested) || workers_requested < 1L)
  stop("N_WORKERS must be a positive integer")
workers <- min(workers_requested, parallel::detectCores())

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
cells <- as.character(wb_cells(bundle, "all"))
metadata <- as.data.table(wb_meta(bundle, cells = cells))
stopifnot(length(cells) == 26723L,
  identical(as.character(metadata$cell_id), cells),
  all(is.finite(metadata$nFeature_RNA)))
pt_all <- as.numeric(metadata$pt)
pt_defined <- is.finite(pt_all)
pt <- pt_all[pt_defined]
pt_cells <- cells[pt_defined]
stopifnot(length(pt_cells) == 12642L, sum(!pt_defined) == 14081L)
if (complexity_residualization) {
  complexity_design <- cbind(1, splines::ns(log1p(as.numeric(
    metadata$nFeature_RNA[pt_defined])), df = 4L))
  complexity_qr <- qr(complexity_design)
  stopifnot(complexity_qr$rank == ncol(complexity_design))
  log_event("Enabled cell-level RNA-feature residualization:",
    "UCell ~ ns(log1p(nFeature_RNA), df=4)")
}

width_fraction <- 0.05
step_fraction <- 0.01
min_cells <- 20L
q_cutoff <- 0.01
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

rank_path <- file.path(stage_dir, "trajectory_cell_rankings_maxrank1500.rds")
rank_cache_source <- Sys.getenv("GO_UCELL_RANKINGS_PATH", "")
if (nzchar(rank_cache_source)) {
  rank_cache_source <- normalizePath(rank_cache_source, mustWork = TRUE)
  log_event("Reusing validated trajectory-cell UCell rankings from",
    rank_cache_source)
  rankings <- readRDS(rank_cache_source)
} else if (file.exists(rank_path)) {
  log_event("Resuming saved trajectory-cell UCell rankings")
  rankings <- readRDS(rank_path)
} else {
  all_features <- as.character(wb_feature_ids(bundle, "RNA_counts"))
  rank_chunk_size <- 400L
  rank_chunk_dir <- file.path(stage_dir, "ranking_chunks")
  dir.create(rank_chunk_dir, showWarnings = FALSE)
  chunk_starts <- seq.int(1L, length(pt_cells), by = rank_chunk_size)
  rank_chunks <- vector("list", length(chunk_starts))
  log_event("Computing reusable sparse UCell rankings in",
    length(chunk_starts), "cell chunks")
  for (chunk_i in seq_along(chunk_starts)) {
    start <- chunk_starts[[chunk_i]]
    idx <- seq.int(start, min(length(pt_cells),
      start + rank_chunk_size - 1L))
    chunk_path <- file.path(rank_chunk_dir,
      sprintf("rank_chunk_%03d.rds", chunk_i))
    if (file.exists(chunk_path)) {
      rank_chunks[[chunk_i]] <- readRDS(chunk_path)
    } else {
      chunk_cells <- pt_cells[idx]
      counts <- methods::as(wb_slice(bundle, matrix_id = "RNA_counts",
        features = all_features, cells = chunk_cells), "dgCMatrix")
      stopifnot(identical(rownames(counts), all_features),
        identical(colnames(counts), chunk_cells))
      rank_chunks[[chunk_i]] <- UCell::StoreRankings_UCell(counts,
        maxRank = 1500L, chunk.size = rank_chunk_size,
        ncores = min(workers, length(chunk_cells)), ties.method = "average",
        force.gc = TRUE)
      saveRDS(rank_chunks[[chunk_i]], chunk_path, compress = FALSE)
      rm(counts); invisible(gc())
    }
    log_event("ranked", max(idx), "of", length(pt_cells), "cells")
  }
  rankings <- do.call(cbind, rank_chunks)
  rm(rank_chunks); invisible(gc())
  stopifnot(identical(colnames(rankings), pt_cells))
  saveRDS(rankings, rank_path, compress = FALSE)
}
stopifnot(identical(colnames(rankings), pt_cells))

# UCell 2.12's rankings2Uscore() infers maxRank as max(stored rank)+1.
# With average ties that can become 1500.5 rather than the requested 1500.
# Apply UCell's own U-statistic algebra directly to the sparse stored ranks so
# maxRank remains exactly 1500, matching scWorkbench's direct scoring route.
rank_presence <- rankings
rank_presence@x[] <- 1
ranking_gene_index <- setNames(seq_len(nrow(rankings)), rownames(rankings))
score_from_rankings <- function(signatures) {
  set_names <- names(signatures)
  sizes <- lengths(signatures)
  stopifnot(length(set_names) > 0L, all(sizes > 0L),
    all(unlist(signatures, use.names = FALSE) %in% rownames(rankings)))
  membership <- Matrix::sparseMatrix(
    i = rep(seq_along(signatures), sizes),
    j = unname(ranking_gene_index[unlist(signatures, use.names = FALSE)]),
    x = 1, dims = c(length(signatures), nrow(rankings)))
  stored_sum <- as.matrix(membership %*% rankings)
  stored_n <- as.matrix(membership %*% rank_presence)
  full_rank_sum <- stored_sum + (sizes - stored_n) * 1500
  rank_sum_min <- sizes * (sizes + 1) / 2
  denominator <- sizes * 1500 - rank_sum_min
  score <- 1 - sweep(sweep(full_rank_sum, 1L, rank_sum_min, "-"),
    1L, denominator, "/")
  score[score < 0 & score > -1e-12] <- 0
  stopifnot(all(is.finite(score)), all(score >= 0), all(score <= 1))
  rownames(score) <- set_names
  colnames(score) <- pt_cells
  score
}

# Validate the optimized precomputed-rank route against all retained native
# scWorkbench cytokine scores before using it for GO.
validation_path <- file.path(stage_dir, "native_ucell_validation.tsv")
if (!file.exists(validation_path)) {
  cyto_members <- members_all[analysis_group == "cytokines"]
  cyto_sets <- split(cyto_members$feature_id, cyto_members$feature_set)
  cyto_sets <- lapply(cyto_sets, unique)
  direct <- score_from_rankings(cyto_sets)
  stopifnot(identical(rownames(direct), names(cyto_sets)),
    identical(colnames(direct), pt_cells))
  native <- read_tsv_gz(file.path(cytokine_report,
    "ucell_cell_scores.tsv.gz"))[cell_id %chin% pt_cells]
  native_wide <- dcast(native, cell_id ~ key, value.var = "value")
  setkey(native_wide, cell_id)
  native_wide <- native_wide[pt_cells]
  cyto_names <- names(cyto_sets)
  native_matrix <- as.matrix(native_wide[, ..cyto_names])
  differences <- abs(t(direct) - native_matrix)
  validation <- data.table(feature_set = names(cyto_sets),
    n_cells = ncol(direct),
    maximum_absolute_difference = apply(differences, 2L, max),
    mean_absolute_difference = colMeans(differences))
  write_tsv(validation, validation_path)
  stopifnot(max(validation$maximum_absolute_difference) < 1e-12)
  rm(direct, native, native_wide, native_matrix, differences)
  invisible(gc())
  log_event("Validated optimized rankings against all 38 native cytokine sets")
} else {
  validation <- fread(validation_path)
  stopifnot(nrow(validation) == 38L,
    max(validation$maximum_absolute_difference) < 1e-12)
}

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

analyze_scores <- function(scores, fs, seed) {
  baseline <- median(scores)
  observed <- vapply(membership_index,
    function(idx) median(scores[idx]), numeric(1))
  effect <- observed - baseline
  observed_range <- max(observed) - min(observed)
  peak_index <- local_maxima_plateau(observed)
  peak_index <- peak_index[effect[peak_index] > .Machine$double.eps]
  curve <- data.table(feature_set = fs, windows,
    ucell_median = observed, trajectory_median = baseline,
    ucell_effect = effect)
  null_max <- numeric(n_permutations)
  null_range <- numeric(n_permutations)
  batch_size <- 100L
  set.seed(seed)
  for (start in seq.int(1L, n_permutations, by = batch_size)) {
    take <- min(batch_size, n_permutations - start + 1L)
    permuted <- replicate(take, sample(scores, length(scores), replace = FALSE))
    if (is.null(dim(permuted))) permuted <- matrix(permuted, ncol = 1L)
    batch_min <- rep(Inf, take)
    batch_max <- rep(-Inf, take)
    for (idx in membership_index) {
      window_median <- matrixStats::colMedians(permuted[idx, , drop = FALSE])
      batch_min <- pmin(batch_min, window_median)
      batch_max <- pmax(batch_max, window_median - baseline)
    }
    pos <- start:(start + take - 1L)
    null_max[pos] <- batch_max
    null_range[pos] <- (batch_max + baseline) - batch_min
  }
  null_sd <- sd(null_max)
  candidates <- if (length(peak_index)) data.table(
    feature_set = fs, peak_window = peak_index,
    peak_pseudotime = windows$window_center[peak_index],
    region_left = windows$window_left[peak_index],
    region_right = windows$window_right[peak_index],
    n_cells = windows$n_cells[peak_index],
    peak_ucell_median = observed[peak_index], trajectory_median = baseline,
    peak_ucell_effect = effect[peak_index],
    null_max_mean = mean(null_max), null_max_sd = null_sd,
    permutation_z = if (is.finite(null_sd) && null_sd > 0)
      (effect[peak_index] - mean(null_max)) / null_sd else NA_real_,
    exceedances = vapply(effect[peak_index],
      function(x) sum(null_max >= x), integer(1)),
    permutations_used = n_permutations) else NULL
  if (!is.null(candidates)) candidates[, p_region :=
    (1 + exceedances) / (1 + permutations_used)]
  list(curve = curve, candidates = candidates,
    permutation = data.table(feature_set = fs,
      candidate_peaks = length(peak_index),
      permutations_used = n_permutations,
      null_max_mean = mean(null_max), null_max_sd = null_sd,
      null_max_q95 = unname(quantile(null_max, .95)),
      null_max_q99 = unname(quantile(null_max, .99)),
      observed_dynamic_range = observed_range,
      null_range_mean = mean(null_range), null_range_sd = sd(null_range),
      null_range_q95 = unname(quantile(null_range, .95)),
      null_range_q99 = unname(quantile(null_range, .99)),
      dynamic_exceedances = sum(null_range >= observed_range),
      p_dynamic = (1 + sum(null_range >= observed_range)) /
        (1 + n_permutations), stop_reason = "complete"))
}

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
    representative <- order(g$q_region, g$p_region,
      -g$peak_ucell_effect)[[1L]]
    out <- copy(g[representative])
    out[, `:=`(region_left = min(g$region_left),
      region_right = max(g$region_right),
      n_constituent_peaks = nrow(g),
      constituent_peak_pts = paste(signif(g$peak_pseudotime, 7), collapse = ";"),
      constituent_peak_ids = paste(g$peak_id, collapse = ";"))]
    out
  }))
}

group_seed_offset <- c(bp = 0L, mf = 100000L, cc = 200000L)
score_batch_size <- 100L
group_completions <- list()
for (g in groups) {
  group_dir <- file.path(stage_dir, g)
  if (file.exists(file.path(group_dir, "completion.json"))) {
    log_event("Reusing completed group", g)
    group_completions[[g]] <- read_json(file.path(group_dir,
      "completion.json"), simplifyVector = TRUE)
    next
  }
  group_stage <- paste0(group_dir, ".incomplete")
  if (dir.exists(group_stage)) unlink(group_stage, recursive = TRUE)
  dir.create(group_stage, recursive = TRUE)
  inv <- copy(inventory_all[analysis_group == g])
  setorder(inv, feature_set)
  mem <- members_all[analysis_group == g]
  signatures <- split(mem$feature_id, mem$feature_set)
  signatures <- lapply(signatures[inv$feature_set], unique)
  stopifnot(!anyNA(names(signatures)),
    identical(names(signatures), inv$feature_set))
  curve_parts <- list(); candidate_parts <- list(); permutation_parts <- list()
  adjustment_qc_parts <- list()
  cursor <- 0L
  log_event(g, "scoring and testing", nrow(inv), "GO sets")
  for (start in seq.int(1L, nrow(inv), by = score_batch_size)) {
    idx <- seq.int(start, min(nrow(inv), start + score_batch_size - 1L))
    batch_names <- inv$feature_set[idx]
    score_matrix <- score_from_rankings(signatures[batch_names])
    stopifnot(identical(rownames(score_matrix), batch_names),
      identical(colnames(score_matrix), pt_cells))
    if (complexity_residualization) {
      raw_score_matrix <- score_matrix
      score_matrix <- t(qr.resid(complexity_qr, t(raw_score_matrix)))
      rownames(score_matrix) <- rownames(raw_score_matrix)
      colnames(score_matrix) <- colnames(raw_score_matrix)
      raw_centered <- raw_score_matrix - rowMeans(raw_score_matrix)
      raw_ss <- rowSums(raw_centered^2)
      residual_ss <- rowSums(score_matrix^2)
      zero_median_windows <- integer(nrow(raw_score_matrix))
      for (window_idx in membership_index) {
        zero_median_windows <- zero_median_windows +
          as.integer(abs(matrixStats::rowMedians(raw_score_matrix[, window_idx,
            drop = FALSE])) <= .Machine$double.eps)
      }
      adjustment_qc_parts[[cursor + 1L]] <- data.table(
        analysis_group = g, feature_set = batch_names,
        model = "lm(UCell ~ ns(log1p(nFeature_RNA), df=4))",
        model_r_squared = fifelse(raw_ss > .Machine$double.eps,
          1 - residual_ss / raw_ss, NA_real_),
        raw_zero_cell_fraction = rowMeans(raw_score_matrix <=
          .Machine$double.eps),
        raw_zero_median_window_fraction = zero_median_windows /
          length(membership_index),
        raw_cell_sd = matrixStats::rowSds(raw_score_matrix),
        residual_cell_sd = matrixStats::rowSds(score_matrix))
      rm(raw_score_matrix, raw_centered, raw_ss, residual_ss,
        zero_median_windows)
    }
    batch_results <- parallel::mclapply(seq_along(batch_names), function(jj) {
      global_index <- idx[[jj]]
      analyze_scores(as.numeric(score_matrix[jj, ]), batch_names[[jj]],
        base_seed + group_seed_offset[[g]] + global_index - 1L)
    }, mc.cores = min(workers, length(batch_names)), mc.preschedule = TRUE)
    if (any(vapply(batch_results, inherits, logical(1), "try-error")))
      stop("Numeric analysis failed for ", g, " batch starting ", start)
    cursor <- cursor + 1L
    curve_parts[[cursor]] <- rbindlist(lapply(batch_results, `[[`, "curve"))
    cand <- Filter(Negate(is.null), lapply(batch_results, `[[`, "candidates"))
    candidate_parts[[cursor]] <- if (length(cand)) rbindlist(cand) else NULL
    permutation_parts[[cursor]] <- rbindlist(lapply(batch_results, `[[`,
      "permutation"))
    rm(score_matrix, batch_results); invisible(gc())
    log_event(g, "completed", max(idx), "of", nrow(inv), "sets")
  }
  curves <- rbindlist(curve_parts)
  candidates <- rbindlist(Filter(Negate(is.null), candidate_parts))
  permutation_summary <- rbindlist(permutation_parts)
  adjustment_qc <- if (complexity_residualization)
    rbindlist(adjustment_qc_parts) else NULL
  curves[, analysis_group := g]
  setcolorder(curves, c("analysis_group", setdiff(names(curves),
    "analysis_group")))
  candidates[, analysis_group := g]
  if (inference_mode == "dynamic_range") {
    permutation_summary[, q_dynamic := p.adjust(p_dynamic, method = "BH")]
    permutation_summary[, dynamic_significant := q_dynamic < q_cutoff]
    candidates[permutation_summary[, .(feature_set, p_dynamic, q_dynamic,
      dynamic_significant)], on = "feature_set",
      `:=`(p_dynamic = i.p_dynamic, q_dynamic = i.q_dynamic,
        dynamic_significant = i.dynamic_significant)]
    candidates[, `:=`(q_region = q_dynamic,
      significant = dynamic_significant & p_region < region_fwer_cutoff)]
  } else {
    candidates[, q_region := p.adjust(p_region, method = "BH")]
    candidates[, significant := q_region < q_cutoff]
  }
  setorder(candidates, q_region, p_region, -peak_ucell_effect,
    feature_set, peak_pseudotime)
  candidates[, peak_id := sprintf(paste0("ucell_", g, "_peak_%06d"),
    seq_len(.N))]
  significant <- candidates[significant == TRUE]
  if (nrow(significant)) {
    final_regions <- rbindlist(lapply(split(seq_len(nrow(significant)),
      significant$feature_set), function(ii) merge_peaks(significant[ii])))
    setorder(final_regions, q_region, p_region, -peak_ucell_effect,
      feature_set, region_left)
    final_regions[, region_id := sprintf(paste0("ucell_", g,
      "_region_%06d"), seq_len(.N))]
  } else {
    final_regions <- data.table()
  }
  if (inference_mode == "dynamic_range") {
    region_counts <- if (nrow(final_regions)) final_regions[, .(
      n_regions = .N, best_region_fwer_p = min(p_region),
      strongest_ucell_effect = max(peak_ucell_effect),
      peak_pseudotimes = paste(signif(peak_pseudotime[order(region_left)], 6),
        collapse = ";")), by = feature_set] else data.table(
          feature_set = character(), n_regions = integer(),
          best_region_fwer_p = numeric(), strongest_ucell_effect = numeric(),
          peak_pseudotimes = character())
    significant_sets <- merge(permutation_summary[
      dynamic_significant == TRUE, .(feature_set, observed_dynamic_range,
        p_dynamic, q_dynamic)], region_counts,
      by = "feature_set", all.x = TRUE, sort = FALSE)
    significant_sets[, set_name := inv$set_name[
      match(feature_set, inv$feature_set)]]
    significant_sets[is.na(n_regions), n_regions := 0L]
    setorder(significant_sets, q_dynamic, p_dynamic,
      -observed_dynamic_range, feature_set)
  } else if (nrow(final_regions)) {
    significant_sets <- final_regions[, .(
      set_name = inv$set_name[match(feature_set, inv$feature_set)],
      n_regions = .N, best_q = min(q_region), best_p = min(p_region),
      strongest_permutation_z = suppressWarnings(max(permutation_z,
        na.rm = TRUE)),
      strongest_ucell_effect = max(peak_ucell_effect),
      peak_pseudotimes = paste(signif(peak_pseudotime[order(region_left)], 6),
        collapse = ";")), by = feature_set]
    setorder(significant_sets, best_q, best_p, -strongest_permutation_z,
      feature_set)
  } else significant_sets <- data.table()
  permutation_summary[, analysis_group := g]
  setcolorder(permutation_summary, c("analysis_group",
    setdiff(names(permutation_summary), "analysis_group")))
  write_tsv_gz(curves, file.path(group_stage, "window_curves.tsv.gz"))
  write_tsv_gz(candidates,
    file.path(group_stage, "all_candidate_peaks.tsv.gz"))
  write_tsv(significant,
    file.path(group_stage, "significant_individual_peaks.tsv"))
  write_tsv(final_regions,
    file.path(group_stage, "significant_positive_regions.tsv"))
  write_tsv(significant_sets,
    file.path(group_stage, "significant_gene_sets.tsv"))
  write_tsv(permutation_summary,
    file.path(group_stage, "permutation_summary.tsv"))
  if (complexity_residualization) write_tsv_gz(adjustment_qc,
    file.path(group_stage, "complexity_adjustment_qc.tsv.gz"))
  if (inference_mode == "dynamic_range") write_tsv(permutation_summary,
    file.path(group_stage, "gene_set_dynamic_tests.tsv"))
  counts <- data.table(analysis_group = g, tested_gene_sets = nrow(inv),
    candidate_peaks = nrow(candidates),
    significant_individual_peaks = nrow(significant),
    passing_gene_sets = nrow(significant_sets),
    nonpassing_gene_sets = nrow(inv) - nrow(significant_sets),
    final_merged_regions = nrow(final_regions))
  write_tsv(counts, file.path(group_stage, "filter_counts.tsv"))
  group_completion <- list(status = "complete", analysis_group = g,
    inference_mode = inference_mode,
    tested_gene_sets = nrow(inv), candidate_peaks = nrow(candidates),
    significant_individual_peaks = nrow(significant),
    passing_gene_sets = nrow(significant_sets),
    nonpassing_gene_sets = nrow(inv) - nrow(significant_sets),
    significant_positive_regions = nrow(final_regions),
    score = if (complexity_residualization)
      paste("residual of native-UCell-compatible RNA_counts max_rank=1500",
        "from ns(log1p(nFeature_RNA), df=4)") else
      "native-UCell-compatible RNA_counts max_rank=1500",
    complexity_residualization_enabled = complexity_residualization,
    complexity_residualization_model = if (complexity_residualization)
      "lm(UCell ~ ns(log1p(nFeature_RNA), df=4))" else NULL,
    raw_zero_median_warning_sets = if (complexity_residualization)
      adjustment_qc[raw_zero_median_window_fraction > 0.5, .N] else NULL,
    window_statistic = if (complexity_residualization)
      "median RNA-feature-residualized UCell" else "median UCell",
    test_statistic = if (inference_mode == "dynamic_range")
      "maximum window median minus minimum window median" else
      "window median minus full-trajectory median",
    permutation_null = if (inference_mode == "dynamic_range")
      "peak-to-trough range across trajectory" else
      "maximum positive window deviation across trajectory",
    permutations = n_permutations,
    multiple_testing = if (inference_mode == "dynamic_range")
      paste("BH across all", g, "gene sets") else
      paste("BH across all", g, "candidate peaks"),
    q_cutoff = q_cutoff, magnitude_or_prominence_filter = FALSE,
    positive_region_fwer_p_cutoff = if (inference_mode == "dynamic_range")
      region_fwer_cutoff else NULL,
    retained_windows = nrow(windows), minimum_window_cells = min_cells,
    merge_gap_fraction_strictly_less_than = merge_gap_fraction)
  write_json_file(group_completion,
    file.path(group_stage, "completion.json"))
  stopifnot(file.rename(group_stage, group_dir))
  group_completions[[g]] <- group_completion
  log_event(g, "complete:", nrow(significant_sets), "passing sets and",
    nrow(final_regions), "merged regions")
  rm(curves, candidates, significant, final_regions, significant_sets,
    permutation_summary, curve_parts, candidate_parts, permutation_parts,
    adjustment_qc, adjustment_qc_parts)
  invisible(gc())
}

summary_table <- rbindlist(lapply(groups, function(g)
  fread(file.path(stage_dir, g, "filter_counts.tsv"))))
write_tsv(summary_table, file.path(stage_dir, "group_summary.tsv"))
writeLines(c(if (inference_mode == "dynamic_range")
    "# GO UCell trajectory-dynamics permutation analysis" else
    "# GO UCell median-permutation analysis", "",
  "Native-UCell-compatible scores were calculated from one reusable RNA-count rank matrix and reduced immediately to trajectory-window summaries; the full cell-by-GO long table was not materialized.", "",
  "| Group | Tested | Passing | Non-passing | Significant peaks | Merged regions |",
  "|---|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(summary_table)), function(i) paste0("| ",
    toupper(summary_table$analysis_group[[i]]), " | ",
    summary_table$tested_gene_sets[[i]], " | ",
    summary_table$passing_gene_sets[[i]], " | ",
    summary_table$nonpassing_gene_sets[[i]], " | ",
    summary_table$significant_individual_peaks[[i]], " | ",
    summary_table$final_merged_regions[[i]], " |"), character(1)), "",
  "- Exact Bdev metadata `pt`: 12,642 defined and 14,081 missing entries preserved.",
  if (complexity_residualization)
    "- Score: native-compatible UCell on RNA counts, max_rank=1500, then cell-level residualization with `UCell ~ ns(log1p(nFeature_RNA), df=4)`; no pseudotime in the adjustment model."
  else "- UCell: RNA counts, max_rank=1500, no smoothing.",
  "- Native validation: all 38 retained cytokine cell scores matched to <1e-12.",
  "- Windows: 5% width, 1% step, minimum 20 cells; 97 retained.",
  if (inference_mode == "dynamic_range")
    "- Dynamic inference: peak-to-trough window-median range, including endpoints; 1,000 permutations and BH q<0.01 across gene sets separately within BP, MF, and CC."
  else "- Positive-only maximum-window permutation inference: 1,000 permutations and BH q<0.01 separately within BP, MF, and CC.",
  if (inference_mode == "dynamic_range")
    "- Positive-region localization: every positive local maximum is compared with the permutation maximum positive deviation; within-set trajectory-wide p<0.01; no candidate-score cutoff."
  else NULL,
  "- No raw-UCell magnitude or relative-prominence cutoff.",
  "- Multiple peaks retained; uncovered gaps strictly less than one window are merged.",
  if (complexity_residualization)
    "- Sparse-score audit tables report raw zero-cell fractions and the fraction of windows with zero raw-UCell median; these are warnings and do not alter inference."
  else NULL),
  file.path(stage_dir, "REPORT.md"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
write_json_file(list(status = "complete", groups = group_completions,
  inference_mode = inference_mode,
  score = if (complexity_residualization)
    paste("residual of native-UCell-compatible RNA_counts max_rank=1500",
      "from ns(log1p(nFeature_RNA), df=4)") else
    "native-UCell-compatible RNA_counts max_rank=1500",
  complexity_residualization_enabled = complexity_residualization,
  complexity_residualization_model = if (complexity_residualization)
    "lm(UCell ~ ns(log1p(nFeature_RNA), df=4))" else NULL,
  optimized_strategy = paste("precompute ranks once; score GO sets in batches;",
    "retain window summaries rather than full long cell-by-set table"),
  native_validation_sets = 38L,
  native_validation_maximum_absolute_difference =
    max(fread(validation_path)$maximum_absolute_difference),
  defined_pt_cells = length(pt_cells), missing_pt_cells = sum(!pt_defined),
  windows = nrow(windows), permutations = n_permutations,
  groupwise_multiple_testing = TRUE,
  multiple_testing_units = if (inference_mode == "dynamic_range")
    "gene sets" else "candidate peaks", q_cutoff = q_cutoff,
  positive_region_fwer_p_cutoff = if (inference_mode == "dynamic_range")
    region_fwer_cutoff else NULL,
  magnitude_or_prominence_filter = FALSE, random_seed = base_seed,
  workers = workers,
  reused_rankings_path = if (nzchar(rank_cache_source)) rank_cache_source else NULL,
  bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = fingerprint_after, bundle_unchanged = TRUE,
  promoted = FALSE), file.path(stage_dir, "completion.json"))

wb_close_bundle(bundle)
stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed GO UCell permutation report at", output_dir)
