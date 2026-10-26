#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 93_ebf1_depth_complexity_matched_da.R BUNDLE",
  "SOURCE_DA_REPORT OUTPUT_DIR SESSION_ROOT EXPECTED_FINGERPRINT",
  "CALIPER_SD SCRIPT_VERSION"
), call. = FALSE)

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
source_dir <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
session_root <- normalizePath(args[[4L]], mustWork = TRUE)
expected_fp <- args[[5L]]
caliper_sd <- as.numeric(args[[6L]])
script_version <- args[[7L]]
stopifnot(is.finite(caliper_sd), caliper_sd > 0)

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    paste(..., collapse = " "), "\n")
  flush.console()
}
smd <- function(x, group) {
  a <- x[group == "WT"]
  b <- x[group == "KO"]
  denom <- sqrt((stats::var(a) + stats::var(b)) / 2)
  if (is.finite(denom) && denom > 0) (mean(a) - mean(b)) / denom else NA_real_
}

# Maximum-cardinality, minimum-distance ordered matching for a one-dimensional
# complexity score. Both inputs must be sorted; matches beyond the caliper are
# forbidden. This avoids a greedy-order artifact without external packages.
ordered_caliper_match <- function(a, b, caliper) {
  n <- length(a)
  m <- length(b)
  if (!n || !m) return(matrix(integer(), ncol = 2L,
    dimnames = list(NULL, c("a", "b"))))
  nm <- matrix(0L, n + 1L, m + 1L)
  cost <- matrix(0, n + 1L, m + 1L)
  action <- matrix(0L, n + 1L, m + 1L)
  for (i in seq_len(n)) for (j in seq_len(m)) {
    candidates_n <- c(nm[i, j + 1L], nm[i + 1L, j])
    candidates_cost <- c(cost[i, j + 1L], cost[i + 1L, j])
    candidates_action <- c(1L, 2L) # skip a, skip b
    d <- abs(a[[i]] - b[[j]])
    if (d <= caliper) {
      candidates_n <- c(candidates_n, nm[i, j] + 1L)
      candidates_cost <- c(candidates_cost, cost[i, j] + d)
      candidates_action <- c(candidates_action, 3L)
    }
    best_n <- max(candidates_n)
    eligible <- which(candidates_n == best_n)
    k <- eligible[which.min(candidates_cost[eligible])]
    nm[i + 1L, j + 1L] <- candidates_n[[k]]
    cost[i + 1L, j + 1L] <- candidates_cost[[k]]
    action[i + 1L, j + 1L] <- candidates_action[[k]]
  }
  i <- n
  j <- m
  pairs <- list()
  while (i > 0L && j > 0L) {
    act <- action[i + 1L, j + 1L]
    if (act == 3L) {
      pairs[[length(pairs) + 1L]] <- c(i, j)
      i <- i - 1L
      j <- j - 1L
    } else if (act == 1L) {
      i <- i - 1L
    } else {
      j <- j - 1L
    }
  }
  if (!length(pairs)) return(matrix(integer(), ncol = 2L,
    dimnames = list(NULL, c("a", "b"))))
  out <- do.call(rbind, rev(pairs))
  colnames(out) <- c("a", "b")
  out
}

assignment_path <- file.path(source_dir,
  "shared_region_cell_assignments.tsv.gz")
original_da_path <- file.path(source_dir,
  "all_targets_ebf1_pm2mb_shared_min10_wt_vs_ko_da.tsv")
stopifnot(file.exists(assignment_path), file.exists(original_da_path))
assignments <- fread(cmd = paste("zcat", shQuote(assignment_path)))
original_da <- fread(original_da_path)
stopifnot(!anyDuplicated(assignments$cell_id), nrow(original_da) == 502L,
  uniqueN(original_da$feature_id) == 251L)

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$current_state$fingerprint, expected_fp),
  all(assignments$cell_id %in% wb_cells(bundle, "set:ko_query")))

# Compute total counts and detected peaks over the entire 97,788-peak matrix in
# bounded blocks. These are matching covariates, not Ebf1-window outcomes.
features <- as.character(wb_feature_ids(bundle, "ATAC_bdev_mm10_peaks"))
cells <- as.character(assignments$cell_id)
blocks <- split(features, ceiling(seq_along(features) / 512L))
peak_counts <- numeric(length(cells))
detected_peaks <- numeric(length(cells))
log_event("Computing full-peak depth and complexity in", length(blocks),
  "bounded blocks")
for (bi in seq_along(blocks)) {
  m <- as.matrix(wb_matrix_read_features(bundle, "ATAC_bdev_mm10_peaks",
    blocks[[bi]], cells, missing = "NA"))
  stopifnot(!anyNA(m), identical(colnames(m), cells))
  peak_counts <- peak_counts + colSums(m)
  detected_peaks <- detected_peaks + colSums(m > 0)
  if (bi %% 25L == 0L) log_event("Depth block", bi, "of", length(blocks))
}
assignments[, `:=`(peak_counts = peak_counts,
  detected_peaks = detected_peaks,
  log1p_peak_counts = log1p(peak_counts),
  log1p_detected_peaks = log1p(detected_peaks))]

pair_tables <- list()
balance_tables <- list()
retention_tables <- list()
for (target_value in c("Flt3l", "Il7")) {
  z <- copy(assignments[knockout_target == target_value])
  x <- scale(as.matrix(z[, .(log1p_peak_counts, log1p_detected_peaks)]))
  pc <- stats::prcomp(x, center = FALSE, scale. = FALSE)
  score <- pc$x[, 1L]
  if (cor(score, z$log1p_peak_counts) < 0) score <- -score
  z[, complexity_score := score]
  target_sd <- stats::sd(z$complexity_score)
  caliper <- caliper_sd * target_sd

  target_pairs <- list()
  target_retention <- list()
  for (bin_value in sort(unique(z$bin))) {
    b <- z[bin == bin_value]
    wt <- b[genotype == "WT"][order(complexity_score, cell_id)]
    ko <- b[genotype == "KO"][order(complexity_score, cell_id)]
    common_left <- max(min(wt$complexity_score), min(ko$complexity_score))
    common_right <- min(max(wt$complexity_score), max(ko$complexity_score))
    wt_overlap <- wt[complexity_score >= common_left &
      complexity_score <= common_right]
    ko_overlap <- ko[complexity_score >= common_left &
      complexity_score <= common_right]
    match_index <- ordered_caliper_match(wt_overlap$complexity_score,
      ko_overlap$complexity_score, caliper)
    matched <- if (nrow(match_index)) data.table(
      knockout_target = target_value,
      bin = bin_value,
      pair_in_bin = seq_len(nrow(match_index)),
      wt_cell_id = wt_overlap$cell_id[match_index[, "a"]],
      ko_cell_id = ko_overlap$cell_id[match_index[, "b"]],
      wt_peak_counts = wt_overlap$peak_counts[match_index[, "a"]],
      ko_peak_counts = ko_overlap$peak_counts[match_index[, "b"]],
      wt_detected_peaks = wt_overlap$detected_peaks[match_index[, "a"]],
      ko_detected_peaks = ko_overlap$detected_peaks[match_index[, "b"]],
      wt_complexity_score = wt_overlap$complexity_score[match_index[, "a"]],
      ko_complexity_score = ko_overlap$complexity_score[match_index[, "b"]],
      score_distance = abs(wt_overlap$complexity_score[match_index[, "a"]] -
        ko_overlap$complexity_score[match_index[, "b"]]),
      caliper_score_units = caliper) else data.table()
    target_pairs[[as.character(bin_value)]] <- matched
    target_retention[[as.character(bin_value)]] <- data.table(
      knockout_target = target_value, bin = bin_value,
      original_wt = nrow(wt), original_ko = nrow(ko),
      common_support_wt = nrow(wt_overlap),
      common_support_ko = nrow(ko_overlap),
      matched_pairs = nrow(matched), caliper_score_units = caliper)
  }
  pairs <- rbindlist(target_pairs, fill = TRUE)
  pairs[, pair_id := sprintf("%s_bin%02d_pair%04d", tolower(target_value),
    bin, pair_in_bin)]
  pair_tables[[target_value]] <- pairs
  retention_tables[[target_value]] <- rbindlist(target_retention)

  matched_ids <- c(pairs$wt_cell_id, pairs$ko_cell_id)
  z[, matched := cell_id %in% matched_ids]
  balance_tables[[target_value]] <- rbindlist(lapply(c(FALSE, TRUE),
    function(is_matched) {
      q <- if (is_matched) z[matched == TRUE] else z
      data.table(knockout_target = target_value,
        population = if (is_matched) "matched" else "original",
        cells_per_genotype = min(q[, .N, by = genotype]$N),
        wt_cells = q[genotype == "WT", .N],
        ko_cells = q[genotype == "KO", .N],
        peak_counts_smd = smd(q$log1p_peak_counts, q$genotype),
        detected_peaks_smd = smd(q$log1p_detected_peaks, q$genotype),
        wt_median_peak_counts = median(q[genotype == "WT", peak_counts]),
        ko_median_peak_counts = median(q[genotype == "KO", peak_counts]),
        wt_median_detected_peaks = median(q[genotype == "WT", detected_peaks]),
        ko_median_detected_peaks = median(q[genotype == "KO", detected_peaks]))
    }))
}

pairs <- rbindlist(pair_tables)
retention <- rbindlist(retention_tables)
balance <- rbindlist(balance_tables)
stopifnot(!anyDuplicated(pairs$wt_cell_id), !anyDuplicated(pairs$ko_cell_id),
  all(pairs$score_distance <= pairs$caliper_score_units + 1e-12))
write_tsv(pairs, file.path(stage_dir, "matched_cell_pairs.tsv"))
write_tsv(retention, file.path(stage_dir, "matching_retention_by_bin.tsv"))
write_tsv(balance, file.path(stage_dir, "depth_complexity_balance.tsv"))
fwrite(assignments[, matched := cell_id %in%
  c(pairs$wt_cell_id, pairs$ko_cell_id)],
  file.path(stage_dir, "cell_depth_complexity_and_match_status.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip")

# Reuse the native full-matrix normalization staged by analysis 91 and rerun
# native Wilcoxon DA on the matched cells only.
session_id <- "ebf1_pm2mb_shared_support_da_v1"
workspace <- file.path(session_root, session_id)
stopifnot(dir.exists(workspace))
session <- wb_session_resume(bundle, workspace)
stopifnot(is.null(session$durable$promotion_receipt))
find_run <- function(session, marker, operation) {
  ids <- names(session$runs)[vapply(session$runs, function(x) {
    note <- if (is.null(x$notes)) "" else as.character(x$notes)
    identical(x$op_name, operation) && startsWith(note, paste0("[", marker, "]"))
  }, logical(1L))]
  stopifnot(length(ids) <= 1L)
  if (length(ids)) session$runs[[ids[[1L]]]] else NULL
}
normal_run <- find_run(session, "normalize_full_atac", "matrix_log_normalize")
stopifnot(!is.null(normal_run), "matrix" %in% names(normal_run$outputs))
normalized_id <- "ATAC_bdev_mm10_peaks_log1p_cp10k_ebf1_pm2mb_shared_v1"
matched_de <- list()
for (target_value in c("Flt3l", "Il7")) {
  p <- pairs[knockout_target == target_value]
  wt_cells <- as.character(p$wt_cell_id)
  ko_cells <- as.character(p$ko_cell_id)
  marker <- paste0("de_", tolower(target_value),
    "_depth_complexity_matched_caliper", gsub("\\.", "p", caliper_sd))
  de_run <- find_run(session, marker, "de_wilcox")
  if (is.null(de_run)) {
    log_event("Native matched DA for", target_value, length(wt_cells),
      "WT and", length(ko_cells), "KO cells")
    result <- wb_session_run_op(session, bundle, "de_wilcox",
      inputs = list(matrix = normalized_id, cell_set = "set:ko_query",
        min_cells = 10L, group_a = wt_cells, group_b = ko_cells,
        features = unique(original_da[target == target_value, feature_id]),
        min_pct = 0, p_adjust_method = "BH"),
      notes = paste0("[", marker, "] Group A=WT; group B=KO; exact-bin ",
        "matching on full-peak count/detection complexity PC1; common support; ",
        "caliper=", caliper_sd, " target SD; Ebf1 TSS +/-2 Mb; native ",
        "full-matrix log1p-CP10K; descriptive cell-level sensitivity; no promotion."))
    session <- result$session
    de_run <- result$run
  }
  d <- as.data.table(de_run$outputs$de_table$table)
  source <- original_da[target == target_value]
  d <- merge(d, source[, .(feature_id, chromosome, peak_start_1based,
    peak_end_1based, peak_midpoint_1based, midpoint_relative_to_tss_bp,
    original_pct_wt = pct_a, original_pct_ko = pct_b,
    original_avg_log1p_fc = avg_log1p_fc,
    original_p_val_adj = p_val_adj,
    original_direction = direction,
    original_significant_bh05 = significant_bh05)],
    by = "feature_id", all.x = TRUE)
  d[, `:=`(target = target_value,
    matched_pairs = nrow(p),
    matched_detection_difference_wt_minus_ko = pct_a - pct_b,
    matched_direction = fifelse(avg_log1p_fc > 0, "WT_higher",
      fifelse(avg_log1p_fc < 0, "KO_higher", "equal")),
    matched_significant_bh05 = is.finite(p_val_adj) & p_val_adj < 0.05)]
  d <- d[order(p_val_adj, -abs(avg_log1p_fc), feature_id)]
  matched_de[[target_value]] <- d
  write_tsv(d, file.path(stage_dir, paste0(tolower(target_value),
    "_depth_complexity_matched_ebf1_pm2mb_da.tsv")))
}
combined <- rbindlist(matched_de, use.names = TRUE)
write_tsv(combined, file.path(stage_dir,
  "all_targets_depth_complexity_matched_ebf1_pm2mb_da.tsv"))

summary <- combined[, .(
  matched_pairs = unique(matched_pairs),
  tested_peaks = .N,
  matched_significant_bh05 = sum(matched_significant_bh05),
  matched_significant_wt_higher = sum(matched_significant_bh05 &
    matched_direction == "WT_higher"),
  matched_significant_ko_higher = sum(matched_significant_bh05 &
    matched_direction == "KO_higher"),
  original_significant_bh05 = sum(original_significant_bh05),
  original_significant_retained_bh05 = sum(original_significant_bh05 &
    matched_significant_bh05),
  original_significant_same_direction = sum(original_significant_bh05 &
    original_direction == matched_direction),
  effect_spearman = cor(original_avg_log1p_fc, avg_log1p_fc,
    method = "spearman")), by = target]
write_tsv(summary, file.path(stage_dir, "matched_da_summary.tsv"))

plot_data <- combined[, .(target, feature_id,
  original = original_avg_log1p_fc, matched = avg_log1p_fc,
  original_significant_bh05, matched_significant_bh05)]
p <- ggplot(plot_data, aes(original, matched,
  color = matched_significant_bh05)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
    color = "grey45") +
  geom_point(alpha = 0.75, size = 1.5) +
  facet_wrap(~target, scales = "free") +
  scale_color_manual(values = c(`TRUE` = "#7A3E9D", `FALSE` = "grey70"),
    name = "Matched BH < 0.05") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom") +
  labs(x = "Original pooled avg_log1p_fc (WT - KO)",
    y = "Depth/complexity-matched avg_log1p_fc (WT - KO)",
    title = "Ebf1-window DA sensitivity to ATAC depth/complexity matching")
ggsave(file.path(stage_dir, "figures",
  "original_vs_depth_complexity_matched_effects.png"), p,
  width = 9, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "original_vs_depth_complexity_matched_effects.pdf"), p,
  width = 9, height = 5.5, device = cairo_pdf, bg = "white")

write_json(list(status = "complete_unpromoted_report_only",
  script_version = script_version,
  source_report = source_dir,
  bundle_fingerprint = expected_fp,
  matching = list(exact_strata = "knockout_target x original pseudotime bin",
    covariates = c("log1p total counts across 97,788 Bdev peaks",
      "log1p number of detected Bdev peaks"),
    score = "oriented PC1 of standardized covariates within target",
    common_support = TRUE, caliper_target_sd = caliper_sd,
    algorithm = "maximum-cardinality minimum-distance ordered matching without replacement"),
  balance = split(balance, balance$knockout_target),
  results = split(summary, summary$target),
  biological_replicates_per_genotype = 1L,
  inference = "descriptive sensitivity analysis",
  bundle_mutated = FALSE, promotion = FALSE),
  file.path(stage_dir, "diagnostics.json"))

balance_lines <- apply(balance, 1L, function(x) paste0("- `",
  x[["knockout_target"]], " ", x[["population"]], "`: ",
  x[["wt_cells"]], " WT / ", x[["ko_cells"]], " KO; count SMD ",
  sprintf("%.3f", as.numeric(x[["peak_counts_smd"]])),
  "; detected-peak SMD ",
  sprintf("%.3f", as.numeric(x[["detected_peaks_smd"]])), "."))
summary_lines <- apply(summary, 1L, function(x) paste0("- `",
  x[["target"]], "`: ", x[["matched_pairs"]], " pairs; ",
  x[["matched_significant_bh05"]], " matched BH-significant peaks (",
  x[["matched_significant_wt_higher"]], " WT-higher; ",
  x[["matched_significant_ko_higher"]], " KO-higher); ",
  x[["original_significant_retained_bh05"]], "/",
  x[["original_significant_bh05"]], " original calls remain significant."))
writeLines(c(
  "# Ebf1 ±2 Mb depth/complexity-matched DA sensitivity", "",
  "Cells were matched without replacement within every original pseudotime bin using a common-support, caliper-constrained ATAC complexity score based on total counts and detected peaks across all 97,788 Bdev peaks.", "",
  "## Balance", "", balance_lines, "",
  "Absolute standardized mean differences below 0.1 are conventionally considered small.", "",
  "## Matched differential accessibility", "", summary_lines, "",
  "The matched cell sets were retested with native scWorkbench full-matrix log1p-CP10K normalization and `de_wilcox`; BH adjustment is within each 251-peak comparison.", "",
  "Because genotype is represented by one biological sample, matching cannot separate genotype biology from sample batch and the result remains descriptive.", "",
  "- [Summary](matched_da_summary.tsv)",
  "- [Balance](depth_complexity_balance.tsv)",
  "- [Matched pairs](matched_cell_pairs.tsv)",
  "- [Complete matched DA](all_targets_depth_complexity_matched_ebf1_pm2mb_da.tsv)",
  "- [Effect comparison](figures/original_vs_depth_complexity_matched_effects.png)"
), file.path(stage_dir, "REPORT.md"))

capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed Ebf1 depth/complexity-matched DA at", output_dir, "\n")
