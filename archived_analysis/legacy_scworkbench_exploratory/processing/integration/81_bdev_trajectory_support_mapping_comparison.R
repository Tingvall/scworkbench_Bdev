#!/usr/bin/env Rscript

# Calibrate trajectory membership exclusively from held-out Bdev cells, then
# compare the locked Bdev-derived gates with the current frozen-Bdev query
# mapping. Report-only; no source or staged bundle is modified or promoted.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) stop(paste(
  "Usage: 81_bdev_trajectory_support_mapping_comparison.R MAPPING_BUNDLE",
  "FROZEN_STAGE_BUNDLE BDEV_BUNDLE CURRENT_PROJECTION QUERY_NEIGHBORS",
  "SIGNAL_SCORES REGION_TABLE MAPPING_RECEIPT OUTPUT_DIR SCRIPT_VERSION"
), call. = FALSE)

mapping_bundle_path <- normalizePath(args[[1L]], mustWork = TRUE)
frozen_stage_path <- normalizePath(args[[2L]], mustWork = TRUE)
bdev_bundle_path <- normalizePath(args[[3L]], mustWork = TRUE)
projection_path <- normalizePath(args[[4L]], mustWork = TRUE)
query_neighbor_path <- normalizePath(args[[5L]], mustWork = TRUE)
signal_path <- normalizePath(args[[6L]], mustWork = TRUE)
region_path <- normalizePath(args[[7L]], mustWork = TRUE)
mapping_receipt_path <- normalizePath(args[[8L]], mustWork = TRUE)
output_dir <- normalizePath(args[[9L]], mustWork = FALSE)
script_version <- args[[10L]]

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output already exists: ", stage_dir, call. = FALSE)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage_dir, "figures"), showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(RANN)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
})

seed <- 1729L
k <- 30L
n_bins <- 100L
specificity_target <- 0.90
set.seed(seed)

write_json <- function(x, path) jsonlite::write_json(
  x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA,
  null = "null", na = "null")
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
weighted_on_fraction <- function(index, distance, on_label, k = 30L) {
  index <- index[, seq_len(k), drop = FALSE]
  distance <- distance[, seq_len(k), drop = FALSE]
  labs <- matrix(on_label[index], nrow = nrow(index))
  sigma <- distance[, k]
  sigma[sigma <= 0 | !is.finite(sigma)] <- 1
  weights <- exp(-(distance / sigma)^2)
  rowSums(weights * labs) / rowSums(weights)
}
auc_rank <- function(score, truth) {
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[truth]) - sum(seq_len(sum(truth)))) /
    (sum(truth) * sum(!truth))
}
operating_metrics <- function(score, truth, cutoff) list(
  cutoff = cutoff,
  sensitivity = mean(score[truth] >= cutoff),
  specificity = mean(score[!truth] < cutoff),
  balanced_accuracy = (mean(score[truth] >= cutoff) +
    mean(score[!truth] < cutoff)) / 2,
  accepted_fraction = mean(score >= cutoff))

mapping <- wb_open_bundle(mapping_bundle_path, read_only = TRUE,
  validate = FALSE)
frozen <- wb_open_bundle(frozen_stage_path, read_only = TRUE,
  validate = FALSE)
bdev <- wb_open_bundle(bdev_bundle_path, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(mapping), silent = TRUE), add = TRUE)
on.exit(try(wb_close_bundle(frozen), silent = TRUE), add = TRUE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)

mapping_manifest <- wb_manifest(mapping)
frozen_manifest <- wb_manifest(frozen)
bdev_manifest <- wb_manifest(bdev)
mapping_fingerprint_before <- mapping_manifest$current_state$fingerprint
bdev_fingerprint_before <- bdev_manifest$current_state$fingerprint
stopifnot(identical(mapping_manifest$bundle_id,
    "bdev_integrated_ko_mapping_v1"),
  identical(bdev_manifest$bundle_id, "bdev_reference"))

reference_ids <- as.character(wb_cells(mapping, "set:bdev_reference"))
query_ids <- as.character(wb_cells(mapping, "set:ko_query"))
reference_md <- as.data.table(wb_meta(mapping, cells = reference_ids))
query_md <- as.data.table(wb_meta(mapping, cells = query_ids))
reference_latent <- as.matrix(wb_space(frozen,
  "peakvi.full.reference", cells = reference_ids))
stopifnot(length(reference_ids) == 26723L, length(query_ids) == 23252L,
  identical(rownames(reference_latent), reference_ids),
  identical(dim(reference_latent), c(26723L, 20L)),
  all(is.finite(reference_latent)))
reference_md[, trajectory_truth := is.finite(reference_pt)]
stopifnot(sum(reference_md$trajectory_truth) == 12642L,
  sum(!reference_md$trajectory_truth) == 14081L)

mapping_receipt <- read_json(mapping_receipt_path, simplifyVector = TRUE)
stopifnot(identical(mapping_receipt$scope,
    "Full Bdev-trained PeakVI fit and frozen-model KO projection"),
  isTRUE(mapping_receipt$checks$projection_without_refit),
  identical(mapping_receipt$parameters$n_cells, 26723L))

# Leave one original Bdev sample out of the reference-neighbor pool. The PeakVI
# representation itself is the already frozen Bdev-only model; no query cell is
# used anywhere in cutoff fitting or selection.
folds <- sort(unique(as.character(reference_md$sample_id)))
cv <- rbindlist(lapply(folds, function(fold) {
  test <- which(as.character(reference_md$sample_id) == fold)
  train <- which(as.character(reference_md$sample_id) != fold)
  nn <- RANN::nn2(reference_latent[train, , drop = FALSE],
    reference_latent[test, , drop = FALSE], k = 50L,
    searchtype = "standard", eps = 0)
  data.table(
    cell_id = reference_ids[test], held_out_sample = fold,
    trajectory_truth = reference_md$trajectory_truth[test],
    trajectory_support = weighted_on_fraction(nn$nn.idx, nn$nn.dists,
      reference_md$trajectory_truth[train], k),
    nearest_distance = nn$nn.dists[, 1L],
    median30_distance = apply(nn$nn.dists[, seq_len(k), drop = FALSE],
      1L, median),
    kth30_distance = nn$nn.dists[, k])
}))
stopifnot(nrow(cv) == length(reference_ids), !anyDuplicated(cv$cell_id),
  all(is.finite(cv$trajectory_support)),
  all(cv$trajectory_support >= 0 & cv$trajectory_support <= 1))

ord <- order(cv$trajectory_support, decreasing = TRUE)
truth_ord <- cv$trajectory_truth[ord]
roc <- data.table(
  threshold = c(Inf, cv$trajectory_support[ord]),
  sensitivity = c(0, cumsum(truth_ord) / sum(truth_ord)),
  false_positive_rate = c(0, cumsum(!truth_ord) / sum(!truth_ord)))
roc[, specificity := 1 - false_positive_rate]
roc[, youden_j := sensitivity + specificity - 1]
balanced_i <- which.max(roc$youden_j)
balanced_cutoff <- roc$threshold[[balanced_i]]
conservative_cutoff <- as.numeric(quantile(
  cv[trajectory_truth == FALSE, trajectory_support], specificity_target,
  type = 1))
thresholds <- rbindlist(list(
  data.table(policy = "Current general-support gate", cutoff = NA_real_,
    target = "No Bdev trajectory-membership cutoff"),
  data.table(policy = "Bdev balanced cutoff", cutoff = balanced_cutoff,
    target = "Maximum pooled held-out Youden J"),
  data.table(policy = "Bdev 90%-specific cutoff",
    cutoff = conservative_cutoff,
    target = "90% pooled specificity for held-out off-trajectory Bdev cells")))
thresholds[, `:=`(
  sensitivity = c(NA_real_,
    mean(cv[trajectory_truth == TRUE, trajectory_support] >= balanced_cutoff),
    mean(cv[trajectory_truth == TRUE, trajectory_support] >=
      conservative_cutoff)),
  specificity = c(NA_real_,
    mean(cv[trajectory_truth == FALSE, trajectory_support] < balanced_cutoff),
    mean(cv[trajectory_truth == FALSE, trajectory_support] <
      conservative_cutoff)))]
thresholds[, balanced_accuracy := (sensitivity + specificity) / 2]

cv_by_sample <- rbindlist(lapply(thresholds[is.finite(cutoff)]$policy,
  function(policy_i) {
    cutoff_i <- thresholds[policy == policy_i, cutoff]
    cv[, .(
      n_cells = .N,
      n_on_trajectory = sum(trajectory_truth),
      n_off_trajectory = sum(!trajectory_truth),
      sensitivity = mean(trajectory_support[trajectory_truth] >= cutoff_i),
      specificity = mean(trajectory_support[!trajectory_truth] < cutoff_i)),
      by = held_out_sample][, policy := policy_i]
  }))
setcolorder(cv_by_sample, c("policy", "held_out_sample",
  setdiff(names(cv_by_sample), c("policy", "held_out_sample"))))

# Apply the locked Bdev-only scores to the already projected query. Scores use
# the saved Bdev-restricted neighbors from the frozen reference model.
neighbors <- fread(cmd = paste("zcat", shQuote(query_neighbor_path)))[rank <= k]
reference_key <- reference_md[, .(
  reference_cell_id = cell_id, trajectory_truth)]
neighbors <- reference_key[neighbors, on = "reference_cell_id"]
stopifnot(nrow(neighbors) == length(query_ids) * k,
  !anyNA(neighbors$trajectory_truth))
neighbors[, sigma := max(distance), by = query_cell_id]
neighbors[sigma <= 0 | !is.finite(sigma), sigma := 1]
neighbors[, weight := exp(-(distance / sigma)^2)]
query_support <- neighbors[, .(
  trajectory_support = sum(weight * trajectory_truth) / sum(weight),
  trajectory_neighbor_fraction_unweighted = mean(trajectory_truth),
  all_neighbor_median30_distance = median(distance)), by = query_cell_id]

projection <- fread(cmd = paste("zcat", shQuote(projection_path)))
stopifnot(nrow(projection) == length(query_ids),
  setequal(projection$cell_id, query_ids))
query <- merge(projection, query_support, by.x = "cell_id",
  by.y = "query_cell_id", all.x = TRUE, sort = FALSE)
query[, current_accepted := support_status %in% c("supported", "borderline") &
  is.finite(mapped_pt) & is.finite(projected_x) & is.finite(projected_y)]
query[, balanced_accepted := current_accepted &
  trajectory_support >= balanced_cutoff]
query[, conservative_accepted := current_accepted &
  trajectory_support >= conservative_cutoff]

policy_levels <- thresholds$policy
policy_query <- rbindlist(list(
  query[, .(cell_id, source_cell_id, sample_id, genotype, knockout_target,
    mapped_pt, projected_x, projected_y, trajectory_support,
    accepted = current_accepted,
    policy = "Current general-support gate")],
  query[, .(cell_id, source_cell_id, sample_id, genotype, knockout_target,
    mapped_pt, projected_x, projected_y, trajectory_support,
    accepted = balanced_accepted,
    policy = "Bdev balanced cutoff")],
  query[, .(cell_id, source_cell_id, sample_id, genotype, knockout_target,
    mapped_pt, projected_x, projected_y, trajectory_support,
    accepted = conservative_accepted,
    policy = "Bdev 90%-specific cutoff")]))
policy_query[, policy := factor(policy, levels = policy_levels)]

acceptance <- policy_query[, .(
  n_query_cells = .N,
  n_currently_mapped = sum(is.finite(mapped_pt) & is.finite(projected_x)),
  n_accepted = sum(accepted),
  accepted_fraction_all_query = mean(accepted),
  accepted_fraction_currently_mapped = sum(accepted) /
    sum(is.finite(mapped_pt) & is.finite(projected_x)),
  trajectory_support_q25 = quantile(trajectory_support, .25),
  trajectory_support_median = median(trajectory_support),
  trajectory_support_q75 = quantile(trajectory_support, .75)),
  by = .(policy, sample_id, genotype, knockout_target)]
setorder(acceptance, policy, knockout_target, genotype)

# Common Bdev pseudotime bins are fixed before evaluating query distributions.
bdev_source_ids <- as.character(wb_cells(bdev, "all"))
bdev_md <- as.data.table(wb_meta(bdev, cells = bdev_source_ids))
bdev_umap <- as.data.table(wb_space(bdev, "umap.atac", cells = bdev_source_ids))
stopifnot(identical(as.character(bdev_md$cell_id), bdev_source_ids),
  identical(as.character(bdev_umap$cell_id), bdev_source_ids),
  identical(as.numeric(bdev_md$pt),
    as.numeric(reference_md$reference_pt)))
pt_range <- range(bdev_md$pt, na.rm = TRUE)
breaks <- seq(pt_range[[1L]], pt_range[[2L]], length.out = n_bins + 1L)
bin_grid <- data.table(bin = seq_len(n_bins),
  bin_left = head(breaks, -1L), bin_right = tail(breaks, -1L))
bin_grid[, bin_midpoint := (bin_left + bin_right) / 2]
assign_bin <- function(pt) pmin(n_bins, pmax(1L,
  findInterval(pt, breaks, all.inside = TRUE)))
bdev_counts <- data.table(bin = assign_bin(bdev_md$pt[is.finite(bdev_md$pt)]))[
  , .(bdev_n_cells = .N), by = bin]
bdev_density <- merge(bin_grid, bdev_counts, by = "bin", all.x = TRUE)
bdev_density[is.na(bdev_n_cells), bdev_n_cells := 0L]
bdev_density[, bdev_fraction := bdev_n_cells / sum(bdev_n_cells)]
bdev_density[, bdev_rolling_fraction := frollmean(bdev_fraction, 3L,
  align = "center", fill = NA_real_)]
current_totals <- policy_query[policy == "Current general-support gate" &
  accepted & is.finite(mapped_pt), .(fixed_denominator = .N),
  by = .(knockout_target, genotype)]

density <- rbindlist(lapply(policy_levels, function(policy_i) {
  rbindlist(lapply(c("Flt3l", "Il7"), function(target_i) {
    rbindlist(lapply(c("WT", "KO"), function(genotype_i) {
      cells <- policy_query[policy == policy_i & knockout_target == target_i &
        genotype == genotype_i & accepted & is.finite(mapped_pt)]
      counts <- cells[, .(n_cells = .N), by = .(bin = assign_bin(mapped_pt))]
      ans <- merge(copy(bdev_density), counts, by = "bin", all.x = TRUE)
      ans[is.na(n_cells), n_cells := 0L]
      fixed_denominator_i <- current_totals[
        knockout_target == target_i & genotype == genotype_i,
        fixed_denominator]
      stopifnot(length(fixed_denominator_i) == 1L,
        fixed_denominator_i > 0L)
      ans[, `:=`(policy = policy_i, target = target_i,
        genotype = genotype_i, total_accepted = nrow(cells),
        fixed_denominator = fixed_denominator_i,
        conditional_cell_fraction = n_cells / nrow(cells),
        fixed_denominator_cell_fraction = n_cells / fixed_denominator_i)]
      ans[, `:=`(
        conditional_rolling_fraction = frollmean(conditional_cell_fraction,
          3L, align = "center", fill = NA_real_),
        fixed_denominator_rolling_fraction = frollmean(
          fixed_denominator_cell_fraction, 3L, align = "center",
          fill = NA_real_))]
      ans[, `:=`(
        conditional_bdev_normalized_rolling_fraction =
          conditional_rolling_fraction / bdev_rolling_fraction,
        fixed_denominator_bdev_normalized_rolling_fraction =
          fixed_denominator_rolling_fraction / bdev_rolling_fraction)]
      ans[!is.finite(conditional_bdev_normalized_rolling_fraction),
        conditional_bdev_normalized_rolling_fraction := NA_real_]
      ans[!is.finite(fixed_denominator_bdev_normalized_rolling_fraction),
        fixed_denominator_bdev_normalized_rolling_fraction := NA_real_]
      ans
    }))
  }))
}))
density[, policy := factor(policy, levels = policy_levels)]

# Detect WT enrichment regions separately under each locked policy using the
# pre-existing rule: maximal contiguous bins with query/Bdev ratio > 1.
wt_above <- density[genotype == "WT" &
  is.finite(fixed_denominator_bdev_normalized_rolling_fraction) &
  fixed_denominator_bdev_normalized_rolling_fraction > 1]
setorder(wt_above, policy, target, bin)
wt_above[, run := cumsum(seq_len(.N) == 1L |
  policy != shift(policy, fill = policy[[1L]]) |
  target != shift(target, fill = target[[1L]]) |
  bin != shift(bin, fill = bin[[1L]] - 1L) + 1L)]
enrichment_regions <- wt_above[, .(
  first_bin = min(bin), last_bin = max(bin), n_bins = .N,
  pt_left = min(bin_left), pt_right = max(bin_right),
  maximum_wt_ratio = max(
    fixed_denominator_bdev_normalized_rolling_fraction)),
  by = .(policy, target, run)]
setorder(enrichment_regions, policy, target, pt_left)
enrichment_regions[, peak_number := seq_len(.N), by = .(policy, target)]
enrichment_regions[, region_id := paste0(target, "_peak", peak_number)]

signal <- fread(cmd = paste("zcat", shQuote(signal_path)))[
  panel_label %in% c("Flt3L UP", "IL-7 UP")]
signal[, expected_target := fifelse(panel_label == "Flt3L UP", "Flt3l", "Il7")]
signal <- signal[knockout_target == expected_target,
  .(source_cell_id = cell_id, panel_label, aggregate_score)]
signal_policy <- merge(policy_query[accepted & is.finite(mapped_pt)], signal,
  by = "source_cell_id")
signal_bins <- signal_policy[, .(
  n_cells = .N, score_sum = sum(aggregate_score),
  bin_mean = mean(aggregate_score)),
  by = .(policy, knockout_target, genotype,
    bin = assign_bin(mapped_pt), panel_label)]
signal_trajectory <- merge(
  CJ(policy = factor(policy_levels, levels = policy_levels),
    knockout_target = c("Flt3l", "Il7"), genotype = c("WT", "KO"),
    bin = seq_len(n_bins), unique = TRUE),
  signal_bins, by = c("policy", "knockout_target", "genotype", "bin"),
  all.x = TRUE)
signal_trajectory <- merge(signal_trajectory,
  bin_grid, by = "bin", all.x = TRUE)
signal_trajectory[is.na(n_cells), `:=`(n_cells = 0L, score_sum = 0)]
setorder(signal_trajectory, policy, knockout_target, genotype, bin)
signal_trajectory[, `:=`(
  rolling_n_cells = frollsum(n_cells, 5L, align = "center", fill = NA_real_),
  rolling_score_sum = frollsum(score_sum, 5L, align = "center",
    fill = NA_real_)), by = .(policy, knockout_target, genotype)]
signal_trajectory[, rolling_mean := fifelse(rolling_n_cells >= 20,
  rolling_score_sum / rolling_n_cells, NA_real_)]
signal_trajectory[, panel_label := fifelse(knockout_target == "Flt3l",
  "Flt3L UP", "IL-7 UP")]

original_regions <- fread(region_path)[target %in% c("Flt3l", "Il7")]
region_retention <- rbindlist(lapply(seq_len(nrow(original_regions)), function(i) {
  reg <- original_regions[i]
  policy_query[knockout_target == reg$target & is.finite(mapped_pt) &
    mapped_pt >= reg$pt_left & mapped_pt < reg$pt_right,
    .(n_original_region_cells = .N, n_accepted = sum(accepted),
      accepted_fraction = mean(accepted)), by = .(policy, genotype)][,
      `:=`(region_id = reg$region_id, target = reg$target,
        pt_left = reg$pt_left, pt_right = reg$pt_right)]
}))
setcolorder(region_retention, c("region_id", "target", "policy", "genotype",
  "pt_left", "pt_right", "n_original_region_cells", "n_accepted",
  "accepted_fraction"))
setorder(region_retention, region_id, policy, genotype)

write_tsv(cv, file.path(stage_dir,
  "bdev_leave_one_sample_out_trajectory_support.tsv.gz"))
write_tsv(roc, file.path(stage_dir, "bdev_trajectory_support_roc.tsv"))
write_tsv(thresholds, file.path(stage_dir,
  "bdev_derived_trajectory_support_cutoffs.tsv"))
write_tsv(cv_by_sample, file.path(stage_dir,
  "bdev_cutoff_performance_by_held_out_sample.tsv"))
write_tsv(query[, .(cell_id, source_cell_id, sample_id, genotype,
  knockout_target, support_status, mapped_pt, trajectory_support,
  trajectory_neighbor_fraction_unweighted, current_accepted,
  balanced_accepted, conservative_accepted)],
  file.path(stage_dir, "query_trajectory_support_and_acceptance.tsv.gz"))
write_tsv(acceptance, file.path(stage_dir, "query_acceptance_summary.tsv"))
write_tsv(density, file.path(stage_dir,
  "query_bdev_normalized_density_by_policy.tsv"))
write_tsv(enrichment_regions[, !"run"], file.path(stage_dir,
  "wt_enrichment_regions_by_policy.tsv"))
write_tsv(signal_trajectory, file.path(stage_dir,
  "respective_signal_score_trajectory_by_policy.tsv"))
write_tsv(region_retention, file.path(stage_dir,
  "original_region_retention_by_policy.tsv"))

genotype_colours <- c("WT" = "#0072B2", "KO" = "#D55E00")
policy_labels <- c(
  "Current general-support gate" = "Current\ngeneral support",
  "Bdev balanced cutoff" = "Bdev-only\nbalanced cutoff",
  "Bdev 90%-specific cutoff" = "Bdev-only\n90% specificity")
target_labels <- c("Flt3l" = "Flt3L", "Il7" = "IL-7")
base_theme <- theme_bw(base_size = 10) + theme(
  panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
  legend.position = "bottom")

cv_plot <- copy(cv)
cv_plot[, truth_label := factor(ifelse(trajectory_truth,
  "On trajectory (finite pt)", "Off trajectory (missing pt)"),
  levels = c("Off trajectory (missing pt)", "On trajectory (finite pt)"))]
p_score <- ggplot(cv_plot,
    aes(trajectory_support, fill = truth_label, colour = truth_label)) +
  geom_density(alpha = .22, linewidth = .7) +
  geom_vline(xintercept = balanced_cutoff, colour = "#009E73",
    linetype = "dashed", linewidth = .7) +
  geom_vline(xintercept = conservative_cutoff, colour = "#CC79A7",
    linetype = "dotdash", linewidth = .7) +
  facet_wrap(~held_out_sample, nrow = 1) +
  scale_fill_manual(values = c("Off trajectory (missing pt)" = "grey55",
    "On trajectory (finite pt)" = "#56B4E9")) +
  scale_colour_manual(values = c("Off trajectory (missing pt)" = "grey35",
    "On trajectory (finite pt)" = "#0072B2")) +
  labs(x = "Weighted fraction of Bdev neighbors with defined pt",
    y = "Density", fill = "Known Bdev status", colour = "Known Bdev status",
    title = "A. Bdev-only leave-one-sample-out calibration",
    subtitle = "Vertical lines: balanced (green) and pooled 90%-specific (magenta) cutoffs") +
  base_theme
p_roc <- ggplot(roc, aes(false_positive_rate, sensitivity)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
    colour = "grey50") +
  geom_path(linewidth = .9, colour = "#3B77A8") +
  annotate("point", x = 1 - thresholds[policy == "Bdev balanced cutoff",
    specificity], y = thresholds[policy == "Bdev balanced cutoff", sensitivity],
    colour = "#009E73", size = 2.5) +
  annotate("point", x = 1 - thresholds[policy == "Bdev 90%-specific cutoff",
    specificity], y = thresholds[policy == "Bdev 90%-specific cutoff",
    sensitivity], colour = "#CC79A7", size = 2.5) +
  coord_equal() + labs(x = "False-positive rate", y = "Sensitivity",
    title = sprintf("B. Held-out ROC (AUC %.3f)",
      auc_rank(cv$trajectory_support, cv$trajectory_truth))) + base_theme +
  theme(legend.position = "none")
p_perf <- ggplot(cv_by_sample,
    aes(held_out_sample, sensitivity, fill = policy)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  geom_point(aes(y = specificity, colour = policy),
    position = position_dodge(width = .8), shape = 21, fill = "white",
    size = 2.2, stroke = .8) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = c("Bdev balanced cutoff" = "#009E73",
    "Bdev 90%-specific cutoff" = "#CC79A7"), labels = policy_labels) +
  scale_colour_manual(values = c("Bdev balanced cutoff" = "#00794F",
    "Bdev 90%-specific cutoff" = "#8E4A72"), labels = policy_labels) +
  labs(x = "Held-out Bdev sample", y = "Rate",
    fill = "Cutoff", colour = "Specificity points",
    title = "C. Cross-sample sensitivity and specificity",
    subtitle = "Bars: sensitivity; open points: specificity") + base_theme
calibration_figure <- p_score / (p_roc | p_perf) +
  plot_annotation(title = "Trajectory-membership threshold learned only from Bdev") +
  plot_layout(heights = c(1.05, 1))
ggsave(file.path(stage_dir, "figures",
  "bdev_trajectory_support_calibration.pdf"), calibration_figure,
  width = 12, height = 8.5, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures",
  "bdev_trajectory_support_calibration.png"), calibration_figure,
  width = 12, height = 8.5, dpi = 220, bg = "white")

reference_plot <- data.table(x = bdev_umap$x, y = bdev_umap$y)
umap_query <- copy(policy_query[is.finite(projected_x) & is.finite(projected_y)])
umap_query[, display_status := fifelse(accepted, as.character(genotype),
  "Rejected")]
umap_query[, policy_display := factor(policy_labels[as.character(policy)],
  levels = unname(policy_labels[policy_levels]))]
umap_query[, sample_display := factor(sample_id,
  levels = c("FLT3LKO_BM_WT", "FLT3LKO_BM_KO", "IL7KO_BM_WT",
    "IL7KO_BM_KO"))]
umap_query[, draw_order := display_status != "Rejected"]
setorder(umap_query, draw_order)
p_umap <- ggplot() +
  geom_point(data = reference_plot, aes(x, y), colour = "grey91", size = .07) +
  geom_point(data = umap_query,
    aes(projected_x, projected_y, colour = display_status),
    size = .18, alpha = .58) +
  facet_grid(sample_display ~ policy_display) + coord_equal() +
  scale_colour_manual(values = c("Rejected" = "grey72", genotype_colours),
    breaks = c("Rejected", "WT", "KO")) +
  labs(x = "Bdev umap.atac 1", y = "Bdev umap.atac 2",
    colour = "Query status",
    title = "Frozen-Bdev query projection before and after Bdev-only trajectory gates",
    subtitle = "Grey background: unchanged Bdev UMAP; grey query points fail the indicated gate") +
  base_theme
ggsave(file.path(stage_dir, "figures",
  "query_umap_mapping_policy_comparison.pdf"), p_umap,
  width = 12, height = 13, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures",
  "query_umap_mapping_policy_comparison.png"), p_umap,
  width = 12, height = 13, dpi = 220, bg = "white")

density[, `:=`(
  target_display = factor(target_labels[target], levels = unname(target_labels)),
  policy_display = factor(policy_labels[as.character(policy)],
    levels = unname(policy_labels[policy_levels])))]
p_density <- ggplot(density,
    aes(bin_midpoint, fixed_denominator_bdev_normalized_rolling_fraction,
      colour = genotype, group = genotype)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey30") +
  geom_line(linewidth = .8, na.rm = TRUE) +
  facet_grid(target_display ~ policy_display, scales = "free_y") +
  scale_colour_manual(values = genotype_colours) +
  coord_cartesian(xlim = pt_range) +
  labs(x = "Mapped Bdev pseudotime",
    y = "Retained query / Bdev rolling fraction", colour = "Genotype",
    title = "A. Enrichment peaks under current and Bdev-only trajectory gates",
    subtitle = "All policies use the current-mapped sample size as a fixed denominator; ratio > 1 indicates retained-cell enrichment") +
  base_theme
signal_trajectory[, `:=`(
  target_display = factor(target_labels[knockout_target],
    levels = unname(target_labels)),
  policy_display = factor(policy_labels[as.character(policy)],
    levels = unname(policy_labels[policy_levels])))]
p_signal <- ggplot(signal_trajectory,
    aes(bin_midpoint, rolling_mean, colour = genotype, group = genotype)) +
  geom_line(linewidth = .8, na.rm = TRUE) +
  facet_grid(target_display ~ policy_display, scales = "free_y") +
  scale_colour_manual(values = genotype_colours) +
  coord_cartesian(xlim = pt_range) +
  labs(x = "Mapped Bdev pseudotime",
    y = "5-bin cell-weighted mean signal score", colour = "Genotype",
    title = "B. Respective cytokine-response signal along the retained trajectory",
    subtitle = "Flt3L samples: Flt3L UP score; IL-7 samples: IL-7 UP score; windows require at least 20 cells") +
  base_theme
trajectory_figure <- p_density / p_signal +
  plot_annotation(
    title = "Effect of Bdev-derived trajectory support on mapped enrichment and signal") +
  plot_layout(heights = c(1, 1))
ggsave(file.path(stage_dir, "figures",
  "enrichment_and_signal_policy_comparison.pdf"), trajectory_figure,
  width = 13, height = 10, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures",
  "enrichment_and_signal_policy_comparison.png"), trajectory_figure,
  width = 13, height = 10, dpi = 220, bg = "white")

early_retention <- region_retention[region_id == "Il7_peak1"]
peak2_retention <- region_retention[region_id == "Il7_peak2"]
report <- c(
  "# Frozen-Bdev trajectory-support mapping comparison", "",
  "## What was compared", "",
  paste(
    "The existing mapping already uses a PeakVI model fitted on Bdev alone",
    "(26,723 reference cells) and projects WT/KO with `transform_without_refit`.",
    "Therefore this report compares the current general Bdev-support gate with",
    "two trajectory-membership cutoffs learned exclusively from Bdev."), "",
  "Bdev truth was defined before examining WT/KO: finite preserved `pt` is on trajectory; missing `pt` is off trajectory. Each original Bdev sample was mapped to the other three samples in the frozen 20-dimensional reference space.", "",
  "## Bdev calibration", "",
  sprintf("- Weighted on-trajectory neighbor fraction: held-out AUC %.3f.",
    auc_rank(cv$trajectory_support, cv$trajectory_truth)),
  sprintf("- Balanced cutoff %.4f: sensitivity %.1f%%, specificity %.1f%%.",
    balanced_cutoff,
    100 * thresholds[policy == "Bdev balanced cutoff", sensitivity],
    100 * thresholds[policy == "Bdev balanced cutoff", specificity]),
  sprintf("- Pooled 90%%-specific cutoff %.4f: sensitivity %.1f%%, specificity %.1f%%.",
    conservative_cutoff,
    100 * thresholds[policy == "Bdev 90%-specific cutoff", sensitivity],
    100 * thresholds[policy == "Bdev 90%-specific cutoff", specificity]),
  "- Performance varies across held-out Bdev samples, so the conservative gate is included as sensitivity analysis rather than assumed universally optimal.", "",
  "## Query comparison", "",
  "The UMAP figure shows which query populations survive each gate. The trajectory figure compares Bdev-normalized mapped-cell enrichment and the corresponding Flt3L-UP or IL-7-UP aggregate signal score.", "",
  "For a fair comparison, density curves retain the current-mapped sample size as the denominator under every policy. Filtered subsets are not renormalized to 100%; conditional-on-retained densities are exported separately in the density table.", "",
  sprintf("For the original IL-7 early interval, the balanced gate retains WT %.1f%% and KO %.1f%%; the 90%%-specific gate retains WT %.1f%% and KO %.1f%%.",
    100 * early_retention[policy == "Bdev balanced cutoff" & genotype == "WT",
      accepted_fraction],
    100 * early_retention[policy == "Bdev balanced cutoff" & genotype == "KO",
      accepted_fraction],
    100 * early_retention[policy == "Bdev 90%-specific cutoff" & genotype == "WT",
      accepted_fraction],
    100 * early_retention[policy == "Bdev 90%-specific cutoff" & genotype == "KO",
      accepted_fraction]),
  sprintf("For the original IL-7 peak-2 interval, the corresponding retention is WT %.1f%% and KO %.1f%% under the balanced gate, and WT %.1f%% and KO %.1f%% under the 90%%-specific gate.",
    100 * peak2_retention[policy == "Bdev balanced cutoff" & genotype == "WT",
      accepted_fraction],
    100 * peak2_retention[policy == "Bdev balanced cutoff" & genotype == "KO",
      accepted_fraction],
    100 * peak2_retention[policy == "Bdev 90%-specific cutoff" & genotype == "WT",
      accepted_fraction],
    100 * peak2_retention[policy == "Bdev 90%-specific cutoff" & genotype == "KO",
      accepted_fraction]), "",
  sprintf("Across complete samples, the balanced gate retains %.1f%%–%.1f%% of currently mapped query cells and the 90%%-specific gate retains %.1f%%–%.1f%%.",
    100 * min(acceptance[policy == "Bdev balanced cutoff",
      accepted_fraction_currently_mapped]),
    100 * max(acceptance[policy == "Bdev balanced cutoff",
      accepted_fraction_currently_mapped]),
    100 * min(acceptance[policy == "Bdev 90%-specific cutoff",
      accepted_fraction_currently_mapped]),
    100 * max(acceptance[policy == "Bdev 90%-specific cutoff",
      accepted_fraction_currently_mapped])), "",
  "## Interpretation limits", "",
  paste(
    "The cutoff is independent of WT/KO composition, but it is specific to the",
    "frozen Bdev reference, its feature space, and its preserved definition of",
    "trajectory membership. The latent model was fitted to all Bdev cells before",
    "the sample-held-out neighbor calibration; no WT/KO cell influenced model",
    "fitting, score definition, or cutoff selection."), "",
  "All outputs are report-only and unpromoted. Source bundles and Bdev `pt` were unchanged.", "",
  "## Figures", "",
  "- `figures/bdev_trajectory_support_calibration.pdf`", 
  "- `figures/query_umap_mapping_policy_comparison.pdf`",
  "- `figures/enrichment_and_signal_policy_comparison.pdf`", "",
  sprintf("Script version: `%s`.", script_version))
writeLines(report, file.path(stage_dir, "REPORT.md"))

write_tsv(data.table(role = c("mapping_bundle", "frozen_stage_bundle",
  "bdev_bundle", "current_projection", "query_neighbors", "signal_scores",
  "original_regions", "mapping_receipt"), path = c(mapping_bundle_path,
  frozen_stage_path, bdev_bundle_path, projection_path, query_neighbor_path,
  signal_path, region_path, mapping_receipt_path)),
  file.path(stage_dir, "input_manifest.tsv"))
write_json(list(
  status = "complete_unpromoted_report_only",
  script_version = script_version,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  reference_cells = length(reference_ids),
  reference_on_trajectory = sum(reference_md$trajectory_truth),
  reference_off_trajectory = sum(!reference_md$trajectory_truth),
  query_cells = length(query_ids), k = k,
  held_out_unit = "Bdev original sample",
  auc = auc_rank(cv$trajectory_support, cv$trajectory_truth),
  thresholds = thresholds,
  mapping_bundle_fingerprint_before = mapping_fingerprint_before,
  mapping_bundle_fingerprint_after =
    wb_manifest(mapping)$current_state$fingerprint,
  bdev_bundle_fingerprint_before = bdev_fingerprint_before,
  bdev_bundle_fingerprint_after = wb_manifest(bdev)$current_state$fingerprint,
  reference_pt_preserved = identical(as.numeric(bdev_md$pt),
    as.numeric(reference_md$reference_pt)), promoted = FALSE),
  file.path(stage_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

stopifnot(identical(mapping_fingerprint_before,
    wb_manifest(mapping)$current_state$fingerprint),
  identical(bdev_fingerprint_before,
    wb_manifest(bdev)$current_state$fingerprint))
if (!file.rename(stage_dir, output_dir))
  stop("Could not atomically finalize output directory", call. = FALSE)
cat("Comparison report:", file.path(output_dir, "REPORT.md"), "\n")
