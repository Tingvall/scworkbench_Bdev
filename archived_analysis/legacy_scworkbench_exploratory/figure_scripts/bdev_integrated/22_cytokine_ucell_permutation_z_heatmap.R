#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 22_cytokine_ucell_permutation_z_heatmap.R",
  "<native_ucell_report> <permutation_report> <row_membership.tsv>",
  "<bb.tsv> <output_dir>"))
source_dir <- normalizePath(args[[1L]], mustWork = TRUE)
permutation_dir <- normalizePath(args[[2L]], mustWork = TRUE)
membership_path <- normalizePath(args[[3L]], mustWork = TRUE)
bb_path <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(matrixStats)
})

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

completion <- read_json(file.path(permutation_dir, "completion.json"),
  simplifyVector = TRUE)
stopifnot(identical(completion$status, "complete"),
  completion$tested_gene_sets == 38L,
  completion$permutations == 1000L,
  completion$q_cutoff == 0.01,
  completion$retained_windows == 97L,
  isTRUE(completion$native_window_medians_reproduced))
n_permutations <- as.integer(completion$permutations)
base_seed <- as.integer(completion$random_seed)

curves <- read_tsv_gz(file.path(permutation_dir, "window_curves.tsv.gz"))
significant_sets <- fread(file.path(permutation_dir,
  "significant_gene_sets.tsv"))$feature_set
membership <- fread(membership_path)
setorder(membership, row_position_top_to_bottom)
feature_sets <- sort(unique(curves$feature_set))
stopifnot(length(feature_sets) == 38L,
  length(significant_sets) == 15L,
  setequal(feature_sets, membership$feature_set),
  all(significant_sets %chin% feature_sets),
  uniqueN(curves$window_index) == 97L,
  nrow(curves) == 38L * 97L)

pt_table <- fread(file.path(source_dir, "source_pseudotime.tsv"))
setnames(pt_table, c("cell_id", "pt")); setorder(pt_table, cell_id)
stopifnot(nrow(pt_table) == completion$defined_pt_cells,
  all(is.finite(pt_table$pt)))
pt <- pt_table$pt
scores_long <- read_tsv_gz(file.path(source_dir,
  "ucell_cell_scores.tsv.gz"))[, .(key, cell_id, value)]
scores_long <- scores_long[cell_id %chin% pt_table$cell_id]
stopifnot(nrow(scores_long) == 38L * nrow(pt_table),
  !anyDuplicated(scores_long, by = c("key", "cell_id")))
score_wide <- dcast(scores_long, cell_id ~ key, value.var = "value")
setorder(score_wide, cell_id)
stopifnot(identical(score_wide$cell_id, pt_table$cell_id),
  setequal(setdiff(names(score_wide), "cell_id"), feature_sets))
score_matrix <- t(as.matrix(score_wide[, ..feature_sets]))
rownames(score_matrix) <- feature_sets
rm(scores_long, score_wide); invisible(gc())

window_map <- unique(curves[, .(grid_window_index, window_index,
  window_center, window_left, window_right, n_cells)])
setorder(window_map, window_index)
membership_index <- lapply(seq_len(nrow(window_map)), function(i)
  which(pt >= window_map$window_left[[i]] & pt <=
    window_map$window_right[[i]]))
stopifnot(identical(lengths(membership_index), window_map$n_cells))

calculate_one <- function(ii) {
  fs <- feature_sets[[ii]]
  scores <- as.numeric(score_matrix[ii, ])
  baseline <- median(scores)
  observed <- curves[feature_set == fs][order(window_index)]$ucell_effect
  null_sum <- numeric(length(membership_index))
  null_sum_sq <- numeric(length(membership_index))
  null_less <- integer(length(membership_index))
  null_equal <- integer(length(membership_index))
  batch_size <- 100L
  set.seed(base_seed + ii - 1L)
  for (start in seq.int(1L, n_permutations, by = batch_size)) {
    take <- min(batch_size, n_permutations - start + 1L)
    permuted <- replicate(take, sample(scores, length(scores), replace = FALSE))
    if (is.null(dim(permuted))) permuted <- matrix(permuted, ncol = 1L)
    null_curves <- vapply(membership_index, function(idx)
      matrixStats::colMedians(permuted[idx, , drop = FALSE]) - baseline,
      numeric(take))
    null_curves <- t(null_curves)
    null_sum <- null_sum + rowSums(null_curves)
    null_sum_sq <- null_sum_sq + rowSums(null_curves^2)
    delta <- sweep(null_curves, 1L, observed, "-")
    null_less <- null_less + rowSums(delta < -1e-12)
    null_equal <- null_equal + rowSums(abs(delta) <= 1e-12)
  }
  null_mean <- null_sum / n_permutations
  null_var <- (null_sum_sq - n_permutations * null_mean^2) /
    (n_permutations - 1L)
  null_var[null_var < 0 & null_var > -1e-15] <- 0
  null_sd <- sqrt(null_var)
  # Convert the mid-rank empirical null CDF to a standard-normal score. This
  # remains finite and maps an all-tied zero null to z=0, unlike division by a
  # zero/near-zero permutation SD for sparse UCell medians.
  empirical_cdf <- (null_less + 0.5 * null_equal + 0.5) /
    (n_permutations + 1)
  z <- qnorm(empirical_cdf)
  data.table(feature_set = fs, window_map,
    ucell_effect = observed, permutation_null_mean = null_mean,
    permutation_null_sd = null_sd,
    permutations_below_observed = null_less,
    permutations_tied_observed = null_equal,
    empirical_null_cdf_midrank = empirical_cdf,
    ucell_permutation_z = z)
}

workers_requested <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", "8")))
if (!is.finite(workers_requested) || workers_requested < 1L)
  stop("N_WORKERS must be positive")
workers <- min(workers_requested, parallel::detectCores(), length(feature_sets))
log_event("Calculating window-specific UCell null moments with", workers,
  "workers")
parts <- parallel::mclapply(seq_along(feature_sets), calculate_one,
  mc.cores = workers, mc.preschedule = TRUE)
if (any(vapply(parts, inherits, logical(1), "try-error")))
  stop("At least one window-z worker failed")
z_long <- rbindlist(parts)
stopifnot(nrow(z_long) == 38L * 97L)

prior_order <- membership$feature_set
passing_order <- prior_order[prior_order %chin% significant_sets]
nonpassing_order <- prior_order[!prior_order %chin% significant_sets]
display_order <- c(passing_order, nonpassing_order)
z_long[, status := factor(ifelse(feature_set %chin% significant_sets,
  "Passing: ≥1 region at BH q < 0.01", "Not passing"),
  levels = c("Passing: ≥1 region at BH q < 0.01", "Not passing"))]
# Factor levels run bottom-to-top; reversing each section preserves the earlier
# UCell heatmap's top-to-bottom row order within both status groups.
z_long[, feature_label := factor(short_label(feature_set),
  levels = rev(short_label(display_order)))]

finite_z <- z_long$ucell_permutation_z[is.finite(z_long$ucell_permutation_z)]
color_limit <- max(3, as.numeric(quantile(abs(finite_z), .995)))
bb <- fread(bb_path)
stopifnot(nrow(bb) == 2L,
  identical(as.character(bb$label), c("up1", "up2")))
bb[, boundary_label := c("BB1", "BB2")]

p <- ggplot(z_long, aes(window_center, feature_label,
    fill = ucell_permutation_z)) +
  geom_raster() +
  geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
    linewidth = .45, linetype = "dashed", color = "#374151") +
  geom_text(data = bb, aes(x = pt_median, y = Inf, label = boundary_label),
    inherit.aes = FALSE, vjust = 1.25, size = 3, fontface = "bold",
    color = "#374151") +
  facet_grid(status ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
    midpoint = 0, limits = c(-color_limit, color_limit),
    oob = scales::squish, na.value = "#D1D5DB",
    name = "UCell empirical\npermutation z") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(drop = TRUE) +
  labs(title = "Cytokine UCell empirical permutation z-scores",
    subtitle = paste0("Passing and non-passing sets use the same scale and ",
      "the previous UCell row order; 1,000 permutations per set"),
    x = "Existing Bdev pseudotime (pt)", y = "Cytokine gene set") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),
    strip.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"))

ggsave(file.path(stage_dir,
  "cytokine_ucell_permutation_z_passing_vs_not_passing.png"), p,
  width = 12, height = 11, dpi = 240, bg = "white")
ggsave(file.path(stage_dir,
  "cytokine_ucell_permutation_z_passing_vs_not_passing.pdf"), p,
  width = 12, height = 11, device = cairo_pdf, bg = "white")
fwrite(z_long, file.path(stage_dir,
  "cytokine_ucell_window_permutation_z.tsv.gz"), sep = "\t", quote = FALSE,
  na = "NA", compress = "gzip")
fwrite(data.table(feature_set = display_order,
  short_name = short_label(display_order),
  status = ifelse(display_order %chin% significant_sets, "passing",
    "not_passing"),
  row_position_within_status = c(seq_along(passing_order),
    seq_along(nonpassing_order)),
  prior_ucell_row_position = membership$row_position_top_to_bottom[
    match(display_order, membership$feature_set)]),
  file.path(stage_dir, "heatmap_row_order.tsv"), sep = "\t", quote = FALSE)

writeLines(c("# Cytokine UCell permutation-z heatmap", "",
  paste0("- Passing sets: ", length(passing_order)),
  paste0("- Non-passing sets: ", length(nonpassing_order)),
  paste0("- Shared symmetric display limit: ±", signif(color_limit, 5)),
  "- Each window z-score is the standard-normal transform of its mid-rank empirical permutation-null percentile.",
  "- This empirical definition maps an all-tied zero null to z=0 and avoids unstable division by zero or near-zero null SD.",
  "- Passing status comes from the separate maximum-window permutation test with BH q < 0.01; the heatmap z-score itself was not used to redefine passing status.",
  "- No UCell magnitude or prominence cutoff was applied.",
  "- Rows preserve the prior UCell heatmap order within passing status."),
  file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", gene_sets = 38L,
  passing_sets = length(passing_order),
  nonpassing_sets = length(nonpassing_order), windows = 97L,
  permutations = n_permutations,
  heatmap_value = paste("window-specific mid-rank empirical permutation",
    "normal-score z of median UCell deviation"),
  passing_definition = "at least one final region at BH q < 0.01",
  shared_color_scale = TRUE, symmetric_color_limit = color_limit,
  row_order = "prior UCell heatmap order within passing status",
  source_permutation_report = permutation_dir, random_seed = base_seed,
  workers = workers), file.path(stage_dir, "completion.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)

stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed UCell permutation-z heatmap at", output_dir)
