#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 61_refine_cc_z20_permutations.R",
  "<bundle_dir> <inventory.tsv.gz> <memberships.tsv.gz>",
  "<source_cc_dir> <output_dir> <seed> <permutations>"))

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
inventory_path <- normalizePath(args[[2L]], mustWork = TRUE)
membership_path <- normalizePath(args[[3L]], mustWork = TRUE)
source_dir <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
base_seed <- as.integer(args[[6L]])
n_permutations <- as.integer(args[[7L]])
stopifnot(is.finite(base_seed), is.finite(n_permutations),
  n_permutations >= 10000L, n_permutations %% 10000L == 0L)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(jsonlite)
  library(splines)
})
read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ..., "\n")
  flush.console()
}

group <- "cc"
z_cutoff <- 20
checkpoint <- 10000L
source_regions <- fread(file.path(source_dir,
  "significant_positive_regions.tsv"))
source_candidates <- read_tsv_gz(file.path(source_dir,
  "all_candidate_peaks.tsv.gz"))
selected_regions <- source_regions[peak_z > z_cutoff]
selected_sets <- unique(selected_regions$feature_set)
stopifnot(nrow(selected_regions) == 88L, length(selected_sets) == 88L,
  !anyDuplicated(selected_regions$peak_id),
  all(selected_regions$peak_id %chin% source_candidates$peak_id))

inventory_all <- read_tsv_gz(inventory_path)[analysis_group == group &
  eligible == TRUE]
inventory_all[, original_inventory_index := .I]
inventory <- inventory_all[feature_set %chin% selected_sets]
members <- read_tsv_gz(membership_path)[analysis_group == group &
  feature_set %chin% selected_sets]
setorder(inventory, feature_set); setorder(members, feature_set, feature_id)
stopifnot(nrow(inventory) == length(selected_sets),
  setequal(inventory$feature_set, selected_sets),
  setequal(inventory$feature_set, unique(members$feature_set)))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
cells <- as.character(wb_cells(bundle, "all"))
metadata <- as.data.table(wb_meta(bundle, cells = cells))
stopifnot(length(cells) == 26723L,
  identical(as.character(metadata$cell_id), cells),
  all(is.finite(metadata$nCount_RNA)), all(metadata$nCount_RNA > 0))
pt_all <- as.numeric(metadata$pt); pt_defined <- is.finite(pt_all)
pt <- pt_all[pt_defined]; pt_cells <- cells[pt_defined]
stopifnot(length(pt) == 12642L, sum(!pt_defined) == 14081L)

log_depth <- log1p(as.numeric(metadata$nCount_RNA)[pt_defined])
depth_z <- (log_depth - mean(log_depth)) / sd(log_depth)
depth_design <- cbind(1, ns(pt, df = 8L), depth_z)
depth_qr <- qr(depth_design)

pt_limits <- range(pt); pt_span <- diff(pt_limits)
centers <- pt_limits[[1L]] + seq(0, 1, by = 0.01) * pt_span
window_width <- 0.05 * pt_span
windows <- data.table(grid_window_index = seq_along(centers),
  window_center = centers,
  window_left = pmax(pt_limits[[1L]], centers - window_width / 2),
  window_right = pmin(pt_limits[[2L]], centers + window_width / 2))
membership_index <- lapply(seq_len(nrow(windows)), function(i)
  which(pt >= windows$window_left[[i]] & pt <= windows$window_right[[i]]))
windows[, n_cells := lengths(membership_index)]
windows <- windows[n_cells >= 30L]
membership_index <- lapply(seq_len(nrow(windows)), function(i)
  which(pt >= windows$window_left[[i]] & pt <= windows$window_right[[i]]))
weight_matrix <- matrix(0, nrow(windows), length(pt))
for (i in seq_len(nrow(windows)))
  weight_matrix[i, membership_index[[i]]] <- 1 / length(membership_index[[i]])

raw_score_matrix <- matrix(NA_real_, nrow(inventory), length(pt),
  dimnames = list(inventory$feature_set, pt_cells))
log_event("Reading scores for", nrow(inventory), "selected CC gene sets")
for (ii in seq_len(nrow(inventory))) {
  fs <- inventory$feature_set[[ii]]
  genes <- members[feature_set == fs, feature_id]
  expression <- as.matrix(wb_matrix_read_features(bundle,
    "RNA_log1p_cp10k", genes, pt_cells, missing = "NA"))
  raw_score_matrix[ii, ] <- colMeans(expression, na.rm = TRUE)
}

analyze_one <- function(ii) {
  fs <- inventory$feature_set[[ii]]
  target <- selected_regions[feature_set == fs]
  raw_scores <- raw_score_matrix[ii, ]
  beta <- unname(qr.coef(depth_qr, raw_scores)[ncol(depth_design)])
  scores <- raw_scores - beta * depth_z
  null_mean <- mean(scores)
  null_sd <- sd(scores) * sqrt((1 - windows$n_cells / length(scores)) /
    windows$n_cells)
  observed_curve <- as.numeric(weight_matrix %*% scores)
  observed_z <- (observed_curve - null_mean) / null_sd
  stopifnot(all(target$peak_window >= 1L),
    all(target$peak_window <= length(observed_z)))
  reproduced_z <- observed_z[target$peak_window]
  if (!isTRUE(all.equal(reproduced_z, target$peak_z,
    tolerance = 1e-10, check.attributes = FALSE)))
    stop("Observed z-score mismatch for ", fs)

  exceed_checkpoint <- integer(nrow(target))
  exceed_final <- integer(nrow(target))
  set.seed(base_seed + inventory$original_inventory_index[[ii]] - 1L)
  for (start in seq.int(1L, n_permutations, by = 100L)) {
    take <- min(100L, n_permutations - start + 1L)
    permuted <- replicate(take,
      sample(scores, length(scores), replace = FALSE))
    if (is.null(dim(permuted))) permuted <- matrix(permuted, ncol = 1L)
    null_curves <- weight_matrix %*% permuted
    null_z <- sweep(sweep(null_curves, 1L, null_mean, "-"), 1L,
      null_sd, "/")
    null_max <- apply(null_z, 2L, max)
    increments <- vapply(reproduced_z,
      function(v) sum(null_max >= v), integer(1))
    exceed_final <- exceed_final + increments
    if (start + take - 1L <= checkpoint)
      exceed_checkpoint <- exceed_checkpoint + increments
  }
  data.table(analysis_group = group, feature_set = fs,
    set_name = target$set_name, peak_id = target$peak_id,
    region_id = target$region_id, peak_window = target$peak_window,
    peak_pseudotime = target$peak_pseudotime, peak_z = reproduced_z,
    original_exceedances_10k = target$exceedances,
    replay_exceedances_10k = exceed_checkpoint,
    replay_matches_original = exceed_checkpoint == target$exceedances,
    refined_exceedances = exceed_final,
    refined_permutations = n_permutations,
    p_region_refined = (1 + exceed_final) / (1 + n_permutations))
}

workers_requested <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", "8")))
if (!is.finite(workers_requested) || workers_requested < 1L)
  stop("N_WORKERS must be a positive integer")
workers <- min(workers_requested, parallel::detectCores(), nrow(inventory))
log_event("Starting", n_permutations, "permutations with", workers, "workers")
started <- Sys.time()
results <- parallel::mclapply(seq_len(nrow(inventory)), analyze_one,
  mc.cores = workers, mc.preschedule = TRUE)
if (any(vapply(results, inherits, logical(1), "try-error")))
  stop("At least one parallel refinement failed")
refined <- rbindlist(results)
elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
stopifnot(nrow(refined) == nrow(selected_regions),
  all(refined$replay_matches_original))

# Preserve the original CC family for BH adjustment: only replace the p-values
# of the selected peak IDs, then adjust all 3,012 original candidate peaks.
all_p <- source_candidates[, .(peak_id, p_for_adjustment = p_region)]
all_p[refined, on = "peak_id", p_for_adjustment := i.p_region_refined]
all_p[, q_region_refined_full_cc_family := p.adjust(p_for_adjustment,
  method = "BH")]
refined[all_p, on = "peak_id",
  q_region_refined_full_cc_family := i.q_region_refined_full_cc_family]
setorder(refined, -peak_z, feature_set)
fwrite(refined, file.path(stage_dir, "cc_z20_regions_100k.tsv"), sep = "\t",
  quote = FALSE, na = "NA")

summary <- data.table(
  selected_gene_sets = uniqueN(refined$feature_set),
  selected_regions = nrow(refined), permutations = n_permutations,
  replay_matches = sum(refined$replay_matches_original),
  zero_exceedance_regions = sum(refined$refined_exceedances == 0L),
  nonzero_exceedance_regions = sum(refined$refined_exceedances > 0L),
  minimum_p = min(refined$p_region_refined),
  maximum_p = max(refined$p_region_refined),
  minimum_q_full_cc_family = min(refined$q_region_refined_full_cc_family),
  maximum_q_full_cc_family = max(refined$q_region_refined_full_cc_family),
  permutation_elapsed_seconds = elapsed_seconds)
fwrite(summary, file.path(stage_dir, "summary.tsv"), sep = "\t",
  quote = FALSE, na = "NA")

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
completion <- list(status = "complete", analysis_group = group,
  scope = "final regions with peak z > 20 in the source 30-cell analysis",
  selected_gene_sets = uniqueN(refined$feature_set),
  selected_regions = nrow(refined), permutations = n_permutations,
  null_statistic = "maximum positive window z-score across the trajectory",
  checkpoint_replay = paste0(checkpoint,
    " original-seed permutations reproduced all source exceedance counts"),
  multiple_testing = paste0("selected p-values substituted into all ",
    nrow(source_candidates), " original CC candidate peaks before BH"),
  random_seed = base_seed, workers = workers,
  permutation_elapsed_seconds = elapsed_seconds,
  pseudotime_source = "exact Bdev metadata pt", minimum_window_cells = 30L,
  bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = fingerprint_after,
  bundle_unchanged = identical(fingerprint_before, fingerprint_after))
write_json(completion, file.path(stage_dir, "completion.json"))

report <- c(
  "# CC z > 20 targeted 100k-permutation pilot", "",
  paste0("- Source: `", source_dir, "`"),
  paste0("- Selected gene sets/regions: ", uniqueN(refined$feature_set),
    "/", nrow(refined)),
  paste0("- Permutations per selected gene set: ", n_permutations),
  "- Null statistic: maximum positive window z-score across the full trajectory",
  paste0("- Original 10k replay matches: ",
    sum(refined$replay_matches_original), "/", nrow(refined)),
  paste0("- Regions with zero exceedances at 100k: ",
    sum(refined$refined_exceedances == 0L), "/", nrow(refined)),
  paste0("- Empirical p-value range: ",
    format(min(refined$p_region_refined), scientific = TRUE), " to ",
    format(max(refined$p_region_refined), scientific = TRUE)),
  paste0("- BH q-value range using the full original CC family: ",
    format(min(refined$q_region_refined_full_cc_family), scientific = TRUE),
    " to ",
    format(max(refined$q_region_refined_full_cc_family), scientific = TRUE)),
  paste0("- Permutation runtime: ", round(elapsed_seconds / 60, 2), " minutes"),
  paste0("- Bundle unchanged: ", identical(fingerprint_before,
    fingerprint_after)))
writeLines(report, file.path(stage_dir, "REPORT.md"))
stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed targeted CC refinement at", output_dir)

