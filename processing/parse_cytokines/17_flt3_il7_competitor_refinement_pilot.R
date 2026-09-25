#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 2L)) {
  stop(paste("Usage: 17_flt3_il7_competitor_refinement_pilot.R",
    "<output_directory> [candidate_top_n]"))
}

output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(UCell)
  library(Matrix)
  library(jsonlite)
  library(ggplot2)
  library(patchwork)
})
source("src/bundle_config.R")
source("src/bundle_validation.R")

candidate_top_n <- if (length(args) == 2L) as.integer(args[[2L]]) else 100L
cache_top_n <- 100L
stopifnot(is.finite(candidate_top_n), candidate_top_n >= 1L,
  candidate_top_n <= cache_top_n)
analysis_id <- paste0("flt3_il7_competitor_refinement_top", candidate_top_n,
  "_v1")
namespace <- paste0("parse_cytokine_rest_only_refinement_pilot_top",
  candidate_top_n, "_v1")
targets <- c("Flt3L", "IL-7")
directions <- c("up", "down")
folds <- c("1", "2", "3")
minimum_size <- 20L
log2fc_minimum <- 0.5
fdr_maximum <- 0.05
max_rank <- 1500L
rank_outer_chunk <- 400L
rank_inner_chunk <- 100L
ncores <- as.integer(Sys.getenv("PARSE_THREADS", "4"))

read_tsv <- function(path) read.delim(path, check.names = FALSE,
  na.strings = c("NA", ""), stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
slug <- function(x) gsub("(^_+|_+$)", "", tolower(gsub("[^A-Za-z0-9]+", "_", x)))
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}

registry <- load_bundle_registry()
parent_cfg <- read_project_yaml("config/parse_cytokines.yml")
cfg <- read_project_yaml("config/parse_cytokines_rest_only.yml")
ref_cfg <- read_project_yaml("config/bdev_reference.yml")
runtime <- registry$bundles$parse_cytokines$runtime
stopifnot(identical(as.character(packageVersion("scworkbench")),
  runtime$scworkbench_version), identical(as.character(packageVersion("UCell")),
  "2.12.0"), is.finite(ncores), ncores >= 1L)

source_dir <- parent_cfg$signature_analysis$discovery$report
input_paths <- c(
  fits = file.path(source_dir, "limma_voom_fits.rds"),
  parent_config = "config/parse_cytokines.yml",
  review_config = "config/parse_cytokines_rest_only.yml",
  reference_config = "config/bdev_reference.yml",
  registry = "config/bundles.yml",
  methods = "processing/parse_cytokines/COMPETITOR_REFINEMENT_PILOT_WORKFLOW.md",
  script = "processing/parse_cytokines/17_flt3_il7_competitor_refinement_pilot.R",
  parse_manifest = file.path(cfg$bundle_dir, "manifest.json"),
  bdev_manifest = file.path(registry$bundles$bdev_reference$bundle_dir,
    "manifest.json")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) stop("Missing inputs: ", paste(missing_inputs,
  collapse = ", "))

for (id in c("parse_cytokines", "parse_cytokines_rest_only", "bdev_reference")) {
  checked <- validate_registered_bundle(registry, id,
    run_recipient_validation = FALSE)
  stopifnot(isTRUE(checked$identity_ok), isTRUE(checked$fingerprint_ok))
  validation_path <- if (id == "bdev_reference") ref_cfg$validation$report else
    registry$bundles[[id]]$validation_report
  stopifnot(isTRUE(jsonlite::read_json(validation_path,
    simplifyVector = TRUE)$ok))
}

# Bdev supplies the SCT_counts feature universe and the final display matrix.
bdev <- wb_open_bundle(registry$bundles$bdev_reference$bundle_dir,
  read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
bdev_manifest <- wb_manifest(bdev)
bdev_fingerprint <- bdev_manifest$current_state$fingerprint
bdev_sct_count_features <- as.character(wb_feature_ids(bdev, "SCT_counts"))
bdev_sct_data_features <- as.character(wb_feature_ids(bdev, "SCT_data"))
stopifnot(length(bdev_sct_count_features) == 20533L,
  !anyDuplicated(bdev_sct_count_features), !anyDuplicated(bdev_sct_data_features))

fits <- readRDS(input_paths[["fits"]])
features <- as.character(fits$filtered_genes)
stopifnot(identical(rownames(fits$full$coefficients), features),
  all(vapply(fits$folds[folds], function(x) identical(
    rownames(x$coefficients), features), logical(1))))
bdev_eligible <- features %in% bdev_sct_count_features

rank_genes <- function(stats, target, direction) {
  sign <- if (direction == "up") 1 else -1
  column <- paste0("specificity::", target)
  coefficient <- as.numeric(stats$coefficients[, column])
  statistic <- as.numeric(stats$t[, column])
  p_value <- as.numeric(stats$p_value[, column])
  fdr <- as.numeric(stats$adjusted_p_value[, column])
  oriented_log2fc <- sign * coefficient
  oriented_t <- sign * statistic
  keep <- bdev_eligible & is.finite(oriented_log2fc) & is.finite(oriented_t) &
    oriented_log2fc > log2fc_minimum
  index <- which(keep)
  strict <- is.finite(fdr[index]) & fdr[index] < fdr_maximum
  fdr_order <- fdr[index]
  fdr_order[!is.finite(fdr_order)] <- Inf
  index <- index[order(-as.integer(strict), -oriented_t[index],
    -oriented_log2fc[index], fdr_order, features[index])]
  data.frame(feature_id = features[index], rank = seq_along(index),
    specificity_log2fc = coefficient[index], oriented_log2fc = oriented_log2fc[index],
    specificity_t = statistic[index], oriented_t = oriented_t[index],
    specificity_p_value = p_value[index], specificity_fdr = fdr[index],
    strict_significant = is.finite(fdr[index]) & fdr[index] < fdr_maximum,
    stringsAsFactors = FALSE)
}

direction_rows <- list()
fold_rankings <- list()
candidate_sets <- list()
candidate_members <- list()
pool_rows <- list()
row_counter <- 0L
for (target in targets) for (direction in directions) {
  row_counter <- row_counter + 1L
  signature_id <- paste(namespace, slug(target), direction, sep = ".")
  local_folds <- setNames(lapply(folds, function(fold) {
    rank_genes(fits$folds[[fold]], target, direction)
  }), folds)
  fold_rankings[[signature_id]] <- local_folds
  top50 <- lapply(local_folds, function(x) head(x$feature_id, candidate_top_n))
  top100 <- lapply(local_folds, function(x) head(x$feature_id, cache_top_n))
  candidate <- Reduce(intersect, top50)
  fold_rank_matrix <- vapply(folds, function(fold) match(candidate,
    local_folds[[fold]]$feature_id), integer(length(candidate)))
  candidate_order <- data.frame(feature_id = candidate,
    mean_fold_rank = rowMeans(fold_rank_matrix),
    max_fold_rank = apply(fold_rank_matrix, 1L, max),
    min_fold_rank = apply(fold_rank_matrix, 1L, min),
    stringsAsFactors = FALSE)
  candidate_order <- candidate_order[order(candidate_order$mean_fold_rank,
    candidate_order$max_fold_rank, candidate_order$feature_id), ]
  candidate <- candidate_order$feature_id
  construction_ok <- length(candidate) >= minimum_size
  if (construction_ok) candidate_sets[[signature_id]] <- candidate
  direction_rows[[row_counter]] <- data.frame(signature_id = signature_id,
    cytokine = target, direction = direction,
    fold1_eligible = nrow(local_folds[["1"]]),
    fold2_eligible = nrow(local_folds[["2"]]),
    fold3_eligible = nrow(local_folds[["3"]]),
    fold1_top100_size = length(top100[["1"]]),
    fold2_top100_size = length(top100[["2"]]),
    fold3_top100_size = length(top100[["3"]]),
    candidate_size = length(candidate), construction_ok = construction_ok,
    construction_status = if (construction_ok) "PASS" else
      "FAIL_FEWER_THAN_20_IN_THREE_FOLD_TOP50_INTERSECTION",
    stringsAsFactors = FALSE)
  if (length(candidate)) {
    local <- candidate_order
    local$signature_id <- signature_id
    local$cytokine <- target
    local$direction <- direction
    local$candidate_rank <- seq_len(nrow(local))
    candidate_members[[length(candidate_members) + 1L]] <- local[, c(
      "signature_id", "cytokine", "direction", "candidate_rank", "feature_id",
      "mean_fold_rank", "max_fold_rank", "min_fold_rank")]
  }
  for (fold in folds) {
    x <- local_folds[[fold]]
    use <- seq_len(min(cache_top_n, nrow(x)))
    if (length(use)) {
      local <- x[use, ]
      local$signature_id <- signature_id
      local$cytokine <- target
      local$direction <- direction
      local$fold <- fold
      pool_rows[[length(pool_rows) + 1L]] <- local[, c("signature_id",
        "cytokine", "direction", "fold", "rank", "feature_id",
        "strict_significant", "oriented_log2fc", "oriented_t",
        "specificity_fdr")]
    }
  }
}
direction_inventory <- do.call(rbind, direction_rows)
candidate_membership <- do.call(rbind, candidate_members)
top100_membership <- do.call(rbind, pool_rows)
pool_genes <- sort(unique(top100_membership$feature_id))
stopifnot(nrow(direction_inventory) == 4L, length(pool_genes) > 0L,
  all(unlist(candidate_sets, use.names = FALSE) %in% pool_genes),
  all(pool_genes %in% bdev_sct_count_features))
write_tsv(direction_inventory, file.path(stage_dir, "candidate_inventory.tsv"))
write_tsv(candidate_membership, file.path(stage_dir, "candidate_membership.tsv"))
write_tsv(top100_membership, file.path(stage_dir,
  "top100_fold_membership.tsv"))
write_tsv(data.frame(pool_rank = seq_along(pool_genes), feature_id = pool_genes),
  file.path(stage_dir, "top100_union_genes.tsv"))
log_event("Constructed", length(candidate_sets), "of 4 top", candidate_top_n,
  "-intersection candidates;",
  length(pool_genes), "genes in reusable top100 union")

# Parse metadata and exact sample-level UCell capped-rank contribution cache.
parse <- wb_open_bundle(cfg$bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(parse), silent = TRUE), add = TRUE)
parse_fingerprint <- wb_manifest(parse)$current_state$fingerprint
cells <- as.character(wb_cells(parse, "all"))
md <- wb_meta(parse, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "sample_id", "condition", "replicate",
  "treatment_category"))
stopifnot(length(cells) == 154209L, identical(as.character(md$cell_id), cells),
  !anyNA(md))
sample_md <- unique(md[, c("sample_id", "condition", "replicate",
  "treatment_category")])
sample_md <- sample_md[match(unique(md$sample_id), sample_md$sample_id), ]
stopifnot(nrow(sample_md) == 81L, all(table(sample_md$condition) == 3L),
  all(table(sample_md$replicate) == 27L))
sample_factor <- factor(md$sample_id, levels = sample_md$sample_id)
sample_cell_counts <- tabulate(as.integer(sample_factor), nbins = nrow(sample_md))
parse_features <- as.character(wb_feature_ids(parse, "RNA_counts"))
stopifnot(all(pool_genes %in% parse_features))

cache_rds <- file.path(stage_dir, "top100_union_sample_rank_cache.rds")
progress_rds <- file.path(stage_dir, "rank_cache_progress.rds")
if (file.exists(cache_rds)) {
  rank_cache <- readRDS(cache_rds)
  stopifnot(identical(rank_cache$genes, pool_genes),
    identical(rank_cache$samples, sample_md$sample_id),
    identical(rank_cache$sample_cell_counts, sample_cell_counts))
  log_event("Reusing completed sample-rank contribution cache")
} else {
  if (file.exists(progress_rds)) {
    progress <- readRDS(progress_rds)
    stopifnot(identical(progress$genes, pool_genes),
      identical(progress$samples, sample_md$sample_id),
      identical(progress$sample_cell_counts, sample_cell_counts))
    rank_sums <- progress$rank_sums
    next_cell <- progress$next_cell
    log_event("Resuming contribution cache at cell", next_cell)
  } else {
    rank_sums <- matrix(0, nrow = length(pool_genes), ncol = nrow(sample_md),
      dimnames = list(pool_genes, sample_md$sample_id))
    next_cell <- 1L
  }
  starts <- if (next_cell <= length(cells)) seq.int(next_cell, length(cells),
    by = rank_outer_chunk) else integer()
  for (chunk_index in seq_along(starts)) {
    start <- starts[[chunk_index]]
    end <- min(start + rank_outer_chunk - 1L, length(cells))
    chunk_cells <- cells[start:end]
    x <- wb_slice(parse, matrix_id = "RNA_counts", features = parse_features,
      cells = chunk_cells)
    x <- methods::as(x, "dgCMatrix")
    ranks <- UCell::StoreRankings_UCell(x, maxRank = max_rank,
      chunk.size = rank_inner_chunk, ncores = ncores, ties.method = "average")
    local <- as.matrix(ranks[pool_genes, , drop = FALSE])
    local[local == 0] <- max_rank
    stopifnot(identical(rownames(local), pool_genes),
      identical(colnames(local), chunk_cells), all(is.finite(local)),
      all(local >= 1 & local <= max_rank))
    observed_sums <- rowsum(t(local), md$sample_id[start:end], reorder = FALSE)
    local_sums <- matrix(0, nrow = length(pool_genes), ncol = nrow(sample_md),
      dimnames = list(pool_genes, sample_md$sample_id))
    local_sums[, match(rownames(observed_sums), sample_md$sample_id)] <-
      t(observed_sums)
    rank_sums <- rank_sums + local_sums
    next_value <- end + 1L
    if (chunk_index %% 10L == 0L || end == length(cells)) {
      saveRDS(list(genes = pool_genes, samples = sample_md$sample_id,
        sample_cell_counts = sample_cell_counts, rank_sums = rank_sums,
        next_cell = next_value), progress_rds, compress = FALSE)
      log_event("Rank contributions", end, "/", length(cells), "cells")
    }
    rm(x, ranks, local)
    if (chunk_index %% 10L == 0L) gc()
  }
  rank_means <- sweep(rank_sums, 2L, sample_cell_counts, "/")
  rank_cache <- list(genes = pool_genes, samples = sample_md$sample_id,
    sample_cell_counts = sample_cell_counts, mean_capped_rank = rank_means,
    max_rank = max_rank, ties_method = "average", ucell_version =
      as.character(packageVersion("UCell")))
  saveRDS(rank_cache, cache_rds, compress = "xz")
}
rank_means <- rank_cache$mean_capped_rank
rank_cache_table <- data.frame(feature_id = rep(rownames(rank_means),
  times = ncol(rank_means)), sample_id = rep(colnames(rank_means),
  each = nrow(rank_means)), mean_capped_rank = as.vector(rank_means),
  stringsAsFactors = FALSE)
write_tsv_gz(rank_cache_table, file.path(stage_dir,
  "top100_union_sample_rank_contributions.tsv.gz"))

score_samples_algebraic <- function(genes) {
  n <- length(genes)
  stopifnot(n >= 1L, all(genes %in% rownames(rank_means)))
  rank_sum <- colSums(rank_means[genes, , drop = FALSE])
  minimum <- n * (n + 1) / 2
  score <- 1 - (rank_sum - minimum) / (n * max_rank - minimum)
  names(score) <- colnames(rank_means)
  score
}

condition_metric <- function(sample_values, target, direction) {
  stopifnot(identical(names(sample_values), sample_md$sample_id))
  oriented <- if (direction == "up") sample_values else -sample_values
  medians <- vapply(unique(sample_md$condition), function(condition) {
    median(oriented[sample_md$condition == condition])
  }, numeric(1))
  target_value <- medians[[target]]
  competitors <- medians[names(medians) != target]
  competitor_value <- max(competitors)
  tied <- sort(names(competitors)[competitors == competitor_value])
  spread <- sd(medians)
  gap <- if (is.finite(spread) && spread > 0) {
    (target_value - competitor_value) / spread
  } else 0
  list(gap = gap, target = target_value, competitor = competitor_value,
    strongest = paste(tied, collapse = ";"), spread = spread,
    rank = 1L + sum(competitors > target_value), pass = target_value > competitor_value,
    condition_values = medians)
}

optimization_rows <- list()
final_sets <- list()
opt_counter <- 0L
for (i in seq_len(nrow(direction_inventory))) {
  info <- direction_inventory[i, ]
  id <- info$signature_id
  if (!isTRUE(info$construction_ok)) next
  genes <- candidate_sets[[id]]
  original_order <- setNames(seq_along(genes), genes)
  iteration <- 0L
  repeat {
    current_scores <- score_samples_algebraic(genes)
    current <- condition_metric(current_scores, info$cytokine, info$direction)
    if (length(genes) <= minimum_size) {
      opt_counter <- opt_counter + 1L
      optimization_rows[[opt_counter]] <- data.frame(signature_id = id,
        cytokine = info$cytokine, direction = info$direction,
        iteration = iteration, set_size_before = length(genes),
        gap_before = current$gap, tested_gene = NA_character_,
        deletion_gap = NA_real_, gap_improvement = NA_real_,
        selected_for_removal = FALSE, stop_reason = "minimum_size_reached",
        strongest_competitor_before = current$strongest,
        deletion_strongest_competitor = NA_character_,
        stringsAsFactors = FALSE)
      break
    }
    trials <- lapply(genes, function(gene) {
      trial <- condition_metric(score_samples_algebraic(setdiff(genes, gene)),
        info$cytokine, info$direction)
      data.frame(tested_gene = gene, deletion_gap = trial$gap,
        gap_improvement = trial$gap - current$gap,
        original_candidate_rank = original_order[[gene]],
        deletion_strongest_competitor = trial$strongest,
        stringsAsFactors = FALSE)
    })
    trials <- do.call(rbind, trials)
    best_improvement <- max(trials$gap_improvement)
    selectable <- which(trials$gap_improvement == best_improvement)
    selected <- selectable[which.min(trials$original_candidate_rank[selectable])]
    remove_gene <- trials$tested_gene[[selected]]
    improves <- is.finite(best_improvement) && best_improvement > 0
    for (j in seq_len(nrow(trials))) {
      opt_counter <- opt_counter + 1L
      optimization_rows[[opt_counter]] <- data.frame(signature_id = id,
        cytokine = info$cytokine, direction = info$direction,
        iteration = iteration, set_size_before = length(genes),
        gap_before = current$gap, tested_gene = trials$tested_gene[[j]],
        deletion_gap = trials$deletion_gap[[j]],
        gap_improvement = trials$gap_improvement[[j]],
        selected_for_removal = improves && j == selected,
        stop_reason = if (!improves && j == selected) "no_positive_improvement" else NA_character_,
        strongest_competitor_before = current$strongest,
        deletion_strongest_competitor = trials$deletion_strongest_competitor[[j]],
        stringsAsFactors = FALSE)
    }
    if (!improves) break
    genes <- genes[genes != remove_gene]
    iteration <- iteration + 1L
  }
  final_sets[[id]] <- genes
}
optimization_trace <- if (length(optimization_rows)) do.call(rbind,
  optimization_rows) else data.frame()
write_tsv(optimization_trace, file.path(stage_dir, "optimization_trace.tsv"))

final_membership <- do.call(rbind, lapply(names(final_sets), function(id) {
  info <- direction_inventory[direction_inventory$signature_id == id, ]
  genes <- final_sets[[id]]
  candidate <- candidate_membership[candidate_membership$signature_id == id, ]
  data.frame(signature_id = id, cytokine = info$cytokine,
    direction = info$direction, final_rank = seq_along(genes),
    candidate_rank = candidate$candidate_rank[match(genes, candidate$feature_id)],
    feature_id = genes, present_in_bdev_sct_counts = TRUE,
    present_in_bdev_sct_data = genes %in% bdev_sct_data_features,
    stringsAsFactors = FALSE)
}))
if (is.null(final_membership)) final_membership <- data.frame()
write_tsv(final_membership, file.path(stage_dir, "final_gene_membership.tsv"))

# Native scWorkbench UCell check for all constructed initial/final sets.
native_signatures <- list()
native_map_rows <- list()
for (id in names(candidate_sets)) {
  info <- direction_inventory[direction_inventory$signature_id == id, ]
  for (stage in c("initial", "final")) {
    key <- paste(id, stage, sep = ".")
    genes <- if (stage == "initial") candidate_sets[[id]] else final_sets[[id]]
    native_signatures[[key]] <- genes
    native_map_rows[[length(native_map_rows) + 1L]] <- data.frame(score_key = key,
      signature_id = id, cytokine = info$cytokine, direction = info$direction,
      set_stage = stage, set_size = length(genes), stringsAsFactors = FALSE)
  }
}
native_map <- do.call(rbind, native_map_rows)
native_session_workspace <- NA_character_
native_summary <- data.frame()
native_sample_scores <- data.frame()
native_condition_results <- data.frame()
if (length(native_signatures)) {
  session_id <- paste0("flt3_il7_refinement_top", candidate_top_n,
    "_native_v1")
  session_root <- file.path(cfg$storage_root, "sessions", cfg$run_id)
  workspace <- file.path(session_root, session_id)
  session <- if (dir.exists(workspace)) wb_session_resume(parse, workspace) else
    wb_session_create(parse, workspace_root = session_root, session_id = session_id)
  if (!is.null(session$durable$promotion_receipt)) stop("Native review session was promoted")
  score_runs <- names(session$runs)[vapply(session$runs, function(x)
    identical(x$op_name, "ucell"), logical(1))]
  matching <- score_runs[vapply(score_runs, function(id) {
    tab <- session$runs[[id]]$outputs$signature_summary$table
    setequal(as.character(tab$signature), names(native_signatures))
  }, logical(1))]
  if (!length(matching)) {
    if (length(score_runs)) stop("Unexpected UCell run in native session")
    log_event("Running native scWorkbench UCell for", length(native_signatures),
      "initial/final sets")
    result <- wb_session_run_op(session, parse, "ucell", inputs = list(
      matrix = "RNA_counts", cell_set = "all", cell_ids = cells,
      min_cells = 10L, features = native_signatures, name = session_id,
      max_rank = max_rank, ncores = ncores, parallel_backend = "auto",
      chunk_size = rank_outer_chunk, missing_features_policy = "error"),
      notes = "Native verification and responder scoring for unpromoted Flt3L/IL-7 refinement pilot")
    session <- result$session
    run <- result$run
  } else {
    stopifnot(length(matching) == 1L)
    run <- session$runs[[matching[[1L]]]]
    log_event("Reusing completed native UCell run")
  }
  native_session_workspace <- workspace
  native_summary <- run$outputs$signature_summary$table
  score <- run$outputs$signature_scores
  stopifnot(length(score$value) == length(cells) * length(native_signatures),
    all(is.finite(score$value)), setequal(unique(score$key), names(native_signatures)))
  sample_rows <- list()
  condition_rows <- list()
  key_rle <- rle(as.character(score$key))
  ends <- cumsum(key_rle$lengths)
  starts <- ends - key_rle$lengths + 1L
  stopifnot(all(key_rle$lengths == length(cells)))
  for (j in seq_along(key_rle$values)) {
    key <- key_rle$values[[j]]
    map <- native_map[native_map$score_key == key, ]
    values <- as.numeric(score$value[starts[[j]]:ends[[j]]])
    score_cells <- as.character(score$cell_id[starts[[j]]:ends[[j]]])
    stopifnot(identical(score_cells, cells), nrow(map) == 1L)
    sample_mean <- as.numeric(rowsum(values, sample_factor,
      reorder = FALSE)) / sample_cell_counts
    thresholds <- setNames(vapply(folds, function(replicate) {
      use <- md$condition == "control" & as.character(md$replicate) == replicate
      probability <- if (map$direction == "up") 0.95 else 0.05
      as.numeric(quantile(values[use], probability, names = FALSE, type = 7))
    }, numeric(1)), folds)
    responding <- if (map$direction == "up") {
      values > thresholds[as.character(md$replicate)]
    } else {
      values < thresholds[as.character(md$replicate)]
    }
    responder_fraction <- as.numeric(rowsum(as.numeric(responding), sample_factor,
      reorder = FALSE)) / sample_cell_counts
    algebraic <- score_samples_algebraic(native_signatures[[key]])
    max_error <- max(abs(sample_mean - algebraic[sample_md$sample_id]))
    stopifnot(max_error < 1e-10)
    local_samples <- sample_md
    local_samples$score_key <- key
    local_samples$signature_id <- map$signature_id
    local_samples$cytokine <- map$cytokine
    local_samples$direction <- map$direction
    local_samples$set_stage <- map$set_stage
    local_samples$set_size <- map$set_size
    local_samples$native_sample_mean <- sample_mean
    local_samples$algebraic_sample_mean <- algebraic[sample_md$sample_id]
    local_samples$absolute_error <- abs(local_samples$native_sample_mean -
      local_samples$algebraic_sample_mean)
    local_samples$responder_fraction <- responder_fraction
    local_samples$responder_threshold <- thresholds[
      as.character(local_samples$replicate)]
    sample_rows[[j]] <- local_samples
    mean_metric <- condition_metric(setNames(sample_mean, sample_md$sample_id),
      map$cytokine, map$direction)
    responder_values <- vapply(unique(sample_md$condition), function(condition) {
      median(responder_fraction[sample_md$condition == condition])
    }, numeric(1))
    target_value <- responder_values[[map$cytokine]]
    competitors <- responder_values[names(responder_values) != map$cytokine]
    competitor_value <- max(competitors)
    spread <- sd(responder_values)
    responder_gap <- if (is.finite(spread) && spread > 0) {
      (target_value - competitor_value) / spread
    } else 0
    condition_rows[[j]] <- data.frame(score_key = key,
      signature_id = map$signature_id, cytokine = map$cytokine,
      direction = map$direction, set_stage = map$set_stage,
      set_size = map$set_size, native_algebraic_max_abs_sample_mean_error = max_error,
      mean_target_response = mean_metric$target,
      mean_strongest_competitor_response = mean_metric$competitor,
      mean_strongest_competitor = mean_metric$strongest,
      mean_response_sd_across_27_conditions = mean_metric$spread,
      mean_standardized_gap = mean_metric$gap,
      mean_target_rank = mean_metric$rank, mean_pass = mean_metric$pass,
      responder_target_response = target_value,
      responder_strongest_competitor_response = competitor_value,
      responder_strongest_competitor = paste(sort(names(competitors)[
        competitors == competitor_value]), collapse = ";"),
      responder_response_sd_across_27_conditions = spread,
      responder_standardized_gap = responder_gap,
      responder_target_rank = 1L + sum(competitors > target_value),
      responder_pass = target_value > competitor_value,
      stringsAsFactors = FALSE)
  }
  native_sample_scores <- do.call(rbind, sample_rows)
  native_condition_results <- do.call(rbind, condition_rows)
  native_condition_results$overall_pass <- with(native_condition_results,
    mean_pass | responder_pass)
}
write_tsv(native_summary, file.path(stage_dir, "native_ucell_summary.tsv"))
write_tsv(native_sample_scores, file.path(stage_dir, "sample_metrics.tsv"))
write_tsv(native_condition_results, file.path(stage_dir,
  "initial_final_metrics.tsv"))

# One-row-per-direction result table, including construction failures.
summary_rows <- lapply(seq_len(nrow(direction_inventory)), function(i) {
  info <- direction_inventory[i, ]
  id <- info$signature_id
  if (!isTRUE(info$construction_ok)) {
    return(data.frame(signature_id = id, cytokine = info$cytokine,
      direction = info$direction, candidate_size = info$candidate_size,
      construction_ok = FALSE, final_size = NA_integer_, genes_removed = NA_integer_,
      initial_mean_gap = NA_real_, final_mean_gap = NA_real_,
      mean_gap_improvement = NA_real_, final_mean_strongest_competitor = NA_character_,
      final_mean_rank = NA_integer_, final_mean_pass = NA,
      final_responder_gap = NA_real_, final_responder_strongest_competitor = NA_character_,
      final_responder_rank = NA_integer_, final_responder_pass = NA,
      final_overall_pass = FALSE, final_genes = NA_character_,
      removed_genes = NA_character_, stringsAsFactors = FALSE))
  }
  initial <- native_condition_results[native_condition_results$signature_id == id &
    native_condition_results$set_stage == "initial", ]
  final <- native_condition_results[native_condition_results$signature_id == id &
    native_condition_results$set_stage == "final", ]
  initial_genes <- candidate_sets[[id]]
  final_genes <- final_sets[[id]]
  removed <- initial_genes[!initial_genes %in% final_genes]
  stopifnot(nrow(initial) == 1L, nrow(final) == 1L)
  data.frame(signature_id = id, cytokine = info$cytokine,
    direction = info$direction, candidate_size = length(initial_genes),
    construction_ok = TRUE, final_size = length(final_genes),
    genes_removed = length(removed), initial_mean_gap = initial$mean_standardized_gap,
    final_mean_gap = final$mean_standardized_gap,
    mean_gap_improvement = final$mean_standardized_gap - initial$mean_standardized_gap,
    final_mean_strongest_competitor = final$mean_strongest_competitor,
    final_mean_rank = final$mean_target_rank, final_mean_pass = final$mean_pass,
    final_responder_gap = final$responder_standardized_gap,
    final_responder_strongest_competitor = final$responder_strongest_competitor,
    final_responder_rank = final$responder_target_rank,
    final_responder_pass = final$responder_pass,
    final_overall_pass = final$overall_pass,
    final_genes = paste(final_genes, collapse = ";"),
    removed_genes = paste(removed, collapse = ";"), stringsAsFactors = FALSE)
})
signature_summary <- do.call(rbind, summary_rows)
write_tsv(signature_summary, file.path(stage_dir, "signature_summary.tsv"))

# Bdev SCT_data expression on the existing ATAC UMAP.
figure_dir <- file.path(stage_dir, "figures")
dir.create(figure_dir, showWarnings = FALSE)
bdev_cells <- as.character(wb_cells(bdev, "all"))
space <- wb_space(bdev, "umap.atac", cells = bdev_cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), bdev_cells))
  coordinates <- data.frame(cell_id = bdev_cells, x = as.numeric(space$x),
    y = as.numeric(space$y), stringsAsFactors = FALSE)
} else {
  stopifnot(identical(rownames(space), bdev_cells), ncol(space) == 2L)
  coordinates <- data.frame(cell_id = bdev_cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]), stringsAsFactors = FALSE)
}
plot_members <- final_membership
available_genes <- unique(plot_members$feature_id[
  plot_members$present_in_bdev_sct_data])
expression <- if (length(available_genes)) {
  as.matrix(wb_matrix_read_features(bdev, "SCT_data", available_genes,
    bdev_cells, missing = "NA"))
} else matrix(numeric(), nrow = 0L, ncol = length(bdev_cells))
stopifnot(identical(rownames(expression), available_genes),
  identical(colnames(expression), bdev_cells))

expression_rows <- list()
plot_values <- list()
for (i in seq_len(nrow(plot_members))) {
  gene <- plot_members$feature_id[[i]]
  available <- gene %in% available_genes
  values <- if (available) as.numeric(expression[gene, ]) else
    rep(NA_real_, length(bdev_cells))
  finite <- is.finite(values)
  detected <- finite & values > 0
  positive <- values[detected]
  cap <- if (length(positive)) as.numeric(quantile(positive, 0.99,
    names = FALSE, type = 7)) else NA_real_
  scaled <- rep(NA_real_, length(values))
  if (is.finite(cap) && cap > 0) scaled[finite] <- pmin(pmax(values[finite], 0), cap) / cap
  key <- paste(plot_members$signature_id[[i]], gene, sep = "::")
  plot_values[[key]] <- scaled
  expression_rows[[i]] <- data.frame(signature_id = plot_members$signature_id[[i]],
    cytokine = plot_members$cytokine[[i]], direction = plot_members$direction[[i]],
    final_rank = plot_members$final_rank[[i]], feature_id = gene,
    present_in_sct_data = available, cells = length(values),
    detected_cells = sum(detected), detected_percent = 100 * mean(detected),
    sct_data_min = if (any(finite)) min(values[finite]) else NA_real_,
    sct_data_median = if (any(finite)) median(values[finite]) else NA_real_,
    sct_data_max = if (any(finite)) max(values[finite]) else NA_real_,
    positive_p99_display_cap = cap, stringsAsFactors = FALSE)
}
expression_summary <- if (length(expression_rows)) do.call(rbind,
  expression_rows) else data.frame()
write_tsv(expression_summary, file.path(stage_dir,
  "bdev_sct_data_expression_summary.tsv"))

make_group_plot <- function(cytokine, direction, columns = 5L) {
  members <- plot_members[plot_members$cytokine == cytokine &
    plot_members$direction == direction, ]
  if (!nrow(members)) return(ggplot() + theme_void() + labs(title = paste(cytokine,
    direction, "— construction failed")))
  summaries <- expression_summary[match(paste(members$signature_id, members$feature_id),
    paste(expression_summary$signature_id, expression_summary$feature_id)), ]
  labels <- ifelse(summaries$present_in_sct_data,
    paste0(summaries$feature_id, "\n", sprintf("%.1f%% detected",
      summaries$detected_percent)), paste0(summaries$feature_id,
        "\nnot in SCT_data"))
  rows <- lapply(seq_len(nrow(members)), function(i) {
    key <- paste(members$signature_id[[i]], members$feature_id[[i]], sep = "::")
    data.frame(x = coordinates$x, y = coordinates$y,
      relative_expression = plot_values[[key]], gene_panel = labels[[i]],
      panel_rank = i, stringsAsFactors = FALSE)
  })
  plot_data <- do.call(rbind, rows)
  plot_data <- plot_data[order(plot_data$panel_rank,
    ifelse(is.na(plot_data$relative_expression), -Inf,
      plot_data$relative_expression)), ]
  plot_data$gene_panel <- factor(plot_data$gene_panel, levels = labels)
  ggplot(plot_data, aes(x, y, color = relative_expression)) +
    geom_point(size = 0.10, alpha = 0.9, stroke = 0) +
    facet_wrap(~gene_panel, ncol = columns) +
    scale_color_viridis_c(option = "magma", limits = c(0, 1),
      oob = scales::squish, na.value = "#E5E7EB", breaks = c(0, 0.5, 1),
      labels = c("low", "mid", "p99+"), name = "Within-gene\nSCT_data") +
    coord_equal() + labs(title = paste(cytokine, direction, "final refined set")) +
    theme_void(base_size = 8) + theme(strip.text = element_text(size = 7.3,
      face = "bold", lineheight = 0.95), plot.title = element_text(size = 12,
      face = "bold"), legend.position = "bottom")
}

figure_files <- character()
for (direction in directions) {
  p1 <- make_group_plot("Flt3L", direction)
  p2 <- make_group_plot("IL-7", direction)
  combined <- (p1 / p2) + plot_layout(guides = "collect") +
    plot_annotation(title = paste0("Flt3L and IL-7 refined genes — ", direction),
      subtitle = paste("Bdev SCT_data on existing ATAC-derived UMAP;",
        "panels follow final set order"),
      caption = paste("Values are read directly from SCT_data and scaled within",
        "each gene to the positive-cell 99th percentile. Gray denotes zero,",
        "non-positive, non-finite, or unavailable values."),
      theme = theme(plot.title = element_text(size = 17, face = "bold"),
        plot.subtitle = element_text(size = 11),
        plot.caption = element_text(size = 8, hjust = 0))) &
    theme(legend.position = "bottom")
  counts <- table(plot_members$cytokine[plot_members$direction == direction])
  max_panels <- if (length(counts)) max(counts) else 1L
  height <- max(12, 4 + 2.35 * ceiling(max_panels / 5))
  path <- file.path(figure_dir, paste0("flt3_il7_refined_", direction,
    "_sct_data_atac_umap.png"))
  ggsave(path, combined, width = 15, height = height, dpi = 220, units = "in",
    bg = "white", limitsize = FALSE)
  figure_files[direction] <- path
}

stopifnot(identical(parse_fingerprint,
  wb_manifest(parse)$current_state$fingerprint), identical(bdev_fingerprint,
  wb_manifest(bdev)$current_state$fingerprint),
  all(file.exists(figure_files)))

report <- c("# Flt3L/IL-7 competitor-refinement pilot", "",
  paste0("This is a separate, unpromoted top-", candidate_top_n,
    " pilot. Per-gene rank contributions were cached for the union of all fold-specific top-100 lists."), "",
  paste0("- Reusable top-100 union: **", length(pool_genes), " genes**."),
  paste0("- Constructed top-", candidate_top_n,
    "-intersection candidates: **", length(candidate_sets),
    "/4**."),
  paste0("- Final overall passes: **", sum(signature_summary$final_overall_pass),
    "/4**."), "",
  "| Cytokine | Direction | Candidate | Final | Removed | Initial mean gap | Final mean gap | Mean competitor | Mean | Responder | Overall |",
  "|---|---|---:|---:|---:|---:|---:|---|---|---|---|")
for (i in seq_len(nrow(signature_summary))) {
  x <- signature_summary[i, ]
  report <- c(report, paste0("| ", x$cytokine, " | ", x$direction, " | ",
    x$candidate_size, " | ", ifelse(is.na(x$final_size), "-", x$final_size),
    " | ", ifelse(is.na(x$genes_removed), "-", x$genes_removed), " | ",
    ifelse(is.na(x$initial_mean_gap), "-", sprintf("%.4f", x$initial_mean_gap)),
    " | ", ifelse(is.na(x$final_mean_gap), "-", sprintf("%.4f", x$final_mean_gap)),
    " | ", ifelse(is.na(x$final_mean_strongest_competitor), "-",
      x$final_mean_strongest_competitor), " | ",
    ifelse(isTRUE(x$final_mean_pass), "PASS", "FAIL"), " | ",
    ifelse(isTRUE(x$final_responder_pass), "PASS", "FAIL"), " | **",
    ifelse(isTRUE(x$final_overall_pass), "PASS", "FAIL"), "** |"))
}
report <- c(report, "", "## Outputs", "",
  "- `signature_summary.tsv`: one-row-per-direction final metrics.",
  "- `initial_final_metrics.tsv`: native initial and final mean/responder metrics.",
  "- `candidate_membership.tsv` and `final_gene_membership.tsv`: exact memberships.",
  "- `optimization_trace.tsv`: every tested deletion and selected removal.",
  "- `top100_union_sample_rank_cache.rds`: reusable sample-level capped-rank cache.",
  "- `top100_union_sample_rank_contributions.tsv.gz`: tabular cache export.",
  "- `figures/`: Bdev SCT_data expression on the existing ATAC UMAP.", "",
  "Refinement and final testing use the same complete dataset, so pass status is an internally optimized discrimination result.")
writeLines(report, file.path(stage_dir, "REPORT.md"))
file.copy(input_paths[["methods"]], file.path(stage_dir, "METHODS.md"),
  overwrite = TRUE)

input_manifest <- data.frame(path = normalizePath(input_paths, mustWork = TRUE),
  sha256 = vapply(input_paths, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))
completion <- list(status = "complete", promoted = FALSE,
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  analysis_id = analysis_id, namespace = namespace,
  parse_bundle = cfg$bundle_dir, parse_bundle_fingerprint = parse_fingerprint,
  parse_bundle_unchanged = TRUE,
  bdev_bundle = registry$bundles$bdev_reference$bundle_dir,
  bdev_bundle_fingerprint = bdev_fingerprint, bdev_bundle_unchanged = TRUE,
  bdev_feature_universe = "SCT_counts", bdev_expression_matrix = "SCT_data",
  targets = targets, directions = directions, candidate_top_n = candidate_top_n,
  cached_top_n = cache_top_n, top100_union_genes = length(pool_genes),
  minimum_final_size = minimum_size, constructed_candidates = length(candidate_sets),
  final_passes = sum(signature_summary$final_overall_pass),
  native_session_workspace = native_session_workspace,
  native_max_abs_sample_mean_error = if (nrow(native_condition_results))
    max(native_condition_results$native_algebraic_max_abs_sample_mean_error) else NA,
  figures = as.list(sub(paste0("^", stage_dir, "/"), "", figure_files)),
  result_promotions = 0L, runtime = runtime)
write_json(completion, file.path(stage_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

output_files <- list.files(stage_dir, recursive = TRUE, full.names = TRUE)
output_files <- output_files[basename(output_files) != "output_manifest.tsv"]
output_manifest <- data.frame(file = substring(output_files,
  nchar(stage_dir) + 2L), bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE)
write_tsv(output_manifest[order(output_manifest$file), ],
  file.path(stage_dir, "output_manifest.tsv"))

if (!file.rename(stage_dir, output_dir)) stop("Failed to finalize output directory")
log_event("Completed unpromoted Flt3L/IL-7 refinement pilot at", output_dir)
