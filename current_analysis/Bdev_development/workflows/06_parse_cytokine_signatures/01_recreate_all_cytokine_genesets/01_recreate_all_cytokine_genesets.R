#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)

parse_args <- function(args) {
  if (length(args) %% 2L != 0L) stop("Arguments must be --name value pairs")
  out <- list()
  for (i in seq.int(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[[i]])
    if (identical(key, args[[i]])) stop("Expected --name, got: ", args[[i]])
    out[[gsub("-", "_", key)]] <- args[[i + 1L]]
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("fits", "bdev_sct_count_features", "bdev_sct_data_features",
  "rank_cache", "sample_metrics", "reference_report", "output")
missing_args <- required_args[!required_args %in% names(args)]
if (length(missing_args)) stop("Missing arguments: ", paste(missing_args, collapse = ", "))

input_paths <- c(
  fits = args$fits,
  bdev_sct_count_features = args$bdev_sct_count_features,
  bdev_sct_data_features = args$bdev_sct_data_features,
  rank_cache = args$rank_cache,
  sample_metrics = args$sample_metrics,
  reference_candidate_inventory = file.path(args$reference_report, "candidate_inventory.tsv"),
  reference_candidate_membership = file.path(args$reference_report, "candidate_membership.tsv"),
  reference_final_membership = file.path(args$reference_report, "final_gene_membership.tsv"),
  reference_signature_summary = file.path(args$reference_report, "signature_summary.tsv")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) stop("Missing inputs: ", paste(missing_inputs, collapse = ", "))

output_dir <- normalizePath(args$output, mustWork = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(output_dir) || dir.exists(output_dir)) stop("Output already exists: ", output_dir)
if (file.exists(stage_dir) || dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(file.path(stage_dir, "01_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage_dir, "03_reports"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage_dir, "04_logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage_dir, "06_manifests"), recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) read.delim(path, check.names = FALSE,
  na.strings = "NA", stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_tsv_gz <- function(x, path) {
  connection <- gzfile(path, "wt")
  on.exit(close(connection), add = TRUE)
  write.table(x, connection, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
sha256 <- function(path) {
  value <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(value, "status")
  if (!length(value) || (length(status) && status != 0L)) stop("sha256sum failed: ", path)
  sub("[[:space:]].*$", "", value[[1L]])
}
slug <- function(x) gsub("(^_+|_+$)", "", tolower(gsub("[^A-Za-z0-9]+", "_", x)))
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}

namespace <- "parse_cytokine_rest_only_refinement_all_top100_v1"
directions <- c("up", "down")
folds <- c("1", "2", "3")
candidate_top_n <- 100L
minimum_size <- 20L
log2fc_minimum <- 0.5
fdr_maximum <- 0.05
max_rank <- 1500L

log_event("Loading frozen discovery statistics and feature universes")
fits <- readRDS(input_paths[["fits"]])
features <- as.character(fits$filtered_genes)
bdev_sct_count_features <- readLines(input_paths[["bdev_sct_count_features"]], warn = FALSE)
bdev_sct_data_features <- readLines(input_paths[["bdev_sct_data_features"]], warn = FALSE)
stopifnot(length(bdev_sct_count_features) == 20533L,
  !anyDuplicated(bdev_sct_count_features), !anyDuplicated(bdev_sct_data_features),
  identical(rownames(fits$full$coefficients), features),
  all(vapply(fits$folds[folds], function(x) identical(
    rownames(x$coefficients), features), logical(1))))
bdev_eligible <- features %in% bdev_sct_count_features
targets <- sub("^specificity::", "", grep("^specificity::",
  colnames(fits$full$coefficients), value = TRUE))
stopifnot(length(targets) == 26L, !anyDuplicated(targets), !"control" %in% targets)

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

log_event("Reconstructing three-fold top-100 intersection candidates")
direction_rows <- list()
candidate_sets <- list()
candidate_members <- list()
row_counter <- 0L
for (target in targets) for (direction in directions) {
  row_counter <- row_counter + 1L
  signature_id <- paste(namespace, slug(target), direction, sep = ".")
  local_folds <- setNames(lapply(folds, function(fold) {
    rank_genes(fits$folds[[fold]], target, direction)
  }), folds)
  top100 <- lapply(local_folds, function(x) head(x$feature_id, candidate_top_n))
  candidate <- Reduce(intersect, top100)
  fold_rank_matrix <- vapply(folds, function(fold) match(candidate,
    local_folds[[fold]]$feature_id), integer(length(candidate)))
  candidate_order <- data.frame(feature_id = candidate,
    mean_fold_rank = rowMeans(fold_rank_matrix),
    max_fold_rank = apply(fold_rank_matrix, 1L, max),
    min_fold_rank = apply(fold_rank_matrix, 1L, min), stringsAsFactors = FALSE)
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
      "FAIL_FEWER_THAN_20_IN_THREE_FOLD_TOP100_INTERSECTION",
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
}
direction_inventory <- do.call(rbind, direction_rows)
candidate_membership <- do.call(rbind, candidate_members)
stopifnot(nrow(direction_inventory) == 52L, length(candidate_sets) == 43L)

reference_inventory <- read_tsv(input_paths[["reference_candidate_inventory"]])
reference_candidates <- read_tsv(input_paths[["reference_candidate_membership"]])
rownames(direction_inventory) <- NULL
rownames(candidate_membership) <- NULL
rownames(reference_inventory) <- NULL
rownames(reference_candidates) <- NULL
inventory_columns <- intersect(names(reference_inventory), names(direction_inventory))
candidate_columns <- intersect(names(reference_candidates), names(candidate_membership))
inventory_equal <- isTRUE(all.equal(direction_inventory[, inventory_columns],
  reference_inventory[, inventory_columns], check.attributes = FALSE))
candidate_equal <- isTRUE(all.equal(candidate_membership[, candidate_columns],
  reference_candidates[, candidate_columns], check.attributes = FALSE))
if (!inventory_equal || !candidate_equal) {
  write_tsv(direction_inventory, file.path(stage_dir, "03_reports",
    "candidate_inventory_mismatch.tsv"))
  write_tsv(candidate_membership, file.path(stage_dir, "03_reports",
    "candidate_membership_mismatch.tsv"))
  stop("Candidate reconstruction mismatch: inventory=", inventory_equal,
    "; membership=", candidate_equal, "; membership comparison: ",
    paste(all.equal(candidate_membership[, candidate_columns],
      reference_candidates[, candidate_columns]), collapse = "; "))
}

log_event("Loading frozen full-precision sample-level capped-rank cache")
rank_cache <- readRDS(input_paths[["rank_cache"]])
required_cache_fields <- c("genes", "samples", "sample_cell_counts",
  "mean_capped_rank", "max_rank", "ties_method")
stopifnot(all(required_cache_fields %in% names(rank_cache)),
  identical(as.integer(rank_cache$max_rank), max_rank),
  identical(rank_cache$ties_method, "average"),
  identical(rownames(rank_cache$mean_capped_rank), rank_cache$genes),
  identical(colnames(rank_cache$mean_capped_rank), rank_cache$samples),
  all(is.finite(rank_cache$mean_capped_rank)))
rank_means <- rank_cache$mean_capped_rank
stopifnot(all(unlist(candidate_sets, use.names = FALSE) %in% rownames(rank_means)))
rm(rank_cache)

sample_table <- read_tsv(input_paths[["sample_metrics"]])
sample_columns <- c("sample_id", "condition", "replicate", "treatment_category")
stopifnot(all(sample_columns %in% names(sample_table)))
sample_md <- unique(sample_table[, sample_columns])
sample_md <- sample_md[match(colnames(rank_means), sample_md$sample_id), ]
stopifnot(nrow(sample_md) == 81L, !anyNA(sample_md),
  identical(sample_md$sample_id, colnames(rank_means)),
  all(table(sample_md$condition) == 3L))
rm(sample_table)

score_samples_algebraic <- function(genes) {
  n <- length(genes)
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
    rank = 1L + sum(competitors > target_value),
    pass = target_value > competitor_value)
}

log_event("Running greedy strongest-competitor refinement")
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
    current <- condition_metric(score_samples_algebraic(genes),
      info$cytokine, info$direction)
    if (length(genes) <= minimum_size) {
      opt_counter <- opt_counter + 1L
      optimization_rows[[opt_counter]] <- data.frame(signature_id = id,
        cytokine = info$cytokine, direction = info$direction,
        iteration = iteration, set_size_before = length(genes),
        gap_before = current$gap, tested_gene = NA_character_,
        deletion_gap = NA_real_, gap_improvement = NA_real_,
        selected_for_removal = FALSE, stop_reason = "minimum_size_reached",
        strongest_competitor_before = current$strongest,
        deletion_strongest_competitor = NA_character_, stringsAsFactors = FALSE)
      break
    }
    trials <- do.call(rbind, lapply(genes, function(gene) {
      trial <- condition_metric(score_samples_algebraic(setdiff(genes, gene)),
        info$cytokine, info$direction)
      data.frame(tested_gene = gene, deletion_gap = trial$gap,
        gap_improvement = trial$gap - current$gap,
        original_candidate_rank = original_order[[gene]],
        deletion_strongest_competitor = trial$strongest,
        stringsAsFactors = FALSE)
    }))
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
optimization_trace <- do.call(rbind, optimization_rows)

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

summary_rows <- lapply(seq_len(nrow(direction_inventory)), function(i) {
  info <- direction_inventory[i, ]
  id <- info$signature_id
  if (!isTRUE(info$construction_ok)) {
    return(data.frame(signature_id = id, cytokine = info$cytokine,
      direction = info$direction, candidate_size = info$candidate_size,
      construction_ok = FALSE, construction_status = info$construction_status,
      final_size = NA_integer_, genes_removed = NA_integer_,
      initial_mean_gap = NA_real_, final_mean_gap = NA_real_,
      mean_gap_improvement = NA_real_, final_mean_strongest_competitor = NA_character_,
      final_mean_rank = NA_integer_, final_mean_pass = FALSE,
      selected_for_collection = FALSE, final_genes = NA_character_,
      removed_genes = NA_character_, stringsAsFactors = FALSE))
  }
  initial_genes <- candidate_sets[[id]]
  final_genes <- final_sets[[id]]
  initial <- condition_metric(score_samples_algebraic(initial_genes),
    info$cytokine, info$direction)
  final <- condition_metric(score_samples_algebraic(final_genes),
    info$cytokine, info$direction)
  removed <- initial_genes[!initial_genes %in% final_genes]
  data.frame(signature_id = id, cytokine = info$cytokine,
    direction = info$direction, candidate_size = length(initial_genes),
    construction_ok = TRUE, construction_status = info$construction_status,
    final_size = length(final_genes), genes_removed = length(removed),
    initial_mean_gap = initial$gap, final_mean_gap = final$gap,
    mean_gap_improvement = final$gap - initial$gap,
    final_mean_strongest_competitor = final$strongest,
    final_mean_rank = final$rank, final_mean_pass = final$pass,
    selected_for_collection = final$pass,
    final_genes = paste(final_genes, collapse = ";"),
    removed_genes = paste(removed, collapse = ";"), stringsAsFactors = FALSE)
})
signature_summary <- do.call(rbind, summary_rows)
passing_ids <- signature_summary$signature_id[signature_summary$selected_for_collection]
passing_membership <- final_membership[final_membership$signature_id %in% passing_ids, ]
stopifnot(length(passing_ids) == 38L, length(unique(passing_membership$signature_id)) == 38L)

log_event("Comparing reconstructed collection with archived reference")
reference_final <- read_tsv(input_paths[["reference_final_membership"]])
reference_summary <- read_tsv(input_paths[["reference_signature_summary"]])
rownames(final_membership) <- NULL
rownames(reference_final) <- NULL
membership_columns <- c("signature_id", "cytokine", "direction", "final_rank",
  "candidate_rank", "feature_id", "present_in_bdev_sct_counts",
  "present_in_bdev_sct_data")
membership_comparison <- all.equal(final_membership[, membership_columns],
  reference_final[, membership_columns], check.attributes = FALSE)
if (!isTRUE(membership_comparison)) {
  write_tsv(final_membership, file.path(stage_dir, "03_reports",
    "final_membership_mismatch.tsv"))
  stop("Final membership mismatch: ", paste(membership_comparison, collapse = "; "))
}
summary_match <- match(signature_summary$signature_id, reference_summary$signature_id)
stopifnot(!anyNA(summary_match))
reference_summary <- reference_summary[summary_match, ]
comparison <- data.frame(signature_id = signature_summary$signature_id,
  construction_ok_match = signature_summary$construction_ok == reference_summary$construction_ok,
  final_size_match = (is.na(signature_summary$final_size) & is.na(reference_summary$final_size)) |
    signature_summary$final_size == reference_summary$final_size,
  final_genes_match = (is.na(signature_summary$final_genes) & is.na(reference_summary$final_genes)) |
    signature_summary$final_genes == reference_summary$final_genes,
  mean_gap_absolute_error = abs(signature_summary$final_mean_gap -
    reference_summary$final_mean_gap),
  mean_pass_match = ifelse(!signature_summary$construction_ok &
      !reference_summary$construction_ok, TRUE,
    signature_summary$final_mean_pass == reference_summary$final_mean_pass),
  selected_match = signature_summary$selected_for_collection ==
    reference_summary$final_overall_pass,
  reference_responder_pass = reference_summary$final_responder_pass,
  stringsAsFactors = FALSE)
comparison$mean_gap_absolute_error[is.na(comparison$mean_gap_absolute_error)] <- 0
stopifnot(all(comparison$construction_ok_match), all(comparison$final_size_match),
  all(comparison$final_genes_match), max(comparison$mean_gap_absolute_error) < 1e-12,
  all(comparison$mean_pass_match), all(comparison$selected_match),
  all(reference_summary$final_mean_pass[reference_summary$final_overall_pass]))

tables_dir <- file.path(stage_dir, "01_tables")
write_tsv(direction_inventory, file.path(tables_dir, "01_candidate_inventory.tsv"))
write_tsv(candidate_membership, file.path(tables_dir, "02_candidate_membership.tsv"))
write_tsv(final_membership, file.path(tables_dir, "03_final_gene_membership.tsv"))
write_tsv(passing_membership, file.path(tables_dir, "04_passing_gene_membership.tsv"))
writeLines(vapply(passing_ids, function(id) {
  local <- passing_membership[passing_membership$signature_id == id, ]
  description <- paste("Parse cytokine", local$cytokine[[1L]],
    local$direction[[1L]], "standalone top100 competitor refinement")
  paste(c(id, description, local$feature_id), collapse = "\t")
}, character(1)), file.path(tables_dir, "05_passing_gene_sets.gmt"))
write_tsv(signature_summary, file.path(tables_dir, "06_signature_summary.tsv"))
write_tsv_gz(optimization_trace, file.path(tables_dir, "07_optimization_trace.tsv.gz"))
write_tsv(comparison, file.path(tables_dir, "08_archived_equivalence.tsv"))

report <- data.frame(metric = c("attempted_directions", "constructed_candidates",
  "construction_failures", "selected_gene_sets", "distinct_selected_genes",
  "archived_membership_exact", "archived_selection_exact",
  "maximum_mean_gap_absolute_error", "scworkbench_used"),
  value = c(nrow(direction_inventory), length(candidate_sets),
    nrow(direction_inventory) - length(candidate_sets), length(passing_ids),
    length(unique(passing_membership$feature_id)), TRUE, TRUE,
    format(max(comparison$mean_gap_absolute_error), scientific = TRUE), FALSE),
  stringsAsFactors = FALSE)
write_tsv(report, file.path(stage_dir, "03_reports", "01_reconstruction_summary.tsv"))
capture.output(sessionInfo(), file = file.path(stage_dir, "04_logs", "01_session_info.txt"))

input_manifest <- data.frame(input_id = names(input_paths),
  path = unname(input_paths), sha256 = vapply(input_paths, sha256, character(1)),
  stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(stage_dir, "06_manifests", "01_input_manifest.tsv"))
output_files <- list.files(stage_dir, recursive = TRUE, full.names = TRUE)
output_files <- output_files[basename(output_files) != "02_output_manifest.tsv"]
output_manifest <- data.frame(
  path = substring(output_files, nchar(stage_dir) + 2L),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, sha256, character(1)), stringsAsFactors = FALSE)
write_tsv(output_manifest[order(output_manifest$path), ],
  file.path(stage_dir, "06_manifests", "02_output_manifest.tsv"))

if (!file.rename(stage_dir, output_dir)) stop("Failed to finalize output directory")
log_event("Standalone reconstruction complete:", output_dir)
