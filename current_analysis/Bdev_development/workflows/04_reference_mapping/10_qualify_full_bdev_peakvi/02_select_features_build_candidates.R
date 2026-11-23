#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 02_select_features_build_candidates.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(Matrix)
  library(yaml)
  library(jsonlite)
  library(digest)
})

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_feature_selection.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), ..., "\n")
  flush.console()
}
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

source_manifest <- file.path(cfg$inputs$source_bundle$path, "manifest.json")
stopifnot(identical(sha256(source_manifest),
    cfg$inputs$source_bundle$manifest_sha256),
  identical(sha256(cfg$inputs$authenticated_l4$pseudotime_path),
    cfg$inputs$authenticated_l4$pseudotime_sha256),
  identical(sha256(cfg$inputs$authenticated_l4$weights_path),
    cfg$inputs$authenticated_l4$weights_sha256),
  identical(cfg$inputs$authenticated_l4$orientation, "native"))

source <- wb_open_bundle(cfg$inputs$source_bundle$path, validate = FALSE,
  read_only = TRUE)
on.exit(try(wb_close_bundle(source), silent = TRUE), add = TRUE)
cells <- as.character(wb_cells(source, "all"))
features <- as.character(wb_feature_ids(source, cfg$inputs$count_matrix))
stopifnot(length(cells) == cfg$inputs$source_bundle$expected_cells,
  length(features) == cfg$inputs$expected_features,
  !anyDuplicated(cells), !anyDuplicated(features))
metadata <- wb_meta(source, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$cluster_field,
  cfg$inputs$sample_field, cfg$inputs$depth_field))
stopifnot(identical(as.character(metadata$cell_id), cells))
cluster <- as.character(metadata[[cfg$inputs$cluster_field]])
sample_id <- as.character(metadata[[cfg$inputs$sample_field]])
depth <- as.numeric(metadata[[cfg$inputs$depth_field]])
stopifnot(!anyNA(cluster), !anyNA(sample_id), all(is.finite(depth)), all(depth > 0))
cluster_levels <- sort(unique(cluster), method = "radix")
cluster_index <- match(cluster, cluster_levels)
cluster_sizes <- tabulate(cluster_index, nbins = length(cluster_levels))

pt_family <- readRDS(cfg$inputs$authenticated_l4$pseudotime_path)
pt_selected <- as.character(pt_family$key) == cfg$inputs$authenticated_l4$key
pt <- stats::setNames(as.numeric(pt_family$value[pt_selected]),
  as.character(pt_family$cell_id[pt_selected]))
weight_family <- readRDS(cfg$inputs$authenticated_l4$weights_path)
weight_selected <- as.character(weight_family$key) == cfg$inputs$authenticated_l4$key
l4_weight <- stats::setNames(as.numeric(weight_family$value[weight_selected]),
  as.character(weight_family$cell_id[weight_selected]))
stopifnot(!anyDuplicated(names(l4_weight)), all(names(l4_weight) %in% cells))
metadata$l4_weight <- 0
metadata$l4_weight[match(names(l4_weight), cells)] <- unname(l4_weight)
metadata$native_l4 <- unname(pt[cells])
metadata$l4_weight_axis_status <- ifelse(cells %in% names(l4_weight),
  "authenticated_weight", "excluded_from_fitted_lineage_family")

quantiles <- unique(as.numeric(stats::quantile(log1p(depth),
  probs = seq(0, 1, length.out = cfg$feature_selection$depth_bins + 1L),
  names = FALSE, type = 7)))
stopifnot(length(quantiles) >= 3L)
depth_bin <- as.integer(cut(log1p(depth), breaks = quantiles,
  include.lowest = TRUE, labels = FALSE))
n_depth <- max(depth_bin)

split_half <- integer(length(cells))
qualification_partition <- rep("reference", length(cells))
for (k in seq_along(cluster_levels)) {
  idx <- which(cluster_index == k)
  idx <- idx[order(idx, method = "radix")]
  split_half[idx] <- 1L + ((seq_along(idx) - 1L) %% 2L)
  qualification_partition[idx[seq_along(idx) %% 5L == 0L]] <- "held_out"
}
reference_cells <- cells[qualification_partition == "reference"]
held_out_cells <- cells[qualification_partition == "held_out"]
stopifnot(length(reference_cells) + length(held_out_cells) == length(cells),
  !length(intersect(reference_cells, held_out_cells)))

make_design <- function(which_half = NULL) {
  keep <- if (is.null(which_half)) rep(TRUE, length(cells)) else split_half == which_half
  j <- (depth_bin - 1L) * length(cluster_levels) + cluster_index
  Matrix::sparseMatrix(i = which(keep), j = j[keep], x = 1,
    dims = c(length(cells), n_depth * length(cluster_levels)))
}
designs <- list(full = make_design(), half1 = make_design(1L), half2 = make_design(2L))
design_totals <- lapply(designs, Matrix::colSums)

xlogp <- function(x, p) {
  out <- numeric(length(x)); use <- x > 0
  out[use] <- x[use] * log(p[use]); out
}
binom_ll <- function(y, n) {
  p <- y / n
  matrix(xlogp(as.vector(y), as.vector(p)) +
      xlogp(as.vector(n - y), as.vector(1 - p)), nrow = nrow(y))
}
collapse_depth <- function(group_counts) {
  out <- matrix(0, nrow(group_counts), n_depth)
  for (b in seq_len(n_depth)) {
    cols <- ((b - 1L) * length(cluster_levels) + 1L):(b * length(cluster_levels))
    out[, b] <- rowSums(group_counts[, cols, drop = FALSE])
  }
  out
}
score_chunk <- function(binary) {
  grouped <- lapply(designs, function(d) as.matrix(binary %*% d))
  overall <- vector("list", length(grouped)); names(overall) <- names(grouped)
  cluster_scores <- vector("list", length(grouped)); names(cluster_scores) <- names(grouped)
  for (route in names(grouped)) {
    g <- grouped[[route]]
    n_group <- matrix(design_totals[[route]], nrow(g), ncol(g), byrow = TRUE)
    y_depth <- collapse_depth(g)
    n_depth_route <- collapse_depth(n_group)
    ll_null <- rowSums(binom_ll(y_depth, n_depth_route))
    overall[[route]] <- -2 * ll_null
    scores <- matrix(0, nrow(g), length(cluster_levels))
    for (k in seq_along(cluster_levels)) {
      cols <- (seq_len(n_depth) - 1L) * length(cluster_levels) + k
      y_in <- g[, cols, drop = FALSE]
      n_in <- n_group[, cols, drop = FALSE]
      y_out <- y_depth - y_in
      n_out <- n_depth_route - n_in
      scores[, k] <- pmax(0, 2 * (rowSums(binom_ll(y_in, n_in)) +
        rowSums(binom_ll(y_out, n_out)) - ll_null))
    }
    cluster_scores[[route]] <- scores
  }
  n_half <- c(sum(split_half == 1L), sum(split_half == 2L))
  stable_cluster <- sqrt(pmax(cluster_scores$half1 / n_half[[1L]], 0) *
    pmax(cluster_scores$half2 / n_half[[2L]], 0))
  stable_overall <- sqrt(pmax(overall[["half1"]] / n_half[[1L]], 0) *
    pmax(overall[["half2"]] / n_half[[2L]], 0))
  list(prevalence = Matrix::rowSums(binary), cluster = stable_cluster,
    overall = stable_overall)
}

log_event("Scoring", length(features), "Bdev peaks in chunks of",
  cfg$feature_selection$feature_chunk_size)
cluster_score_matrix <- matrix(NA_real_, length(features), length(cluster_levels),
  dimnames = list(features, cluster_levels))
overall_score <- prevalence <- numeric(length(features))
chunk_size <- as.integer(cfg$feature_selection$feature_chunk_size)
for (start in seq.int(1L, length(features), by = chunk_size)) {
  idx <- start:min(length(features), start + chunk_size - 1L)
  values <- wb_slice(source, cfg$inputs$count_matrix, features = features[idx],
    cells = cells)
  stopifnot(inherits(values, "sparseMatrix"), identical(rownames(values), features[idx]))
  values@x[] <- 1
  scored <- score_chunk(values)
  prevalence[idx] <- scored$prevalence
  overall_score[idx] <- scored$overall
  cluster_score_matrix[idx, ] <- scored$cluster
  log_event("Scored features", min(idx), "through", max(idx))
  rm(values, scored); invisible(gc())
}

eligible <- prevalence >= cfg$feature_selection$minimum_accessible_cells &
  prevalence <= cfg$feature_selection$maximum_accessible_fraction * length(cells)
cluster_target_index <- max.col(cluster_score_matrix, ties.method = "first")
cluster_target <- cluster_levels[cluster_target_index]
cluster_score <- cluster_score_matrix[cbind(seq_along(features), cluster_target_index)]

allocate_quotas <- function(total) {
  weights <- sqrt(cluster_sizes)
  raw <- total * weights / sum(weights)
  quota <- floor(raw)
  remainder <- total - sum(quota)
  if (remainder > 0L) {
    add <- order(raw - quota, decreasing = TRUE, method = "radix")[seq_len(remainder)]
    quota[add] <- quota[add] + 1L
  }
  quota
}
select_cluster <- function(total) {
  quota <- allocate_quotas(total)
  selected <- integer()
  for (k in seq_along(cluster_levels)) {
    candidates <- which(eligible & cluster_target_index == k)
    candidates <- candidates[order(-cluster_score[candidates], candidates,
      method = "radix")]
    selected <- c(selected, head(candidates, quota[[k]]))
  }
  selected <- unique(selected)
  if (length(selected) < total) {
    fill <- which(eligible & !seq_along(features) %in% selected)
    fill <- fill[order(-cluster_score[fill], fill, method = "radix")]
    selected <- c(selected, head(fill, total - length(selected)))
  }
  stopifnot(length(selected) == total)
  list(index = selected, quota = quota)
}

cluster10 <- select_cluster(cfg$feature_selection$cluster_features)
cluster20 <- select_cluster(cfg$feature_selection$target_features)
overall_order <- which(eligible)
overall_order <- overall_order[order(-overall_score[overall_order], overall_order,
  method = "radix")]
overall20 <- head(overall_order, cfg$feature_selection$target_features)
overall10_fill <- overall_order[!overall_order %in% cluster10$index]
overall10_fill <- head(overall10_fill, cfg$feature_selection$overall_features)
hybrid20 <- c(cluster10$index, overall10_fill)
stopifnot(length(unique(hybrid20)) == cfg$feature_selection$target_features,
  length(cluster20$index) == cfg$feature_selection$target_features,
  length(overall20) == cfg$feature_selection$target_features)

cluster_rank <- rank(-cluster_score, ties.method = "first")
overall_rank <- rank(-overall_score, ties.method = "first")
selection_route <- rep("not_selected", length(features))
selection_route[cluster10$index] <- "cluster_10000"
selection_route[overall10_fill] <- "overall_10000"
score_table <- data.frame(feature_id = features, source_order = seq_along(features),
  accessible_cells = prevalence, accessible_fraction = prevalence / length(cells),
  eligible = eligible, cluster_target = cluster_target,
  cluster_score = cluster_score, cluster_rank = cluster_rank,
  overall_score = overall_score, overall_rank = overall_rank,
  hybrid_selection_route = selection_route,
  in_hybrid20k = seq_along(features) %in% hybrid20,
  in_cluster20k = seq_along(features) %in% cluster20$index,
  in_overall20k = seq_along(features) %in% overall20,
  stringsAsFactors = FALSE)
score_path <- gzfile(file.path(dirs[["tables"]], "01_peak_scores.tsv.gz"), "wt")
write.table(score_table, score_path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA"); close(score_path)

selected_table <- score_table[score_table$in_hybrid20k, , drop = FALSE]
selected_table <- selected_table[order(selected_table$source_order), , drop = FALSE]
write_tsv(selected_table, file.path(dirs[["tables"]], "02_hybrid20k_peak_manifest.tsv"))
quota_table <- data.frame(cluster = cluster_levels, cells = cluster_sizes,
  cluster10k_quota = cluster10$quota,
  cluster20k_quota = cluster20$quota, stringsAsFactors = FALSE)
write_tsv(quota_table, file.path(dirs[["tables"]], "03_cluster_quotas.tsv"))

metadata$partition <- qualification_partition
metadata$split_half <- split_half
metadata$depth_bin <- depth_bin
write_tsv(metadata, file.path(dirs[["tables"]], "04_bdev_cell_design.tsv"))

pool_size <- as.integer(cfg$qualification$pseudobulk_pool_size)
pools <- list(); pool_cursor <- 0L
for (k in seq_along(cluster_levels)) {
  ids <- held_out_cells[cluster[match(held_out_cells, cells)] == cluster_levels[[k]]]
  n_pool <- min(length(ids) %/% pool_size,
    as.integer(cfg$qualification$pseudobulks_per_cluster))
  if (!n_pool) next
  for (p in seq_len(n_pool)) {
    pool_cursor <- pool_cursor + 1L
    members <- ids[((p - 1L) * pool_size + 1L):(p * pool_size)]
    pools[[pool_cursor]] <- list(id = sprintf("bdev_pseudobulk__cluster_%s__%03d",
      cluster_levels[[k]], p), cluster = cluster_levels[[k]], members = members)
  }
}
pseudobulk_ids <- vapply(pools, `[[`, "", "id")
query_cells <- c(held_out_cells, pseudobulk_ids)
all_candidate_cells <- c(reference_cells, query_cells)

if (dir.exists(cfg$storage$candidate_bundle))
  stop("Candidate bundle already exists; preserve it and use a new run identity")
bundle <- wb_init_bundle(cfg$storage$candidate_bundle, overwrite = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
wb_define_cells(bundle, all_candidate_cells)

synthetic_meta <- data.frame(cell_id = pseudobulk_ids,
  stringsAsFactors = FALSE, check.names = FALSE)
synthetic_meta[[cfg$inputs$cluster_field]] <- vapply(pools, `[[`, "", "cluster")
synthetic_meta[[cfg$inputs$sample_field]] <- "Bdev_heldout_pool"
synthetic_meta[[cfg$inputs$depth_field]] <- vapply(pools, function(x)
  sum(metadata[[cfg$inputs$depth_field]][match(x$members, metadata$cell_id)]), numeric(1))
synthetic_meta$l4_weight <- vapply(pools, function(x)
  mean(metadata$l4_weight[match(x$members, metadata$cell_id)]), numeric(1))
synthetic_meta$native_l4 <- vapply(pools, function(x) {
  value <- metadata$native_l4[match(x$members, metadata$cell_id)]
  if (any(is.finite(value))) stats::median(value[is.finite(value)]) else NA_real_
}, numeric(1))
synthetic_meta$l4_weight_axis_status <- "derived_Bdev_pseudobulk_truth"
synthetic_meta$partition <- "pseudobulk"
synthetic_meta$split_half <- NA_integer_
synthetic_meta$depth_bin <- NA_integer_
candidate_meta <- rbind(metadata[match(c(reference_cells, held_out_cells), metadata$cell_id), ],
  synthetic_meta)
candidate_meta <- candidate_meta[match(all_candidate_cells, candidate_meta$cell_id), ]
stopifnot(identical(as.character(candidate_meta$cell_id), all_candidate_cells))
wb_attach_metadata(bundle, candidate_meta, cell_ids = all_candidate_cells,
  source = "Bdev-only full-reference feature qualification design",
  created_by = "scworkbench_Bdev development")

candidate_axes <- list(
  hybrid20k = sort(hybrid20),
  cluster20k = sort(cluster20$index),
  overall20k = sort(overall20),
  full191433 = seq_along(features))
candidate_inventory <- list()
for (candidate in names(candidate_axes)) {
  feature_index <- candidate_axes[[candidate]]
  axis <- features[feature_index]
  log_event("Materializing candidate", candidate, "with", length(axis), "features")
  values <- wb_slice(source, cfg$inputs$count_matrix, features = axis, cells = cells)
  reference_matrix <- values[, reference_cells, drop = FALSE]
  held_out_matrix <- values[, held_out_cells, drop = FALSE]
  pseudo_columns <- lapply(pools, function(pool) {
    aggregate <- Matrix::rowSums(values[, pool$members, drop = FALSE])
    nz <- which(aggregate != 0)
    Matrix::sparseMatrix(i = nz, j = rep.int(1L, length(nz)), x = aggregate[nz],
      dims = c(length(axis), 1L), dimnames = list(axis, pool$id))
  })
  pseudobulk_matrix <- if (length(pseudo_columns)) {
    out <- do.call(cbind, pseudo_columns)
    rownames(out) <- axis; colnames(out) <- pseudobulk_ids; out
  } else Matrix::Matrix(0, nrow = length(axis), ncol = 0L, sparse = TRUE,
    dimnames = list(axis, character()))
  query_matrix <- cbind(held_out_matrix, pseudobulk_matrix)
  stopifnot(identical(colnames(query_matrix), query_cells))
  wb_attach_matrix(bundle, paste0(candidate, "_reference_counts"), reference_matrix,
    layer = "counts", feature_type = "peak", cell_ids = reference_cells,
    description = paste(candidate, "Bdev-only PeakVI qualification reference"),
    source = "active_03c_Bdev_ATAC_counts_Bdev_only_feature_selection",
    created_by = "scworkbench_Bdev development")
  wb_attach_matrix(bundle, paste0(candidate, "_query_counts"), query_matrix,
    layer = "counts", feature_type = "peak", cell_ids = query_cells,
    description = paste(candidate, "held-out Bdev cells and Bdev pseudobulks"),
    source = "active_03c_Bdev_ATAC_counts_Bdev_only_feature_selection",
    created_by = "scworkbench_Bdev development")
  candidate_inventory[[candidate]] <- list(features = length(axis),
    feature_sha256 = digest::digest(axis, algo = "sha256", serialize = TRUE),
    reference_cells = length(reference_cells), query_cells = length(query_cells))
  rm(values, reference_matrix, held_out_matrix, pseudobulk_matrix, query_matrix,
    pseudo_columns); invisible(gc())
}
wb_close_bundle(bundle); bundle <- NULL
wb_close_bundle(source); source <- NULL

selection_summary <- list(
  schema_version = "bdev.full_peakvi_feature_selection.v1",
  status = "complete",
  run_identity = cfg$run_identity,
  source_cells = length(cells), source_features = length(features),
  cluster_count = length(cluster_levels), reference_cells = length(reference_cells),
  held_out_cells = length(held_out_cells), pseudobulks = length(pseudobulk_ids),
  hybrid_cluster_features = sum(selection_route == "cluster_10000"),
  hybrid_overall_features = sum(selection_route == "overall_10000"),
  external_query_use = FALSE, real_immgen_use = FALSE,
  candidates = candidate_inventory,
  source_manifest_unchanged = identical(sha256(source_manifest),
    cfg$inputs$source_bundle$manifest_sha256))
stopifnot(selection_summary$hybrid_cluster_features == 10000L,
  selection_summary$hybrid_overall_features == 10000L,
  selection_summary$source_manifest_unchanged)
write_json(selection_summary,
  file.path(dirs[["reports"]], "01_feature_selection_summary.json"))
write_json(list(schema_version = "bdev.full_peakvi_feature_execution.v1",
  run_identity = cfg$run_identity, config_path = cfg_path,
  config_sha256 = sha256(cfg_path), source_manifest_sha256 = sha256(source_manifest),
  external_query_use = FALSE, real_immgen_use = FALSE),
  file.path(dirs[["manifests"]], "02_feature_selection_receipt.json"))
log_event("Feature selection and candidate bundle construction complete")
