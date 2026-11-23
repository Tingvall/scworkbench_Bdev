#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 01_run_sparsity_adapter.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(digest)
  library(Matrix)
  library(ggplot2)
  library(ggrepel)
})

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
promoted_ref <- function(promoted, output, kind) {
  take <- promoted$promoted$status == "ok" & promoted$promoted$output == output
  stopifnot(sum(take) == 1L)
  wb_ref_artifact(promoted$promoted$artifact_id[take][[1L]], kind)
}
weighted_quantile <- function(value, weight, probabilities) {
  keep <- is.finite(value) & is.finite(weight) & weight > 0
  value <- value[keep]; weight <- weight[keep]
  index <- order(value, seq_along(value), method = "radix")
  value <- value[index]; weight <- weight[index]
  cumulative <- cumsum(weight) / sum(weight)
  vapply(probabilities, function(p) value[[which(cumulative >= p)[[1L]]]],
    numeric(1))
}
space_2d <- function(x) {
  if (is.data.frame(x) && all(c("cell_id", "x", "y") %in% names(x)))
    return(data.frame(cell_id = as.character(x$cell_id), UMAP_1 = x$x,
      UMAP_2 = x$y, stringsAsFactors = FALSE))
  x <- as.matrix(x)
  stopifnot(ncol(x) == 2L, !is.null(rownames(x)))
  data.frame(cell_id = rownames(x), UMAP_1 = x[, 1L], UMAP_2 = x[, 2L],
    stringsAsFactors = FALSE)
}
quantile_row <- function(x, prefix) {
  q <- stats::quantile(x, c(0, .25, .5, .75, 1), names = FALSE, na.rm = TRUE)
  stats::setNames(as.list(as.numeric(q)), paste0(prefix,
    c("_min", "_q25", "_median", "_q75", "_max")))
}

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
.libPaths(c(cfg$runtime$run_library, .libPaths()))
suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbLatent)
  library(scwbRepresentation)
})

root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_sparsity_adapter.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

source_manifest <- file.path(cfg$inputs$source_mapping_bundle$path, "manifest.json")
umap_manifest <- file.path(cfg$inputs$approved_atac_umap$bundle_path, "manifest.json")
stopifnot(
  identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256),
  identical(sha256(source_manifest),
    cfg$inputs$source_mapping_bundle$manifest_sha256),
  identical(sha256(umap_manifest),
    cfg$inputs$approved_atac_umap$manifest_sha256),
  identical(sha256(cfg$inputs$authenticated_l4$path),
    cfg$inputs$authenticated_l4$sha256),
  identical(cfg$inputs$authenticated_l4$orientation, "native"),
  identical(sha256(cfg$inputs$audit_only_comparators$full_depth_fixed_umap$path),
    cfg$inputs$audit_only_comparators$full_depth_fixed_umap$sha256),
  identical(sha256(cfg$inputs$audit_only_comparators$depth_only_design$path),
    cfg$inputs$audit_only_comparators$depth_only_design$sha256),
  as.character(packageVersion("scworkbench")) == cfg$runtime$scworkbench$version,
  as.character(packageVersion("scwbLatent")) == cfg$runtime$scwbLatent$version,
  as.character(packageVersion("scwbRepresentation")) ==
    cfg$runtime$scwbRepresentation$version)
for (package in c("scwbLatent", "scwbRepresentation")) {
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  stopifnot(identical(sha256(manifest), cfg$runtime[[package]]$manifest_sha256))
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
}
stopifnot(all(c("scwblatent::peakvi_fit", "scwblatent::peakvi_project",
  "scwblatent::model_bound_latent_union", "scwbrepresentation::knn_neighbors") %in%
  wb_list_ops()))

source_manifest_before <- sha256(source_manifest)
umap_manifest_before <- sha256(umap_manifest)
source <- wb_open_bundle(cfg$inputs$source_mapping_bundle$path, validate = FALSE,
  read_only = TRUE)
reference_space_existing <- wb_space(source,
  cfg$inputs$source_mapping_bundle$reference_space)
full_query_latent <- as.matrix(wb_space(source,
  cfg$inputs$source_mapping_bundle$full_query_space))
reference_cells <- rownames(reference_space_existing)
query_cells <- rownames(full_query_latent)
feature_ids <- wb_feature_ids(source,
  cfg$inputs$source_mapping_bundle$reference_matrix)
reference_counts <- wb_slice(source,
  cfg$inputs$source_mapping_bundle$reference_matrix,
  features = feature_ids, cells = reference_cells)
bulk_counts <- wb_slice(source, cfg$inputs$source_mapping_bundle$bulk_matrix,
  features = feature_ids, cells = query_cells)
wb_close_bundle(source)
stopifnot(length(feature_ids) == cfg$inputs$source_mapping_bundle$expected_features,
  length(reference_cells) == cfg$inputs$source_mapping_bundle$expected_reference_cells,
  length(query_cells) == cfg$inputs$source_mapping_bundle$expected_libraries,
  identical(colnames(reference_counts), reference_cells),
  identical(colnames(bulk_counts), query_cells))

# The only Bdev information used to construct queries is global matrix depth and
# sparsity. No metadata are loaded here.
reference_depth <- as.numeric(Matrix::colSums(reference_counts))
reference_nonzero <- as.integer(Matrix::colSums(reference_counts > 0))
target_depth <- as.integer(round(stats::median(reference_depth)))
window_fraction <- as.numeric(cfg$adapter$sparsity_depth_window_fraction)
in_sparsity_pool <- abs(reference_depth - target_depth) <=
  target_depth * window_fraction
sparsity_pool <- reference_nonzero[in_sparsity_pool]
stopifnot(target_depth == as.integer(cfg$adapter$expected_target_fragment_count),
  sum(in_sparsity_pool) == as.integer(cfg$adapter$expected_sparsity_pool_cells),
  all(sparsity_pool > 0L), all(sparsity_pool < target_depth))

structure_summary <- data.frame(
  statistic = c("reference_cells", "reference_features", "target_fragment_count",
    "sparsity_depth_window_fraction", "sparsity_pool_cells",
    "reference_nonzero_q25", "reference_nonzero_median", "reference_nonzero_q75",
    "pool_nonzero_q25", "pool_nonzero_median", "pool_nonzero_q75"),
  value = c(length(reference_cells), length(feature_ids), target_depth,
    window_fraction, sum(in_sparsity_pool),
    stats::quantile(reference_nonzero, .25), stats::median(reference_nonzero),
    stats::quantile(reference_nonzero, .75),
    stats::quantile(sparsity_pool, .25), stats::median(sparsity_pool),
    stats::quantile(sparsity_pool, .75)),
  construction_role = c(rep("global_reference_structure", 11L)),
  stringsAsFactors = FALSE)
write_tsv(structure_summary,
  file.path(dirs[["tables"]], "01_global_reference_structure.tsv"))
sparsity_frequency <- as.data.frame(table(target_nonzero_peaks = sparsity_pool),
  stringsAsFactors = FALSE)
sparsity_frequency$target_nonzero_peaks <- as.integer(
  as.character(sparsity_frequency$target_nonzero_peaks))
sparsity_frequency$reference_cell_count <- as.integer(sparsity_frequency$Freq)
sparsity_frequency$Freq <- NULL
write_tsv(sparsity_frequency,
  file.path(dirs[["tables"]], "02_sparsity_pool_distribution.tsv"))

n_replicates <- as.integer(cfg$adapter$pseudocells_per_library)
rank_power <- as.numeric(cfg$adapter$rank_power)
oversampling <- as.numeric(cfg$adapter$draw_oversampling_factor)
set.seed(as.integer(cfg$adapter$seed))
n_pseudo <- length(query_cells) * n_replicates
i_parts <- vector("list", n_pseudo)
x_parts <- vector("list", n_pseudo)
design_rows <- vector("list", n_pseudo)
rank_audit_rows <- vector("list", length(query_cells))
cursor <- 0L

for (library_index in seq_along(query_cells)) {
  gsm <- query_cells[[library_index]]
  values <- bulk_counts[, gsm, drop = FALSE]
  nz <- which(values[, 1L] > 0)
  observed <- as.numeric(values[nz, 1L])
  count_rank <- rank(observed, ties.method = "average")
  rank_percentile <- count_rank / length(count_rank)
  selection_weight <- count_rank ^ rank_power
  target_k <- sample(sparsity_pool, n_replicates, replace = TRUE)
  stopifnot(max(target_k) <= length(nz), max(target_k) < target_depth)
  draw_budget <- as.integer(ceiling(max(target_k) * oversampling))
  draw_stream <- matrix(sample.int(length(nz), n_replicates * draw_budget,
    replace = TRUE, prob = selection_weight), nrow = draw_budget,
    ncol = n_replicates)
  selection_count <- integer(length(nz))
  fragment_total <- numeric(length(nz))

  for (replicate in seq_len(n_replicates)) {
    k <- as.integer(target_k[[replicate]])
    selected <- unique(draw_stream[, replicate])
    while (length(selected) < k) {
      need <- k - length(selected)
      extra <- sample.int(length(nz), max(128L, as.integer(ceiling(need * 1.5))),
        replace = TRUE, prob = selection_weight)
      selected <- unique(c(selected, extra))
    }
    selected <- selected[seq_len(k)]
    remaining <- target_depth - k
    extra_counts <- as.numeric(stats::rmultinom(1L, size = remaining,
      prob = count_rank[selected]))
    adapted_counts <- 1 + extra_counts
    cursor <- cursor + 1L
    pseudo_id <- sprintf("%s__sparse_r%03d", gsm, replicate)
    i_parts[[cursor]] <- nz[selected]
    x_parts[[cursor]] <- adapted_counts
    selection_count <- selection_count + tabulate(selected, nbins = length(nz))
    fragment_total[selected] <- fragment_total[selected] + adapted_counts
    design_rows[[cursor]] <- data.frame(query_cell_id = pseudo_id, gsm = gsm,
      replicate = replicate, target_matrix_depth = target_depth,
      observed_matrix_depth = sum(adapted_counts), target_nonzero_peaks = k,
      observed_nonzero_peaks = length(selected),
      source_matrix_depth = as.numeric(Matrix::colSums(values)),
      source_nonzero_peaks = length(nz),
      mean_selected_rank_percentile = mean(rank_percentile[selected]),
      fraction_selected_from_bulk_top_decile = mean(
        rank_percentile[selected] >= .9), stringsAsFactors = FALSE)
  }
  rank_audit_rows[[library_index]] <- data.frame(gsm = gsm,
    source_nonzero_peaks = length(nz),
    accessibility_frequency_vs_bulk_rank_spearman = suppressWarnings(
      cor(selection_count / n_replicates, count_rank, method = "spearman")),
    fragment_total_vs_bulk_rank_spearman = suppressWarnings(
      cor(fragment_total, count_rank, method = "spearman")),
    top_decile_mean_accessibility_frequency = mean(
      selection_count[rank_percentile >= .9] / n_replicates),
    bottom_decile_mean_accessibility_frequency = mean(
      selection_count[rank_percentile <= .1] / n_replicates),
    stringsAsFactors = FALSE)
}
stopifnot(cursor == n_pseudo)
design <- do.call(rbind, design_rows); rownames(design) <- NULL
rank_audit <- do.call(rbind, rank_audit_rows); rownames(rank_audit) <- NULL
column_lengths <- lengths(i_parts)
pseudo_counts <- Matrix::sparseMatrix(
  i = unlist(i_parts, use.names = FALSE),
  j = rep.int(seq_len(n_pseudo), column_lengths),
  x = unlist(x_parts, use.names = FALSE),
  dims = c(length(feature_ids), n_pseudo),
  dimnames = list(feature_ids, design$query_cell_id))
stopifnot(ncol(pseudo_counts) == cfg$acceptance$expected_pseudocells,
  identical(colnames(pseudo_counts), design$query_cell_id),
  all(Matrix::colSums(pseudo_counts) == target_depth),
  all(Matrix::colSums(pseudo_counts > 0) == design$target_nonzero_peaks),
  all(design$observed_matrix_depth == target_depth),
  all(design$observed_nonzero_peaks == design$target_nonzero_peaks))
write_tsv(design, file.path(dirs[["tables"]], "03_pseudocell_design.tsv"))
write_tsv(rank_audit,
  file.path(dirs[["tables"]], "04_bulk_rank_preservation_audit.tsv"))
rm(i_parts, x_parts, design_rows, draw_stream); invisible(gc())

if (dir.exists(cfg$storage$diagnostic_bundle))
  stop("Refusing to overwrite existing diagnostic bundle: ",
    cfg$storage$diagnostic_bundle)
bundle <- wb_init_bundle(cfg$storage$diagnostic_bundle, overwrite = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
all_cells <- c(reference_cells, design$query_cell_id)
wb_define_cells(bundle, all_cells)
wb_attach_matrix(bundle, "Bdev_ATAC_reference_counts", reference_counts,
  layer = "counts", feature_type = "peak", cell_ids = reference_cells,
  description = "Qualified Bdev-only PeakVI reference counts",
  source = cfg$inputs$source_mapping_bundle$path,
  created_by = "scworkbench_Bdev label-blind sparsity adapter diagnostic")
wb_attach_matrix(bundle, "ImmGen_ATAC_label_blind_sparse_pseudocells",
  pseudo_counts, layer = "counts", feature_type = "peak",
  cell_ids = design$query_cell_id,
  description = "Label-blind rank-weighted ImmGen pseudo-cells with global reference-matched depth and sparsity",
  source = cfg$inputs$source_mapping_bundle$path,
  created_by = "scworkbench_Bdev label-blind sparsity adapter diagnostic")
rm(reference_counts, bulk_counts, pseudo_counts); invisible(gc())

# The current runtime has no authenticated cross-bundle model-import primitive.
# Reconstruct the already qualified deterministic model recipe, then require the
# exact qualified reference-latent checksum before projecting any query.
fit_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_peakvi_fit"))
fit <- wb_session_run_op(fit_session, bundle, "scwblatent::peakvi_fit",
  inputs = list(matrix = "Bdev_ATAC_reference_counts",
    n_latent = as.integer(cfg$peakvi$n_latent),
    n_layers = as.integer(cfg$peakvi$n_layers),
    n_hidden = as.integer(cfg$peakvi$n_hidden),
    dropout = as.numeric(cfg$peakvi$dropout),
    max_epochs = as.integer(cfg$peakvi$max_epochs),
    lr = as.numeric(cfg$peakvi$learning_rate),
    early_stopping = isTRUE(cfg$peakvi$early_stopping),
    seed = as.integer(cfg$peakvi$seed),
    chunk_size = as.integer(cfg$peakvi$chunk_size), batch_key = NULL,
    use_gpu = FALSE, cell_set = "all", min_cells = 1L),
  notes = "Deterministic reconstruction of the qualified Bdev-only PeakVI reference")
reference_space_id <- as.character(fit$run$outputs$space$space_id %||%
  fit$run$outputs$space$allocated_id)[[1L]]
fit_summary <- fit$run$outputs$fit_summary$table
stopifnot(fit_summary$latent_sha256[[1L]] ==
  cfg$peakvi$expected_reference_latent_sha256)
fit$session <- wb_session_resume(bundle, fit$session$durable$workspace)
fit_promoted <- wb_session_promote_current(bundle, fit$session,
  outputs = c("model", "fit_summary"), name_prefix = "sparsity_adapter_peakvi_reference",
  description = "Checksum-authenticated deterministic Bdev-only PeakVI reference",
  created_by = "scworkbench_Bdev label-blind sparsity adapter diagnostic",
  fail_fast = TRUE)
model_ref <- promoted_ref(fit_promoted, "model", "model")
fit_promoted$session <- wb_session_discard(fit_promoted$session)

project_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_peakvi_project"))
project <- wb_session_run_op(project_session, bundle, "scwblatent::peakvi_project",
  inputs = list(model = model_ref,
    query_matrix = "ImmGen_ATAC_label_blind_sparse_pseudocells",
    space_id_out = "sparsity_adapter_peakvi_queries",
    chunk_size = as.integer(cfg$peakvi$chunk_size), batch_key = NULL,
    use_gpu = FALSE, pickle_execution_trust = cfg$peakvi$pickle_execution_trust,
    cell_set = "all", min_cells = 1L),
  notes = "Projection of label-blind depth-and-sparsity-matched ImmGen pseudo-cells without query retraining")
stopifnot(project$run$outputs$projection_summary$table$transformation[[1L]] ==
  "transfer_without_refit")
project$session <- wb_session_resume(bundle, project$session$durable$workspace)
project_promoted <- wb_session_promote_current(bundle, project$session,
  outputs = "projection_summary", name_prefix = "sparsity_adapter_projection",
  description = "Label-blind sparsity-adapted ImmGen PeakVI projection",
  created_by = "scworkbench_Bdev label-blind sparsity adapter diagnostic",
  fail_fast = TRUE)
projection_ref <- promoted_ref(project_promoted, "projection_summary", "obs_table")
project_promoted$session <- wb_session_discard(project_promoted$session)

union_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_latent_union"))
union <- wb_session_run_op(union_session, bundle,
  "scwblatent::model_bound_latent_union", inputs = list(
    reference_space = reference_space_id,
    query_space = "sparsity_adapter_peakvi_queries", model = model_ref,
    projection_summary = projection_ref,
    space_id_out = "sparsity_adapter_reference_and_queries",
    chunk_size = as.integer(cfg$mapping$query_chunk_size)),
  notes = "Authenticated Bdev/sparsity-adapted-query latent union")
union$session <- wb_session_resume(bundle, union$session$durable$workspace)
union_promoted <- wb_session_promote_current(bundle, union$session,
  outputs = "union_summary", name_prefix = "sparsity_adapter_union",
  description = "Authenticated model-bound sparsity-adapter latent union",
  created_by = "scworkbench_Bdev label-blind sparsity adapter diagnostic",
  fail_fast = TRUE)
union_promoted$session <- wb_session_discard(union_promoted$session)

knn_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_knn"))
knn <- wb_session_run_op(knn_session, bundle,
  "scwbrepresentation::knn_neighbors", inputs = list(
    space = "sparsity_adapter_reference_and_queries",
    query_cells = design$query_cell_id, candidate_cells = reference_cells,
    k = as.integer(cfg$mapping$neighbor_k), include_query = FALSE,
    backend = cfg$mapping$backend,
    space_chunk_size = as.integer(cfg$mapping$space_chunk_size),
    query_chunk_size = as.integer(cfg$mapping$query_chunk_size)),
  notes = "Reference neighbors for label-blind sparsity-adapted queries")
neighbors_all <- knn$run$outputs$neighbors$table
stopifnot(nrow(neighbors_all) == nrow(design) * cfg$mapping$neighbor_k)
knn$session <- wb_session_resume(bundle, knn$session$durable$workspace)
knn$session <- wb_session_discard(knn$session)

query_latent <- as.matrix(wb_space(bundle, "sparsity_adapter_peakvi_queries"))
query_latent <- query_latent[design$query_cell_id, , drop = FALSE]
latent_table <- data.frame(query_cell_id = rownames(query_latent),
  gsm = design$gsm, replicate = design$replicate, query_latent,
  check.names = FALSE, stringsAsFactors = FALSE)
write_tsv(latent_table,
  file.path(dirs[["tables"]], "05_pseudocell_peakvi_latent.tsv"))

bdev <- wb_open_bundle(cfg$inputs$approved_atac_umap$bundle_path,
  validate = FALSE, read_only = TRUE)
atac_umap <- space_2d(wb_space(bdev, cfg$inputs$approved_atac_umap$space))
wb_close_bundle(bdev)
stopifnot(nrow(atac_umap) == cfg$inputs$approved_atac_umap$expected_cells,
  all(reference_cells %in% atac_umap$cell_id))
reference_umap <- atac_umap[match(reference_cells, atac_umap$cell_id), ]
rownames(reference_umap) <- reference_umap$cell_id

family <- readRDS(cfg$inputs$authenticated_l4$path)
selected_l4 <- as.character(family$key) == cfg$inputs$authenticated_l4$key &
  as.character(family$cell_id) %in% reference_cells &
  is.finite(as.numeric(family$value))
l4 <- stats::setNames(as.numeric(family$value[selected_l4]),
  as.character(family$cell_id[selected_l4]))
stopifnot(length(l4) == length(reference_cells))
trajectory <- merge(data.frame(cell_id = names(l4), native_l4 = unname(l4),
  stringsAsFactors = FALSE), atac_umap, by = "cell_id", sort = FALSE)

neighbors <- neighbors_all[neighbors_all$neighbor_rank <=
  as.integer(cfg$mapping$primary_k), , drop = FALSE]
neighbors$gsm <- design$gsm[match(neighbors$query_cell_id,
  design$query_cell_id)]
neighbors$replicate <- design$replicate[match(neighbors$query_cell_id,
  design$query_cell_id)]
neighbors$kernel_weight <- NA_real_
neighbors$normalized_weight <- NA_real_
position_rows <- vector("list", nrow(design))
for (i in seq_len(nrow(design))) {
  cell <- design$query_cell_id[[i]]
  index <- which(neighbors$query_cell_id == cell)
  index <- index[order(neighbors$neighbor_rank[index])]
  stopifnot(length(index) == as.integer(cfg$mapping$primary_k))
  distance <- neighbors$distance[index]
  kth <- distance[[length(distance)]]
  weight <- if (kth == 0) as.numeric(distance == 0) else
    exp(-(distance / kth)^2)
  normalized <- weight / sum(weight)
  neighbors$kernel_weight[index] <- weight
  neighbors$normalized_weight[index] <- normalized
  xy <- as.matrix(reference_umap[neighbors$cell_id[index],
    c("UMAP_1", "UMAP_2")])
  center <- colSums(xy * normalized)
  spread <- sqrt(sum(normalized * ((xy[, 1L] - center[[1L]])^2 +
    (xy[, 2L] - center[[2L]])^2)))
  l4_q <- weighted_quantile(unname(l4[neighbors$cell_id[index]]), weight,
    c(.25, .5, .75))
  position_rows[[i]] <- data.frame(query_cell_id = cell,
    gsm = design$gsm[[i]], replicate = design$replicate[[i]],
    UMAP_1 = center[[1L]], UMAP_2 = center[[2L]],
    neighbor_spread = spread, nearest_distance = min(distance),
    median_neighbor_distance = stats::median(distance), kth_distance = kth,
    native_l4_q25 = l4_q[[1L]], native_l4_median = l4_q[[2L]],
    native_l4_q75 = l4_q[[3L]], stringsAsFactors = FALSE)
}
positions <- do.call(rbind, position_rows); rownames(positions) <- NULL
write_tsv(neighbors,
  file.path(dirs[["tables"]], "06_pseudocell_neighbor_ledger.tsv"))
write_tsv(positions,
  file.path(dirs[["tables"]], "07_pseudocell_fixed_umap_positions.tsv"))

library_rows <- lapply(query_cells, function(gsm) {
  p <- positions[positions$gsm == gsm, , drop = FALSE]
  z <- query_latent[p$query_cell_id, , drop = FALSE]
  center <- c(stats::median(p$UMAP_1), stats::median(p$UMAP_2))
  spatial_distance <- sqrt((p$UMAP_1 - center[[1L]])^2 +
    (p$UMAP_2 - center[[2L]])^2)
  data.frame(gsm = gsm, n_pseudocells = nrow(p),
    UMAP_1_median = center[[1L]], UMAP_1_q05 = quantile(p$UMAP_1, .05),
    UMAP_1_q95 = quantile(p$UMAP_1, .95),
    UMAP_2_median = center[[2L]], UMAP_2_q05 = quantile(p$UMAP_2, .05),
    UMAP_2_q95 = quantile(p$UMAP_2, .95),
    median_radial_position_spread = median(spatial_distance),
    native_l4_q25 = quantile(p$native_l4_median, .25),
    native_l4_median = median(p$native_l4_median),
    native_l4_q75 = quantile(p$native_l4_median, .75),
    median_neighbor_distance = median(p$median_neighbor_distance),
    median_neighbor_spread = median(p$neighbor_spread),
    latent_3_median = median(z[, 3L]), stringsAsFactors = FALSE)
})
library_summary <- do.call(rbind, library_rows); rownames(library_summary) <- NULL
library_summary <- merge(library_summary, rank_audit, by = "gsm", sort = FALSE)
library_summary <- library_summary[match(query_cells, library_summary$gsm), ]

# Parameter-convergence audit uses only prefixes of the fixed seeded replicates.
prefixes <- unique(pmin(c(25L, 50L, 75L, 100L), n_replicates))
convergence_rows <- list(); convergence_cursor <- 0L
for (gsm in query_cells) {
  p_all <- positions[positions$gsm == gsm, , drop = FALSE]
  p_all <- p_all[order(p_all$replicate), ]
  final_center <- c(median(p_all$UMAP_1), median(p_all$UMAP_2))
  final_l4 <- median(p_all$native_l4_median)
  for (n_used in prefixes) {
    p <- p_all[seq_len(n_used), ]
    center <- c(median(p$UMAP_1), median(p$UMAP_2))
    convergence_cursor <- convergence_cursor + 1L
    convergence_rows[[convergence_cursor]] <- data.frame(gsm = gsm,
      pseudocells_used = n_used, UMAP_1_median = center[[1L]],
      UMAP_2_median = center[[2L]], native_l4_median = median(p$native_l4_median),
      distance_from_100_pseudocell_center = sqrt(sum((center - final_center)^2)),
      absolute_l4_difference_from_100 = abs(median(p$native_l4_median) - final_l4),
      stringsAsFactors = FALSE)
  }
}
convergence <- do.call(rbind, convergence_rows); rownames(convergence) <- NULL
write_tsv(convergence,
  file.path(dirs[["tables"]], "08_pseudocell_number_convergence.tsv"))

# Consensus neighbors summarize the ensemble; no annotation or expected label is
# involved in ranking them.
consensus_rows <- list(); consensus_cursor <- 0L
for (gsm in query_cells) {
  x <- neighbors[neighbors$gsm == gsm, , drop = FALSE]
  selected_count <- table(x$cell_id)
  total_weight <- tapply(x$normalized_weight, x$cell_id, sum)
  cells <- names(total_weight)
  ranked <- data.frame(gsm = gsm, cell_id = cells,
    selection_fraction = as.numeric(selected_count[cells]) / n_replicates,
    mean_normalized_weight = as.numeric(total_weight[cells]) / n_replicates,
    stringsAsFactors = FALSE)
  ranked <- ranked[order(-ranked$mean_normalized_weight,
    -ranked$selection_fraction, ranked$cell_id, method = "radix"), ]
  ranked <- ranked[seq_len(min(as.integer(cfg$mapping$primary_k), nrow(ranked))), ]
  ranked$consensus_rank <- seq_len(nrow(ranked))
  consensus_cursor <- consensus_cursor + 1L
  consensus_rows[[consensus_cursor]] <- ranked
}
consensus <- do.call(rbind, consensus_rows); rownames(consensus) <- NULL
write_tsv(consensus,
  file.path(dirs[["tables"]], "09_library_consensus_neighbors.tsv"))
reuse <- as.data.frame(table(cell_id = consensus$cell_id), stringsAsFactors = FALSE)
reuse$cell_id <- as.character(reuse$cell_id)
names(reuse)[[2L]] <- "libraries_selecting_cell"
reuse <- reuse[order(-reuse$libraries_selecting_cell, reuse$cell_id), ]
write_tsv(reuse,
  file.path(dirs[["tables"]], "10_consensus_neighbor_reuse.tsv"))

# These comparisons are deliberately post-projection and audit-only. They do not
# enter pseudo-cell construction or parameter selection.
full_position <- read.delim(
  cfg$inputs$audit_only_comparators$full_depth_fixed_umap$path,
  check.names = FALSE, quote = "\"", na.strings = "NA")
full_position <- full_position[, c("gsm", "UMAP_1", "UMAP_2",
  "accepted_native_l4")]
names(full_position) <- c("gsm", "full_UMAP_1", "full_UMAP_2", "full_native_l4")
library_summary <- merge(library_summary, full_position, by = "gsm", sort = FALSE)
library_summary <- library_summary[match(query_cells, library_summary$gsm), ]
library_summary$fixed_umap_shift_from_full <- sqrt(
  (library_summary$UMAP_1_median - library_summary$full_UMAP_1)^2 +
  (library_summary$UMAP_2_median - library_summary$full_UMAP_2)^2)
library_summary$native_l4_shift_from_full <- library_summary$native_l4_median -
  library_summary$full_native_l4
write_tsv(library_summary,
  file.path(dirs[["tables"]], "11_library_adapter_summary.tsv"))

latent_rows <- list(); latent_cursor <- 0L
for (gsm in query_cells) {
  z <- query_latent[design$query_cell_id[design$gsm == gsm], , drop = FALSE]
  adapted <- apply(z, 2L, median)
  full <- full_query_latent[gsm, ]
  for (dimension in seq_along(adapted)) {
    latent_cursor <- latent_cursor + 1L
    latent_rows[[latent_cursor]] <- data.frame(gsm = gsm,
      dimension = paste0("latent_", dimension), full_bulk = full[[dimension]],
      adapted_median = adapted[[dimension]],
      adapted_minus_full = adapted[[dimension]] - full[[dimension]],
      stringsAsFactors = FALSE)
  }
}
latent_comparison <- do.call(rbind, latent_rows); rownames(latent_comparison) <- NULL
write_tsv(latent_comparison,
  file.path(dirs[["tables"]], "12_full_vs_adapter_latent.tsv"))

depth_only_design <- read.delim(
  cfg$inputs$audit_only_comparators$depth_only_design$path,
  check.names = FALSE, quote = "\"", na.strings = "NA")
sparsity_plot_data <- rbind(
  data.frame(method = "Reference cells near median depth", nonzero_peaks = sparsity_pool),
  data.frame(method = "Depth-only pseudo-cells", nonzero_peaks = depth_only_design$nonzero_peaks),
  data.frame(method = "Rank + sparsity adapter", nonzero_peaks = design$observed_nonzero_peaks))
sparsity_plot_data$method <- factor(sparsity_plot_data$method,
  levels = c("Reference cells near median depth", "Depth-only pseudo-cells",
    "Rank + sparsity adapter"))
p_sparsity <- ggplot(sparsity_plot_data, aes(nonzero_peaks, fill = method,
  color = method)) +
  geom_density(alpha = .22, linewidth = .8) +
  scale_fill_manual(values = c("#009E73", "#D55E00", "#0072B2")) +
  scale_color_manual(values = c("#009E73", "#D55E00", "#0072B2")) +
  labs(title = "The adapter independently matches single-cell sparsity",
    subtitle = paste0("All adapted pseudo-cells have ", target_depth,
      " fragments; construction used no labels or expected positions"),
    x = "Nonzero peaks per cell or pseudo-cell", y = "Density",
    fill = NULL, color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank())
ggsave(file.path(dirs[["figures"]], "01_sparsity_matching.png"), p_sparsity,
  width = 9, height = 6, dpi = 220, bg = "white")
ggsave(file.path(dirs[["figures"]], "01_sparsity_matching.pdf"), p_sparsity,
  width = 9, height = 6, device = cairo_pdf, bg = "white")

p_umap <- ggplot() +
  geom_point(data = atac_umap, aes(UMAP_1, UMAP_2), color = "#D0D0D0",
    size = .24, alpha = .32) +
  geom_point(data = trajectory, aes(UMAP_1, UMAP_2, color = native_l4),
    size = .48, alpha = .50) +
  scale_color_viridis_c(option = "inferno", name = "Native L4") +
  geom_point(data = positions, aes(UMAP_1, UMAP_2), color = "#0072B2",
    size = .55, alpha = .10) +
  geom_segment(data = library_summary,
    aes(x = full_UMAP_1, y = full_UMAP_2, xend = UMAP_1_median,
      yend = UMAP_2_median), inherit.aes = FALSE,
    arrow = arrow(length = grid::unit(.08, "inches")), color = "#555555",
    linewidth = .45, alpha = .75) +
  geom_point(data = library_summary, aes(full_UMAP_1, full_UMAP_2),
    inherit.aes = FALSE, shape = 4, color = "#333333", size = 2.5,
    stroke = .8) +
  geom_point(data = library_summary, aes(UMAP_1_median, UMAP_2_median),
    inherit.aes = FALSE, shape = 21, fill = "#56B4E9", color = "#111111",
    size = 3.0, stroke = .6) +
  geom_label_repel(data = library_summary,
    aes(UMAP_1_median, UMAP_2_median, label = gsm), inherit.aes = FALSE,
    seed = as.integer(cfg$adapter$seed), max.overlaps = Inf, size = 2.25,
    box.padding = .25, point.padding = .25, min.segment.length = 0,
    label.size = .12, fill = grDevices::adjustcolor("white", .88),
    max.time = 15) +
  labs(title = "Label-blind sparsity adapter on the fixed Bdev ATAC UMAP",
    subtitle = "Blue clouds are 100 pseudo-cells; crosses are full-bulk positions; arrows end at adapted medians",
    x = "ATAC UMAP 1", y = "ATAC UMAP 2",
    caption = "UMAP coordinates are unchanged; pseudo-cells use PeakVI-neighbor-weighted barycentric placement.") +
  coord_equal(clip = "off") + theme_classic(base_size = 11) +
  theme(plot.margin = margin(8, 35, 8, 8))
ggsave(file.path(dirs[["figures"]], "02_fixed_umap_adapter_positions.png"),
  p_umap, width = 14, height = 10, dpi = 220, bg = "white")
ggsave(file.path(dirs[["figures"]], "02_fixed_umap_adapter_positions.pdf"),
  p_umap, width = 14, height = 10, device = cairo_pdf, bg = "white")

positions$gsm <- factor(positions$gsm, levels = query_cells)
p_l4 <- ggplot(positions, aes(gsm, native_l4_median)) +
  geom_violin(fill = "#56B4E9", color = "#333333", linewidth = .3,
    scale = "width") +
  geom_boxplot(width = .13, outlier.shape = NA, fill = "white", linewidth = .3) +
  geom_point(data = library_summary, aes(gsm, full_native_l4),
    inherit.aes = FALSE, shape = 4, color = "#D55E00", size = 2.2,
    stroke = .9) +
  labs(title = "Native-L4 distributions after label-blind adaptation",
    subtitle = "Orange crosses show the original full-bulk values; GSM order is not biologically arranged",
    x = NULL, y = "Neighbor-weighted native L4") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 70,
    hjust = 1, vjust = 1), panel.grid.minor = element_blank())
ggsave(file.path(dirs[["figures"]], "03_library_l4_distributions.png"), p_l4,
  width = 14, height = 7, dpi = 220, bg = "white")
ggsave(file.path(dirs[["figures"]], "03_library_l4_distributions.pdf"), p_l4,
  width = 14, height = 7, device = cairo_pdf, bg = "white")

p_convergence <- ggplot(convergence, aes(pseudocells_used,
  distance_from_100_pseudocell_center, group = gsm)) +
  geom_line(alpha = .35, color = "#0072B2") +
  geom_point(alpha = .65, color = "#0072B2", size = 1.2) +
  stat_summary(aes(group = 1), fun = median, geom = "line",
    color = "#D55E00", linewidth = 1.1) +
  scale_x_continuous(breaks = prefixes) +
  labs(title = "Adapter median-position convergence",
    subtitle = "Blue lines are libraries; orange is the across-library median",
    x = "Seeded pseudo-cells used", y = "Distance from the 100-pseudo-cell median") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave(file.path(dirs[["figures"]], "04_pseudocell_number_convergence.png"),
  p_convergence, width = 8, height = 6, dpi = 220, bg = "white")
ggsave(file.path(dirs[["figures"]], "04_pseudocell_number_convergence.pdf"),
  p_convergence, width = 8, height = 6, device = cairo_pdf, bg = "white")

source_manifest_after <- sha256(source_manifest)
umap_manifest_after <- sha256(umap_manifest)
stopifnot(identical(source_manifest_before, source_manifest_after),
  identical(umap_manifest_before, umap_manifest_after))
conv50 <- convergence[convergence$pseudocells_used == min(50L, n_replicates), ]
summary <- list(
  schema_version = "bdev.immgen_label_blind_sparsity_adapter_summary.v1",
  status = "complete_diagnostic", run_identity = cfg$run_identity,
  libraries = length(query_cells), pseudocells = nrow(design),
  target_fragment_count = target_depth,
  reference_sparsity_pool_cells = length(sparsity_pool),
  sparsity_pool_nonzero_peaks = quantile_row(sparsity_pool, "pool"),
  adapted_nonzero_peaks = quantile_row(design$observed_nonzero_peaks, "adapted"),
  depth_only_nonzero_peaks = quantile_row(depth_only_design$nonzero_peaks,
    "depth_only"),
  median_rank_accessibility_spearman = median(
    rank_audit$accessibility_frequency_vs_bulk_rank_spearman),
  median_rank_fragment_spearman = median(
    rank_audit$fragment_total_vs_bulk_rank_spearman),
  median_fixed_umap_shift_from_full = median(
    library_summary$fixed_umap_shift_from_full),
  median_absolute_l4_shift_from_full = median(abs(
    library_summary$native_l4_shift_from_full), na.rm = TRUE),
  median_50_to_100_position_difference = median(
    conv50$distance_from_100_pseudocell_center),
  consensus_unique_reference_neighbors = nrow(reuse),
  maximum_consensus_library_reuse = max(reuse$libraries_selecting_cell),
  adapter_construction = list(
    bdev_inputs = c("global reference median fragment count",
      "global reference sparsity among cells within 10% of that count"),
    bulk_input = "within-library peak count rank",
    labels_used = FALSE, expected_positions_used = FALSE,
    metadata_loaded_for_construction = FALSE),
  fixed_map_semantics = "PeakVI-neighbor-weighted barycenter on unchanged approved ATAC UMAP; not umap_transform",
  model_semantics = "qualified deterministic recipe reconstructed; exact reference latent checksum required; query projected without retraining",
  source_manifests_unchanged = TRUE, scientific_status = "diagnostic_only",
  mapping_accepted = FALSE, paper_analysis = FALSE,
  implementation_unlocked = FALSE, sif_built = FALSE)
write_json(summary,
  file.path(dirs[["reports"]], "01_sparsity_adapter_summary.json"))
write_json(c(summary, list(config_path = cfg_path,
  config_sha256 = sha256(cfg_path), source_manifest_sha256 = source_manifest_after,
  approved_umap_manifest_sha256 = umap_manifest_after,
  reference_latent_sha256 = fit_summary$latent_sha256[[1L]])),
  file.path(dirs[["manifests"]], "01_execution_receipt.json"))
wb_close_bundle(bundle)
cat("Completed label-blind sparsity adapter for", nrow(design),
  "pseudo-cells across", length(query_cells), "libraries\n")
