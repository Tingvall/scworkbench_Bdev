#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 01_run_median_depth_test.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(digest)
  library(Matrix)
  library(ggplot2)
})

read_config <- function(path) {
  current <- yaml::read_yaml(path)
  base_name <- current$base_config
  current$base_config <- NULL
  if (is.null(base_name)) return(current)
  utils::modifyList(read_config(normalizePath(file.path(dirname(path), base_name))),
    current, keep.null = TRUE)
}
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
promoted_ref <- function(promoted, output, kind) {
  selected <- promoted$promoted$status == "ok" & promoted$promoted$output == output
  stopifnot(sum(selected) == 1L)
  wb_ref_artifact(promoted$promoted$artifact_id[selected][[1L]], kind)
}
weighted_quantile <- function(value, weight, probabilities) {
  keep <- is.finite(value) & is.finite(weight) & weight > 0
  value <- value[keep]; weight <- weight[keep]
  index <- order(value, seq_along(value), method = "radix")
  value <- value[index]; weight <- weight[index]
  cumulative <- cumsum(weight) / sum(weight)
  vapply(probabilities, function(p) value[[which(cumulative >= p)[[1L]]]], numeric(1))
}
pairwise_jaccard <- function(groups) {
  if (length(groups) < 2L) return(NA_real_)
  values <- numeric(); cursor <- 0L
  for (i in seq_len(length(groups) - 1L)) for (j in (i + 1L):length(groups)) {
    cursor <- cursor + 1L
    values[[cursor]] <- length(intersect(groups[[i]], groups[[j]])) /
      length(union(groups[[i]], groups[[j]]))
  }
  mean(values)
}
reuse_stats <- function(ledger, label) {
  occurrence <- table(ledger$cell_id)
  data.frame(projection = label, assignments = nrow(ledger),
    unique_reference_cells = length(occurrence),
    cells_in_2plus_libraries = sum(occurrence >= 2L),
    cells_in_5plus_libraries = sum(occurrence >= 5L),
    cells_in_10plus_libraries = sum(occurrence >= 10L),
    maximum_library_occurrence = max(occurrence), stringsAsFactors = FALSE)
}

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- read_config(cfg_path)
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
sink(file.path(dirs[["logs"]], "01_median_depth_test.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256),
  identical(sha256(file.path(cfg$inputs$source_mapping_bundle$path, "manifest.json")),
    cfg$inputs$source_mapping_bundle$manifest_sha256),
  identical(sha256(cfg$inputs$full_depth_neighbors$path),
    cfg$inputs$full_depth_neighbors$sha256),
  identical(sha256(cfg$inputs$full_depth_mapping$path),
    cfg$inputs$full_depth_mapping$sha256),
  identical(sha256(cfg$inputs$authenticated_l4$path),
    cfg$inputs$authenticated_l4$sha256),
  identical(cfg$inputs$authenticated_l4$orientation, "native"))
for (package in c("scwbRepresentation", "scwbLatent")) {
  stopifnot(as.character(packageVersion(package)) == cfg$runtime[[package]]$version)
  manifest <- system.file("scworkbench", "operation-pack.json", package = package)
  stopifnot(identical(sha256(manifest), cfg$runtime[[package]]$manifest_sha256))
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
}
stopifnot(all(c("scwblatent::peakvi_fit", "scwblatent::peakvi_project",
  "scwblatent::model_bound_latent_union", "scwbrepresentation::knn_neighbors") %in%
  wb_list_ops()))

source_manifest_before <- sha256(file.path(cfg$inputs$source_mapping_bundle$path,
  "manifest.json"))
mapping <- read.delim(cfg$inputs$full_depth_mapping$path, check.names = FALSE,
  quote = "\"", na.strings = "NA")
full_neighbors <- read.delim(cfg$inputs$full_depth_neighbors$path,
  check.names = FALSE, quote = "\"", na.strings = "NA")
full_neighbors <- full_neighbors[full_neighbors$neighbor_rank <=
  as.integer(cfg$mapping$primary_k), , drop = FALSE]

source <- wb_open_bundle(cfg$inputs$source_mapping_bundle$path, validate = FALSE,
  read_only = TRUE)
reference_space <- wb_space(source, cfg$inputs$source_mapping_bundle$reference_space)
reference_cells <- rownames(reference_space)
query_cells <- as.character(mapping$gsm)
feature_ids <- wb_feature_ids(source, cfg$inputs$source_mapping_bundle$reference_matrix)
reference_counts <- wb_slice(source, cfg$inputs$source_mapping_bundle$reference_matrix,
  features = feature_ids, cells = reference_cells)
bulk_counts <- wb_slice(source, cfg$inputs$source_mapping_bundle$bulk_matrix,
  features = feature_ids, cells = query_cells)
reference_meta <- wb_meta(source, cells = reference_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$sample_field))
wb_close_bundle(source)
stopifnot(ncol(bulk_counts) == cfg$acceptance$expected_libraries,
  nrow(bulk_counts) == cfg$inputs$expected_features,
  all(Matrix::colSums(bulk_counts) >= cfg$downsampling$target_matrix_depth))

target <- as.integer(cfg$downsampling$target_matrix_depth)
n_replicates <- as.integer(cfg$downsampling$replicates_per_library)
set.seed(as.integer(cfg$downsampling$seed))
pseudo_columns <- vector("list", length(query_cells) * n_replicates)
design_rows <- vector("list", length(pseudo_columns)); cursor <- 0L
for (gsm in query_cells) {
  values <- bulk_counts[, gsm, drop = FALSE]
  nz <- which(values[, 1L] > 0)
  probability <- as.numeric(values[nz, 1L])
  probability <- probability / sum(probability)
  for (replicate in seq_len(n_replicates)) {
    cursor <- cursor + 1L
    pseudo_id <- sprintf("%s__median_r%02d", gsm, replicate)
    sampled <- as.numeric(stats::rmultinom(1L, size = target, prob = probability))
    keep <- which(sampled > 0)
    pseudo_columns[[cursor]] <- Matrix::sparseMatrix(i = nz[keep],
      j = rep.int(1L, length(keep)),
      x = sampled[keep], dims = c(length(feature_ids), 1L),
      dimnames = list(feature_ids, pseudo_id))
    design_rows[[cursor]] <- data.frame(query_cell_id = pseudo_id, gsm = gsm,
      population = mapping$population[match(gsm, mapping$gsm)],
      replicate = replicate, seed_stream_position = cursor,
      target_matrix_depth = target, observed_matrix_depth = sum(sampled),
      nonzero_peaks = length(keep), source_matrix_depth = Matrix::colSums(values),
      source_nonzero_peaks = length(nz), stringsAsFactors = FALSE)
  }
}
pseudo_counts <- do.call(cbind, pseudo_columns)
design <- do.call(rbind, design_rows); rownames(design) <- NULL
stopifnot(ncol(pseudo_counts) == cfg$acceptance$expected_pseudocells,
  identical(colnames(pseudo_counts), design$query_cell_id),
  all(Matrix::colSums(pseudo_counts) == target),
  all(design$observed_matrix_depth == target))
write_tsv(design, file.path(dirs[["tables"]], "01_pseudocell_design.tsv"))

if (dir.exists(cfg$storage$test_bundle))
  stop("Refusing to overwrite existing test bundle: ", cfg$storage$test_bundle)
bundle <- wb_init_bundle(cfg$storage$test_bundle, overwrite = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
all_cells <- c(reference_cells, design$query_cell_id)
wb_define_cells(bundle, all_cells)
wb_attach_matrix(bundle, "Bdev_ATAC_reference_counts", reference_counts,
  layer = "counts", feature_type = "peak", cell_ids = reference_cells,
  description = "Qualified Bdev-only PeakVI reference counts",
  source = cfg$inputs$source_mapping_bundle$path,
  created_by = "scworkbench_Bdev median-depth ImmGen diagnostic")
wb_attach_matrix(bundle, "ImmGen_ATAC_median_depth_pseudocells", pseudo_counts,
  layer = "counts", feature_type = "peak", cell_ids = design$query_cell_id,
  description = "Thirty 18,952-count multinomial pseudo-cells per ImmGen library",
  source = cfg$inputs$source_mapping_bundle$path,
  created_by = "scworkbench_Bdev median-depth ImmGen diagnostic")
pseudo_meta <- data.frame(cell_id = design$query_cell_id, stringsAsFactors = FALSE)
pseudo_meta[[cfg$inputs$sample_field]] <- design$gsm
metadata <- rbind(reference_meta, pseudo_meta)
metadata <- metadata[match(all_cells, metadata$cell_id), , drop = FALSE]
wb_attach_metadata(bundle, metadata, cell_ids = all_cells,
  source = "qualified Bdev metadata plus ImmGen pseudo-cell design",
  created_by = "scworkbench_Bdev median-depth ImmGen diagnostic")
rm(reference_counts, bulk_counts, pseudo_counts, pseudo_columns); invisible(gc())

fit_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_peakvi_fit"))
fit <- wb_session_run_op(fit_session, bundle, "scwblatent::peakvi_fit",
  inputs = list(matrix = "Bdev_ATAC_reference_counts",
    n_latent = as.integer(cfg$peakvi$n_latent),
    n_layers = as.integer(cfg$peakvi$n_layers),
    n_hidden = as.integer(cfg$peakvi$n_hidden), dropout = as.numeric(cfg$peakvi$dropout),
    max_epochs = as.integer(cfg$peakvi$max_epochs),
    lr = as.numeric(cfg$peakvi$learning_rate),
    early_stopping = isTRUE(cfg$peakvi$early_stopping),
    seed = as.integer(cfg$peakvi$seed), chunk_size = as.integer(cfg$peakvi$chunk_size),
    batch_key = NULL, use_gpu = FALSE, cell_set = "all", min_cells = 1L),
  notes = "Deterministic qualified Bdev-only reference reconstruction")
reference_space_id <- as.character(fit$run$outputs$space$space_id %||%
  fit$run$outputs$space$allocated_id)[[1L]]
fit_summary <- fit$run$outputs$fit_summary$table
stopifnot(fit_summary$latent_sha256[[1L]] ==
  cfg$peakvi$expected_reference_latent_sha256)
fit$session <- wb_session_resume(bundle, fit$session$durable$workspace)
fit_promoted <- wb_session_promote_current(bundle, fit$session,
  outputs = c("model", "fit_summary"), name_prefix = "median_depth_peakvi_reference",
  description = "Checksum-authenticated Bdev-only PeakVI reference",
  created_by = "scworkbench_Bdev median-depth ImmGen diagnostic", fail_fast = TRUE)
model_ref <- promoted_ref(fit_promoted, "model", "model")
fit_promoted$session <- wb_session_discard(fit_promoted$session)

project_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_peakvi_project"))
project <- wb_session_run_op(project_session, bundle, "scwblatent::peakvi_project",
  inputs = list(model = model_ref,
    query_matrix = "ImmGen_ATAC_median_depth_pseudocells",
    space_id_out = "median_depth_immgen_peakvi_queries",
    chunk_size = as.integer(cfg$peakvi$chunk_size), batch_key = NULL,
    use_gpu = FALSE, pickle_execution_trust = cfg$peakvi$pickle_execution_trust,
    cell_set = "all", min_cells = 1L),
  notes = "Frozen-recipe projection of 840 median-depth ImmGen pseudo-cells")
stopifnot(project$run$outputs$projection_summary$table$transformation[[1L]] ==
  "transfer_without_refit")
project$session <- wb_session_resume(bundle, project$session$durable$workspace)
project_promoted <- wb_session_promote_current(bundle, project$session,
  outputs = "projection_summary", name_prefix = "median_depth_immgen_projection",
  description = "Median-depth ImmGen pseudo-cell PeakVI projection",
  created_by = "scworkbench_Bdev median-depth ImmGen diagnostic", fail_fast = TRUE)
projection_ref <- promoted_ref(project_promoted, "projection_summary", "obs_table")
project_promoted$session <- wb_session_discard(project_promoted$session)

union_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_latent_union"))
union <- wb_session_run_op(union_session, bundle,
  "scwblatent::model_bound_latent_union", inputs = list(
    reference_space = reference_space_id,
    query_space = "median_depth_immgen_peakvi_queries", model = model_ref,
    projection_summary = projection_ref,
    space_id_out = "median_depth_reference_and_queries",
    chunk_size = as.integer(cfg$mapping$query_chunk_size)),
  notes = "Authenticated Bdev/median-depth-ImmGen latent union")
union$session <- wb_session_resume(bundle, union$session$durable$workspace)
union_promoted <- wb_session_promote_current(bundle, union$session,
  outputs = "union_summary", name_prefix = "median_depth_immgen_union",
  description = "Authenticated median-depth model-bound latent union",
  created_by = "scworkbench_Bdev median-depth ImmGen diagnostic", fail_fast = TRUE)
union_promoted$session <- wb_session_discard(union_promoted$session)

knn_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_knn"))
knn <- wb_session_run_op(knn_session, bundle,
  "scwbrepresentation::knn_neighbors", inputs = list(
    space = "median_depth_reference_and_queries",
    query_cells = design$query_cell_id, candidate_cells = reference_cells,
    k = as.integer(cfg$mapping$neighbor_k), include_query = FALSE,
    backend = cfg$mapping$backend,
    space_chunk_size = as.integer(cfg$mapping$space_chunk_size),
    query_chunk_size = as.integer(cfg$mapping$query_chunk_size)),
  notes = "Median-depth ImmGen pseudo-cell neighbors in qualified PeakVI space")
neighbors <- knn$run$outputs$neighbors$table
neighbors <- neighbors[neighbors$neighbor_rank <= as.integer(cfg$mapping$primary_k), ]
knn$session <- wb_session_resume(bundle, knn$session$durable$workspace)
knn$session <- wb_session_discard(knn$session)

family <- readRDS(cfg$inputs$authenticated_l4$path)
selected <- as.character(family$key) == cfg$inputs$authenticated_l4$key &
  as.character(family$cell_id) %in% reference_cells &
  is.finite(as.numeric(family$value))
l4 <- stats::setNames(as.numeric(family$value[selected]),
  as.character(family$cell_id[selected]))
stopifnot(length(l4) == length(reference_cells),
  all(neighbors$cell_id %in% names(l4)))
supported_threshold <- unique(mapping$supported_threshold)
out_threshold <- unique(mapping$out_of_reference_threshold)
stopifnot(length(supported_threshold) == 1L, length(out_threshold) == 1L)

prediction_rows <- vector("list", nrow(design))
neighbors$kernel_weight <- NA_real_; neighbors$normalized_weight <- NA_real_
for (i in seq_len(nrow(design))) {
  cell <- design$query_cell_id[[i]]
  index <- which(neighbors$query_cell_id == cell)
  index <- index[order(neighbors$neighbor_rank[index])]
  distance <- neighbors$distance[index]
  statistic <- median(distance)
  support <- if (statistic <= supported_threshold) "supported" else if (
    statistic <= out_threshold) "borderline" else "out_of_reference"
  kth <- distance[[length(distance)]]
  weight <- if (kth == 0) as.numeric(distance == 0) else exp(-(distance / kth)^2)
  neighbors$kernel_weight[index] <- weight
  neighbors$normalized_weight[index] <- weight / sum(weight)
  quantiles <- weighted_quantile(unname(l4[neighbors$cell_id[index]]), weight,
    c(.25, .5, .75))
  prediction_rows[[i]] <- data.frame(query_cell_id = cell, gsm = design$gsm[[i]],
    population = design$population[[i]], replicate = design$replicate[[i]],
    support_call = support, distance_statistic = statistic,
    weighted_q25 = quantiles[[1L]], mapped_native_l4 = if (
      support == "out_of_reference") NA_real_ else quantiles[[2L]],
    weighted_q75 = quantiles[[3L]], stringsAsFactors = FALSE)
}
prediction <- do.call(rbind, prediction_rows)
write_tsv(neighbors, file.path(dirs[["tables"]], "02_pseudocell_neighbor_ledger.tsv"))
write_tsv(prediction, file.path(dirs[["tables"]], "03_pseudocell_l4_predictions.tsv"))

library_rows <- lapply(query_cells, function(gsm) {
  p <- prediction[prediction$gsm == gsm, , drop = FALSE]
  ledger <- neighbors[neighbors$query_cell_id %in% p$query_cell_id, , drop = FALSE]
  replicate_sets <- split(as.character(ledger$cell_id), ledger$query_cell_id)
  finite <- is.finite(p$mapped_native_l4)
  data.frame(gsm = gsm, population = p$population[[1L]],
    n_replicates = nrow(p), supported_fraction = mean(p$support_call == "supported"),
    borderline_fraction = mean(p$support_call == "borderline"),
    out_of_reference_fraction = mean(p$support_call == "out_of_reference"),
    median_distance_statistic = median(p$distance_statistic),
    downsample_l4_q25 = if (any(finite)) quantile(p$mapped_native_l4[finite], .25) else NA,
    downsample_l4_median = if (any(finite)) median(p$mapped_native_l4[finite]) else NA,
    downsample_l4_q75 = if (any(finite)) quantile(p$mapped_native_l4[finite], .75) else NA,
    downsample_l4_min = if (any(finite)) min(p$mapped_native_l4[finite]) else NA,
    downsample_l4_max = if (any(finite)) max(p$mapped_native_l4[finite]) else NA,
    replicate_neighbor_mean_jaccard = pairwise_jaccard(replicate_sets),
    unique_reference_neighbors = length(unique(ledger$cell_id)),
    full_depth_support_call = mapping$support_call[match(gsm, mapping$gsm)],
    full_depth_l4 = mapping$accepted_native_l4[match(gsm, mapping$gsm)],
    stringsAsFactors = FALSE)
})
library_summary <- do.call(rbind, library_rows)
library_summary$l4_shift <- library_summary$downsample_l4_median -
  library_summary$full_depth_l4

consensus_rows <- list(); cursor <- 0L
for (gsm in query_cells) {
  query_ids <- design$query_cell_id[design$gsm == gsm]
  x <- neighbors[neighbors$query_cell_id %in% query_ids, , drop = FALSE]
  cells <- unique(as.character(x$cell_id))
  rows <- lapply(cells, function(cell) {
    take <- x$cell_id == cell
    data.frame(gsm = gsm,
      population = mapping$population[match(gsm, mapping$gsm)], cell_id = cell,
      selection_fraction = sum(take) / n_replicates,
      mean_normalized_weight = sum(x$normalized_weight[take]) / n_replicates,
      mean_distance_when_selected = mean(x$distance[take]), stringsAsFactors = FALSE)
  })
  ranked <- do.call(rbind, rows)
  ranked <- ranked[order(-ranked$selection_fraction,
    -ranked$mean_normalized_weight, ranked$cell_id, method = "radix"), ]
  ranked <- ranked[seq_len(as.integer(cfg$mapping$primary_k)), ]
  ranked$consensus_rank <- seq_len(nrow(ranked))
  cursor <- cursor + 1L; consensus_rows[[cursor]] <- ranked
}
consensus <- do.call(rbind, consensus_rows); rownames(consensus) <- NULL
write_tsv(consensus, file.path(dirs[["tables"]], "04_library_consensus_neighbors.tsv"))

library_summary$full_consensus_neighbor_jaccard <- vapply(query_cells, function(gsm) {
  a <- full_neighbors$cell_id[full_neighbors$query_cell_id == gsm]
  b <- consensus$cell_id[consensus$gsm == gsm]
  length(intersect(a, b)) / length(union(a, b))
}, numeric(1))
write_tsv(library_summary, file.path(dirs[["tables"]],
  "05_library_downsample_comparison.tsv"))

full_reuse <- full_neighbors[, c("query_cell_id", "cell_id")]
consensus_reuse <- consensus[, c("gsm", "cell_id")]
names(consensus_reuse)[[1L]] <- "query_cell_id"
reuse <- rbind(reuse_stats(full_reuse, "full_depth"),
  reuse_stats(consensus_reuse, "median_depth_consensus"))
write_tsv(reuse, file.path(dirs[["tables"]], "06_neighbor_reuse_comparison.tsv"))

prediction$population <- factor(prediction$population,
  levels = unique(mapping$population))
mapping$population <- factor(mapping$population, levels = levels(prediction$population))
p_l4 <- ggplot(prediction[is.finite(prediction$mapped_native_l4), ],
  aes(mapped_native_l4)) +
  geom_density(aes(fill = population), alpha = 0.65, color = "#333333") +
  geom_vline(data = mapping, aes(xintercept = accepted_native_l4,
    color = support_call), linewidth = 0.65) +
  facet_wrap(~ population, ncol = 4, scales = "free_y") +
  scale_fill_discrete(guide = "none") +
  scale_color_manual(values = c(supported = "#0072B2", borderline = "#D55E00")) +
  labs(title = "Median-depth ImmGen pseudo-cell L4 distributions",
    subtitle = "Thirty 18,952-count projections per library; vertical lines show full-depth L4",
    x = "Mapped native L4", y = "Density", color = "Full-depth call",
    caption = "Diagnostic only; the original and downsampled mappings are not accepted downstream.") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 8.5),
    panel.grid.minor = element_blank())
ggsave(file.path(dirs[["figures"]], "01_population_downsample_l4.png"), p_l4,
  width = 14, height = 11, dpi = 220, bg = "white")
ggsave(file.path(dirs[["figures"]], "01_population_downsample_l4.pdf"), p_l4,
  width = 14, height = 11, device = cairo_pdf, bg = "white")

p_compare <- ggplot(library_summary, aes(full_depth_l4, downsample_l4_median)) +
  geom_abline(slope = 1, intercept = 0, color = "#777777", linetype = 2) +
  geom_linerange(aes(ymin = downsample_l4_q25, ymax = downsample_l4_q75,
    color = supported_fraction), linewidth = 0.7) +
  geom_point(aes(shape = full_depth_support_call, color = supported_fraction),
    size = 2.8, stroke = 0.8) +
  geom_text(aes(label = population), size = 2.4, nudge_x = 0.15,
    check_overlap = TRUE) +
  scale_color_viridis_c(option = "C", limits = c(0, 1),
    name = "Downsampled\nsupported fraction") +
  scale_shape_manual(values = c(supported = 16, borderline = 17),
    name = "Full-depth call") +
  coord_equal() +
  labs(title = "Full-depth versus median-depth ImmGen L4",
    subtitle = "Points are medians and vertical bars are Q25–Q75 across 30 pseudo-cells",
    x = "Full-depth mapped L4", y = "Median-depth mapped L4") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(dirs[["figures"]], "02_full_vs_downsample_l4.png"), p_compare,
  width = 9, height = 8, dpi = 220, bg = "white")
ggsave(file.path(dirs[["figures"]], "02_full_vs_downsample_l4.pdf"), p_compare,
  width = 9, height = 8, device = cairo_pdf, bg = "white")

reuse_plot_rows <- rbind(
  data.frame(projection = "full depth", occurrence = as.integer(table(
    factor(full_reuse$cell_id, levels = unique(full_reuse$cell_id))))),
  data.frame(projection = "median-depth consensus", occurrence = as.integer(table(
    factor(consensus_reuse$cell_id, levels = unique(consensus_reuse$cell_id))))))
p_reuse <- ggplot(reuse_plot_rows, aes(occurrence, fill = projection)) +
  geom_histogram(binwidth = 1, boundary = 0.5, position = "dodge", color = "white") +
  scale_x_continuous(breaks = seq_len(max(reuse_plot_rows$occurrence))) +
  labs(title = "Reference-cell reuse across the 28 ImmGen libraries",
    subtitle = "Median-depth uses each library's top-30 consensus across 30 pseudo-cells",
    x = "Number of libraries selecting a reference cell", y = "Reference cells",
    fill = "Projection") + theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(dirs[["figures"]], "03_neighbor_reuse_comparison.png"), p_reuse,
  width = 9, height = 6, dpi = 220, bg = "white")
ggsave(file.path(dirs[["figures"]], "03_neighbor_reuse_comparison.pdf"), p_reuse,
  width = 9, height = 6, device = cairo_pdf, bg = "white")

source_manifest_after <- sha256(file.path(cfg$inputs$source_mapping_bundle$path,
  "manifest.json"))
stopifnot(identical(source_manifest_before, source_manifest_after))
summary <- list(schema_version = "bdev.immgen_median_depth_test_summary.v1",
  status = "complete_diagnostic", run_identity = cfg$run_identity,
  target_matrix_depth = target, replicates_per_library = n_replicates,
  libraries = length(query_cells), pseudocells = nrow(design),
  pseudo_nonzero_peaks = as.list(stats::setNames(as.numeric(quantile(
    design$nonzero_peaks, c(.25, .5, .75))), c("q25", "median", "q75"))),
  full_depth_reuse = as.list(reuse[reuse$projection == "full_depth", ]),
  median_depth_consensus_reuse = as.list(reuse[
    reuse$projection == "median_depth_consensus", ]),
  median_library_neighbor_stability = median(
    library_summary$replicate_neighbor_mean_jaccard),
  median_full_consensus_jaccard = median(
    library_summary$full_consensus_neighbor_jaccard),
  median_absolute_l4_shift = median(abs(library_summary$l4_shift), na.rm = TRUE),
  libraries_all_replicates_supported = sum(library_summary$supported_fraction == 1),
  source_manifest_unchanged = TRUE, scientific_status = "diagnostic_only",
  paper_analysis = FALSE, implementation_unlocked = FALSE, sif_built = FALSE)
write_json(summary, file.path(dirs[["reports"]], "01_median_depth_test_summary.json"))
write_json(c(summary, list(config_path = cfg_path, config_sha256 = sha256(cfg_path),
  source_manifest_sha256 = source_manifest_after)),
  file.path(dirs[["manifests"]], "01_execution_receipt.json"))
wb_close_bundle(bundle)
cat("Completed median-depth test for", nrow(design), "pseudo-cells\n")
