#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 01_build_diagnostic.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(digest)
  library(ggplot2)
})

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
weighted_quantile <- function(value, weight, probability) {
  keep <- is.finite(value) & is.finite(weight) & weight > 0
  value <- value[keep]; weight <- weight[keep]
  index <- order(value, seq_along(value), method = "radix")
  value <- value[index]; weight <- weight[index]
  value[[which(cumsum(weight) / sum(weight) >= probability)[[1L]]]]
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

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
.libPaths(c(cfg$runtime$run_library, .libPaths()))
suppressPackageStartupMessages(library(scworkbench))

root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_build_diagnostic.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(
  identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256),
  identical(sha256(file.path(cfg$inputs$bdev_bundle$path, "manifest.json")),
    cfg$inputs$bdev_bundle$manifest_sha256),
  identical(sha256(cfg$inputs$mapping_table$path), cfg$inputs$mapping_table$sha256),
  identical(sha256(cfg$inputs$authenticated_l4$path), cfg$inputs$authenticated_l4$sha256),
  identical(sha256(cfg$inputs$neighbor_ledger$path), cfg$inputs$neighbor_ledger$sha256),
  identical(sha256(cfg$inputs$fixed_umap_placement$path),
    cfg$inputs$fixed_umap_placement$sha256),
  identical(cfg$inputs$authenticated_l4$orientation, "native"),
  as.character(packageVersion("scworkbench")) == cfg$runtime$scworkbench_version)
source_manifest_before <- sha256(file.path(cfg$inputs$bdev_bundle$path, "manifest.json"))

mapping <- read.delim(cfg$inputs$mapping_table$path, check.names = FALSE,
  quote = "\"", na.strings = "NA")
placement <- read.delim(cfg$inputs$fixed_umap_placement$path, check.names = FALSE,
  quote = "\"", na.strings = "NA")
neighbors <- read.delim(cfg$inputs$neighbor_ledger$path, check.names = FALSE,
  quote = "\"", na.strings = "NA")
k <- as.integer(cfg$diagnostic$primary_k)
neighbors <- neighbors[neighbors$neighbor_rank <= k, , drop = FALSE]
stopifnot(nrow(mapping) == 28L, length(unique(mapping$population)) == 16L,
  nrow(placement) == 28L, nrow(neighbors) == 28L * k,
  setequal(mapping$gsm, placement$gsm),
  setequal(mapping$gsm, unique(neighbors$query_cell_id)))

family <- readRDS(cfg$inputs$authenticated_l4$path)
take <- as.character(family$key) == cfg$inputs$authenticated_l4$key &
  is.finite(as.numeric(family$value))
reference_l4 <- stats::setNames(as.numeric(family$value[take]),
  as.character(family$cell_id[take]))
stopifnot(!anyDuplicated(names(reference_l4)),
  all(neighbors$cell_id %in% names(reference_l4)))

bdev <- wb_open_bundle(cfg$inputs$bdev_bundle$path, validate = FALSE, read_only = TRUE)
atac_umap <- space_2d(wb_space(bdev, cfg$inputs$bdev_bundle$atac_umap_space))
driver_cells <- unique(as.character(neighbors$cell_id))
driver_meta <- wb_meta(bdev, cells = driver_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), cfg$inputs$bdev_bundle$sample_field,
  cfg$inputs$bdev_bundle$annotation_field))
wb_close_bundle(bdev)
stopifnot(nrow(atac_umap) == cfg$inputs$bdev_bundle$expected_cells,
  all(driver_cells %in% atac_umap$cell_id),
  identical(as.character(driver_meta$cell_id), driver_cells))

neighbors$population <- mapping$population[match(neighbors$query_cell_id, mapping$gsm)]
neighbors$support_call <- mapping$support_call[match(neighbors$query_cell_id, mapping$gsm)]
neighbors$query_accepted_native_l4 <- mapping$accepted_native_l4[
  match(neighbors$query_cell_id, mapping$gsm)]
neighbors$reference_native_l4 <- unname(reference_l4[neighbors$cell_id])
neighbors$reference_UMAP_1 <- atac_umap$UMAP_1[match(neighbors$cell_id, atac_umap$cell_id)]
neighbors$reference_UMAP_2 <- atac_umap$UMAP_2[match(neighbors$cell_id, atac_umap$cell_id)]
neighbors$query_UMAP_1 <- placement$UMAP_1[match(neighbors$query_cell_id, placement$gsm)]
neighbors$query_UMAP_2 <- placement$UMAP_2[match(neighbors$query_cell_id, placement$gsm)]
neighbors$reference_sample <- driver_meta[[cfg$inputs$bdev_bundle$sample_field]][
  match(neighbors$cell_id, driver_meta$cell_id)]
neighbors$reference_cluster <- driver_meta[[cfg$inputs$bdev_bundle$annotation_field]][
  match(neighbors$cell_id, driver_meta$cell_id)]

query_names <- unique(as.character(neighbors$query_cell_id))
neighbors$kernel_weight <- NA_real_
neighbors$normalized_weight <- NA_real_
for (query in query_names) {
  index <- which(neighbors$query_cell_id == query)
  index <- index[order(neighbors$neighbor_rank[index])]
  distance <- neighbors$distance[index]
  kth <- distance[[length(distance)]]
  weight <- if (kth == 0) as.numeric(distance == 0) else exp(-(distance / kth)^2)
  neighbors$kernel_weight[index] <- weight
  neighbors$normalized_weight[index] <- weight / sum(weight)
}
write_tsv(neighbors, file.path(dirs[["tables"]], "01_neighbor_l4_umap_driver_ledger.tsv"))

library_rows <- lapply(query_names, function(query) {
  x <- neighbors[neighbors$query_cell_id == query, , drop = FALSE]
  w <- x$normalized_weight
  data.frame(gsm = query, population = x$population[[1L]],
    support_call = x$support_call[[1L]], n_neighbors = nrow(x),
    effective_neighbors = 1 / sum(w^2),
    neighbor_l4_min = min(x$reference_native_l4),
    neighbor_l4_q25 = weighted_quantile(x$reference_native_l4, w, .25),
    neighbor_l4_median = weighted_quantile(x$reference_native_l4, w, .5),
    neighbor_l4_q75 = weighted_quantile(x$reference_native_l4, w, .75),
    neighbor_l4_max = max(x$reference_native_l4),
    neighbor_l4_weighted_mean = weighted.mean(x$reference_native_l4, w),
    recorded_accepted_native_l4 = x$query_accepted_native_l4[[1L]],
    reconstructed_UMAP_1 = weighted.mean(x$reference_UMAP_1, w),
    reconstructed_UMAP_2 = weighted.mean(x$reference_UMAP_2, w),
    recorded_UMAP_1 = x$query_UMAP_1[[1L]],
    recorded_UMAP_2 = x$query_UMAP_2[[1L]], stringsAsFactors = FALSE)
})
library_summary <- do.call(rbind, library_rows)
library_summary$l4_reconstruction_error <- library_summary$neighbor_l4_median -
  library_summary$recorded_accepted_native_l4
library_summary$coordinate_reconstruction_error <- sqrt(
  (library_summary$reconstructed_UMAP_1 - library_summary$recorded_UMAP_1)^2 +
  (library_summary$reconstructed_UMAP_2 - library_summary$recorded_UMAP_2)^2)
write_tsv(library_summary, file.path(dirs[["tables"]],
  "02_library_l4_umap_reconstruction.tsv"))
stopifnot(max(abs(library_summary$l4_reconstruction_error)) <=
    cfg$diagnostic$l4_tolerance,
  max(library_summary$coordinate_reconstruction_error) <=
    cfg$diagnostic$coordinate_tolerance)

population_names <- unique(as.character(mapping$population))
neighbors$population_weight <- NA_real_
for (population in population_names) {
  libraries <- unique(neighbors$query_cell_id[neighbors$population == population])
  index <- which(neighbors$population == population)
  neighbors$population_weight[index] <- neighbors$normalized_weight[index] /
    length(libraries)
}
population_rows <- lapply(population_names, function(population) {
  x <- neighbors[neighbors$population == population, , drop = FALSE]
  w <- x$population_weight
  library_values <- unique(x[, c("query_cell_id", "query_accepted_native_l4")])
  data.frame(population = population,
    n_libraries = nrow(library_values), n_neighbor_assignments = nrow(x),
    n_unique_reference_cells = length(unique(x$cell_id)),
    effective_neighbor_assignments = 1 / sum(w^2),
    neighbor_l4_min = min(x$reference_native_l4),
    neighbor_l4_q25 = weighted_quantile(x$reference_native_l4, w, .25),
    neighbor_l4_median = weighted_quantile(x$reference_native_l4, w, .5),
    neighbor_l4_q75 = weighted_quantile(x$reference_native_l4, w, .75),
    neighbor_l4_max = max(x$reference_native_l4),
    neighbor_l4_weighted_mean = weighted.mean(x$reference_native_l4, w),
    mapped_library_l4_min = min(library_values$query_accepted_native_l4),
    mapped_library_l4_median = median(library_values$query_accepted_native_l4),
    mapped_library_l4_max = max(library_values$query_accepted_native_l4),
    mapped_library_l4_range = diff(range(library_values$query_accepted_native_l4)),
    stringsAsFactors = FALSE)
})
population_summary <- do.call(rbind, population_rows)
population_summary <- population_summary[order(population_summary$neighbor_l4_median,
  population_summary$population), , drop = FALSE]
write_tsv(population_summary, file.path(dirs[["tables"]],
  "03_population_l4_distribution_summary.tsv"))

cluster_key <- paste(neighbors$population, neighbors$reference_cluster, sep = "\r")
cluster_groups <- split(seq_len(nrow(neighbors)), cluster_key)
cluster_rows <- lapply(cluster_groups, function(index) {
  x <- neighbors[index, , drop = FALSE]
  data.frame(population = x$population[[1L]],
    reference_cluster = x$reference_cluster[[1L]],
    neighbor_assignments = nrow(x), unique_reference_cells = length(unique(x$cell_id)),
    total_population_weight = sum(x$population_weight),
    weighted_mean_l4 = weighted.mean(x$reference_native_l4, x$population_weight),
    stringsAsFactors = FALSE)
})
cluster_summary <- do.call(rbind, cluster_rows)
cluster_summary <- cluster_summary[order(cluster_summary$population,
  -cluster_summary$total_population_weight, cluster_summary$reference_cluster), ]
write_tsv(cluster_summary, file.path(dirs[["tables"]],
  "04_population_reference_cluster_drivers.tsv"))

population_levels <- population_summary$population
neighbors$population <- factor(neighbors$population, levels = population_levels)
mapping$population <- factor(mapping$population, levels = population_levels)
p_l4 <- ggplot(neighbors, aes(reference_native_l4, weight = population_weight)) +
  geom_density(aes(fill = population), adjust = cfg$diagnostic$l4_density_adjust,
    color = "#363636", linewidth = 0.35, alpha = 0.72) +
  geom_vline(data = mapping, aes(xintercept = accepted_native_l4,
    color = support_call), linewidth = 0.65, alpha = 0.8) +
  facet_wrap(~ population, ncol = 4, scales = "free_y") +
  scale_fill_discrete(guide = "none") +
  scale_color_manual(values = c(supported = "#0072B2", borderline = "#D55E00"),
    name = "Mapped library L4") +
  labs(title = "Native L4 distributions of the Bdev neighbors driving each ImmGen population",
    subtitle = "Density pools each library equally; vertical lines are its weighted-median mapped L4",
    x = "Native fitted Slingshot L4 of primary-30 Bdev neighbor", y = "Weighted density",
    caption = "Diagnostic only: the underlying ImmGen mapping was not scientifically accepted.") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 8.5),
    legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(dirs[["figures"]], "01_population_neighbor_l4_distributions.png"),
  p_l4, width = 14, height = 11, dpi = 220, bg = "white")
ggsave(file.path(dirs[["figures"]], "01_population_neighbor_l4_distributions.pdf"),
  p_l4, width = 14, height = 11, device = cairo_pdf, bg = "white")

top_n <- as.integer(cfg$diagnostic$top_driver_segments_per_library)
segment_rows <- do.call(rbind, lapply(query_names, function(query) {
  x <- neighbors[neighbors$query_cell_id == query, , drop = FALSE]
  x[order(-x$normalized_weight, x$neighbor_rank), ][seq_len(top_n), ]
}))
p_umap <- ggplot() +
  geom_point(data = atac_umap, aes(UMAP_1, UMAP_2), color = "#D0D0D0",
    size = 0.12, alpha = 0.28) +
  geom_segment(data = segment_rows, aes(x = reference_UMAP_1, y = reference_UMAP_2,
    xend = query_UMAP_1, yend = query_UMAP_2, alpha = normalized_weight),
    color = "#555555", linewidth = 0.24, show.legend = FALSE) +
  geom_point(data = neighbors, aes(reference_UMAP_1, reference_UMAP_2,
    color = reference_native_l4, size = normalized_weight), alpha = 0.72) +
  geom_point(data = placement, aes(UMAP_1, UMAP_2, shape = support_call),
    fill = "white", color = "black", stroke = 0.9, size = 2.8) +
  facet_wrap(~ population, ncol = 4) +
  scale_color_viridis_c(option = "inferno", name = "Neighbor L4") +
  scale_size_continuous(range = c(0.25, 2.1), name = "Normalized\nPeakVI weight") +
  scale_shape_manual(values = c(supported = 21, borderline = 24),
    name = "Library support") +
  coord_equal() +
  labs(title = "What drives each ImmGen population's location on the fixed ATAC UMAP",
    subtitle = "Colored points are the primary-30 PeakVI neighbors; lines show the ten strongest contributors per library",
    x = "Fixed ATAC UMAP 1", y = "Fixed ATAC UMAP 2",
    caption = "White symbols are recorded weighted-barycenter placements; UMAP was not refit or transformed.") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 8.5),
    panel.grid = element_blank(), legend.position = "bottom")
ggsave(file.path(dirs[["figures"]], "02_population_fixed_umap_drivers.png"),
  p_umap, width = 15, height = 12, dpi = 220, bg = "white")
ggsave(file.path(dirs[["figures"]], "02_population_fixed_umap_drivers.pdf"),
  p_umap, width = 15, height = 12, device = cairo_pdf, bg = "white")

source_manifest_after <- sha256(file.path(cfg$inputs$bdev_bundle$path, "manifest.json"))
stopifnot(identical(source_manifest_before, source_manifest_after))
summary <- list(schema_version = "bdev.immgen_l4_umap_driver_summary.v1",
  status = "complete_diagnostic", run_identity = cfg$run_identity,
  libraries = nrow(mapping), populations = length(population_names), primary_k = k,
  unique_reference_driver_cells = length(driver_cells),
  max_absolute_l4_reconstruction_error = max(abs(library_summary$l4_reconstruction_error)),
  max_coordinate_reconstruction_error = max(library_summary$coordinate_reconstruction_error),
  l4_estimator = "distance-kernel-weighted median of primary-30 Bdev-neighbor native L4",
  fixed_umap_estimator = "distance-kernel-weighted arithmetic mean of primary-30 neighbor fixed umap.atac coordinates",
  weight_kernel = "exp(-(distance / kth_distance)^2)",
  source_manifest_unchanged = TRUE, native_l4_orientation = "native",
  scientific_status = "diagnostic_only_mapping_not_accepted",
  paper_analysis = FALSE, implementation_unlocked = FALSE, sif_built = FALSE)
write_json(summary, file.path(dirs[["reports"]], "01_diagnostic_summary.json"))
write_json(c(summary, list(config_path = cfg_path, config_sha256 = sha256(cfg_path),
  bdev_manifest_sha256 = source_manifest_after)),
  file.path(dirs[["manifests"]], "01_execution_receipt.json"))
cat("Completed L4/fixed-UMAP driver diagnostic for", length(population_names),
  "populations and", nrow(mapping), "libraries\n")

