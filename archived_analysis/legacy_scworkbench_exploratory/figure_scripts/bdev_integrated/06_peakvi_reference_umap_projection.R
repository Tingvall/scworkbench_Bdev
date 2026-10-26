#!/usr/bin/env Rscript

# Display KO query cells on the unchanged native Bdev ATAC UMAP by taking a
# distance-weighted barycenter of their Bdev-restricted PeakVI neighbors.
# This is a report-only projection; it does not create or promote a bundle space.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "Usage: 06_peakvi_reference_umap_projection.R",
    "BDEV_BUNDLE QUERY_MAPPING_TSV_GZ NEIGHBORS_TSV_GZ MAPPING_METRICS_JSON REPORT_DIR"
  ), call. = FALSE)
}

bdev_dir <- normalizePath(args[[1L]], mustWork = TRUE)
query_mapping_path <- normalizePath(args[[2L]], mustWork = TRUE)
neighbors_path <- normalizePath(args[[3L]], mustWork = TRUE)
mapping_metrics_path <- normalizePath(args[[4L]], mustWork = TRUE)
report_dir <- normalizePath(args[[5L]], mustWork = FALSE)
if (file.exists(report_dir) || dir.exists(report_dir)) {
  stop("Report directory already exists; refusing to overwrite: ", report_dir,
       call. = FALSE)
}

options(width = 180, warn = 1)
suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(report_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

write_json <- function(x, path) jsonlite::write_json(
  x, path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null",
  digits = NA
)
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

bdev <- wb_open_bundle(bdev_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
manifest_before <- wb_manifest(bdev)
fingerprint_before <- manifest_before$current_state$fingerprint
stopifnot(identical(manifest_before$bundle_id, "bdev_reference"))

bdev_cells <- wb_cells(bdev, "all")
bdev_md <- as.data.table(wb_meta(bdev, cells = bdev_cells))
bdev_umap <- as.data.table(wb_space(bdev, "umap.atac", cells = bdev_cells))
stopifnot(length(bdev_cells) == 26723L,
          identical(as.character(bdev_md$cell_id), bdev_cells),
          identical(as.character(bdev_umap$cell_id), bdev_cells),
          all(is.finite(bdev_umap$x)), all(is.finite(bdev_umap$y)))

reference <- data.table(
  reference_cell_id = paste0("bdev::", bdev_cells),
  source_cell_id = bdev_cells,
  x = bdev_umap$x,
  y = bdev_umap$y,
  reference_cluster = as.character(bdev_md[["peaks_snn_res.1.2"]]),
  reference_pt = as.numeric(bdev_md$pt)
)
reference_xy_original <- copy(reference[, .(reference_cell_id, x, y)])

query <- fread(cmd = paste("zcat", shQuote(query_mapping_path)))
neighbors <- fread(cmd = paste("zcat", shQuote(neighbors_path)))
required_query <- c("cell_id", "source_cell_id", "sample_id", "genotype",
                    "knockout_target", "support_status", "median30_distance",
                    "predicted_cluster_k30", "mapped_pt")
required_neighbors <- c("query_cell_id", "rank", "reference_cell_id", "distance")
stopifnot(all(required_query %in% names(query)),
          all(required_neighbors %in% names(neighbors)),
          nrow(query) == 23252L, !anyDuplicated(query$cell_id),
          nrow(neighbors) == 23252L * 50L,
          setequal(neighbors$query_cell_id, query$cell_id),
          all(neighbors$rank %in% seq_len(50L)),
          all(is.finite(neighbors$distance)), all(neighbors$distance >= 0),
          all(neighbors$reference_cell_id %in% reference$reference_cell_id))
rank_check <- neighbors[, .(n = .N, min_rank = min(rank), max_rank = max(rank),
                            unique_rank = uniqueN(rank)), by = query_cell_id]
stopifnot(all(rank_check$n == 50L), all(rank_check$min_rank == 1L),
          all(rank_check$max_rank == 50L), all(rank_check$unique_rank == 50L))

neighbors <- reference[neighbors, on = "reference_cell_id"]
stopifnot(nrow(neighbors) == 23252L * 50L,
          all(is.finite(neighbors$x)), all(is.finite(neighbors$y)))

project_k <- function(dt, k) {
  z <- dt[rank <= k]
  z[, sigma := max(distance), by = query_cell_id]
  z[sigma <= 0 | !is.finite(sigma), sigma := 1]
  z[, weight := exp(-(distance / sigma)^2)]
  ans <- z[, {
    sw <- sum(weight)
    px <- sum(weight * x) / sw
    py <- sum(weight * y) / sw
    list(
      projected_x = px,
      projected_y = py,
      unweighted_x = mean(x),
      unweighted_y = mean(y),
      neighbor_umap_dispersion = sqrt(sum(weight * ((x - px)^2 + (y - py)^2)) / sw),
      effective_neighbor_n = sw^2 / sum(weight^2)
    )
  }, by = query_cell_id]
  setnames(ans, setdiff(names(ans), "query_cell_id"),
           paste0(setdiff(names(ans), "query_cell_id"), "_k", k))
  ans
}

projected <- Reduce(function(x, y) merge(x, y, by = "query_cell_id", sort = FALSE),
                    lapply(c(15L, 30L, 50L), function(k) project_k(neighbors, k)))
projected[, sensitivity_max_distance := pmax(
  sqrt((projected_x_k15 - projected_x_k30)^2 + (projected_y_k15 - projected_y_k30)^2),
  sqrt((projected_x_k15 - projected_x_k50)^2 + (projected_y_k15 - projected_y_k50)^2),
  sqrt((projected_x_k30 - projected_x_k50)^2 + (projected_y_k30 - projected_y_k50)^2)
)]
setnames(projected, "query_cell_id", "cell_id")
projected <- query[projected, on = "cell_id"]
stopifnot(nrow(projected) == nrow(query),
          identical(projected$cell_id, query$cell_id),
          all(is.finite(projected$projected_x_k30)),
          all(is.finite(projected$projected_y_k30)))

# Primary coordinates are deliberately withheld for calibrated OOD cells. Raw
# weighted coordinates remain available for auditing, not biological display.
projected[, `:=`(
  projected_x = fifelse(support_status == "out_of_reference", NA_real_, projected_x_k30),
  projected_y = fifelse(support_status == "out_of_reference", NA_real_, projected_y_k30)
)]
stopifnot(sum(is.na(projected$projected_x)) == sum(projected$support_status == "out_of_reference"),
          sum(is.na(projected$projected_y)) == sum(projected$support_status == "out_of_reference"),
          projected[is.finite(projected_x), all(projected_x >= min(reference$x) & projected_x <= max(reference$x))],
          projected[is.finite(projected_y), all(projected_y >= min(reference$y) & projected_y <= max(reference$y))])

fwrite(projected, file.path(report_dir, "query_reference_umap_projection.tsv.gz"), sep = "\t")

summary_sample <- projected[, .(
  n_cells = .N,
  n_displayed = sum(is.finite(projected_x)),
  supported = sum(support_status == "supported"),
  borderline = sum(support_status == "borderline"),
  out_of_reference = sum(support_status == "out_of_reference"),
  median_neighbor_umap_dispersion = median(neighbor_umap_dispersion_k30),
  q95_neighbor_umap_dispersion = quantile(neighbor_umap_dispersion_k30, 0.95),
  median_k_sensitivity = median(sensitivity_max_distance),
  q95_k_sensitivity = quantile(sensitivity_max_distance, 0.95)
), by = .(sample_id, genotype, knockout_target)]
fwrite(summary_sample, file.path(report_dir, "projection_summary_by_sample.tsv"), sep = "\t")

reference_plot <- reference[, .(x, y)]
query_plot <- projected[is.finite(projected_x), .(
  x = projected_x, y = projected_y, sample_id, genotype, knockout_target,
  support_status, predicted_cluster_k30 = as.character(predicted_cluster_k30),
  mapped_pt, neighbor_umap_dispersion_k30, sensitivity_max_distance
)]

base_theme <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = "right")
save_plot <- function(plot, name, width = 10, height = 7) {
  ggsave(file.path(fig_dir, name), plot, width = width, height = height,
         dpi = 180, bg = "white")
}

p_sample <- ggplot() +
  geom_point(data = reference_plot, aes(x, y), color = "grey82", size = 0.12) +
  geom_point(data = query_plot, aes(x, y, color = genotype), size = 0.25, alpha = 0.55) +
  facet_wrap(~sample_id, ncol = 2) +
  scale_color_manual(values = c(KO = "#D55E00", WT = "#0072B2")) +
  coord_equal() + base_theme +
  labs(title = "KO-query projection onto the unchanged Bdev ATAC UMAP",
       subtitle = "Distance-weighted mean of 30 Bdev neighbors in PeakVI latent space; OOD cells withheld",
       x = "Bdev umap.atac 1", y = "Bdev umap.atac 2", color = "Genotype")
save_plot(p_sample, "reference_umap_projection_by_sample.png", 11, 8)

p_support <- ggplot() +
  geom_point(data = reference_plot, aes(x, y), color = "grey85", size = 0.12) +
  geom_point(data = query_plot, aes(x, y, color = support_status), size = 0.25, alpha = 0.55) +
  scale_color_manual(values = c(supported = "#0072B2", borderline = "#E69F00")) +
  coord_equal() + base_theme +
  labs(title = "Reference-anchored query locations by mapping support",
       subtitle = "Out-of-reference cells are not assigned primary display coordinates",
       x = "Bdev umap.atac 1", y = "Bdev umap.atac 2", color = "Support")
save_plot(p_support, "reference_umap_projection_support.png")

cluster_levels <- sort(unique(c(reference$reference_cluster, query_plot$predicted_cluster_k30)))
palette <- setNames(grDevices::hcl.colors(length(cluster_levels), "Dynamic"), cluster_levels)
p_cluster <- ggplot() +
  geom_point(data = reference, aes(x, y, color = reference_cluster), size = 0.16, alpha = 0.35) +
  geom_point(data = query_plot, aes(x, y, color = predicted_cluster_k30), size = 0.22, alpha = 0.55) +
  scale_color_manual(values = palette) + coord_equal() + base_theme +
  labs(title = "Transferred peak-cluster identities on native Bdev geometry",
       subtitle = "Reference: native cluster; query: k=30 PeakVI-neighbor transfer",
       x = "Bdev umap.atac 1", y = "Bdev umap.atac 2", color = "Peak cluster")
save_plot(p_cluster, "reference_umap_projection_cluster.png")

p_uncertainty <- ggplot(query_plot, aes(x, y, color = neighbor_umap_dispersion_k30)) +
  geom_point(size = 0.28, alpha = 0.65) +
  scale_color_viridis_c(option = "magma", direction = -1) + coord_equal() + base_theme +
  labs(title = "Projection uncertainty on the native Bdev UMAP",
       subtitle = "Weighted RMS spread of each cell's 30 reference-neighbor coordinates",
       x = "Bdev umap.atac 1", y = "Bdev umap.atac 2", color = "Neighbor\nspread")
save_plot(p_uncertainty, "reference_umap_projection_uncertainty.png")

p_density <- ggplot() +
  geom_point(data = reference_plot, aes(x, y), color = "grey88", size = 0.1) +
  stat_density_2d(data = query_plot, aes(x, y, color = genotype), bins = 8,
                  linewidth = 0.45) +
  facet_wrap(~sample_id, ncol = 2) +
  scale_color_manual(values = c(KO = "#D55E00", WT = "#0072B2")) +
  coord_equal() + base_theme +
  labs(title = "Projected query density on the unchanged Bdev ATAC UMAP",
       subtitle = "Contours reduce overplotting where many query cells share similar barycentric coordinates",
       x = "Bdev umap.atac 1", y = "Bdev umap.atac 2", color = "Genotype")
save_plot(p_density, "reference_umap_projection_density_by_sample.png", 11, 8)

manifest_after <- wb_manifest(bdev)
stopifnot(identical(fingerprint_before, manifest_after$current_state$fingerprint),
          identical(reference_xy_original, reference[, .(reference_cell_id, x, y)]))

mapping_metrics <- jsonlite::read_json(mapping_metrics_path, simplifyVector = TRUE)
support_counts <- as.list(table(projected$support_status))
dispersion_quantiles <- as.list(setNames(
  as.numeric(quantile(projected$neighbor_umap_dispersion_k30, c(.5, .9, .95, .99))),
  c("q50", "q90", "q95", "q99")
))
sensitivity_quantiles <- as.list(setNames(
  as.numeric(quantile(projected$sensitivity_max_distance, c(.5, .9, .95, .99))),
  c("q50", "q90", "q95", "q99")
))

completion <- list(
  status = "complete_unpromoted_report_only",
  method = list(
    neighbor_space = "20-dimensional PeakVI latent space",
    neighbor_pool = "Bdev reference only",
    primary_k = 30L,
    sensitivity_k = c(15L, 30L, 50L),
    weight = "exp(-(distance / kth_distance)^2)",
    display_space = "unchanged bdev_reference umap.atac",
    primary_ood_policy = "withhold coordinates for out_of_reference cells"
  ),
  cells = list(reference = nrow(reference), query = nrow(projected),
               displayed_query = sum(is.finite(projected$projected_x)),
               support = support_counts),
  diagnostics = list(neighbor_umap_dispersion = dispersion_quantiles,
                     k_sensitivity_max_pairwise_distance = sensitivity_quantiles,
                     upstream_mapping_metrics = mapping_metrics),
  preservation = list(bdev_fingerprint_before = fingerprint_before,
                      bdev_fingerprint_after = manifest_after$current_state$fingerprint,
                      bdev_unchanged = TRUE, original_umap_coordinates_unchanged = TRUE,
                      promotion = FALSE),
  inputs = list(
    bdev_bundle = bdev_dir,
    query_mapping = list(path = query_mapping_path, sha256 = sha256_file(query_mapping_path)),
    neighbors = list(path = neighbors_path, sha256 = sha256_file(neighbors_path)),
    mapping_metrics = list(path = mapping_metrics_path, sha256 = sha256_file(mapping_metrics_path))
  ),
  output = list(table = "query_reference_umap_projection.tsv.gz",
                sample_summary = "projection_summary_by_sample.tsv",
                figures = list.files(fig_dir))
)
write_json(completion, file.path(report_dir, "completion.json"))

report <- c(
  "# PeakVI-neighbor projection onto the native Bdev ATAC UMAP", "",
  "Status: complete, report-only, and unpromoted.", "",
  sprintf("The original `bdev_reference` `umap.atac` coordinates for all %s reference cells are displayed without modification.", format(nrow(reference), big.mark = ",")),
  sprintf("Each of %s query cells received k=15, 30, and 50 distance-weighted barycentric coordinates from its exact Bdev-restricted PeakVI neighbors; k=30 is primary.", format(nrow(projected), big.mark = ",")),
  sprintf("Primary coordinates are shown for %s supported or borderline cells. The %s calibrated out-of-reference cells retain audit coordinates but are withheld from the primary display.", format(sum(is.finite(projected$projected_x)), big.mark = ","), format(sum(projected$support_status == "out_of_reference"), big.mark = ",")),
  "", "## Interpretation", "",
  "This is a reference-anchored visualization, not a refitted UMAP and not evidence that the query coordinates are intrinsic UMAP positions. Weighted averaging can place a cell between disconnected reference regions. Use `neighbor_umap_dispersion_k30` and `sensitivity_max_distance` to identify such ambiguous placements.",
  "", "## Retained outputs", "",
  "- `query_reference_umap_projection.tsv.gz`: per-cell raw k=15/30/50 projections, primary OOD-filtered coordinates, neighbor spread, and k sensitivity.",
  "- `projection_summary_by_sample.tsv`: sample-level support and uncertainty summaries.",
  "- `figures/reference_umap_projection_by_sample.png`: primary sample/genotype view.",
  "- `figures/reference_umap_projection_density_by_sample.png`: density-contour view for overplotting.",
  "- `figures/reference_umap_projection_support.png`: supported versus borderline cells.",
  "- `figures/reference_umap_projection_cluster.png`: native/reference and transferred query peak clusters.",
  "- `figures/reference_umap_projection_uncertainty.png`: neighbor-coordinate dispersion.",
  "", "The Bdev bundle fingerprint and original `umap.atac` coordinates were unchanged; nothing was promoted."
)
writeLines(report, file.path(report_dir, "REPORT.md"))

