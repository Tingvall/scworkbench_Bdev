#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 04_plot_comparison.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(digest)
  library(scworkbench)
  library(slingshot)
  library(ggplot2)
  library(ggrepel)
})

read_config <- function(path) {
  current <- yaml::read_yaml(path)
  base <- current$base_config
  current$base_config <- NULL
  if (is.null(base)) return(current)
  utils::modifyList(read_config(normalizePath(file.path(dirname(path), base))),
    current, keep.null = TRUE)
}
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_tsv <- function(x, path, compress = FALSE) {
  con <- if (compress) gzfile(path, open = "wt") else path
  on.exit(if (compress) close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
space_2d <- function(x) {
  if (is.data.frame(x) && all(c("cell_id", "x", "y") %in% names(x)))
    return(data.frame(cell_id = as.character(x$cell_id), UMAP_1 = x$x,
      UMAP_2 = x$y, stringsAsFactors = FALSE))
  x <- as.matrix(x)
  stopifnot(ncol(x) == 2L, !is.null(rownames(x)))
  data.frame(cell_id = rownames(x), UMAP_1 = x[, 1L], UMAP_2 = x[, 2L],
    stringsAsFactors = FALSE)
}
named_character <- function(x) {
  answer <- unlist(x, use.names = TRUE)
  structure(as.character(answer), names = names(answer))
}

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- read_config(cfg_path)
set.seed(as.integer(cfg$visualization$label_seed))
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "03_plot_comparison.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

mapping_path <- file.path(dirs[["tables"]],
  "09_variance_vs_depth_adjusted_mappings.tsv")
neighbor_path <- file.path(dirs[["tables"]],
  "10_variance_vs_depth_adjusted_neighbors.tsv.gz")
metric_path <- file.path(dirs[["tables"]],
  "11_variance_vs_depth_adjusted_metrics.tsv")
execution_path <- file.path(dirs[["manifests"]],
  "03_external_execution_receipt.json")
input_paths <- c(
  bdev_manifest = file.path(cfg$inputs$source_bundle$path, "manifest.json"),
  l4_model = cfg$inputs$authenticated_l4$model_path,
  l4_pseudotime = cfg$inputs$authenticated_l4$pseudotime_path,
  mappings = mapping_path, neighbors = neighbor_path, metrics = metric_path,
  execution = execution_path)
observed_sha <- vapply(input_paths, sha256, character(1L))
stopifnot(observed_sha[["bdev_manifest"]] ==
    cfg$inputs$source_bundle$manifest_sha256,
  observed_sha[["l4_model"]] == cfg$inputs$authenticated_l4$model_sha256,
  observed_sha[["l4_pseudotime"]] ==
    cfg$inputs$authenticated_l4$pseudotime_sha256,
  cfg$inputs$authenticated_l4$orientation == "native")
execution <- jsonlite::read_json(execution_path, simplifyVector = TRUE)
stopifnot(execution$status == "complete",
  isTRUE(execution$source_manifest_unchanged),
  isTRUE(execution$query_manifest_unchanged),
  isFALSE(execution$current_default_changed))

mappings <- read.delim(mapping_path, check.names = FALSE, quote = "\"",
  na.strings = "NA")
neighbors <- read.delim(gzfile(neighbor_path), check.names = FALSE,
  quote = "\"", na.strings = "NA")
metrics <- read.delim(metric_path, check.names = FALSE, quote = "\"",
  na.strings = "NA")
candidate_order <- unlist(cfg$visualization$candidate_order, use.names = FALSE)
candidate_labels <- named_character(cfg$visualization$candidate_labels)
stopifnot(identical(sort(unique(mappings$candidate)), sort(candidate_order)),
  nrow(mappings) == 2L * cfg$external_evaluation$query_bundle$expected_libraries,
  nrow(neighbors) == 2L * cfg$external_evaluation$query_bundle$expected_libraries *
    cfg$external_evaluation$neighbor_k,
  all(table(mappings$candidate) ==
    cfg$external_evaluation$query_bundle$expected_libraries),
  all(table(interaction(neighbors$candidate, neighbors$gsm, drop = TRUE)) ==
    cfg$external_evaluation$neighbor_k),
  identical(sort(unique(metrics$candidate)), sort(candidate_order)))

bdev <- wb_open_bundle(cfg$inputs$source_bundle$path, validate = FALSE,
  read_only = TRUE)
atac_umap <- space_2d(wb_space(bdev, cfg$visualization$fixed_space))
wb_close_bundle(bdev)
stopifnot(nrow(atac_umap) == cfg$inputs$source_bundle$expected_cells,
  !anyDuplicated(atac_umap$cell_id),
  all(neighbors$reference_cell_id %in% atac_umap$cell_id))
rownames(atac_umap) <- atac_umap$cell_id

l4_family <- readRDS(cfg$inputs$authenticated_l4$pseudotime_path)
l4_keep <- as.character(l4_family$key) == cfg$inputs$authenticated_l4$key &
  is.finite(as.numeric(l4_family$value)) &
  as.character(l4_family$cell_id) %in% atac_umap$cell_id
l4_cells <- data.frame(cell_id = as.character(l4_family$cell_id[l4_keep]),
  native_l4 = as.numeric(l4_family$value[l4_keep]), stringsAsFactors = FALSE)
l4_cells <- merge(l4_cells, atac_umap, by = "cell_id", sort = FALSE)
l4_model <- readRDS(cfg$inputs$authenticated_l4$model_path)
l4_curves <- slingshot::slingCurves(l4_model)
stopifnot(cfg$inputs$authenticated_l4$curve_name %in% names(l4_curves))
raw_curve <- l4_curves[[cfg$inputs$authenticated_l4$curve_name]]
curve_xy <- raw_curve$s[raw_curve$ord, , drop = FALSE]
stopifnot(ncol(curve_xy) == 2L, nrow(curve_xy) > 2L,
  all(is.finite(curve_xy)))
l4_curve <- data.frame(curve_order = seq_len(nrow(curve_xy)),
  UMAP_1 = curve_xy[, 1L], UMAP_2 = curve_xy[, 2L],
  orientation = "native", stringsAsFactors = FALSE)
curve_arrow_points <- tail(l4_curve, 2L)
curve_arrow <- data.frame(x = curve_arrow_points$UMAP_1[[1L]],
  y = curve_arrow_points$UMAP_2[[1L]],
  xend = curve_arrow_points$UMAP_1[[2L]],
  yend = curve_arrow_points$UMAP_2[[2L]])
curve_ends <- rbind(
  transform(l4_curve[1L, ], endpoint = "native L4 start"),
  transform(l4_curve[nrow(l4_curve), ], endpoint = "native L4 end"))

neighbors$UMAP_1 <- atac_umap[neighbors$reference_cell_id, "UMAP_1"]
neighbors$UMAP_2 <- atac_umap[neighbors$reference_cell_id, "UMAP_2"]
neighbor_groups <- split(neighbors, interaction(neighbors$candidate,
  neighbors$gsm, drop = TRUE, lex.order = TRUE))
ellipse_scale <- as.numeric(cfg$visualization$ellipse_scale)
project_group <- function(rows) {
  rows <- rows[order(rows$rank), , drop = FALSE]
  kth <- max(rows$distance)
  weight <- if (kth == 0) as.numeric(rows$distance == 0) else
    exp(-(rows$distance / kth)^2)
  xy <- as.matrix(rows[, c("UMAP_1", "UMAP_2")])
  center <- colSums(xy * weight) / sum(weight)
  centered <- sweep(xy, 2L, center, "-")
  weighted_cov <- crossprod(centered * sqrt(weight / sum(weight)))
  eig <- eigen(weighted_cov, symmetric = TRUE)
  radii <- ellipse_scale * sqrt(pmax(eig$values, 0))
  theta <- seq(0, 2 * pi, length.out = 101L)
  ellipse_xy <- sweep(cbind(cos(theta), sin(theta)) %*%
    diag(radii, nrow = 2L) %*% t(eig$vectors), 2L, center, "+")
  list(placement = data.frame(candidate = rows$candidate[[1L]],
      gsm = rows$gsm[[1L]], UMAP_1 = center[[1L]], UMAP_2 = center[[2L]],
      neighbor_spread = sqrt(sum(weight * rowSums(centered^2)) / sum(weight)),
      nearest_distance = min(rows$distance), kth_distance = kth,
      stringsAsFactors = FALSE),
    ellipse = data.frame(candidate = rows$candidate[[1L]],
      gsm = rows$gsm[[1L]], ellipse_order = seq_along(theta),
      UMAP_1 = ellipse_xy[, 1L], UMAP_2 = ellipse_xy[, 2L],
      stringsAsFactors = FALSE),
    weights = data.frame(candidate = rows$candidate, gsm = rows$gsm,
      rank = rows$rank, reference_cell_id = rows$reference_cell_id,
      kernel_weight = weight, stringsAsFactors = FALSE))
}
projected <- lapply(neighbor_groups, project_group)
placements <- do.call(rbind, lapply(projected, `[[`, "placement"))
ellipses <- do.call(rbind, lapply(projected, `[[`, "ellipse"))
kernel_weights <- do.call(rbind, lapply(projected, `[[`, "weights"))
rownames(placements) <- rownames(ellipses) <- rownames(kernel_weights) <- NULL

map_fields <- c("candidate", "gsm", "slot_id", "population", "lineage",
  "predicted_bdev_cluster", "nearest_distance_over_reference_median",
  "neighbor_l4_weight", "provisional_native_l4", "trajectory_output_status")
placements <- merge(placements, mappings[, map_fields],
  by = c("candidate", "gsm"), sort = FALSE)
ellipses <- merge(ellipses,
  mappings[, c("candidate", "gsm", "population", "lineage")],
  by = c("candidate", "gsm"), sort = FALSE)
placements$candidate <- factor(placements$candidate, levels = candidate_order)
ellipses$candidate <- factor(ellipses$candidate, levels = candidate_order)
placements$candidate_label <- factor(
  candidate_labels[as.character(placements$candidate)],
  levels = candidate_labels[candidate_order])
ellipses$candidate_label <- factor(
  candidate_labels[as.character(ellipses$candidate)],
  levels = candidate_labels[candidate_order])
placements$detail_label <- paste0(placements$population, " · ", placements$gsm)

spoke_rows <- neighbors[neighbors$rank <=
  as.integer(cfg$visualization$detailed_spokes_per_library), ]
spokes <- merge(spoke_rows, placements[, c("candidate", "gsm", "population",
  "UMAP_1", "UMAP_2")], by = c("candidate", "gsm"),
  suffixes = c("_neighbor", "_query"), sort = FALSE)
population_order <- c("LTHSC.34+.BM", "LTHSC.34-.BM", "STHSC.150-.BM",
  "MPP3.48+.BM", "MPP4.135+.BM", "proB.CLP.BM", "proB.FrA.BM",
  "proB.FrBC.BM", "B.FrE.BM", "B.T1.Sp", "B.T2.Sp", "B.T3.Sp",
  "B.Fo.Sp", "B.MZ.Sp", "B.Sp", "B1b.PC")
stopifnot(identical(sort(unique(placements$population)), sort(population_order)))
population_colors <- c(
  "LTHSC.34+.BM" = "#3B4CC0", "LTHSC.34-.BM" = "#5967B5",
  "STHSC.150-.BM" = "#748FFC", "MPP3.48+.BM" = "#3A86C8",
  "MPP4.135+.BM" = "#48BFE3", "proB.CLP.BM" = "#168A83",
  "proB.FrA.BM" = "#43AA8B", "proB.FrBC.BM" = "#79B473",
  "B.FrE.BM" = "#A7C957", "B.T1.Sp" = "#E9C46A",
  "B.T2.Sp" = "#F4A261", "B.T3.Sp" = "#E76F51",
  "B.Fo.Sp" = "#D1495B", "B.MZ.Sp" = "#9B5DE5",
  "B.Sp" = "#7A5195", "B1b.PC" = "#5C4D7D")
placements$population <- factor(placements$population, levels = population_order)
ellipses$population <- factor(ellipses$population, levels = population_order)
spokes$population <- factor(spokes$population, levels = population_order)

padding <- as.numeric(cfg$visualization$fixed_coordinate_padding_fraction)
x_range <- range(atac_umap$UMAP_1); y_range <- range(atac_umap$UMAP_2)
x_limits <- x_range + c(-1, 1) * diff(x_range) * padding
y_limits <- y_range + c(-1, 1) * diff(y_range) * padding
base_map <- function() {
  ggplot() +
    geom_point(data = atac_umap, aes(UMAP_1, UMAP_2), color = "#D5D5D5",
      size = 0.22, alpha = 0.35) +
    geom_point(data = l4_cells, aes(UMAP_1, UMAP_2), color = "#8E8E8E",
      size = 0.32, alpha = 0.42) +
    geom_path(data = l4_curve, aes(UMAP_1, UMAP_2), color = "white",
      linewidth = 2.2, lineend = "round") +
    geom_path(data = l4_curve, aes(UMAP_1, UMAP_2), color = "#151515",
      linewidth = 0.9, lineend = "round") +
    geom_segment(data = curve_arrow,
      aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE,
      color = "#151515", linewidth = 0.9,
      arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed"))
}

overview <- base_map() +
  geom_path(data = ellipses,
    aes(UMAP_1, UMAP_2, group = interaction(candidate, gsm),
      color = population), linewidth = 0.42, alpha = 0.33) +
  geom_point(data = placements,
    aes(UMAP_1, UMAP_2, fill = population, shape = lineage),
    color = "#111111", stroke = 0.45, size = 2.8, alpha = 0.96) +
  geom_label(data = curve_ends, aes(UMAP_1, UMAP_2, label = endpoint),
    inherit.aes = FALSE, size = 2.4, label.size = 0.15,
    fill = grDevices::adjustcolor("white", 0.82), color = "#111111") +
  scale_color_manual(values = population_colors, guide = "none") +
  scale_fill_manual(values = population_colors, name = "ImmGen population") +
  scale_shape_manual(values = c("B" = 21, "Stem&Prog" = 24), name = "Lineage") +
  facet_wrap(~candidate_label, ncol = 2L) +
  coord_equal(xlim = x_limits, ylim = y_limits, clip = "off") +
  labs(title = "Raw-count variance versus depth-adjusted peak selection",
    subtitle = paste("Same frozen ImmGen libraries on one fixed Bdev ATAC UMAP;",
      "points are 30-neighbor kernel barycenters"),
    x = "ATAC UMAP 1", y = "ATAC UMAP 2",
    caption = paste("Black arrow: authenticated native L4. Same 26,723 cells,",
      "PeakVI architecture, coordinate limits, and projection without refit.")) +
  theme_classic(base_size = 10.5) +
  theme(strip.text = element_text(face = "bold", size = 10.5),
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "bottom", legend.box = "vertical",
    legend.text = element_text(size = 7.5),
    plot.margin = margin(8, 20, 8, 8)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
      override.aes = list(size = 3, shape = 21, color = "#111111")),
    shape = guide_legend(nrow = 1))
ggsave(file.path(dirs[["figures"]],
  "01_variance_vs_depth_adjusted_fixed_atac_umap.png"), overview,
  width = 15, height = 8.5, dpi = 240, bg = "white")
ggsave(file.path(dirs[["figures"]],
  "01_variance_vs_depth_adjusted_fixed_atac_umap.pdf"), overview,
  width = 15, height = 8.5, device = cairo_pdf, bg = "white")

detail_files <- character()
for (i in seq_along(candidate_order)) {
  candidate <- candidate_order[[i]]
  points_i <- placements[as.character(placements$candidate) == candidate, ]
  ellipses_i <- ellipses[as.character(ellipses$candidate) == candidate, ]
  spokes_i <- spokes[as.character(spokes$candidate) == candidate, ]
  detail <- base_map() +
    geom_segment(data = spokes_i,
      aes(x = UMAP_1_query, y = UMAP_2_query,
        xend = UMAP_1_neighbor, yend = UMAP_2_neighbor,
        color = population), linewidth = 0.25, alpha = 0.18) +
    geom_path(data = ellipses_i,
      aes(UMAP_1, UMAP_2, group = gsm, color = population),
      linewidth = 0.65, alpha = 0.4) +
    geom_point(data = points_i,
      aes(UMAP_1, UMAP_2, fill = population, shape = lineage),
      color = "#111111", stroke = 0.65, size = 4.0) +
    geom_label_repel(data = points_i,
      aes(UMAP_1, UMAP_2, label = detail_label, color = population),
      seed = as.integer(cfg$visualization$label_seed), max.overlaps = Inf,
      box.padding = 0.38, point.padding = 0.5, force = 2,
      min.segment.length = 0, segment.color = "#555555", segment.size = 0.28,
      size = 2.55, label.size = 0.15,
      fill = grDevices::adjustcolor("white", 0.9), max.time = 20) +
    geom_label(data = curve_ends, aes(UMAP_1, UMAP_2, label = endpoint),
      inherit.aes = FALSE, size = 2.5, label.size = 0.15,
      fill = grDevices::adjustcolor("white", 0.82), color = "#111111") +
    scale_color_manual(values = population_colors, guide = "none") +
    scale_fill_manual(values = population_colors, guide = "none") +
    scale_shape_manual(values = c("B" = 21, "Stem&Prog" = 24), guide = "none") +
    coord_equal(xlim = x_limits, ylim = y_limits, clip = "off") +
    labs(title = paste0(candidate_labels[[candidate]],
        " on the fixed Bdev ATAC UMAP"),
      subtitle = paste("All 28 libraries labeled; ellipses use all 30",
        "neighbors and faint spokes show the nearest five"),
      x = "ATAC UMAP 1", y = "ATAC UMAP 2",
      caption = paste("Frozen-neighbor visualization; native L4 is diagnostic",
        "and does not define the whole-reference mapping.")) +
    theme_classic(base_size = 11.5) +
    theme(plot.title = element_text(face = "bold", size = 15),
      plot.margin = margin(10, 45, 10, 10))
  stem <- sprintf("%02d_%s_fixed_atac_umap_detail", i + 1L, candidate)
  png_path <- file.path(dirs[["figures"]], paste0(stem, ".png"))
  pdf_path <- file.path(dirs[["figures"]], paste0(stem, ".pdf"))
  ggsave(png_path, detail, width = 14, height = 10.5, dpi = 240, bg = "white")
  ggsave(pdf_path, detail, width = 14, height = 10.5,
    device = cairo_pdf, bg = "white")
  detail_files <- c(detail_files, basename(png_path), basename(pdf_path))
}

metric_spec <- data.frame(
  field = c("within_population_to_between_population_latent_distance_ratio",
    "replicate_predicted_cluster_agreement", "unique_reference_neighbors",
    "maximum_reference_hub_fraction",
    "median_nearest_distance_over_reference_median",
    "median_neighbor_sample_entropy",
    "maximum_absolute_depth_latent_spearman",
    "maximum_absolute_saturation_latent_spearman"),
  label = c("Within/between population distance (lower)",
    "Replicate cluster agreement (higher)",
    "Unique reference neighbors (higher)",
    "Maximum reference hub fraction (lower)",
    "Median support-distance ratio (lower)",
    "Neighbor sample entropy (higher)",
    "Maximum |depth–latent rho| (lower)",
    "Maximum |saturation–latent rho| (lower)"),
  stringsAsFactors = FALSE)
metric_long <- do.call(rbind, lapply(seq_len(nrow(metric_spec)), function(i)
  data.frame(candidate = metrics$candidate,
    candidate_label = candidate_labels[metrics$candidate],
    metric = metric_spec$label[[i]],
    value = as.numeric(metrics[[metric_spec$field[[i]]]]),
    stringsAsFactors = FALSE)))
metric_long$candidate_label <- factor(metric_long$candidate_label,
  levels = candidate_labels[candidate_order])
metric_long$metric <- factor(metric_long$metric, levels = rev(metric_spec$label))
metric_plot <- ggplot(metric_long,
    aes(value, candidate_label, color = candidate_label)) +
  geom_line(aes(group = metric), color = "#A8A8A8", linewidth = 0.45) +
  geom_point(size = 3.0) +
  geom_text(aes(label = formatC(value, digits = 3, format = "fg")),
    nudge_y = 0.16, size = 2.7, show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_x", ncol = 2L) +
  scale_color_manual(values = c("#0072B2", "#D55E00"), guide = "none") +
  labs(title = "Frozen bulk-mapping diagnostics",
    subtitle = "Raw-count variance20k compared directly with depth-adjusted overall20k",
    x = "Observed metric", y = NULL,
    caption = "Direction in parentheses is preferred; no metric was used to tune either frozen model.") +
  theme_classic(base_size = 10.5) +
  theme(strip.text = element_text(face = "bold", size = 8.5),
    plot.title = element_text(face = "bold", size = 15),
    panel.spacing = grid::unit(1, "lines"))
ggsave(file.path(dirs[["figures"]], "04_external_metric_comparison.png"),
  metric_plot, width = 13, height = 10, dpi = 240, bg = "white")
ggsave(file.path(dirs[["figures"]], "04_external_metric_comparison.pdf"),
  metric_plot, width = 13, height = 10, device = cairo_pdf, bg = "white")

trajectory <- mappings
trajectory$candidate_label <- factor(candidate_labels[trajectory$candidate],
  levels = candidate_labels[candidate_order])
trajectory$population <- factor(trajectory$population,
  levels = rev(population_order))
l4_range <- range(l4_cells$native_l4, finite = TRUE)
na_position <- l4_range[[1L]] - 0.08 * diff(l4_range)
trajectory$l4_plot_value <- ifelse(is.finite(trajectory$provisional_native_l4),
  trajectory$provisional_native_l4, na_position)
trajectory$l4_call <- ifelse(is.finite(trajectory$provisional_native_l4),
  "finite diagnostic value", "no finite value")
trajectory_plot <- ggplot(trajectory,
    aes(l4_plot_value, population, fill = lineage, shape = l4_call)) +
  geom_vline(xintercept = l4_range[[1L]], color = "#9A9A9A", linetype = 3,
    linewidth = 0.4) +
  geom_point(position = position_jitter(height = 0.12, width = 0),
    size = 2.4, stroke = 0.55, color = "#222222", alpha = 0.88) +
  facet_wrap(~candidate_label, ncol = 2L) +
  scale_shape_manual(values = c("finite diagnostic value" = 21,
    "no finite value" = 4), name = "Post-freeze L4 output") +
  scale_fill_manual(values = c("B" = "#D1495B", "Stem&Prog" = "#3B4CC0"),
    name = "ImmGen lineage") +
  scale_x_continuous(breaks = c(na_position, pretty(l4_range)),
    labels = c("NA", format(pretty(l4_range), trim = TRUE))) +
  labs(title = "Diagnostic native-L4 positions",
    subtitle = "Trajectory assignment is a second-stage diagnostic, not the mapping target",
    x = "Provisional native-L4 pseudotime", y = "ImmGen population",
    caption = "Diagnostic/unaccepted output; native Slingshot orientation is preserved.") +
  theme_classic(base_size = 10.5) +
  theme(strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "bottom")
ggsave(file.path(dirs[["figures"]], "05_diagnostic_native_l4.png"),
  trajectory_plot, width = 14, height = 9, dpi = 240, bg = "white")
ggsave(file.path(dirs[["figures"]], "05_diagnostic_native_l4.pdf"),
  trajectory_plot, width = 14, height = 9, device = cairo_pdf, bg = "white")

placements_out <- placements
placements_out$candidate <- as.character(placements_out$candidate)
placements_out$candidate_label <- as.character(placements_out$candidate_label)
placements_out$population <- as.character(placements_out$population)
ellipses_out <- ellipses
ellipses_out$candidate <- as.character(ellipses_out$candidate)
ellipses_out$candidate_label <- as.character(ellipses_out$candidate_label)
ellipses_out$population <- as.character(ellipses_out$population)
write_tsv(placements_out, file.path(dirs[["tables"]],
  "13_fixed_umap_candidate_placements.tsv"))
write_tsv(ellipses_out, file.path(dirs[["tables"]],
  "14_fixed_umap_neighbor_ellipses.tsv.gz"), compress = TRUE)
write_tsv(l4_curve, file.path(dirs[["tables"]],
  "15_authenticated_native_l4_curve.tsv"))
write_tsv(kernel_weights, file.path(dirs[["tables"]],
  "16_neighbor_kernel_weights.tsv.gz"), compress = TRUE)

input_sha_after <- vapply(input_paths, sha256, character(1L))
stopifnot(identical(observed_sha, input_sha_after),
  nrow(placements) == nrow(mappings))
summary <- list(
  schema_version = "bdev.variance20k_vs_depth_adjusted_figures.v1",
  status = "complete_read_only_visualization",
  run_identity = cfg$run_identity,
  fixed_space = cfg$visualization$fixed_space,
  fixed_reference_cells = nrow(atac_umap),
  candidates = candidate_order,
  libraries_per_candidate = as.list(table(placements$candidate)),
  neighbors_per_library = cfg$external_evaluation$neighbor_k,
  placement_kernel = cfg$visualization$placement_kernel,
  authenticated_l4_orientation = "native",
  trajectory_output_status = "diagnostic_unaccepted",
  model_refit = FALSE,
  benchmark_inputs_mutated = FALSE,
  current_default_changed = FALSE,
  implementation_unlocked = FALSE,
  sif_built = FALSE)
write_json(summary, file.path(dirs[["reports"]],
  "02_variance_vs_depth_adjusted_figure_summary.json"))
write_json(c(summary, list(config_path = cfg_path,
  config_sha256 = sha256(cfg_path), input_sha256 = as.list(input_sha_after),
  detail_figures = as.list(detail_files))),
  file.path(dirs[["manifests"]], "04_figure_execution_receipt.json"))
cat("Completed fixed-UMAP, metric, and diagnostic L4 figures\n")
