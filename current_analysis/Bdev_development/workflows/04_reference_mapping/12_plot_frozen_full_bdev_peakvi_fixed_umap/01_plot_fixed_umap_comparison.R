#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 01_plot_fixed_umap_comparison.R --config <config.yaml>")
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]),
  mustWork = TRUE)

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(digest)
  library(scworkbench)
  library(slingshot)
  library(ggplot2)
  library(ggrepel)
})

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
cfg <- yaml::read_yaml(cfg_path)
.libPaths(unique(c(unlist(cfg$runtime$r_libraries), .libPaths())))
set.seed(as.integer(cfg$runtime$seed))

root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_plot_fixed_umap_comparison.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

input_paths <- c(
  sif = cfg$runtime$sif,
  bdev_manifest = file.path(cfg$inputs$bdev_bundle$path, "manifest.json"),
  l4_model = cfg$inputs$authenticated_l4$model_path,
  l4_pseudotime = cfg$inputs$authenticated_l4$pseudotime_path,
  mappings = cfg$inputs$frozen_benchmark$mappings_path,
  neighbors = cfg$inputs$frozen_benchmark$neighbors_path,
  metrics = cfg$inputs$frozen_benchmark$metrics_path,
  benchmark_summary = cfg$inputs$frozen_benchmark$summary_path,
  benchmark_receipt = cfg$inputs$frozen_benchmark$receipt_path)
expected_sha <- c(
  sif = cfg$runtime$sif_sha256,
  bdev_manifest = cfg$inputs$bdev_bundle$manifest_sha256,
  l4_model = cfg$inputs$authenticated_l4$model_sha256,
  l4_pseudotime = cfg$inputs$authenticated_l4$pseudotime_sha256,
  mappings = cfg$inputs$frozen_benchmark$mappings_sha256,
  neighbors = cfg$inputs$frozen_benchmark$neighbors_sha256,
  metrics = cfg$inputs$frozen_benchmark$metrics_sha256,
  benchmark_summary = cfg$inputs$frozen_benchmark$summary_sha256,
  benchmark_receipt = cfg$inputs$frozen_benchmark$receipt_sha256)
observed_sha <- vapply(input_paths, sha256, character(1L))
stopifnot(identical(unname(observed_sha), unname(expected_sha)),
  cfg$inputs$authenticated_l4$orientation == "native",
  cfg$visualization$trajectory_output_status == "diagnostic_unaccepted")

benchmark_summary <- jsonlite::read_json(cfg$inputs$frozen_benchmark$summary_path,
  simplifyVector = TRUE)
benchmark_receipt <- jsonlite::read_json(cfg$inputs$frozen_benchmark$receipt_path,
  simplifyVector = TRUE)
stopifnot(benchmark_summary$status ==
    "complete_external_comparison_no_automatic_promotion",
  benchmark_receipt$status == "complete",
  isTRUE(benchmark_summary$projection_without_refit),
  isTRUE(benchmark_summary$model_payloads_unchanged))

mappings <- read.delim(cfg$inputs$frozen_benchmark$mappings_path,
  check.names = FALSE, quote = "\"", na.strings = "NA")
neighbors <- read.delim(gzfile(cfg$inputs$frozen_benchmark$neighbors_path),
  check.names = FALSE, quote = "\"", na.strings = "NA")
metrics <- read.delim(cfg$inputs$frozen_benchmark$metrics_path,
  check.names = FALSE, quote = "\"", na.strings = "NA")
candidate_order <- unlist(cfg$visualization$candidate_order, use.names = FALSE)
candidate_labels <- named_character(cfg$visualization$candidate_labels)
stopifnot(identical(sort(unique(mappings$candidate)), sort(candidate_order)),
  nrow(mappings) == cfg$inputs$frozen_benchmark$expected_candidates *
    cfg$inputs$frozen_benchmark$expected_libraries_per_candidate,
  nrow(neighbors) == cfg$inputs$frozen_benchmark$expected_candidates *
    cfg$inputs$frozen_benchmark$expected_libraries_per_candidate *
    cfg$inputs$frozen_benchmark$expected_neighbors_per_library,
  all(table(mappings$candidate) ==
    cfg$inputs$frozen_benchmark$expected_libraries_per_candidate),
  all(table(interaction(neighbors$candidate, neighbors$gsm, drop = TRUE)) ==
    cfg$inputs$frozen_benchmark$expected_neighbors_per_library),
  identical(sort(unique(metrics$candidate)), sort(candidate_order)))

bdev <- wb_open_bundle(cfg$inputs$bdev_bundle$path, validate = FALSE,
  read_only = TRUE)
atac_umap <- space_2d(wb_space(bdev, cfg$inputs$bdev_bundle$fixed_space))
wb_close_bundle(bdev)
stopifnot(nrow(atac_umap) == cfg$inputs$bdev_bundle$expected_cells,
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
curve_arrow <- data.frame(
  x = curve_arrow_points$UMAP_1[[1L]],
  y = curve_arrow_points$UMAP_2[[1L]],
  xend = curve_arrow_points$UMAP_1[[2L]],
  yend = curve_arrow_points$UMAP_2[[2L]])
curve_ends <- rbind(
  transform(l4_curve[1L, ], endpoint = "native L4 start"),
  transform(l4_curve[nrow(l4_curve), ], endpoint = "native L4 end"))

neighbors$UMAP_1 <- atac_umap[neighbors$reference_cell_id, "UMAP_1"]
neighbors$UMAP_2 <- atac_umap[neighbors$reference_cell_id, "UMAP_2"]
group_key <- interaction(neighbors$candidate, neighbors$gsm, drop = TRUE,
  lex.order = TRUE)
neighbor_groups <- split(neighbors, group_key)
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
  unit <- cbind(cos(theta), sin(theta))
  ellipse_xy <- sweep(unit %*% diag(radii, nrow = 2L) %*% t(eig$vectors),
    2L, center, "+")
  spread <- sqrt(sum(weight * rowSums(centered^2)) / sum(weight))
  placement <- data.frame(candidate = rows$candidate[[1L]],
    gsm = rows$gsm[[1L]], UMAP_1 = center[[1L]], UMAP_2 = center[[2L]],
    neighbor_spread = spread, nearest_distance = min(rows$distance),
    kth_distance = kth, stringsAsFactors = FALSE)
  ellipse <- data.frame(candidate = rows$candidate[[1L]],
    gsm = rows$gsm[[1L]], ellipse_order = seq_along(theta),
    UMAP_1 = ellipse_xy[, 1L], UMAP_2 = ellipse_xy[, 2L],
    stringsAsFactors = FALSE)
  list(placement = placement, ellipse = ellipse,
    weights = data.frame(candidate = rows$candidate,
      gsm = rows$gsm, rank = rows$rank,
      reference_cell_id = rows$reference_cell_id, kernel_weight = weight,
      stringsAsFactors = FALSE))
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
placements$candidate_label <- factor(candidate_labels[as.character(placements$candidate)],
  levels = candidate_labels[candidate_order])
ellipses$candidate_label <- factor(candidate_labels[as.character(ellipses$candidate)],
  levels = candidate_labels[candidate_order])
placements$detail_label <- paste0(placements$population, " · ", placements$gsm)

spoke_rows <- neighbors[neighbors$rank <=
    as.integer(cfg$visualization$detailed_spokes_per_library), ]
spokes <- merge(spoke_rows, placements[, c("candidate", "gsm", "population",
  "UMAP_1", "UMAP_2")], by = c("candidate", "gsm"), suffixes = c("_neighbor", "_query"),
  sort = FALSE)

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
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, color = "#151515", linewidth = 0.9,
      arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed"))
}

overview <- base_map() +
  geom_path(data = ellipses,
    aes(UMAP_1, UMAP_2, group = interaction(candidate, gsm), color = population),
    linewidth = 0.42, alpha = 0.33) +
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
  labs(title = "Frozen ImmGen mappings on one fixed Bdev ATAC UMAP",
    subtitle = "Points are 30-neighbor kernel barycenters; ellipses show weighted reference-neighbor spread; black arrow is authenticated native L4",
    x = "ATAC UMAP 1", y = "ATAC UMAP 2",
    caption = "Same 26,723-cell umap.atac coordinates and limits in every panel. Models and neighbors are frozen; no refit or feature tuning.") +
  theme_classic(base_size = 10.5) +
  theme(strip.text = element_text(face = "bold", size = 10.5),
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "bottom", legend.box = "vertical",
    legend.key.width = grid::unit(0.42, "cm"),
    legend.text = element_text(size = 7.5),
    plot.margin = margin(8, 20, 8, 8)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
      override.aes = list(size = 3, shape = 21, color = "#111111")),
    shape = guide_legend(nrow = 1))

ggsave(file.path(dirs[["figures"]], "01_four_model_fixed_atac_umap_overview.png"),
  overview, width = 15, height = 11, dpi = 240, bg = "white")
ggsave(file.path(dirs[["figures"]], "01_four_model_fixed_atac_umap_overview.pdf"),
  overview, width = 15, height = 11, device = cairo_pdf, bg = "white")

detail_files <- character()
for (i in seq_along(candidate_order)) {
  candidate <- candidate_order[[i]]
  points_i <- placements[as.character(placements$candidate) == candidate, ]
  ellipses_i <- ellipses[as.character(ellipses$candidate) == candidate, ]
  spokes_i <- spokes[as.character(spokes$candidate) == candidate, ]
  detail <- base_map() +
    geom_segment(data = spokes_i,
      aes(x = UMAP_1_query, y = UMAP_2_query,
        xend = UMAP_1_neighbor, yend = UMAP_2_neighbor, color = population),
      linewidth = 0.25, alpha = 0.18) +
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
    labs(title = paste0(candidate_labels[[candidate]], " on the fixed Bdev ATAC UMAP"),
      subtitle = "All 28 libraries labeled; ellipses use all 30 frozen neighbors and faint spokes show the nearest five",
      x = "ATAC UMAP 1", y = "ATAC UMAP 2",
      caption = "Black arrow: authenticated fitted Slingshot L4 in native orientation. Placement is a visualization of frozen neighbor ledgers.") +
    theme_classic(base_size = 11.5) +
    theme(plot.title = element_text(face = "bold", size = 15),
      plot.margin = margin(10, 45, 10, 10))
  stem <- sprintf("%02d_%s_fixed_atac_umap_detail", i + 1L, candidate)
  png_path <- file.path(dirs[["figures"]], paste0(stem, ".png"))
  pdf_path <- file.path(dirs[["figures"]], paste0(stem, ".pdf"))
  ggsave(png_path, detail, width = 14, height = 10.5, dpi = 240, bg = "white")
  ggsave(pdf_path, detail, width = 14, height = 10.5, device = cairo_pdf,
    bg = "white")
  detail_files <- c(detail_files, basename(png_path), basename(pdf_path))
}

trajectory <- mappings
trajectory$candidate <- factor(trajectory$candidate, levels = candidate_order)
trajectory$candidate_label <- factor(candidate_labels[as.character(trajectory$candidate)],
  levels = candidate_labels[candidate_order])
trajectory$population <- factor(trajectory$population, levels = rev(population_order))
l4_range <- range(l4_cells$native_l4, finite = TRUE)
na_position <- l4_range[[1L]] - 0.08 * diff(l4_range)
trajectory$l4_plot_value <- ifelse(is.finite(trajectory$provisional_native_l4),
  trajectory$provisional_native_l4, na_position)
trajectory$l4_call <- ifelse(is.finite(trajectory$provisional_native_l4),
  "finite diagnostic value", "no finite value")
finite_medians <- aggregate(provisional_native_l4 ~ candidate + population,
  data = trajectory[is.finite(trajectory$provisional_native_l4), ], FUN = median)
finite_medians$candidate <- factor(finite_medians$candidate, levels = candidate_order)
finite_medians$candidate_label <- factor(candidate_labels[as.character(finite_medians$candidate)],
  levels = candidate_labels[candidate_order])

trajectory_plot <- ggplot(trajectory,
    aes(l4_plot_value, population, fill = lineage, shape = l4_call)) +
  geom_vline(xintercept = l4_range[[1L]], color = "#9A9A9A", linetype = 3,
    linewidth = 0.4) +
  geom_point(position = position_jitter(height = 0.12, width = 0),
    size = 2.4, stroke = 0.55, color = "#222222", alpha = 0.88) +
  geom_point(data = finite_medians,
    aes(provisional_native_l4, population), inherit.aes = FALSE,
    shape = 23, size = 3.0, stroke = 0.7, fill = "#F4D35E", color = "#111111") +
  facet_wrap(~candidate_label, ncol = 2L) +
  scale_shape_manual(values = c("finite diagnostic value" = 21,
    "no finite value" = 4), name = "Post-freeze L4 output") +
  scale_fill_manual(values = c("B" = "#D1495B", "Stem&Prog" = "#3B4CC0"),
    name = "ImmGen lineage") +
  scale_x_continuous(breaks = c(na_position, pretty(l4_range)),
    labels = c("NA", format(pretty(l4_range), trim = TRUE))) +
  labs(title = "Diagnostic native-L4 positions from the four frozen models",
    subtitle = "Each point is one ImmGen library; yellow diamonds are population medians among finite values",
    x = "Provisional native-L4 pseudotime (NA plotted left of the dashed boundary)",
    y = "ImmGen population",
    caption = "Diagnostic/unaccepted trajectory output. Native orientation is preserved; these values are not final trajectory assignments.") +
  theme_classic(base_size = 10.5) +
  theme(strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "bottom", panel.spacing = grid::unit(1.0, "lines"))

ggsave(file.path(dirs[["figures"]], "06_four_model_diagnostic_native_l4.png"),
  trajectory_plot, width = 14, height = 10, dpi = 240, bg = "white")
ggsave(file.path(dirs[["figures"]], "06_four_model_diagnostic_native_l4.pdf"),
  trajectory_plot, width = 14, height = 10, device = cairo_pdf, bg = "white")

placements_out <- placements
placements_out$candidate <- as.character(placements_out$candidate)
placements_out$candidate_label <- as.character(placements_out$candidate_label)
placements_out$population <- as.character(placements_out$population)
ellipses_out <- ellipses
ellipses_out$candidate <- as.character(ellipses_out$candidate)
ellipses_out$candidate_label <- as.character(ellipses_out$candidate_label)
ellipses_out$population <- as.character(ellipses_out$population)
write_tsv(placements_out, file.path(dirs[["tables"]],
  "01_fixed_umap_candidate_placements.tsv"))
write_tsv(ellipses_out, file.path(dirs[["tables"]],
  "02_fixed_umap_neighbor_ellipses.tsv.gz"), compress = TRUE)
write_tsv(l4_curve, file.path(dirs[["tables"]], "03_authenticated_native_l4_curve.tsv"))
write_tsv(kernel_weights, file.path(dirs[["tables"]],
  "04_neighbor_kernel_weights.tsv.gz"), compress = TRUE)

input_sha_after <- vapply(input_paths, sha256, character(1L))
stopifnot(identical(observed_sha, input_sha_after),
  nrow(placements) == nrow(mappings),
  all(table(placements$candidate) ==
    cfg$inputs$frozen_benchmark$expected_libraries_per_candidate),
  identical(levels(placements$candidate), candidate_order))

summary <- list(
  schema_version = "bdev.frozen_peakvi_fixed_umap_comparison_summary.v1",
  status = "complete_read_only_visualization",
  run_identity = cfg$run_identity,
  fixed_space = cfg$inputs$bdev_bundle$fixed_space,
  fixed_reference_cells = nrow(atac_umap),
  authenticated_l4_cells = nrow(l4_cells),
  authenticated_l4_curve_points = nrow(l4_curve),
  native_l4_orientation = "native",
  candidates = candidate_order,
  libraries_per_candidate = as.list(table(placements$candidate)),
  neighbors_per_library = cfg$inputs$frozen_benchmark$expected_neighbors_per_library,
  placement_kernel = cfg$visualization$placement_kernel,
  ellipse_scale = ellipse_scale,
  trajectory_output_status = "diagnostic_unaccepted",
  model_refit = FALSE,
  benchmark_mutated = FALSE,
  source_inputs_unchanged = TRUE,
  automatic_model_promotion = FALSE,
  implementation_unlocked = FALSE,
  sif_built = FALSE)
write_json(summary, file.path(dirs[["reports"]],
  "01_fixed_umap_comparison_summary.json"))
receipt <- c(summary, list(config_path = cfg_path,
  config_sha256 = sha256(cfg_path), input_sha256 = as.list(input_sha_after),
  workflow_path = script_path, workflow_sha256 = sha256(script_path),
  generated_detail_figures = as.list(detail_files)))
write_json(receipt, file.path(dirs[["manifests"]], "01_execution_receipt.json"))
cat("Completed fixed-UMAP and native-L4 figures for", length(candidate_order),
  "frozen models and", nrow(mappings), "candidate-library mappings\n")
