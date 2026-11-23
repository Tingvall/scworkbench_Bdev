#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 01_build_mapping_figures.R --config <config.yaml>", call. = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(slingshot)
  library(yaml)
  library(jsonlite)
  library(digest)
})

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
space_matrix <- function(x) {
  if (is.data.frame(x) && "cell_id" %in% names(x)) {
    ids <- as.character(x$cell_id)
    x <- as.matrix(x[, setdiff(names(x), "cell_id"), drop = FALSE])
    rownames(x) <- ids
  } else x <- as.matrix(x)
  x
}
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), ..., "\n")
  flush.console()
}

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
root <- cfg$storage$results_root
dirs <- file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
log_path <- file.path(root, "04_logs", "01_build_mapping_figures.log")
sink(log_path, split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

log_event("Validating pinned visualization inputs")
stopifnot(
  sha256(cfg$runtime$sif) == cfg$runtime$sif_sha256,
  sha256(cfg$inputs$mapping_table$path) == cfg$inputs$mapping_table$sha256,
  sha256(cfg$inputs$mapping_inventory$path) == cfg$inputs$mapping_inventory$sha256,
  sha256(file.path(cfg$inputs$mapping_bundle$path, "manifest.json")) ==
    cfg$inputs$mapping_bundle$manifest_sha256,
  sha256(file.path(cfg$inputs$fixed_bdev_bundle$path, "manifest.json")) ==
    cfg$inputs$fixed_bdev_bundle$manifest_sha256,
  sha256(cfg$inputs$authenticated_l4$model_path) ==
    cfg$inputs$authenticated_l4$model_sha256,
  sha256(cfg$inputs$authenticated_l4$pseudotime_path) ==
    cfg$inputs$authenticated_l4$pseudotime_sha256,
  sha256(cfg$inputs$query_atac_umap$path) == cfg$inputs$query_atac_umap$sha256,
  sha256(cfg$inputs$query_gex_umap$path) == cfg$inputs$query_gex_umap$sha256)

mapping <- as.data.table(read.delim(gzfile(cfg$inputs$mapping_table$path),
  check.names = FALSE, quote = "\"", na.strings = "NA"))
stopifnot(
  nrow(mapping) == cfg$scope$expected_query_cells,
  uniqueN(mapping$slot_id) == nrow(mapping),
  uniqueN(mapping$cell_id) == nrow(mapping),
  all(c("accepted", "accepted_native_l4") %in% names(mapping)),
  all(mapping$accepted == is.finite(mapping$accepted_native_l4)),
  sum(mapping$accepted) == cfg$scope$expected_accepted_cells)

labels <- c(FLT3LKO_BM_WT = "Flt3l WT", FLT3LKO_BM_KO = "Flt3l KO",
  IL7KO_BM_WT = "Il7 WT", IL7KO_BM_KO = "Il7 KO")
label_levels <- unname(labels)
stopifnot(setequal(unique(mapping$library_id), names(labels)))
mapping[, library_label := factor(labels[library_id], levels = label_levels)]

checkpoint_paths <- sort(list.files(cfg$inputs$bounded_checkpoints$path,
  pattern = "^batch_[0-9]{3}\\.rds$", full.names = TRUE))
stopifnot(length(checkpoint_paths) == cfg$inputs$bounded_checkpoints$expected_files)
log_event("Reading", length(checkpoint_paths), "bounded neighbor checkpoints")
checkpoint_values <- lapply(checkpoint_paths, readRDS)
stopifnot(identical(vapply(checkpoint_values, function(x) x$batch_index,
  integer(1L)), seq_along(checkpoint_values)))
neighbors <- rbindlist(lapply(checkpoint_values, function(x) x$neighbors),
  use.names = TRUE)
rm(checkpoint_values); invisible(gc())
stopifnot(
  nrow(neighbors) == nrow(mapping) * cfg$inputs$bounded_checkpoints$primary_k,
  all(neighbors$neighbor_rank >= 1L &
    neighbors$neighbor_rank <= cfg$inputs$bounded_checkpoints$primary_k),
  setequal(unique(neighbors$query_cell_id), mapping$slot_id))

log_event("Reading unchanged Bdev ATAC UMAP and authenticated native L4")
fixed <- wb_open_bundle(cfg$inputs$fixed_bdev_bundle$path,
  validate = FALSE, read_only = TRUE)
reference_cells <- as.character(wb_cells(fixed, "all"))
reference_umap_matrix <- space_matrix(wb_space(fixed,
  cfg$inputs$fixed_bdev_bundle$fixed_space, cells = reference_cells))
wb_close_bundle(fixed)
stopifnot(length(reference_cells) == cfg$scope$expected_reference_cells,
  identical(rownames(reference_umap_matrix), reference_cells),
  ncol(reference_umap_matrix) == 2L)
reference_umap <- data.table(cell_id = reference_cells,
  UMAP_1 = reference_umap_matrix[, 1L], UMAP_2 = reference_umap_matrix[, 2L])

family <- readRDS(cfg$inputs$authenticated_l4$pseudotime_path)
take <- as.character(family$key) == cfg$inputs$authenticated_l4$key &
  is.finite(as.numeric(family$value)) & family$cell_id %in% reference_cells
trajectory <- data.table(cell_id = as.character(family$cell_id[take]),
  native_l4 = as.numeric(family$value[take]))
trajectory <- merge(trajectory, reference_umap, by = "cell_id", sort = FALSE)
stopifnot(nrow(trajectory) > 0L, uniqueN(trajectory$cell_id) == nrow(trajectory))

slingshot_model <- readRDS(cfg$inputs$authenticated_l4$model_path)
curve_value <- slingCurves(slingshot_model)[[cfg$inputs$authenticated_l4$curve_name]]
stopifnot(!is.null(curve_value))
curve_xy <- curve_value$s[curve_value$ord, , drop = FALSE]
curve <- data.table(UMAP_1 = curve_xy[, 1L], UMAP_2 = curve_xy[, 2L])

log_event("Placing all queries on the unchanged Bdev ATAC UMAP")
k <- as.integer(cfg$inputs$bounded_checkpoints$primary_k)
kth <- neighbors[neighbor_rank == k,
  .(query_cell_id, kth_distance = distance)]
stopifnot(nrow(kth) == nrow(mapping), uniqueN(kth$query_cell_id) == nrow(kth))
neighbors <- merge(neighbors, kth, by = "query_cell_id", sort = FALSE)
neighbors <- merge(neighbors, reference_umap, by = "cell_id", sort = FALSE)
stopifnot(nrow(neighbors) == nrow(mapping) * k)
neighbors[, weight := fifelse(kth_distance == 0,
  as.numeric(distance == 0), exp(-(distance / kth_distance)^2))]
fixed_placement <- neighbors[, {
  xbar <- sum(weight * UMAP_1) / sum(weight)
  ybar <- sum(weight * UMAP_2) / sum(weight)
  list(fixed_atac_UMAP_1 = xbar, fixed_atac_UMAP_2 = ybar,
    neighbor_spread = sqrt(sum(weight * ((UMAP_1 - xbar)^2 +
      (UMAP_2 - ybar)^2)) / sum(weight)),
    nearest_distance = min(distance), kth_distance = max(kth_distance))
}, by = query_cell_id]
rm(neighbors); invisible(gc())

read_native_umap <- function(path, prefix) {
  value <- fread(path)
  stopifnot(identical(names(value), c("Barcode", "UMAP-1", "UMAP-2")),
    uniqueN(value$Barcode) == nrow(value))
  setnames(value, c("Barcode", "UMAP-1", "UMAP-2"),
    c("cell_id", paste0(prefix, "_UMAP_1"), paste0(prefix, "_UMAP_2")))
  value
}
query_atac <- read_native_umap(cfg$inputs$query_atac_umap$path, "query_atac")
query_gex <- read_native_umap(cfg$inputs$query_gex_umap$path, "query_gex")

placement <- merge(mapping, fixed_placement, by.x = "slot_id",
  by.y = "query_cell_id", sort = FALSE)
placement <- merge(placement, query_atac, by = "cell_id", sort = FALSE)
placement <- merge(placement, query_gex, by = "cell_id", sort = FALSE)
placement <- placement[match(mapping$slot_id, slot_id)]
stopifnot(nrow(placement) == nrow(mapping),
  identical(as.character(placement$slot_id), as.character(mapping$slot_id)),
  all(is.finite(placement$fixed_atac_UMAP_1)),
  all(is.finite(placement$query_atac_UMAP_1)),
  all(is.finite(placement$query_gex_UMAP_1)))
fwrite(placement, file.path(root, "01_tables",
  "01_gated_20d_umap_placements.tsv.gz"), sep = "\t", quote = TRUE, na = "NA")

summary_table <- placement[, .(
  total_cells = .N,
  accepted_cells = sum(accepted),
  accepted_percent = 100 * mean(accepted),
  accepted_l4_q25 = quantile(accepted_native_l4[accepted], 0.25, na.rm = TRUE),
  accepted_l4_median = median(accepted_native_l4[accepted], na.rm = TRUE),
  accepted_l4_q75 = quantile(accepted_native_l4[accepted], 0.75, na.rm = TRUE)),
  by = .(library_id, knockout_target, genotype, library_label)]
fwrite(summary_table, file.path(root, "01_tables",
  "02_accepted_native_l4_summary.tsv"), sep = "\t", quote = TRUE, na = "NA")

l4_limits <- range(trajectory$native_l4, finite = TRUE)
background_color <- cfg$visualization$background_color
base_theme <- theme_classic(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"), legend.position = "right")

log_event("Drawing authenticated reference and accepted fixed-UMAP placements")
panel_levels <- c("Bdev reference L4", label_levels)
fixed_background <- rbindlist(lapply(panel_levels, function(panel) {
  ans <- copy(reference_umap)
  ans[, panel_label := panel]
  ans
}))
fixed_background[, panel_label := factor(panel_label, levels = panel_levels)]
trajectory_plot <- copy(trajectory)
trajectory_plot[, panel_label := factor("Bdev reference L4", levels = panel_levels)]
query_fixed_plot <- placement[accepted == TRUE,
  .(UMAP_1 = fixed_atac_UMAP_1, UMAP_2 = fixed_atac_UMAP_2,
    native_l4 = accepted_native_l4,
    panel_label = factor(as.character(library_label), levels = panel_levels))]
curve_plot <- rbindlist(lapply(panel_levels, function(panel) {
  ans <- copy(curve)
  ans[, panel_label := factor(panel, levels = panel_levels)]
  ans
}))

p_fixed <- ggplot() +
  geom_point(data = fixed_background, aes(UMAP_1, UMAP_2),
    color = background_color, size = 0.11, alpha = 0.28) +
  geom_path(data = curve_plot, aes(UMAP_1, UMAP_2), color = "#111111",
    linewidth = 0.5, arrow = grid::arrow(type = "closed",
      length = grid::unit(0.08, "inches"))) +
  geom_point(data = trajectory_plot,
    aes(UMAP_1, UMAP_2, color = native_l4), size = 0.28, alpha = 0.70) +
  geom_point(data = query_fixed_plot,
    aes(UMAP_1, UMAP_2, color = native_l4), size = 0.48, alpha = 0.72) +
  scale_color_viridis_c(option = cfg$visualization$color_scale,
    limits = l4_limits, name = "Native L4 PT") +
  facet_wrap(~panel_label, ncol = 3) + coord_equal() +
  labs(title = "Accepted gated-20D mappings on the fixed Bdev ATAC UMAP",
    subtitle = "Reference L4 at left; query panels show only the 2,974 cells passing the complete acceptance gate",
    x = "Bdev ATAC UMAP 1", y = "Bdev ATAC UMAP 2",
    caption = "30-neighbor distance-kernel placement; Bdev UMAP and native Slingshot L4 were not refit.") +
  base_theme
ggsave(file.path(root, "02_figures",
  "01_accepted_mapping_fixed_bdev_atac_umap.png"), p_fixed,
  width = 15.2, height = 9.6, dpi = 280, bg = "white")
ggsave(file.path(root, "02_figures",
  "01_accepted_mapping_fixed_bdev_atac_umap.pdf"), p_fixed,
  width = 15.2, height = 9.6, device = cairo_pdf, bg = "white")

plot_native_umap <- function(data, x, y, title, subtitle, xlab, ylab,
    file_stem) {
  accepted_data <- data[accepted == TRUE]
  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(color = background_color, size = 0.22, alpha = 0.42) +
    geom_point(data = accepted_data,
      aes(color = accepted_native_l4), size = 0.52, alpha = 0.78) +
    scale_color_viridis_c(option = cfg$visualization$color_scale,
      limits = l4_limits, name = "Accepted native L4 PT") +
    facet_wrap(~library_label, ncol = 2) + coord_equal() +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab,
      caption = "Grey: all QC-selected cells. Color: complete-gate accepted cells only. Source UMAP was not refit.") +
    base_theme
  ggsave(file.path(root, "02_figures", paste0(file_stem, ".png")), p,
    width = 13.4, height = 10.0, dpi = 280, bg = "white")
  ggsave(file.path(root, "02_figures", paste0(file_stem, ".pdf")), p,
    width = 13.4, height = 10.0, device = cairo_pdf, bg = "white")
}

log_event("Drawing accepted pseudotime on source query ATAC and GEX UMAPs")
plot_native_umap(placement, "query_atac_UMAP_1", "query_atac_UMAP_2",
  "Accepted native-L4 mapping on the query ATAC UMAP",
  "Cell Ranger ARC ATAC UMAP; 2,974 complete-gate accepted cells colored by transferred native L4",
  "Query ATAC UMAP 1", "Query ATAC UMAP 2",
  "02_accepted_native_l4_query_atac_umap")
plot_native_umap(placement, "query_gex_UMAP_1", "query_gex_UMAP_2",
  "Accepted native-L4 mapping on the query GEX UMAP",
  "Cell Ranger ARC GEX UMAP; 2,974 complete-gate accepted cells colored by transferred native L4",
  "Query GEX UMAP 1", "Query GEX UMAP 2",
  "03_accepted_native_l4_query_gex_umap")

log_event("Drawing descriptive accepted-pseudotime distributions")
p_distribution <- ggplot(placement[accepted == TRUE],
    aes(accepted_native_l4, fill = library_label)) +
  geom_histogram(bins = 32, boundary = 0, color = "white", linewidth = 0.15,
    alpha = 0.88) +
  facet_wrap(~library_label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3"),
    guide = "none") +
  coord_cartesian(xlim = l4_limits) +
  labs(title = "Native-L4 pseudotime among accepted gated-20D mappings",
    subtitle = "Descriptive distributions only; rejected cells are excluded and no genotype contrast is performed",
    x = "Accepted native L4 PT", y = "Accepted cells") + base_theme
ggsave(file.path(root, "02_figures",
  "04_accepted_native_l4_distributions.png"), p_distribution,
  width = 12.2, height = 8.4, dpi = 280, bg = "white")
ggsave(file.path(root, "02_figures",
  "04_accepted_native_l4_distributions.pdf"), p_distribution,
  width = 12.2, height = 8.4, device = cairo_pdf, bg = "white")

source_hashes_after <- list(
  mapping_table = sha256(cfg$inputs$mapping_table$path),
  mapping_inventory = sha256(cfg$inputs$mapping_inventory$path),
  mapping_bundle_manifest = sha256(file.path(cfg$inputs$mapping_bundle$path,
    "manifest.json")),
  fixed_bdev_bundle_manifest = sha256(file.path(cfg$inputs$fixed_bdev_bundle$path,
    "manifest.json")),
  query_atac_umap = sha256(cfg$inputs$query_atac_umap$path),
  query_gex_umap = sha256(cfg$inputs$query_gex_umap$path))
stopifnot(
  source_hashes_after$mapping_table == cfg$inputs$mapping_table$sha256,
  source_hashes_after$mapping_inventory == cfg$inputs$mapping_inventory$sha256,
  source_hashes_after$mapping_bundle_manifest ==
    cfg$inputs$mapping_bundle$manifest_sha256,
  source_hashes_after$fixed_bdev_bundle_manifest ==
    cfg$inputs$fixed_bdev_bundle$manifest_sha256,
  source_hashes_after$query_atac_umap == cfg$inputs$query_atac_umap$sha256,
  source_hashes_after$query_gex_umap == cfg$inputs$query_gex_umap$sha256)

report <- list(
  schema_version = "bdev.flt3lko_il7ko_gated_20d_mapping_visualization_report.v1",
  status = "complete_passed",
  run_identity = cfg$run_identity,
  reference_cells = length(reference_cells),
  reference_l4_cells = nrow(trajectory),
  query_cells = nrow(placement),
  accepted_query_cells = sum(placement$accepted),
  accepted_percent = 100 * mean(placement$accepted),
  library_summary = lapply(seq_len(nrow(summary_table)), function(index)
    as.list(summary_table[index])),
  acceptance_definition = "accepted == TRUE and finite accepted_native_l4",
  fixed_atac_placement = "30_neighbor_exp_negative_squared_distance_over_kth_distance_barycenter",
  native_query_umaps = c("Cell Ranger ARC ATAC", "Cell Ranger ARC GEX"),
  native_l4_orientation = cfg$inputs$authenticated_l4$orientation,
  model_refit = FALSE, umap_refit = FALSE, slingshot_refit = FALSE,
  genotype_contrast = FALSE, paper_analysis = FALSE,
  implementation_unlocked = FALSE, sif_built = FALSE,
  source_manifests_unchanged = TRUE)
write_json(report, file.path(root, "03_reports",
  "01_gated_20d_mapping_visualization.json"))

receipt <- c(report, list(
  config_path = cfg_path,
  config_sha256 = sha256(cfg_path),
  sif_sha256 = sha256(cfg$runtime$sif),
  scworkbench_version = as.character(packageVersion("scworkbench")),
  source_hashes_after = source_hashes_after))
write_json(receipt, file.path(root, "06_manifests",
  "01_execution_receipt.json"))

log_event("Completed gated 20D mapping visualization")
inventory_paths <- list.files(root, recursive = TRUE, full.names = TRUE)
inventory_paths <- inventory_paths[file.info(inventory_paths)$isdir == FALSE]
inventory_paths <- inventory_paths[basename(inventory_paths) !=
  "02_run_inventory.json"]
inventory_paths <- sort(inventory_paths)
inventory <- list(
  schema_version = "bdev.flt3lko_il7ko_gated_20d_mapping_visualization_inventory.v1",
  status = "complete_passed", run_identity = cfg$run_identity,
  files = lapply(inventory_paths, function(path) list(
    path = substring(path, nchar(root) + 2L),
    bytes = unname(file.info(path)$size), sha256 = sha256(path))),
  figure_pairs = 4L, accepted_query_cells = sum(placement$accepted),
  genotype_contrast_performed = FALSE, implementation_unlocked = FALSE,
  sif_built = FALSE)
write_json(inventory, file.path(root, "06_manifests", "02_run_inventory.json"))
