#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 02_build_all_cells_fixed_umap.R --config <config.yaml>",
    call. = FALSE)

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
invisible(lapply(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(root, "04_logs", "01_build_all_cells_fixed_umap.log"),
  split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

log_event("Validating pinned all-cell visualization inputs")
stopifnot(
  sha256(cfg$runtime$sif) == cfg$runtime$sif_sha256,
  sha256(cfg$inputs$parent_placement_table$path) ==
    cfg$inputs$parent_placement_table$sha256,
  sha256(cfg$inputs$parent_inventory$path) ==
    cfg$inputs$parent_inventory$sha256,
  sha256(file.path(cfg$inputs$fixed_bdev_bundle$path, "manifest.json")) ==
    cfg$inputs$fixed_bdev_bundle$manifest_sha256,
  sha256(cfg$inputs$authenticated_l4$model_path) ==
    cfg$inputs$authenticated_l4$model_sha256,
  sha256(cfg$inputs$authenticated_l4$pseudotime_path) ==
    cfg$inputs$authenticated_l4$pseudotime_sha256)

placement <- as.data.table(read.delim(
  gzfile(cfg$inputs$parent_placement_table$path), check.names = FALSE,
  quote = "\"", na.strings = "NA"))
stopifnot(
  nrow(placement) == cfg$scope$expected_query_cells,
  uniqueN(placement$slot_id) == nrow(placement),
  sum(placement$accepted) == cfg$scope$expected_accepted_cells,
  sum(!placement$accepted) == cfg$scope$expected_rejected_cells,
  all(placement$accepted == is.finite(placement$accepted_native_l4)),
  all(is.finite(placement$fixed_atac_UMAP_1)),
  all(is.finite(placement$fixed_atac_UMAP_2)))

labels <- c(FLT3LKO_BM_WT = "Flt3l WT", FLT3LKO_BM_KO = "Flt3l KO",
  IL7KO_BM_WT = "Il7 WT", IL7KO_BM_KO = "Il7 KO")
label_levels <- unname(labels)
placement[, library_label := factor(labels[library_id], levels = label_levels)]

log_event("Reading unchanged Bdev ATAC UMAP and authenticated native L4")
fixed <- wb_open_bundle(cfg$inputs$fixed_bdev_bundle$path,
  validate = FALSE, read_only = TRUE)
reference_cells <- as.character(wb_cells(fixed, "all"))
reference_matrix <- space_matrix(wb_space(fixed,
  cfg$inputs$fixed_bdev_bundle$fixed_space, cells = reference_cells))
wb_close_bundle(fixed)
stopifnot(length(reference_cells) == cfg$scope$expected_reference_cells,
  identical(rownames(reference_matrix), reference_cells),
  ncol(reference_matrix) == 2L)
reference_umap <- data.table(cell_id = reference_cells,
  UMAP_1 = reference_matrix[, 1L], UMAP_2 = reference_matrix[, 2L])

family <- readRDS(cfg$inputs$authenticated_l4$pseudotime_path)
take <- as.character(family$key) == cfg$inputs$authenticated_l4$key &
  is.finite(as.numeric(family$value)) & family$cell_id %in% reference_cells
trajectory <- data.table(cell_id = as.character(family$cell_id[take]),
  native_l4 = as.numeric(family$value[take]))
trajectory <- merge(trajectory, reference_umap, by = "cell_id", sort = FALSE)

slingshot_model <- readRDS(cfg$inputs$authenticated_l4$model_path)
curve_value <- slingCurves(slingshot_model)[[cfg$inputs$authenticated_l4$curve_name]]
stopifnot(!is.null(curve_value))
curve_xy <- curve_value$s[curve_value$ord, , drop = FALSE]
curve <- data.table(UMAP_1 = curve_xy[, 1L], UMAP_2 = curve_xy[, 2L])

panel_levels <- c("Bdev reference L4", label_levels)
reference_background <- rbindlist(lapply(panel_levels, function(panel) {
  value <- copy(reference_umap)
  value[, panel_label := factor(panel, levels = panel_levels)]
  value
}))
curve_plot <- rbindlist(lapply(panel_levels, function(panel) {
  value <- copy(curve)
  value[, panel_label := factor(panel, levels = panel_levels)]
  value
}))
trajectory[, panel_label := factor("Bdev reference L4", levels = panel_levels)]
query_plot <- placement[, .(
  UMAP_1 = fixed_atac_UMAP_1, UMAP_2 = fixed_atac_UMAP_2,
  accepted, native_l4 = accepted_native_l4,
  panel_label = factor(as.character(library_label), levels = panel_levels))]
rejected_plot <- query_plot[accepted == FALSE]
accepted_plot <- query_plot[accepted == TRUE]
l4_limits <- range(trajectory$native_l4, finite = TRUE)

log_event("Drawing all query cells on the fixed Bdev ATAC UMAP")
p <- ggplot() +
  geom_point(data = reference_background, aes(UMAP_1, UMAP_2),
    color = cfg$visualization$reference_background_color,
    size = 0.11, alpha = 0.25) +
  geom_point(data = rejected_plot, aes(UMAP_1, UMAP_2),
    color = cfg$visualization$rejected_query_color,
    size = 0.32, alpha = 0.38) +
  geom_path(data = curve_plot, aes(UMAP_1, UMAP_2), color = "#111111",
    linewidth = 0.5, arrow = grid::arrow(type = "closed",
      length = grid::unit(0.08, "inches"))) +
  geom_point(data = trajectory, aes(UMAP_1, UMAP_2, color = native_l4),
    size = 0.28, alpha = 0.70) +
  geom_point(data = accepted_plot, aes(UMAP_1, UMAP_2, color = native_l4),
    size = 0.50, alpha = 0.80) +
  scale_color_viridis_c(option = cfg$visualization$color_scale,
    limits = l4_limits, name = "Native L4 PT") +
  facet_wrap(~panel_label, ncol = 3) + coord_equal() +
  labs(
    title = "All gated-20D query mappings on the fixed Bdev ATAC UMAP",
    subtitle = "Dark grey: 20,278 cells failing the complete L4 acceptance gate; color: 2,974 accepted cells",
    x = "Bdev ATAC UMAP 1", y = "Bdev ATAC UMAP 2",
    caption = "All query cells use the same 30-neighbor fixed-UMAP placement; Bdev UMAP and native Slingshot L4 were not refit.") +
  theme_classic(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"), legend.position = "right")
ggsave(file.path(root, "02_figures",
  "01_all_cells_fixed_bdev_atac_umap.png"), p,
  width = 15.2, height = 9.6, dpi = 280, bg = "white")
ggsave(file.path(root, "02_figures",
  "01_all_cells_fixed_bdev_atac_umap.pdf"), p,
  width = 15.2, height = 9.6, device = cairo_pdf, bg = "white")

summary_table <- placement[, .(
  total_cells = .N, accepted_cells = sum(accepted),
  rejected_cells = sum(!accepted), accepted_percent = 100 * mean(accepted)),
  by = .(library_id, knockout_target, genotype, library_label)]
fwrite(summary_table, file.path(root, "01_tables",
  "01_all_cells_fixed_umap_summary.tsv"), sep = "\t", quote = TRUE, na = "NA")

source_hashes_after <- list(
  parent_placement_table = sha256(cfg$inputs$parent_placement_table$path),
  parent_inventory = sha256(cfg$inputs$parent_inventory$path),
  fixed_bdev_bundle_manifest = sha256(file.path(cfg$inputs$fixed_bdev_bundle$path,
    "manifest.json")))
stopifnot(
  source_hashes_after$parent_placement_table ==
    cfg$inputs$parent_placement_table$sha256,
  source_hashes_after$parent_inventory == cfg$inputs$parent_inventory$sha256,
  source_hashes_after$fixed_bdev_bundle_manifest ==
    cfg$inputs$fixed_bdev_bundle$manifest_sha256)

report <- list(
  schema_version = "bdev.flt3lko_il7ko_gated_20d_all_cells_fixed_umap_report.v1",
  status = "complete_passed", run_identity = cfg$run_identity,
  reference_cells = length(reference_cells), reference_l4_cells = nrow(trajectory),
  query_cells = nrow(placement), accepted_query_cells = sum(placement$accepted),
  rejected_query_cells = sum(!placement$accepted),
  rejected_color = cfg$visualization$rejected_query_color,
  accepted_color_field = "accepted_native_l4",
  library_summary = lapply(seq_len(nrow(summary_table)), function(index)
    as.list(summary_table[index])),
  fixed_atac_placement = "inherited_30_neighbor_distance_kernel_barycenter",
  native_l4_orientation = cfg$inputs$authenticated_l4$orientation,
  model_refit = FALSE, umap_refit = FALSE, slingshot_refit = FALSE,
  genotype_contrast = FALSE, paper_analysis = FALSE,
  implementation_unlocked = FALSE, sif_built = FALSE,
  source_manifests_unchanged = TRUE)
write_json(report, file.path(root, "03_reports",
  "01_all_cells_fixed_umap_visualization.json"))
write_json(c(report, list(config_path = cfg_path,
  config_sha256 = sha256(cfg_path), sif_sha256 = sha256(cfg$runtime$sif),
  scworkbench_version = as.character(packageVersion("scworkbench")),
  source_hashes_after = source_hashes_after)),
  file.path(root, "06_manifests", "01_execution_receipt.json"))

log_event("Completed all-cell fixed-UMAP visualization")
inventory_paths <- list.files(root, recursive = TRUE, full.names = TRUE)
inventory_paths <- inventory_paths[file.info(inventory_paths)$isdir == FALSE]
inventory_paths <- inventory_paths[basename(inventory_paths) !=
  "02_run_inventory.json"]
inventory_paths <- sort(inventory_paths)
inventory <- list(
  schema_version = "bdev.flt3lko_il7ko_gated_20d_all_cells_fixed_umap_inventory.v1",
  status = "complete_passed", run_identity = cfg$run_identity,
  files = lapply(inventory_paths, function(path) list(
    path = substring(path, nchar(root) + 2L),
    bytes = unname(file.info(path)$size), sha256 = sha256(path))),
  query_cells = nrow(placement), accepted_query_cells = sum(placement$accepted),
  rejected_query_cells = sum(!placement$accepted),
  genotype_contrast_performed = FALSE, implementation_unlocked = FALSE,
  sif_built = FALSE)
write_json(inventory, file.path(root, "06_manifests", "02_run_inventory.json"))
