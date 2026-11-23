#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 02_build_fixed_umap.R --config <config.yaml>", call. = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbRepresentation)
  library(slingshot)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(yaml)
  library(jsonlite)
  library(digest)
})

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
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

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "05_objects", "06_manifests")),
  c("tables", "figures", "reports", "logs", "objects", "manifests"))
invisible(lapply(c(dirs, sessions = cfg$storage$sessions), dir.create,
  recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "02_build_fixed_umap.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(
  sha256(cfg$runtime$sif) == cfg$runtime$sif_sha256,
  sha256(file.path(cfg$inputs$projection_bundle$path, "manifest.json")) ==
    cfg$inputs$projection_bundle$manifest_sha256,
  sha256(file.path(cfg$inputs$fixed_bdev_bundle$path, "manifest.json")) ==
    cfg$inputs$fixed_bdev_bundle$manifest_sha256,
  sha256(cfg$inputs$authenticated_l4$model_path) ==
    cfg$inputs$authenticated_l4$model_sha256,
  sha256(cfg$inputs$authenticated_l4$pseudotime_path) ==
    cfg$inputs$authenticated_l4$pseudotime_sha256,
  sha256(cfg$inputs$slot_map$path) == cfg$inputs$slot_map$sha256,
  sha256(cfg$inputs$mapping_outcomes$path) ==
    cfg$inputs$mapping_outcomes$sha256,
  cfg$inputs$authenticated_l4$orientation == "native",
  as.character(packageVersion("scworkbench")) == cfg$runtime$scworkbench_version,
  as.character(packageVersion("scwbRepresentation")) ==
    cfg$runtime$scwbRepresentation_version)

representation_manifest <- system.file("scworkbench", "operation-pack.json",
  package = "scwbRepresentation")
stopifnot(nzchar(representation_manifest))
try(wb_operation_pack_unload("scwbrepresentation"), silent = TRUE)
wb_operation_pack_activate(representation_manifest, lib.loc = .libPaths())

bundle <- wb_open_bundle(cfg$storage$execution_bundle, validate = FALSE,
  read_only = FALSE)
all_cells <- as.character(wb_cells(bundle, "all"))
query_cells <- all_cells[grepl("^external_query_[0-9]+$", all_cells)]
reference_cells <- setdiff(all_cells, query_cells)
stopifnot(length(query_cells) == cfg$inputs$projection_bundle$expected_query_cells,
  length(reference_cells) ==
    cfg$inputs$projection_bundle$expected_reference_cells)
knn_session <- wb_session_create(bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_knn_retry1"))
knn <- wb_session_run_op(knn_session, bundle,
  "scwbrepresentation::knn_neighbors", inputs = list(
    space = cfg$inputs$projection_bundle$joint_space,
    query_cells = query_cells, candidate_cells = reference_cells,
    k = as.integer(cfg$placement$neighbor_k), include_query = FALSE,
    space_chunk_size = as.integer(cfg$placement$space_chunk_size)),
  notes = "Exact PeakVI neighbors for bounded FLT3LKO/IL7KO fixed-UMAP display")
neighbors <- as.data.table(knn$run$outputs$neighbors$table)
stopifnot(nrow(neighbors) == length(query_cells) * cfg$placement$neighbor_k,
  identical(unique(as.character(neighbors$query_cell_id)), query_cells))
fwrite(neighbors, file.path(dirs[["tables"]],
  "01_query_to_bdev_neighbors.tsv.gz"), sep = "\t", quote = TRUE, na = "NA")
knn$session <- wb_session_resume(bundle, knn$session$durable$workspace)
knn$session <- wb_session_discard(knn$session)
wb_close_bundle(bundle)

bdev <- wb_open_bundle(cfg$inputs$fixed_bdev_bundle$path, validate = FALSE,
  read_only = TRUE)
atac_umap <- space_2d(wb_space(bdev,
  cfg$inputs$fixed_bdev_bundle$fixed_space))
wb_close_bundle(bdev)
stopifnot(nrow(atac_umap) == cfg$inputs$fixed_bdev_bundle$expected_cells,
  all(reference_cells %in% atac_umap$cell_id))
rownames(atac_umap) <- atac_umap$cell_id

project_one <- function(cell) {
  rows <- neighbors[query_cell_id == cell &
    neighbor_rank <= cfg$placement$projection_k][order(neighbor_rank)]
  stopifnot(nrow(rows) == cfg$placement$projection_k)
  kth <- rows$distance[[nrow(rows)]]
  weight <- if (kth == 0) as.numeric(rows$distance == 0) else
    exp(-(rows$distance / kth)^2)
  xy <- as.matrix(atac_umap[rows$cell_id, c("UMAP_1", "UMAP_2")])
  center <- colSums(xy * weight) / sum(weight)
  spread <- sqrt(sum(weight * ((xy[, 1L] - center[[1L]])^2 +
    (xy[, 2L] - center[[2L]])^2)) / sum(weight))
  data.frame(slot_id = cell, UMAP_1 = center[[1L]], UMAP_2 = center[[2L]],
    neighbor_spread = spread, nearest_distance = min(rows$distance),
    kth_distance = kth, stringsAsFactors = FALSE)
}
placement <- rbindlist(lapply(query_cells, project_one))
slot_map <- fread(cfg$inputs$slot_map$path)
outcomes <- as.data.table(read.delim(gzfile(
  cfg$inputs$mapping_outcomes$path), check.names = FALSE,
  quote = "\"", na.strings = "NA"))
keep <- c("slot_id", "support_call", "accepted_native_l4")
stopifnot(all(keep %in% names(outcomes)))
placement <- merge(placement, slot_map, by = "slot_id", sort = FALSE)
placement <- merge(placement, outcomes[, ..keep], by = "slot_id", sort = FALSE)
placement <- placement[match(query_cells, slot_id)]
stopifnot(nrow(placement) == 28L, !anyNA(placement$library_id),
  all(placement$support_call == "supported"))
placement[, within_library_index := seq_len(.N), by = library_id]
placement[, plot_id := sprintf("%02d", within_library_index)]
library_labels <- c(FLT3LKO_BM_WT = "Flt3l WT", FLT3LKO_BM_KO = "Flt3l KO",
  IL7KO_BM_WT = "Il7 WT", IL7KO_BM_KO = "Il7 KO")
placement[, library_label := factor(library_labels[library_id],
  levels = c("Flt3l WT", "Flt3l KO", "Il7 WT", "Il7 KO"))]
fwrite(placement, file.path(dirs[["tables"]],
  "02_fixed_umap_placements.tsv"), sep = "\t", quote = TRUE, na = "NA")

family <- readRDS(cfg$inputs$authenticated_l4$pseudotime_path)
take <- as.character(family$key) == cfg$inputs$authenticated_l4$key &
  is.finite(as.numeric(family$value)) & family$cell_id %in% atac_umap$cell_id
trajectory <- data.frame(cell_id = as.character(family$cell_id[take]),
  native_l4 = as.numeric(family$value[take]), stringsAsFactors = FALSE)
trajectory <- merge(trajectory, atac_umap, by = "cell_id", sort = FALSE)
model <- readRDS(cfg$inputs$authenticated_l4$model_path)
curve_value <- slingCurves(model)[[cfg$inputs$authenticated_l4$curve_name]]
curve_xy <- curve_value$s[curve_value$ord, , drop = FALSE]
curve <- data.frame(UMAP_1 = curve_xy[, 1L], UMAP_2 = curve_xy[, 2L])
l4_limits <- range(c(trajectory$native_l4,
  placement$accepted_native_l4), finite = TRUE)
x_range <- range(atac_umap$UMAP_1); y_range <- range(atac_umap$UMAP_2)
x_limits <- x_range + c(-1, 1) * diff(x_range) * 0.04
y_limits <- y_range + c(-1, 1) * diff(y_range) * 0.04

p <- ggplot() +
  geom_point(data = atac_umap, aes(UMAP_1, UMAP_2), color = "#D0D0D0",
    size = 0.18, alpha = 0.36) +
  geom_point(data = trajectory, aes(UMAP_1, UMAP_2, color = native_l4),
    size = 0.42, alpha = 0.58) +
  scale_color_viridis_c(option = "inferno", limits = l4_limits,
    name = "Bdev native L4") +
  geom_path(data = curve, aes(UMAP_1, UMAP_2), color = "#111111",
    linewidth = 0.7, alpha = 0.9,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.11, "inches"))) +
  geom_point(data = placement[is.finite(accepted_native_l4)],
    aes(UMAP_1, UMAP_2, fill = accepted_native_l4), shape = 21,
    color = "#00A6D6", stroke = 1.05, size = 4.2) +
  scale_fill_viridis_c(option = "inferno", limits = l4_limits,
    name = "Query native L4") +
  geom_point(data = placement[!is.finite(accepted_native_l4)],
    aes(UMAP_1, UMAP_2), shape = 4, color = "#CC0077",
    stroke = 1.35, size = 4.1) +
  geom_label_repel(data = placement, aes(UMAP_1, UMAP_2, label = plot_id),
    seed = as.integer(cfg$placement$label_seed), max.overlaps = Inf,
    box.padding = 0.25, point.padding = 0.3, force = 1.4,
    min.segment.length = 0, segment.color = "#555555",
    segment.size = 0.25, size = 2.7, label.size = 0.15,
    fill = grDevices::adjustcolor("white", 0.88), max.time = 15) +
  facet_wrap(~library_label, ncol = 2) +
  coord_equal(xlim = x_limits, ylim = y_limits, clip = "off") +
  labs(title = "FLT3LKO and IL7KO fixture on the fixed Bdev ATAC UMAP",
    subtitle = "Seven cells per library; cyan-edged circles have finite native L4; magenta × marks supported mappings without finite L4",
    x = "ATAC UMAP 1", y = "ATAC UMAP 2",
    caption = paste("Points are distance-kernel barycenters of 30 exact PeakVI neighbors.",
      "The 26,723-cell umap.atac and native L4 curve are unchanged; this is not a UMAP transform.")) +
  theme_classic(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 9.5), legend.position = "right",
    plot.margin = margin(8, 20, 8, 8))

png_path <- file.path(dirs[["figures"]],
  "01_flt3lko_il7ko_fixture_fixed_atac_umap.png")
pdf_path <- file.path(dirs[["figures"]],
  "01_flt3lko_il7ko_fixture_fixed_atac_umap.pdf")
ggsave(png_path, p, width = 13.5, height = 10.2, dpi = 240, bg = "white")
ggsave(pdf_path, p, width = 13.5, height = 10.2, device = cairo_pdf,
  bg = "white")

group_summary <- placement[, .(cells = .N,
  finite_native_l4 = sum(is.finite(accepted_native_l4)),
  missing_native_l4 = sum(!is.finite(accepted_native_l4)),
  median_native_l4 = median(accepted_native_l4, na.rm = TRUE),
  median_neighbor_spread = median(neighbor_spread)),
  by = .(knockout_target, genotype, library_id)]
fwrite(group_summary, file.path(dirs[["tables"]],
  "03_fixed_umap_group_summary.tsv"), sep = "\t", quote = TRUE, na = "NA")

source_projection_after <- sha256(file.path(
  cfg$inputs$projection_bundle$path, "manifest.json"))
fixed_bdev_after <- sha256(file.path(cfg$inputs$fixed_bdev_bundle$path,
  "manifest.json"))
stopifnot(source_projection_after ==
    cfg$inputs$projection_bundle$manifest_sha256,
  fixed_bdev_after == cfg$inputs$fixed_bdev_bundle$manifest_sha256)
summary <- list(
  schema_version = "bdev.flt3lko_il7ko_fixture_fixed_umap_summary.v1",
  status = "complete", run_identity = cfg$run_identity,
  fixed_space = cfg$inputs$fixed_bdev_bundle$fixed_space,
  fixed_background_cells = nrow(atac_umap),
  trajectory_cells = nrow(trajectory), query_cells = nrow(placement),
  supported_query_cells = sum(placement$support_call == "supported"),
  finite_native_l4 = sum(is.finite(placement$accepted_native_l4)),
  missing_native_l4 = sum(!is.finite(placement$accepted_native_l4)),
  placement = list(neighbor_space = cfg$inputs$projection_bundle$joint_space,
    neighbor_k = cfg$placement$neighbor_k,
    projection_k = cfg$placement$projection_k,
    weight_kernel = cfg$placement$weight_kernel,
    semantics = "fixed_umap_neighbor_barycenter_not_umap_transform"),
  native_l4_orientation = "native", source_manifests_unchanged = TRUE,
  model_refit = FALSE, umap_refit = FALSE, genotype_contrast = FALSE,
  paper_analysis = FALSE, implementation_unlocked = FALSE)
write_json(summary, file.path(dirs[["reports"]],
  "01_fixed_umap_summary.json"))
write_json(c(summary, list(config_path = cfg_path,
  config_sha256 = sha256(cfg_path), sif_sha256 = sha256(cfg$runtime$sif),
  scworkbench_version = as.character(packageVersion("scworkbench")),
  scwbRepresentation_version = as.character(
    packageVersion("scwbRepresentation")),
  scwbRepresentation_manifest_sha256 = sha256(representation_manifest),
  source_projection_manifest_sha256 = source_projection_after,
  fixed_bdev_manifest_sha256 = fixed_bdev_after)),
  file.path(dirs[["manifests"]], "02_execution_receipt.json"))
cat("Completed fixed-UMAP placement for", nrow(placement), "query cells\n")
