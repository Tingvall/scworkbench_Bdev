#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 02_project_and_plot.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(scwbRepresentation)
  library(ggplot2)
  library(ggrepel)
  library(yaml)
  library(jsonlite)
  library(digest)
})

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
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
sink(file.path(dirs[["logs"]], "02_project_and_plot.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

stopifnot(sha256(cfg$runtime$sif) == cfg$runtime$sif_sha256,
  sha256(file.path(cfg$inputs$mapping_bundle$path, "manifest.json")) ==
    cfg$inputs$mapping_bundle$manifest_sha256,
  sha256(file.path(cfg$inputs$bdev_bundle$path, "manifest.json")) ==
    cfg$inputs$bdev_bundle$manifest_sha256,
  sha256(cfg$inputs$mapping_table$path) == cfg$inputs$mapping_table$sha256,
  sha256(cfg$inputs$authenticated_l4$path) == cfg$inputs$authenticated_l4$sha256,
  cfg$inputs$authenticated_l4$orientation == "native",
  as.character(packageVersion("scworkbench")) == cfg$runtime$scworkbench_version,
  as.character(packageVersion("scwbRepresentation")) == cfg$runtime$scwbRepresentation_version)
manifest <- system.file("scworkbench", "operation-pack.json", package = "scwbRepresentation")
stopifnot(sha256(manifest) == cfg$runtime$scwbRepresentation_manifest_sha256)
wb_operation_pack_activate(manifest, lib.loc = .libPaths())

mapping_bundle <- wb_open_bundle(cfg$storage$execution_bundle, validate = FALSE,
  read_only = FALSE)
reference <- wb_space(mapping_bundle, cfg$inputs$mapping_bundle$reference_space)
query <- wb_space(mapping_bundle, cfg$inputs$mapping_bundle$query_space)
stopifnot(nrow(reference) == cfg$inputs$mapping_bundle$expected_reference_cells,
  nrow(query) == cfg$inputs$mapping_bundle$expected_query_libraries)
reference_cells <- rownames(reference); query_cells <- rownames(query)

knn_session <- wb_session_create(mapping_bundle, cfg$storage$sessions,
  paste0(cfg$run_identity, "_knn"))
knn <- wb_session_run_op(knn_session, mapping_bundle,
  "scwbrepresentation::knn_neighbors", inputs = list(
    space = cfg$inputs$mapping_bundle$joint_space,
    query_cells = query_cells, candidate_cells = reference_cells,
    k = as.integer(cfg$placement$neighbor_k), include_query = FALSE,
    backend = cfg$placement$backend,
    space_chunk_size = as.integer(cfg$placement$space_chunk_size),
    query_chunk_size = as.integer(cfg$placement$query_chunk_size)),
  notes = "Reproduce authenticated ImmGen-to-Bdev PeakVI neighbor ledger for approved ATAC-UMAP placement")
neighbors <- knn$run$outputs$neighbors$table
stopifnot(nrow(neighbors) == length(query_cells) * cfg$placement$neighbor_k,
  identical(unique(as.character(neighbors$query_cell_id)), query_cells))
write_tsv(neighbors, file.path(dirs[["tables"]], "01_immgen_peakvi_neighbor_ledger.tsv"))
knn$session <- wb_session_resume(mapping_bundle, knn$session$durable$workspace)
knn$session <- wb_session_discard(knn$session)
wb_close_bundle(mapping_bundle)

bdev <- wb_open_bundle(cfg$inputs$bdev_bundle$path, validate = FALSE, read_only = TRUE)
atac_umap <- space_2d(wb_space(bdev, cfg$inputs$bdev_bundle$atac_umap_space))
wb_close_bundle(bdev)
stopifnot(nrow(atac_umap) == cfg$inputs$bdev_bundle$expected_cells,
  all(reference_cells %in% atac_umap$cell_id))
reference_umap <- atac_umap[match(reference_cells, atac_umap$cell_id), ]
rownames(reference_umap) <- reference_umap$cell_id

project_one <- function(cell) {
  rows <- neighbors[neighbors$query_cell_id == cell &
    neighbors$neighbor_rank <= cfg$placement$projection_k, ]
  rows <- rows[order(rows$neighbor_rank), ]
  stopifnot(nrow(rows) == cfg$placement$projection_k)
  kth <- rows$distance[[nrow(rows)]]
  weight <- if (kth == 0) as.numeric(rows$distance == 0) else
    exp(-(rows$distance / kth)^2)
  xy <- as.matrix(reference_umap[rows$cell_id, c("UMAP_1", "UMAP_2")])
  center <- colSums(xy * weight) / sum(weight)
  spread <- sqrt(sum(weight * ((xy[, 1] - center[[1L]])^2 +
    (xy[, 2] - center[[2L]])^2)) / sum(weight))
  data.frame(gsm = cell, UMAP_1 = center[[1L]], UMAP_2 = center[[2L]],
    neighbor_spread = spread, nearest_distance = min(rows$distance),
    kth_distance = kth, stringsAsFactors = FALSE)
}
placement <- do.call(rbind, lapply(query_cells, project_one))
mapping <- read.delim(cfg$inputs$mapping_table$path, check.names = FALSE,
  quote = "\"", na.strings = "NA")
placement <- merge(placement, mapping[, c("gsm", "population", "support_call",
  "accepted_native_l4")], by = "gsm", sort = FALSE)
placement <- placement[match(query_cells, placement$gsm), ]
placement$label <- paste0(placement$population, "  (", placement$gsm, ")")
write_tsv(placement, file.path(dirs[["tables"]],
  "02_immgen_approved_atac_umap_placement.tsv"))

family <- readRDS(cfg$inputs$authenticated_l4$path)
keep <- as.character(family$key) == cfg$inputs$authenticated_l4$key &
  is.finite(as.numeric(family$value)) & family$cell_id %in% atac_umap$cell_id
trajectory <- data.frame(cell_id = as.character(family$cell_id[keep]),
  native_l4 = as.numeric(family$value[keep]), stringsAsFactors = FALSE)
trajectory <- merge(trajectory, atac_umap, by = "cell_id", sort = FALSE)
l4_limits <- range(c(trajectory$native_l4, placement$accepted_native_l4), finite = TRUE)

p <- ggplot() +
  geom_point(data = atac_umap, aes(UMAP_1, UMAP_2), color = "#C8C8C8",
    size = 0.28, alpha = 0.42) +
  geom_point(data = trajectory, aes(UMAP_1, UMAP_2, color = native_l4),
    size = 0.72, alpha = 0.8) +
  scale_color_viridis_c(option = "inferno", limits = l4_limits,
    name = "Native fitted L4") +
  geom_point(data = placement[placement$support_call == "supported", ],
    aes(UMAP_1, UMAP_2, fill = accepted_native_l4), shape = 21,
    color = "#111111", stroke = 0.7, size = 4.2) +
  geom_point(data = placement[placement$support_call == "borderline", ],
    aes(UMAP_1, UMAP_2, fill = accepted_native_l4), shape = 24,
    color = "#D55E00", stroke = 1.35, size = 4.6) +
  scale_fill_viridis_c(option = "inferno", limits = l4_limits, guide = "none") +
  geom_label_repel(data = placement, aes(UMAP_1, UMAP_2, label = label),
    seed = as.integer(cfg$placement$label_seed), max.overlaps = Inf,
    box.padding = 0.42, point.padding = 0.48, force = 2.2,
    min.segment.length = 0, segment.color = "#555555", segment.size = 0.28,
    size = 2.8, label.size = 0.18, fill = grDevices::adjustcolor("white", 0.9),
    max.time = 15) +
  labs(title = "ImmGen ATAC libraries on the approved Bdev ATAC UMAP",
    subtitle = "PeakVI k=30 neighbor-weighted placement; all library names shown; orange triangles are borderline support",
    x = "ATAC UMAP 1", y = "ATAC UMAP 2",
    caption = "The established umap.atac coordinates and native Slingshot L4 orientation are unchanged.") +
  coord_equal(clip = "off") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 10.5), plot.margin = margin(10, 30, 10, 10),
    legend.position = "right")
ggsave(file.path(dirs[["figures"]], "01_immgen_on_approved_atac_umap.png"),
  p, width = 14, height = 10, dpi = 220, bg = "white")
ggsave(file.path(dirs[["figures"]], "01_immgen_on_approved_atac_umap.pdf"),
  p, width = 14, height = 10, device = cairo_pdf, bg = "white")

mapping_manifest_after <- sha256(file.path(cfg$inputs$mapping_bundle$path, "manifest.json"))
bdev_manifest_after <- sha256(file.path(cfg$inputs$bdev_bundle$path, "manifest.json"))
stopifnot(mapping_manifest_after == cfg$inputs$mapping_bundle$manifest_sha256,
  bdev_manifest_after == cfg$inputs$bdev_bundle$manifest_sha256,
  nrow(placement) == 28L, all(nzchar(placement$label)))
summary <- list(schema_version = "bdev.immgen_atac_umap_overlay_summary.v1",
  status = "complete", run_identity = cfg$run_identity,
  atac_umap_space = cfg$inputs$bdev_bundle$atac_umap_space,
  background_cells = nrow(atac_umap), trajectory_cells = nrow(trajectory),
  query_libraries = nrow(placement), labels_shown = nrow(placement),
  placement = list(neighbor_space = cfg$inputs$mapping_bundle$joint_space,
    neighbor_k = cfg$placement$neighbor_k, projection_k = cfg$placement$projection_k,
    weight_kernel = cfg$placement$weight_kernel),
  source_manifests_unchanged = TRUE, native_l4_orientation = "native",
  paper_analysis = FALSE, implementation_unlocked = FALSE, sif_built = FALSE)
write_json(summary, file.path(dirs[["reports"]], "01_immgen_atac_umap_summary.json"))
write_json(c(summary, list(config_path = cfg_path, config_sha256 = sha256(cfg_path),
  mapping_manifest_sha256 = mapping_manifest_after,
  bdev_manifest_sha256 = bdev_manifest_after)),
  file.path(dirs[["manifests"]], "02_execution_receipt.json"))
cat("Completed approved ATAC-UMAP placement for", nrow(placement), "labeled libraries\n")
