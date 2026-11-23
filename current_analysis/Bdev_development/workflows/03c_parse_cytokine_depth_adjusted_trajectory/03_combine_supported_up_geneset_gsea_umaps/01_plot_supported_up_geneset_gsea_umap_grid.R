#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config"))
  stop("Usage: 01_plot_supported_up_geneset_gsea_umap_grid.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(scworkbench)
  library(yaml)
  library(jsonlite)
  library(digest)
})

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
root <- cfg$storage$results_root
dirs <- setNames(file.path(root, c("01_tables", "02_figures", "03_reports",
  "04_logs", "06_manifests")), c("tables", "figures", "reports", "logs",
  "manifests"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_supported_up_geneset_gsea_umap_grid.log"),
  split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
for (entry in list(cfg$inputs$parent_inventory, cfg$inputs$figure_index,
    cfg$inputs$cell_assignments))
  stopifnot(identical(sha256(entry$path), entry$sha256))
stopifnot(
  identical(sha256(file.path(cfg$inputs$execution_bundle$path,
    "manifest.json")), cfg$inputs$execution_bundle$manifest_sha256),
  identical(sha256(cfg$runtime$sif), cfg$runtime$sif_sha256),
  identical(as.character(packageVersion("scworkbench")),
    as.character(cfg$runtime$scworkbench_version)))

read_table <- function(path) read.delim(path, check.names = FALSE,
  quote = "\"", na.strings = "NA")
index <- read_table(cfg$inputs$figure_index$path)
assignments <- read_table(cfg$inputs$cell_assignments$path)
selected <- index[index$supported_intervals >=
  cfg$selection$minimum_supported_intervals, , drop = FALSE]
stopifnot(
  nrow(selected) == cfg$acceptance$expected_supported_signatures,
  !anyDuplicated(selected$feature_set),
  !anyDuplicated(assignments[c("feature_set", "cell_id")]),
  setequal(unique(assignments$feature_set), selected$feature_set),
  all(is.finite(assignments$normalized_enrichment_score)),
  all(assignments$cell_id != ""))
selected$panel_number <- seq_len(nrow(selected))
selected <- selected[c("panel_number", "feature_set", "display_label",
  "supported_windows", "supported_intervals", "colored_cells",
  "minimum_nes", "maximum_nes")]
write.table(selected, file.path(dirs[["tables"]],
  "01_supported_up_geneset_grid_index.tsv"), sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")

bundle <- wb_open_bundle(cfg$inputs$execution_bundle$path, validate = FALSE,
  read_only = TRUE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
space <- wb_space(bundle, cfg$inputs$fit_space)
stopifnot(!anyDuplicated(space$cell_id),
  all(assignments$cell_id %in% space$cell_id))

nes_limit <- as.numeric(cfg$selection$shared_nes_limit)
stopifnot(is.finite(nes_limit), nes_limit > 0,
  max(abs(assignments$normalized_enrichment_score)) <= nes_limit)
palette <- grDevices::colorRampPalette(
  c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(256L)
map_color <- function(value) {
  scaled <- (pmax(-nes_limit, pmin(nes_limit, value)) + nes_limit) /
    (2 * nes_limit)
  palette[round(scaled * (length(palette) - 1L)) + 1L]
}

outer_margins <- c(3.8, 3.2, 4.8, 0.4)
draw_panel <- function(row) {
  par(mar = c(0.6, 0.6, 3.2, 0.3))
  plot(space$x, space$y, pch = 16, cex = 0.095, col = "#D2D2D270",
    asp = 1, axes = FALSE, xlab = "", ylab = "", main = "")
  feature_assignments <- assignments[assignments$feature_set ==
    row$feature_set, , drop = FALSE]
  matched <- match(feature_assignments$cell_id, space$cell_id)
  stopifnot(!anyNA(matched), nrow(feature_assignments) == row$colored_cells)
  foreground <- data.frame(x = space$x[matched], y = space$y[matched],
    nes = feature_assignments$normalized_enrichment_score)
  foreground <- foreground[order(abs(foreground$nes)), , drop = FALSE]
  points(foreground$x, foreground$y, pch = 16, cex = 0.39,
    col = map_color(foreground$nes))
  subtitle <- sprintf("%d cells | %d interval%s", row$colored_cells,
    row$supported_intervals, if (row$supported_intervals == 1L) "" else "s")
  mtext(row$display_label, side = 3, line = 1.65, font = 2, cex = 0.86)
  mtext(subtitle, side = 3, line = 0.45, cex = 0.67)
  box(col = "#777777")
}

draw_colorbar <- function() {
  par(mar = c(1.0, 0.4, 3.2, 3.7))
  plot(NA, xlim = c(0, 1), ylim = c(-nes_limit, nes_limit), axes = FALSE,
    xlab = "", ylab = "", type = "n")
  edges <- seq(-nes_limit, nes_limit, length.out = length(palette) + 1L)
  rect(0.16, edges[-length(edges)], 0.66, edges[-1L], col = palette,
    border = NA)
  axis(4, at = seq(-nes_limit, nes_limit, length.out = 5L), las = 1,
    cex.axis = 0.78)
  mtext("Local GSEA NES", side = 4, line = 2.7, cex = 0.85)
  box(col = "#777777")
}

draw_grid <- function() {
  layout_matrix <- matrix(c(1:5, 11, 6:10, 11), nrow = 2L, byrow = TRUE)
  layout(layout_matrix, widths = c(rep(1, 5), 0.23), heights = c(1, 1))
  par(oma = outer_margins)
  for (ii in seq_len(nrow(selected))) draw_panel(selected[ii, , drop = FALSE])
  draw_colorbar()
  mtext("Up-signature supported local-GSEA regions", side = 3, outer = TRUE,
    line = 2.45, font = 2, cex = 1.35)
  mtext("ATAC UMAP 1", side = 1, outer = TRUE, line = 1.55, cex = 0.9)
  mtext("ATAC UMAP 2", side = 2, outer = TRUE, line = 1.45, cex = 0.9)
  mtext("Gray: all Bdev cells | Color: fitted-L4 fixture cells in supported merged local-GSEA intervals",
    side = 1, outer = TRUE, line = 2.65, cex = 0.7)
}

stem <- "01_supported_up_geneset_gsea_regions_umap_grid"
grDevices::png(file.path(dirs[["figures"]], paste0(stem, ".png")),
  width = 18 * 180, height = 8.2 * 180, res = 180, type = "cairo-png")
draw_grid()
grDevices::dev.off()
grDevices::cairo_pdf(file.path(dirs[["figures"]], paste0(stem, ".pdf")),
  width = 18, height = 8.2, onefile = TRUE)
draw_grid()
grDevices::dev.off()

report <- list(
  schema_version = "bdev.parse38_supported_up_gsea_umap_grid_summary.v1",
  status = "complete_development_visualization",
  run_identity = cfg$run_identity,
  parent_run_identity = cfg$scope$parent_run_identity,
  scientific_parent_run_identity = cfg$scope$scientific_parent_run_identity,
  supported_up_signatures = nrow(selected),
  feature_sets = selected$feature_set,
  layout = sprintf("%dx%d", cfg$selection$layout_rows,
    cfg$selection$layout_columns),
  shared_nes_limit = nes_limit,
  background = cfg$selection$background_cells,
  foreground = cfg$selection$foreground_cells,
  parent_runs_mutated = FALSE,
  outputs_promoted = FALSE)
writeLines(jsonlite::toJSON(report, auto_unbox = TRUE, pretty = TRUE),
  file.path(dirs[["reports"]],
    "01_supported_up_geneset_gsea_umap_grid_summary.json"))
receipt <- list(
  schema_version = "bdev.parse38_supported_up_gsea_umap_grid_receipt.v1",
  run_identity = cfg$run_identity,
  config_path = cfg_path,
  config_sha256 = sha256(cfg_path),
  input_parent_inventory_sha256 = sha256(cfg$inputs$parent_inventory$path),
  input_figure_index_sha256 = sha256(cfg$inputs$figure_index$path),
  input_cell_assignments_sha256 = sha256(cfg$inputs$cell_assignments$path),
  figure_png = paste0(stem, ".png"),
  figure_pdf = paste0(stem, ".pdf"))
writeLines(jsonlite::toJSON(receipt, auto_unbox = TRUE, pretty = TRUE),
  file.path(dirs[["manifests"]], "01_plotting_receipt.json"))
cat(sprintf("Created combined %d-panel supported up-signature GSEA UMAP figure in %s\n",
  nrow(selected), dirs[["figures"]]))
