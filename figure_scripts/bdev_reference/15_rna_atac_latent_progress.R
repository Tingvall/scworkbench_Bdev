#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 15_rna_atac_latent_progress.R",
  "<bdev_bundle> <output_dir> <expected_bundle_fingerprint>",
  "<window_size> <window_stride>"))

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
expected_fingerprint <- args[[3L]]
window_size <- as.integer(args[[4L]])
window_stride <- as.integer(args[[5L]])
stopifnot(is.finite(window_size), window_size >= 2L,
  is.finite(window_stride), window_stride >= 1L,
  window_stride <= window_size)
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
})

write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
fingerprint_before <- manifest$current_state$fingerprint
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(fingerprint_before, expected_fingerprint))

all_cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = all_cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "pt", "bb_up"))
stopifnot(identical(as.character(metadata$cell_id), all_cells),
  length(all_cells) == 26723L,
  sum(!is.na(metadata$pt)) == 12642L,
  sum(is.na(metadata$pt)) == 14081L,
  all(is.finite(metadata$pt[!is.na(metadata$pt)])))

# Preserve the authoritative pt values exactly. Cell id is only a deterministic
# tie-breaker for cells having the same pt value.
trajectory_metadata <- metadata[!is.na(metadata$pt), , drop = FALSE]
trajectory_metadata <- trajectory_metadata[
  order(trajectory_metadata$pt, trajectory_metadata$cell_id), , drop = FALSE]
trajectory_cells <- as.character(trajectory_metadata$cell_id)
n_trajectory <- length(trajectory_cells)
stopifnot(n_trajectory >= window_size, !anyDuplicated(trajectory_cells))

window_starts <- seq.int(1L, n_trajectory - window_size + 1L,
  by = window_stride)
final_start <- n_trajectory - window_size + 1L
if (tail(window_starts, 1L) != final_start) {
  window_starts <- c(window_starts, final_start)
}
window_ends <- window_starts + window_size - 1L
n_windows <- length(window_starts)
window_ids <- sprintf("window_%04d", seq_len(n_windows))
stopifnot(window_starts[[1L]] == 1L,
  tail(window_ends, 1L) == n_trajectory,
  all(window_ends - window_starts + 1L == window_size))

log_event("Reading 50-dimensional RNA PCA and ATAC LSI for",
  n_trajectory, "trajectory cells")
rna <- wb_space(bundle, "pca", cells = trajectory_cells)
atac <- wb_space(bundle, "lsi", cells = trajectory_cells)
stopifnot(is.matrix(rna), is.matrix(atac),
  identical(rownames(rna), trajectory_cells),
  identical(rownames(atac), trajectory_cells),
  ncol(rna) == 50L, ncol(atac) == 50L,
  all(is.finite(rna)), all(is.finite(atac)))

centroid_matrix <- function(space) {
  out <- t(vapply(seq_len(n_windows), function(i) {
    colMeans(space[window_starts[[i]]:window_ends[[i]], , drop = FALSE])
  }, numeric(ncol(space))))
  rownames(out) <- window_ids
  colnames(out) <- colnames(space)
  out
}
path_metrics <- function(centroids) {
  delta <- centroids[-1L, , drop = FALSE] -
    centroids[-nrow(centroids), , drop = FALSE]
  step <- c(0, sqrt(rowSums(delta^2)))
  cumulative <- cumsum(step)
  stopifnot(all(is.finite(step)), tail(cumulative, 1L) > 0)
  relative <- cumulative / tail(cumulative, 1L)
  list(step = step, cumulative = cumulative, relative = relative,
    total = tail(cumulative, 1L))
}

rna_centroids <- centroid_matrix(rna)
atac_centroids <- centroid_matrix(atac)
rna_path <- path_metrics(rna_centroids)
atac_path <- path_metrics(atac_centroids)

window_summary <- data.frame(
  window_id = window_ids,
  window_index = seq_len(n_windows),
  start_rank = window_starts,
  end_rank = window_ends,
  n_cells = window_size,
  pt_min = trajectory_metadata$pt[window_starts],
  pt_median = vapply(seq_len(n_windows), function(i)
    median(trajectory_metadata$pt[window_starts[[i]]:window_ends[[i]]]),
    numeric(1)),
  pt_max = trajectory_metadata$pt[window_ends],
  atac_step_distance = atac_path$step,
  atac_cumulative_distance = atac_path$cumulative,
  atac_relative_distance = atac_path$relative,
  rna_step_distance = rna_path$step,
  rna_cumulative_distance = rna_path$cumulative,
  rna_relative_distance = rna_path$relative,
  rna_minus_atac_relative_distance = rna_path$relative - atac_path$relative,
  stringsAsFactors = FALSE)
stopifnot(window_summary$atac_relative_distance[[1L]] == 0,
  tail(window_summary$atac_relative_distance, 1L) == 1,
  window_summary$rna_relative_distance[[1L]] == 0,
  tail(window_summary$rna_relative_distance, 1L) == 1,
  all(diff(window_summary$atac_relative_distance) >= 0),
  all(diff(window_summary$rna_relative_distance) >= 0))

membership <- do.call(rbind, lapply(seq_len(n_windows), function(i) {
  ix <- window_starts[[i]]:window_ends[[i]]
  data.frame(window_id = window_ids[[i]], window_index = i,
    position_in_window = seq_along(ix), trajectory_rank = ix,
    cell_id = trajectory_cells[ix], pt = trajectory_metadata$pt[ix],
    bb_up = trajectory_metadata$bb_up[ix], stringsAsFactors = FALSE)
}))
stopifnot(nrow(membership) == n_windows * window_size,
  all(table(membership$window_id) == window_size))

bb_levels <- c("up1", "up2")
bb_rows <- lapply(seq_along(bb_levels), function(i) {
  level <- bb_levels[[i]]
  values <- metadata$pt[!is.na(metadata$pt) & !is.na(metadata$bb_up) &
    metadata$bb_up == level]
  stopifnot(length(values) > 0L, all(is.finite(values)))
  boundary_pt <- median(values)
  nearest <- which.min(abs(window_summary$pt_median - boundary_pt))
  data.frame(boundary = paste0("BB", i), bb_up = level,
    n_source_cells = length(values), boundary_pt = boundary_pt,
    highlighted_window_id = window_summary$window_id[[nearest]],
    highlighted_window_index = nearest,
    highlighted_window_pt = window_summary$pt_median[[nearest]],
    atac_relative_distance = window_summary$atac_relative_distance[[nearest]],
    rna_relative_distance = window_summary$rna_relative_distance[[nearest]],
    stringsAsFactors = FALSE)
})
bb <- do.call(rbind, bb_rows)
stopifnot(!anyDuplicated(bb$highlighted_window_id))
window_summary$highlight <- NA_character_
window_summary$highlight[match(bb$highlighted_window_id,
  window_summary$window_id)] <- bb$boundary

write_tsv(window_summary, file.path(stage_dir, "window_distance_summary.tsv"))
membership_connection <- gzfile(file.path(stage_dir,
  "window_membership.tsv.gz"), open = "wt")
write_tsv(membership, membership_connection)
close(membership_connection)
write_tsv(bb, file.path(stage_dir, "bb1_bb2_highlights.tsv"))
write_tsv(data.frame(cell_id = trajectory_metadata$cell_id,
  pt = trajectory_metadata$pt, bb_up = trajectory_metadata$bb_up),
  file.path(stage_dir, "trajectory_cells_ordered.tsv"))

plot_data <- window_summary
highlight_data <- plot_data[!is.na(plot_data$highlight), , drop = FALSE]
p <- ggplot(plot_data, aes(atac_relative_distance,
  rna_relative_distance)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
    linewidth = 0.75, color = "#4A4A4A") +
  geom_path(color = "#7B8794", linewidth = 0.45, alpha = 0.75) +
  geom_point(aes(color = pt_median), size = 2.2, alpha = 0.9) +
  geom_point(data = highlight_data, shape = 21, size = 5.2,
    stroke = 1.2, color = "#111111", fill = "#FFD43B") +
  ggrepel::geom_label_repel(data = highlight_data, aes(label = highlight),
    color = "#111111", fill = "#FFF4B8", fontface = "bold", size = 4,
    min.segment.length = 0, seed = 1729, box.padding = 0.45,
    point.padding = 0.5, show.legend = FALSE) +
  scale_color_viridis_c(option = "plasma", direction = -1,
    name = "Window median\npseudotime") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.015, 0.04))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.015, 0.04))) +
  coord_equal(clip = "off") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
    legend.position = "right") +
  labs(
    title = "Relative RNA and ATAC progression along B-cell development",
    subtitle = paste0(n_trajectory, " cells with defined pt; ", window_size,
      " cells/window; ", window_stride, "-cell stride; 50 latent dimensions"),
    x = "Relative ATAC latent path distance (LSI)",
    y = "Relative RNA latent path distance (PCA)",
    caption = paste(
      "Each point is a sliding-window centroid; distances are cumulative",
      "Euclidean steps and are independently scaled from 0 to 1.",
      "The dotted diagonal denotes matched relative progress.", sep = "\n"))

ggsave(file.path(stage_dir, "rna_vs_atac_relative_latent_distance.png"), p,
  width = 9.2, height = 8.2, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "rna_vs_atac_relative_latent_distance.pdf"), p,
  width = 9.2, height = 8.2, device = cairo_pdf, bg = "white")

diagnostics <- list(
  status = "complete",
  bundle_id = manifest$bundle_id,
  bundle_fingerprint = fingerprint_before,
  source_cells = length(all_cells),
  trajectory_definition = "authoritative metadata pt is non-missing",
  trajectory_cells = n_trajectory,
  missing_pt_cells_excluded = sum(is.na(metadata$pt)),
  pseudotime_values_modified = FALSE,
  rna_space = list(id = "pca", dimensions = ncol(rna),
    coordinate_scaling = "unchanged imported coordinates",
    total_cumulative_distance = rna_path$total),
  atac_space = list(id = "lsi", dimensions = ncol(atac),
    coordinate_scaling = "unchanged imported coordinates",
    total_cumulative_distance = atac_path$total),
  windows = list(size = window_size, stride = window_stride,
    count = n_windows, terminal_window_appended =
      (n_trajectory - window_size) %% window_stride != 0L,
    all_windows_full_size = TRUE),
  distance = list(metric = "Euclidean",
    summary = "centroid per window",
    path = "cumulative distance between consecutive window centroids",
    standardization = "divide each modality cumulative path by its endpoint"),
  endpoint_checks = list(
    atac_start = window_summary$atac_relative_distance[[1L]],
    atac_end = tail(window_summary$atac_relative_distance, 1L),
    rna_start = window_summary$rna_relative_distance[[1L]],
    rna_end = tail(window_summary$rna_relative_distance, 1L)),
  bb_highlight_rule = paste(
    "Median source pt for bb_up up1/up2; highlight nearest window median pt"),
  promoted = FALSE,
  evidence_mode = "external question-specific terminal analysis; capsule registration follows")
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))

report <- c(
  "# Bdev RNA-versus-ATAC latent progression",
  "",
  "## Question",
  "",
  paste("Compare relative RNA and ATAC latent-space movement along the",
    "existing Bdev developmental trajectory."),
  "",
  "## Input state",
  "",
  paste0("- Accepted bundle: `", manifest$bundle_id, "`"),
  paste0("- Bundle fingerprint at execution: `", fingerprint_before, "`"),
  paste0("- Cells: ", length(all_cells), " total; ", n_trajectory,
    " with non-missing authoritative `pt`; ", sum(is.na(metadata$pt)),
    " excluded because `pt` is missing."),
  "- RNA representation: imported 50-dimensional `pca` space.",
  "- ATAC representation: imported 50-dimensional `lsi` space.",
  "",
  "## Method",
  "",
  paste0("Cells were ordered by unchanged `pt` (cell id breaks exact ties), then",
    " grouped into full sliding windows of ", window_size, " cells with a ",
    window_stride, "-cell stride. A final full window was appended so that the",
    " endpoint cell is included. For each window and modality, the centroid was",
    " calculated in all 50 imported latent dimensions. Consecutive centroid",
    " distances were Euclidean and accumulated along the ordered windows. Each",
    " modality was independently divided by its final cumulative distance, so",
    " the first window is 0 and the endpoint window is 1."),
  "",
  paste("BB1 and BB2 are the windows whose median `pt` is nearest to the median",
    "source `pt` of `bb_up == up1` and `bb_up == up2`, respectively."),
  "",
  "## Outputs",
  "",
  "- `rna_vs_atac_relative_latent_distance.png` and `.pdf`: requested figure.",
  "- `window_distance_summary.tsv`: raw step, cumulative, and relative distances.",
  "- `window_membership.tsv.gz`: exact membership of every window.",
  "- `bb1_bb2_highlights.tsv`: BB source medians and highlighted windows.",
  "- `trajectory_cells_ordered.tsv`: trajectory cell order with unchanged `pt`.",
  "- `diagnostics.json`: input, method, and endpoint checks.",
  "",
  "## Interpretation boundary",
  "",
  paste("Departure from the dotted diagonal describes relative movement in these",
    "specific imported PCA and LSI representations. It does not establish a",
    "causal ordering between transcription and chromatin accessibility, and it",
    "can depend on latent-space construction, dimensionality, and window size."),
  "",
  "## Provenance",
  "",
  paste("This dataset-specific terminal calculation is registered after execution",
    "as an scWorkbench analysis capsule. scWorkbench binds the input entity",
    "identities and exact script/output bytes; it does not claim engine-observed",
    "execution or replay equivalence."))
writeLines(report, file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

files <- list.files(stage_dir, full.names = TRUE)
output_manifest <- data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size), stringsAsFactors = FALSE)
write_tsv(output_manifest, file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
log_event("Completed Bdev RNA/ATAC latent progression report at", output_dir)
