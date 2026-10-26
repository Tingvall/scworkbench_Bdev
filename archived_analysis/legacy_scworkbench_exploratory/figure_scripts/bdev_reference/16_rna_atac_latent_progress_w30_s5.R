#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 16_rna_atac_latent_progress_w30_s5.R",
  "<bdev_bundle> <output_dir> <expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
expected_fingerprint <- args[[3L]]
window_size <- 30L
window_stride <- 5L
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

trajectory_metadata <- metadata[!is.na(metadata$pt), , drop = FALSE]
trajectory_metadata <- trajectory_metadata[
  order(trajectory_metadata$pt, trajectory_metadata$cell_id), , drop = FALSE]
trajectory_cells <- as.character(trajectory_metadata$cell_id)
n_trajectory <- length(trajectory_cells)

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

log_event("Reading RNA PCA and ATAC LSI for", n_trajectory,
  "trajectory cells and", n_windows, "windows")
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
  list(step = step, cumulative = cumulative,
    relative = cumulative / tail(cumulative, 1L),
    total = tail(cumulative, 1L))
}
centered_change <- function(x) {
  n <- length(x)
  c(x[[2L]] - x[[1L]], x[3:n] - x[1:(n - 2L)],
    x[[n]] - x[[n - 1L]])
}

rna_path <- path_metrics(centroid_matrix(rna))
atac_path <- path_metrics(centroid_matrix(atac))
rna_change <- centered_change(rna_path$relative)
atac_change <- centered_change(atac_path$relative)
stopifnot(all(is.finite(rna_change)), all(is.finite(atac_change)),
  all(rna_change > 0), all(atac_change > 0))
local_slope <- rna_change / atac_change
stopifnot(all(is.finite(local_slope)), all(local_slope > 0))

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
  centered_atac_relative_change = atac_change,
  centered_rna_relative_change = rna_change,
  local_slope_rna_per_atac = local_slope,
  stringsAsFactors = FALSE)
stopifnot(window_summary$atac_relative_distance[[1L]] == 0,
  tail(window_summary$atac_relative_distance, 1L) == 1,
  window_summary$rna_relative_distance[[1L]] == 0,
  tail(window_summary$rna_relative_distance, 1L) == 1,
  all(diff(window_summary$atac_relative_distance) > 0),
  all(diff(window_summary$rna_relative_distance) > 0))

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
bb <- do.call(rbind, lapply(seq_along(bb_levels), function(i) {
  level <- bb_levels[[i]]
  values <- metadata$pt[!is.na(metadata$pt) & !is.na(metadata$bb_up) &
    metadata$bb_up == level]
  stopifnot(length(values) > 0L, all(is.finite(values)))
  boundary_pt <- median(values)
  nearest <- which.min(abs(window_summary$pt_median - boundary_pt))
  data.frame(boundary = paste0("BB", i), bb_up = level,
    n_source_cells = length(values), pt_min = min(values),
    pt_median = boundary_pt, pt_max = max(values),
    highlighted_window_id = window_summary$window_id[[nearest]],
    highlighted_window_index = nearest,
    highlighted_window_pt = window_summary$pt_median[[nearest]],
    atac_relative_distance = window_summary$atac_relative_distance[[nearest]],
    rna_relative_distance = window_summary$rna_relative_distance[[nearest]],
    local_slope_rna_per_atac = window_summary$local_slope_rna_per_atac[[nearest]],
    stringsAsFactors = FALSE)
}))
stopifnot(!anyDuplicated(bb$highlighted_window_id),
  all(bb$pt_min < bb$pt_median), all(bb$pt_median < bb$pt_max))
window_summary$highlight <- NA_character_
window_summary$highlight[match(bb$highlighted_window_id,
  window_summary$window_id)] <- bb$boundary

write_tsv(window_summary, file.path(stage_dir, "window_distance_and_slope.tsv"))
membership_connection <- gzfile(file.path(stage_dir,
  "window_membership.tsv.gz"), open = "wt")
write_tsv(membership, membership_connection)
close(membership_connection)
write_tsv(bb, file.path(stage_dir, "bb1_bb2_regions.tsv"))
write_tsv(data.frame(cell_id = trajectory_metadata$cell_id,
  pt = trajectory_metadata$pt, bb_up = trajectory_metadata$bb_up),
  file.path(stage_dir, "trajectory_cells_ordered.tsv"))

highlight_data <- window_summary[!is.na(window_summary$highlight), , drop = FALSE]
scatter <- ggplot(window_summary, aes(atac_relative_distance,
  rna_relative_distance)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
    linewidth = 0.75, color = "#4A4A4A") +
  geom_path(color = "#7B8794", linewidth = 0.35, alpha = 0.7) +
  geom_point(aes(color = pt_median), size = 1.25, alpha = 0.82) +
  geom_point(data = highlight_data, shape = 21, size = 5.2,
    stroke = 1.2, color = "#111111", fill = "#FFD43B") +
  ggrepel::geom_label_repel(data = highlight_data, aes(label = highlight),
    color = "#111111", fill = "#FFF4B8", fontface = "bold", size = 4,
    min.segment.length = 0, seed = 1731, box.padding = 0.45,
    point.padding = 0.5, show.legend = FALSE) +
  scale_color_viridis_c(option = "plasma", direction = -1,
    name = "Window median\npseudotime") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.015, 0.04))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.015, 0.04))) +
  coord_equal(clip = "off") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right") +
  labs(
    title = "Relative RNA and ATAC progression along B-cell development",
    subtitle = paste0(n_trajectory, " cells with defined pt; 30 cells/window; ",
      "5-cell stride; 50 latent dimensions"),
    x = "Relative ATAC latent path distance (LSI)",
    y = "Relative RNA latent path distance (PCA)",
    caption = paste(
      "Each point is a sliding-window centroid; distances are cumulative",
      "Euclidean steps and are independently scaled from 0 to 1.",
      "The dotted diagonal denotes matched relative progress.", sep = "\n"))

bb_fill <- c(BB1 = "#2C7FB8", BB2 = "#D95F0E")
slope_plot <- ggplot(window_summary,
  aes(pt_median, local_slope_rna_per_atac)) +
  geom_rect(data = bb, aes(xmin = pt_min, xmax = pt_max, fill = boundary),
    ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.16) +
  geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.75,
    color = "#4A4A4A") +
  geom_line(color = "#4B1D82", linewidth = 0.48, alpha = 0.9) +
  geom_vline(data = bb, aes(xintercept = pt_median, color = boundary),
    linetype = "dashed", linewidth = 0.75, show.legend = FALSE) +
  geom_label(data = bb,
    aes(x = pt_median, y = max(window_summary$local_slope_rna_per_atac),
      label = boundary, color = boundary),
    inherit.aes = FALSE, fontface = "bold", fill = "white", size = 3.8,
    linewidth = 0.25, show.legend = FALSE) +
  scale_y_log10() +
  scale_fill_manual(values = bb_fill, name = "Source bb_up region") +
  scale_color_manual(values = bb_fill) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right") +
  labs(
    title = "Local RNA-versus-ATAC progression slope",
    subtitle = paste(
      "Centered finite difference across adjacent 30-cell windows;",
      "BB regions span the source bb_up pseudotime ranges"),
    x = "Bdev pseudotime (window median pt)",
    y = "Local slope: ΔRNA relative distance / ΔATAC relative distance",
    caption = paste(
      "Log10 y-axis. Slope = 1 indicates matched local relative progress;",
      "values below 1 indicate faster relative ATAC movement;",
      "values above 1 indicate faster relative RNA movement.", sep = "\n"))

ggsave(file.path(stage_dir, "rna_vs_atac_relative_latent_distance_w30_s5.png"),
  scatter, width = 9.2, height = 8.2, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "rna_vs_atac_relative_latent_distance_w30_s5.pdf"),
  scatter, width = 9.2, height = 8.2, device = cairo_pdf, bg = "white")
ggsave(file.path(stage_dir, "rna_vs_atac_local_slope_w30_s5.png"),
  slope_plot, width = 11, height = 6.8, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "rna_vs_atac_local_slope_w30_s5.pdf"),
  slope_plot, width = 11, height = 6.8, device = cairo_pdf, bg = "white")

diagnostics <- list(
  status = "complete", bundle_id = manifest$bundle_id,
  bundle_fingerprint = fingerprint_before,
  source_cells = length(all_cells), trajectory_cells = n_trajectory,
  missing_pt_cells_excluded = sum(is.na(metadata$pt)),
  pseudotime_values_modified = FALSE,
  rna_space = list(id = "pca", dimensions = ncol(rna),
    total_cumulative_distance = rna_path$total),
  atac_space = list(id = "lsi", dimensions = ncol(atac),
    total_cumulative_distance = atac_path$total),
  windows = list(size = window_size, stride = window_stride,
    count = n_windows, terminal_window_appended =
      (n_trajectory - window_size) %% window_stride != 0L,
    all_windows_full_size = TRUE),
  distance = list(metric = "Euclidean", summary = "centroid per window",
    path = "cumulative consecutive-window centroid distance",
    standardization = "divide by modality endpoint cumulative distance"),
  slope = list(
    definition = "centered finite difference of RNA relative distance with respect to ATAC relative distance",
    endpoint_rule = "one-sided finite difference",
    display_scale = "log10",
    reference = 1),
  bb_region_rule = "source pt min-to-max for bb_up up1/up2; dashed line at median",
  endpoint_checks = list(atac_start = atac_path$relative[[1L]],
    atac_end = tail(atac_path$relative, 1L),
    rna_start = rna_path$relative[[1L]],
    rna_end = tail(rna_path$relative, 1L)),
  promoted = FALSE,
  evidence_mode = "external question-specific terminal analysis; capsule registration follows")
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))

report <- c(
  "# Bdev RNA-versus-ATAC latent progression: 30-cell windows",
  "",
  "## Input and method",
  "",
  paste0("The analysis uses ", n_trajectory,
    " cells with non-missing authoritative `pt`, 30-cell full sliding windows",
    " with a 5-cell stride, and all 50 dimensions of imported RNA `pca` and",
    " ATAC `lsi`. Pseudotime values and missingness are unchanged."),
  "",
  paste("For each modality, consecutive-window centroid distances are Euclidean,",
    "cumulative, and independently divided by the endpoint. Local slope is the",
    "centered finite difference in RNA relative distance divided by the centered",
    "finite difference in ATAC relative distance. Endpoints use one-sided",
    "differences. The slope plot uses a log10 y-axis and a slope-one reference."),
  "",
  paste("BB1 and BB2 shading spans the full source pseudotime range of",
    "`bb_up == up1` and `bb_up == up2`; dashed lines mark their medians."),
  "",
  "## Outputs",
  "",
  "- `rna_vs_atac_relative_latent_distance_w30_s5.png` and `.pdf`: modality-progress scatter.",
  "- `rna_vs_atac_local_slope_w30_s5.png` and `.pdf`: local slope along pseudotime with BB regions.",
  "- `window_distance_and_slope.tsv`: raw distance and slope values.",
  "- `window_membership.tsv.gz`: exact window memberships.",
  "- `bb1_bb2_regions.tsv`: BB ranges, medians, highlighted windows, and slopes.",
  "- `trajectory_cells_ordered.tsv`: trajectory order with unchanged source metadata.",
  "",
  "## Interpretation boundary",
  "",
  paste("These results describe movement within the imported PCA and LSI",
    "representations. Local slopes can be sensitive to representation geometry",
    "and window parameters and do not establish causal molecular timing."),
  "",
  "## Provenance",
  "",
  paste("The script and principal outputs are registered as an scWorkbench",
    "analysis capsule with content-bound input and output identities."))
writeLines(report, file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
log_event("Completed 30-cell Bdev RNA/ATAC distance and slope report at",
  output_dir)
