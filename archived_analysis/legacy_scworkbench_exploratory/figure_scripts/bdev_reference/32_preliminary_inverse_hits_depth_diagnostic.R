#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 32_preliminary_inverse_hits_depth_diagnostic.R",
  "<bdev_bundle> <window_report> <smoothed_report> <turn_pattern_report>",
  "<gene_screen_report> <output_dir> <expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
window_report <- normalizePath(args[[2L]], mustWork = TRUE)
smoothed_report <- normalizePath(args[[3L]], mustWork = TRUE)
turn_pattern_report <- normalizePath(args[[4L]], mustWork = TRUE)
gene_screen_report <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
expected_fingerprint <- args[[7L]]
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output already exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(Matrix)
  library(matrixStats)
  library(ggplot2)
  library(jsonlite)
  library(splines)
})

read_tsv <- function(path) read.delim(path, check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = "NA")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

trajectory <- read_tsv(file.path(window_report, "trajectory_cells_ordered.tsv"))
windows <- read_tsv(file.path(smoothed_report, "smoothed_slope_windows.tsv"))
membership_connection <- gzfile(file.path(window_report,
  "window_membership.tsv.gz"), open = "rt")
membership <- read_tsv(membership_connection)
close(membership_connection)
all_event_definition <- read_tsv(file.path(turn_pattern_report,
  "turn_event_definition.tsv"))
event_definition <- all_event_definition[match(paste0("T", 1:3),
  all_event_definition$region_id), ]
inverse <- read_tsv(file.path(gene_screen_report,
  "turn_local_balanced_top100_inverse.tsv"))
inverse20 <- head(inverse[order(inverse$inverse_balanced_rank,
  inverse$feature_id), ], 20L)
stopifnot(nrow(trajectory) == 12642L, nrow(windows) == 2524L,
  !anyDuplicated(trajectory$cell_id), all(table(membership$window_id) == 30L),
  identical(event_definition$region_id, paste0("T", 1:3)),
  nrow(inverse20) == 20L)

cell_index <- match(membership$cell_id, trajectory$cell_id)
window_index <- match(membership$window_id, windows$window_id)
stopifnot(!anyNA(cell_index), !anyNA(window_index))
window_design <- sparseMatrix(i = cell_index, j = window_index,
  x = rep(1 / 30, nrow(membership)),
  dims = c(nrow(trajectory), nrow(windows)))
membership_split <- split(cell_index, factor(membership$window_id,
  levels = windows$window_id))

smoothing_fraction <- 0.08
smoothing_half_width <- floor(nrow(windows) * smoothing_fraction / 2)
smoothing_rows <- unlist(lapply(seq_len(nrow(windows)), function(j)
  seq.int(max(1L, j - smoothing_half_width),
    min(nrow(windows), j + smoothing_half_width))), use.names = FALSE)
smoothing_columns <- rep(seq_len(nrow(windows)),
  vapply(seq_len(nrow(windows)), function(j)
    min(nrow(windows), j + smoothing_half_width) -
      max(1L, j - smoothing_half_width) + 1L, numeric(1)))
smoothing_counts <- tabulate(smoothing_columns, nbins = nrow(windows))
smoothing_design <- sparseMatrix(i = smoothing_rows, j = smoothing_columns,
  x = 1 / smoothing_counts[smoothing_columns],
  dims = c(nrow(windows), nrow(windows)))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint))
metadata <- wb_meta(bundle, cells = trajectory$cell_id,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt", "nCount_RNA",
    "nCount_peaks"))
stopifnot(identical(as.character(metadata$cell_id), trajectory$cell_id),
  isTRUE(all.equal(as.numeric(metadata$pt), trajectory$pt,
    check.attributes = FALSE)),
  all(is.finite(metadata$nCount_RNA) & metadata$nCount_RNA > 0),
  all(is.finite(metadata$nCount_peaks) & metadata$nCount_peaks > 0),
  all(inverse20$feature_id %in% wb_feature_ids(bundle, "SCT_data")))

# Log totals stabilize the very different count scales. Means and exact
# cell-level quartiles are first calculated in each original 30-cell window,
# then passed through the identical fixed 8%-window smoother used for genes.
total_cell <- rbind(total_RNA = log1p(as.numeric(metadata$nCount_RNA)),
  total_ATAC = log1p(as.numeric(metadata$nCount_peaks)))
total_raw_mean <- as.matrix(total_cell %*% window_design)
total_raw_q25 <- total_raw_q75 <- matrix(NA_real_, nrow = 2L,
  ncol = nrow(windows), dimnames = list(rownames(total_cell), NULL))
for (j in seq_len(nrow(windows))) {
  q <- rowQuantiles(total_cell[, membership_split[[j]], drop = FALSE],
    probs = c(0.25, 0.75), drop = FALSE)
  total_raw_q25[, j] <- q[, 1L]
  total_raw_q75[, j] <- q[, 2L]
}
total_mean <- as.matrix(total_raw_mean %*% smoothing_design)
total_q25 <- as.matrix(total_raw_q25 %*% smoothing_design)
total_q75 <- as.matrix(total_raw_q75 %*% smoothing_design)
rownames(total_mean) <- rownames(total_q25) <- rownames(total_q75) <-
  rownames(total_cell)

genes <- inverse20$feature_id
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", genes,
  trajectory$cell_id, missing = "NA"))
stopifnot(identical(rownames(expression), genes),
  identical(colnames(expression), trajectory$cell_id),
  all(is.finite(expression)))
gene_raw_mean <- as.matrix(expression %*% window_design)
gene_mean <- as.matrix(gene_raw_mean %*% smoothing_design)
gene_raw_q25 <- gene_raw_q75 <- matrix(NA_real_, nrow = length(genes),
  ncol = nrow(windows), dimnames = list(genes, NULL))
for (j in seq_len(nrow(windows))) {
  q <- rowQuantiles(expression[, membership_split[[j]], drop = FALSE],
    probs = c(0.25, 0.75), drop = FALSE)
  gene_raw_q25[, j] <- q[, 1L]
  gene_raw_q75[, j] <- q[, 2L]
}
gene_q25 <- as.matrix(gene_raw_q25 %*% smoothing_design)
gene_q75 <- as.matrix(gene_raw_q75 %*% smoothing_design)

highlight_regions <- data.frame(region_id = event_definition$region_id,
  xmin = windows$pt_median[event_definition$left_window_index],
  xmax = windows$pt_median[event_definition$right_window_index])
scored_indices <- sort(unique(unlist(lapply(seq_len(nrow(event_definition)),
  function(k) seq.int(event_definition$left_window_index[[k]],
    event_definition$right_window_index[[k]])))))

total_plot_data <- do.call(rbind, lapply(seq_len(2L), function(i)
  data.frame(pt = windows$pt_median,
    signal = factor(c("Total RNA (nCount_RNA)", "Total ATAC (nCount_peaks)")[[i]],
      levels = c("Total RNA (nCount_RNA)", "Total ATAC (nCount_peaks)")),
    mean = total_mean[i, ], q25 = total_q25[i, ], q75 = total_q75[i, ])))
total_plot <- ggplot(total_plot_data, aes(pt)) +
  geom_rect(data = highlight_regions,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = region_id),
    inherit.aes = FALSE, alpha = 0.10, color = NA) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#6A3D9A",
    alpha = 0.18, color = NA) +
  geom_line(aes(y = mean), color = "#542788", linewidth = 0.9) +
  geom_vline(data = all_event_definition, aes(xintercept = pt),
    inherit.aes = FALSE, color = "#555555", linetype = "dashed",
    linewidth = 0.3) +
  facet_wrap(~ signal, ncol = 1L, scales = "free_y") +
  scale_fill_manual(values = c("T1" = "#E41A1C", "T2" = "#377EB8",
    "T3" = "#4DAF4A"), name = "Scored ±100 windows") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
    strip.text = element_text(face = "bold"), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15)) +
  labs(title = "Total RNA and ATAC signal across Bdev pseudotime",
    subtitle = paste("Smoothed 30-cell means with cell-level IQR; full range",
      "shown, shaded ±100-window T1-T3 regions; T4 is unscored context"),
    x = "Bdev pseudotime", y = "log1p total count")
ggsave(file.path(stage_dir, "total_rna_atac_over_pseudotime_with_iqr.png"),
  total_plot, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "total_rna_atac_over_pseudotime_with_iqr.pdf"),
  total_plot, width = 10, height = 7, device = cairo_pdf, bg = "white")

standardize_rows <- function(x) {
  center <- rowMeans(x)
  spread <- sqrt(rowSums((x - center)^2) / (ncol(x) - 1L))
  stopifnot(all(is.finite(spread) & spread > 0))
  list(value = (x - center) / spread, center = center, spread = spread)
}
gene_z <- standardize_rows(gene_mean)
rna_z <- as.numeric(scale(total_mean["total_RNA", ]))
atac_z <- as.numeric(scale(total_mean["total_ATAC", ]))
gene_q25_z <- (gene_q25 - gene_z$center) / gene_z$spread
gene_q75_z <- (gene_q75 - gene_z$center) / gene_z$spread

correlation_rows <- lapply(seq_along(genes), function(i) {
  gene <- genes[[i]]
  g <- gene_mean[gene, ]
  # Residual correlations test depth association left after a flexible common
  # pseudotime trend is removed at the cell level.
  design <- ns(trajectory$pt, df = 8)
  gene_residual <- residuals(lm(expression[gene, ] ~ design))
  rna_residual <- residuals(lm(total_cell["total_RNA", ] ~ design))
  atac_residual <- residuals(lm(total_cell["total_ATAC", ] ~ design))
  data.frame(inverse_rank = inverse20$inverse_balanced_rank[[i]],
    feature_id = gene, ribosomal = grepl("^Rp[sl]", gene),
    full_profile_rho_total_RNA = cor(g, total_mean["total_RNA", ],
      method = "spearman"),
    full_profile_rho_total_ATAC = cor(g, total_mean["total_ATAC", ],
      method = "spearman"),
    scored_profile_rho_total_RNA = cor(g[scored_indices],
      total_mean["total_RNA", scored_indices], method = "spearman"),
    scored_profile_rho_total_ATAC = cor(g[scored_indices],
      total_mean["total_ATAC", scored_indices], method = "spearman"),
    cell_rho_total_RNA = cor(expression[gene, ], total_cell["total_RNA", ],
      method = "spearman"),
    cell_rho_total_ATAC = cor(expression[gene, ], total_cell["total_ATAC", ],
      method = "spearman"),
    pt_residual_cor_total_RNA = cor(gene_residual, rna_residual),
    pt_residual_cor_total_ATAC = cor(gene_residual, atac_residual),
    stringsAsFactors = FALSE)
})
correlations <- do.call(rbind, correlation_rows)
write_tsv(correlations, file.path(stage_dir,
  "inverse_top20_total_signal_correlations.tsv"))

overlay_data <- do.call(rbind, lapply(seq_along(genes), function(i) {
  gene <- genes[[i]]
  stats <- correlations[i, ]
  panel <- sprintf("%02d %s%s | scored rho RNA %.2f ATAC %.2f",
    inverse20$inverse_balanced_rank[[i]], gene,
    if (stats$ribosomal) " [RPL/RPS]" else "",
    stats$scored_profile_rho_total_RNA,
    stats$scored_profile_rho_total_ATAC)
  data.frame(pt = windows$pt_median, gene_panel = panel,
    gene_z = gene_z$value[gene, ], q25_z = gene_q25_z[gene, ],
    q75_z = gene_q75_z[gene, ], total_RNA_z = rna_z,
    total_ATAC_z = atac_z, stringsAsFactors = FALSE)
}))
overlay_data$gene_panel <- factor(overlay_data$gene_panel,
  levels = unique(overlay_data$gene_panel))
overlay_plot <- ggplot(overlay_data, aes(pt)) +
  geom_rect(data = highlight_regions,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = region_id),
    inherit.aes = FALSE, alpha = 0.10, color = NA) +
  geom_hline(yintercept = 0, color = "#999999", linetype = "dotted",
    linewidth = 0.25) +
  geom_ribbon(aes(ymin = q25_z, ymax = q75_z), fill = "#2166AC",
    alpha = 0.14, color = NA) +
  geom_line(aes(y = gene_z, color = "Gene SCT_data"), linewidth = 0.75) +
  geom_line(aes(y = total_RNA_z, color = "Total RNA"), linewidth = 0.55) +
  geom_line(aes(y = total_ATAC_z, color = "Total ATAC"), linewidth = 0.55) +
  geom_vline(data = all_event_definition, aes(xintercept = pt),
    inherit.aes = FALSE, color = "#555555", linetype = "dashed",
    linewidth = 0.25) +
  facet_wrap(~ gene_panel, ncol = 4L) +
  coord_cartesian(ylim = c(-3.5, 3.5)) +
  scale_color_manual(values = c("Gene SCT_data" = "#2166AC",
    "Total RNA" = "#D95F02", "Total ATAC" = "#1B9E77"), name = NULL) +
  scale_fill_manual(values = c("T1" = "#E41A1C", "T2" = "#377EB8",
    "T3" = "#4DAF4A"), name = "Scored ±100 windows") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 7.2),
    panel.spacing = grid::unit(5, "mm"), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15)) +
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  labs(title = "Inverse top 20 compared with total RNA and ATAC signal",
    subtitle = paste("All curves standardized over the full trajectory;",
      "gene ribbon is within-window IQR; labels report scored-region Spearman rho"),
    x = "Bdev pseudotime", y = "Standardized profile")
ggsave(file.path(stage_dir,
  "inverse_top20_with_total_rna_atac_depth.png"), overlay_plot,
  width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(stage_dir,
  "inverse_top20_with_total_rna_atac_depth.pdf"), overlay_plot,
  width = 14, height = 12, device = cairo_pdf, bg = "white")

window_totals <- data.frame(window_id = windows$window_id,
  window_index = windows$window_index, pt_median = windows$pt_median,
  log1p_total_RNA_mean = total_mean["total_RNA", ],
  log1p_total_RNA_q25 = total_q25["total_RNA", ],
  log1p_total_RNA_q75 = total_q75["total_RNA", ],
  log1p_total_ATAC_mean = total_mean["total_ATAC", ],
  log1p_total_ATAC_q25 = total_q25["total_ATAC", ],
  log1p_total_ATAC_q75 = total_q75["total_ATAC", ])
write_tsv(window_totals, file.path(stage_dir,
  "smoothed_window_total_rna_atac.tsv"))

ribosomal <- correlations[correlations$ribosomal, ]
non_ribosomal <- correlations[!correlations$ribosomal, ]
diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  trajectory_cells = nrow(trajectory), windows = nrow(windows),
  window_size = 30L, stride = 5L, smoothing_fraction = smoothing_fraction,
  plotted_genes = as.list(genes), ribosomal_genes = as.list(ribosomal$feature_id),
  total_RNA_field = "nCount_RNA", total_ATAC_field = "nCount_peaks",
  transform = "log1p before window summary and smoothing",
  scored_regions = "T1-T3 +/-100 windows", displayed_range = "entire pseudotime",
  ribosomal_median_abs_scored_rho_RNA = median(abs(
    ribosomal$scored_profile_rho_total_RNA)),
  ribosomal_median_abs_scored_rho_ATAC = median(abs(
    ribosomal$scored_profile_rho_total_ATAC)),
  non_ribosomal_median_abs_scored_rho_RNA = median(abs(
    non_ribosomal$scored_profile_rho_total_RNA)),
  non_ribosomal_median_abs_scored_rho_ATAC = median(abs(
    non_ribosomal$scored_profile_rho_total_ATAC)),
  bundle_mutated = FALSE, capsule_registered = FALSE)
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
report <- c("# Inverse-hit RNA/ATAC depth diagnostic", "",
  paste("The total RNA (`nCount_RNA`) and total ATAC (`nCount_peaks`)",
    "profiles use the same 30-cell windows and fixed 8%-window smoother as",
    "the turn-local gene figures. The complete pseudotime range is shown;",
    "only the fixed +/-100-window T1-T3 regions are shaded."), "",
  paste0("The inverse top 20 contain ", nrow(ribosomal),
    " Rpl/Rps genes: ", paste(ribosomal$feature_id, collapse = ", "), "."),
  paste0("Their median absolute scored-region profile correlation is ",
    sprintf("%.3f", diagnostics$ribosomal_median_abs_scored_rho_RNA),
    " with total RNA and ",
    sprintf("%.3f", diagnostics$ribosomal_median_abs_scored_rho_ATAC),
    " with total ATAC."),
  paste0("For the non-ribosomal inverse hits, the corresponding medians are ",
    sprintf("%.3f", diagnostics$non_ribosomal_median_abs_scored_rho_RNA),
    " and ",
    sprintf("%.3f", diagnostics$non_ribosomal_median_abs_scored_rho_ATAC), "."),
  "", paste("Profile correlation with library size is evidence compatible",
    "with a depth/compositional effect, but is not by itself proof of a",
    "technical artifact. The table also provides cell-level correlations and",
    "correlations after removing an 8-df pseudotime spline."))
writeLines(report, file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
files <- list.files(stage_dir, full.names = TRUE)
write_tsv(data.frame(file = basename(files),
  size_bytes = unname(file.info(files)$size)),
  file.path(stage_dir, "output_manifest.tsv"))

wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output directory")
cat("Completed inverse-hit depth diagnostic at", output_dir, "\n")
