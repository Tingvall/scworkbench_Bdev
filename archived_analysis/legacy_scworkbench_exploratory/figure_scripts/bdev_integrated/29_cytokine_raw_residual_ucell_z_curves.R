#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 29_cytokine_raw_residual_ucell_z_curves.R",
  "<bundle_dir> <native_ucell_report> <dynamic_ucell_report>",
  "<z_cytokine_dir> <bb.tsv> <output_dir>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
native_dir <- normalizePath(args[[2L]], mustWork = TRUE)
dynamic_dir <- normalizePath(args[[3L]], mustWork = TRUE)
z_dir <- normalizePath(args[[4L]], mustWork = TRUE)
bb_path <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
short_label <- function(x) {
  x <- sub("^parse_cytokine_rest_only_refinement_all_top100_v1\\.", "", x)
  x <- gsub("\\.", " ", x)
  toupper(x)
}

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint

pt_table <- fread(file.path(native_dir, "source_pseudotime.tsv"))
setnames(pt_table, c("cell_id", "pt"))
setorder(pt_table, cell_id)
stopifnot(nrow(pt_table) == 12642L, !anyDuplicated(pt_table$cell_id),
  all(is.finite(pt_table$pt)))
metadata <- as.data.table(wb_meta(bundle, cells = pt_table$cell_id))
stopifnot(identical(as.character(metadata$cell_id), pt_table$cell_id),
  all(is.finite(metadata$nFeature_RNA)))

scores_long <- read_tsv_gz(file.path(native_dir,
  "ucell_cell_scores.tsv.gz"))[cell_id %chin% pt_table$cell_id,
    .(feature_set = key, cell_id, raw_ucell = value)]
stopifnot(uniqueN(scores_long$feature_set) == 38L,
  nrow(scores_long) == 38L * nrow(pt_table),
  !anyDuplicated(scores_long, by = c("feature_set", "cell_id")))
scores_wide <- dcast(scores_long, cell_id ~ feature_set,
  value.var = "raw_ucell")
setorder(scores_wide, cell_id)
stopifnot(identical(scores_wide$cell_id, pt_table$cell_id))
feature_sets <- sort(setdiff(names(scores_wide), "cell_id"))
score_matrix <- t(as.matrix(scores_wide[, ..feature_sets]))
rownames(score_matrix) <- feature_sets

pt <- pt_table$pt
complexity_x <- log1p(as.numeric(metadata$nFeature_RNA))
residual_matrix <- matrix(NA_real_, nrow(score_matrix), ncol(score_matrix),
  dimnames = dimnames(score_matrix))
adjustment_qc <- rbindlist(lapply(feature_sets, function(fs) {
  y <- as.numeric(score_matrix[fs, ])
  fit <- lm(y ~ splines::ns(complexity_x, df = 4L))
  residual_matrix[fs, ] <<- residuals(fit)
  data.table(feature_set = fs, short_name = short_label(fs),
    model = "lm(UCell ~ ns(log1p(nFeature_RNA), df=4))",
    model_r_squared = summary(fit)$r.squared,
    raw_spearman_nFeature_RNA = cor(y, metadata$nFeature_RNA,
      method = "spearman"),
    residual_spearman_nFeature_RNA = cor(residuals(fit),
      metadata$nFeature_RNA, method = "spearman"))
}))

pt_limits <- range(pt)
pt_span <- diff(pt_limits)
window_width <- .05 * pt_span
centers <- pt_limits[[1L]] + seq(0, 1, by = .01) * pt_span
windows <- data.table(grid_window_index = seq_along(centers),
  window_center = centers,
  window_left = pmax(pt_limits[[1L]], centers - window_width / 2),
  window_right = pmin(pt_limits[[2L]], centers + window_width / 2))
membership_index <- lapply(seq_len(nrow(windows)), function(i)
  which(pt >= windows$window_left[[i]] & pt <= windows$window_right[[i]]))
windows[, n_cells := lengths(membership_index)]
keep <- windows$n_cells >= 20L
windows <- windows[keep]
membership_index <- membership_index[keep]
windows[, window_index := seq_len(.N)]
stopifnot(nrow(windows) == 97L)

calculated <- rbindlist(lapply(feature_sets, function(fs) {
  data.table(feature_set = fs, windows,
    raw_ucell = vapply(membership_index, function(ii)
      median(score_matrix[fs, ii]), numeric(1)),
    residual_ucell = vapply(membership_index, function(ii)
      median(residual_matrix[fs, ii]), numeric(1)))
}))
stored_raw <- read_tsv_gz(file.path(dynamic_dir, "window_curves.tsv.gz"))
raw_check <- merge(calculated[, .(feature_set, grid_window_index,
    reconstructed_raw = raw_ucell)],
  stored_raw[, .(feature_set, grid_window_index,
    stored_raw = ucell_median)], by = c("feature_set", "grid_window_index"))
max_raw_difference <- max(abs(raw_check$reconstructed_raw -
  raw_check$stored_raw))
stopifnot(nrow(raw_check) == 38L * 97L, max_raw_difference < 1e-12)

z <- read_tsv_gz(file.path(z_dir, "window_curves.tsv.gz"))
stopifnot(setequal(unique(z$feature_set), feature_sets),
  nrow(z) == 38L * 97L)
grid_check <- merge(unique(calculated[, .(grid_window_index,
    ucell_center = window_center)]),
  unique(z[, .(grid_window_index, z_center = window_center)]),
  by = "grid_window_index")
stopifnot(nrow(grid_check) == 97L,
  max(abs(grid_check$ucell_center - grid_check$z_center)) < 1e-10)

curves <- rbindlist(list(
  calculated[, .(feature_set, window_index, window_center,
    method = "Raw UCell", original_value = raw_ucell)],
  calculated[, .(feature_set, window_index, window_center,
    method = "RNA-feature residualized UCell",
    original_value = residual_ucell)],
  z[, .(feature_set, window_index, window_center,
    method = "Original expression z-score", original_value = z_score)]))
curves[, scaled_0_1 := {
  limits <- range(original_value)
  if (diff(limits) > 0) (original_value - limits[[1L]]) / diff(limits) else 0
}, by = .(feature_set, method)]
stopifnot(min(curves$scaled_0_1) >= 0, max(curves$scaled_0_1) <= 1)

# Keep the established cytokine UCell k-means order so related profiles stay
# near one another, reading left-to-right and then top-to-bottom.
membership_path <- file.path(dirname(z_dir), "summary_figures",
  "cytokine_ucell_standalone_comparison",
  "cytokine_ucell_kmeans_membership.tsv")
membership <- fread(membership_path)
setorder(membership, row_position_top_to_bottom)
display_order <- membership$feature_set
stopifnot(setequal(display_order, feature_sets))
curves[, set_label := factor(short_label(feature_set),
  levels = short_label(display_order))]
curves[, method := factor(method, levels = c("Raw UCell",
  "RNA-feature residualized UCell", "Original expression z-score"))]
bb <- fread(bb_path)

p <- ggplot(curves, aes(window_center, scaled_0_1, color = method,
    linetype = method)) +
  geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
    linetype = "dotted", linewidth = .28, color = "grey35") +
  geom_line(linewidth = .45) +
  facet_wrap(~set_label, ncol = 5) +
  scale_color_manual(values = c("Raw UCell" = "#252525",
    "RNA-feature residualized UCell" = "#D7301F",
    "Original expression z-score" = "#2166AC")) +
  scale_linetype_manual(values = c("Raw UCell" = "solid",
    "RNA-feature residualized UCell" = "solid",
    "Original expression z-score" = "longdash")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1),
    expand = expansion(mult = c(.01, .03))) +
  labs(title = "Cytokine signals: raw UCell, RNA-complexity residuals, and expression z-score",
    subtitle = paste("Each cytokine/method curve independently min-max scaled",
      "to 0-1; dotted lines mark BB1 and BB2"),
    x = "Existing Bdev pseudotime (pt)",
    y = "Within-curve scaled value (0-1)", color = NULL,
    linetype = NULL,
    caption = "Scaling is visualization-only; existing permutation calls were not changed.") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
    strip.text = element_text(size = 6.5),
    plot.title = element_text(face = "bold"))

png_path <- file.path(stage_dir,
  "cytokine_raw_residual_ucell_original_z_scaled_0_1_curves.png")
pdf_path <- file.path(stage_dir,
  "cytokine_raw_residual_ucell_original_z_scaled_0_1_curves.pdf")
table_path <- file.path(stage_dir,
  "cytokine_raw_residual_ucell_original_z_scaled_0_1_curves.tsv.gz")
ggsave(png_path, p, width = 17, height = 21, dpi = 240, bg = "white")
ggsave(pdf_path, p, width = 17, height = 21, bg = "white")
fwrite(curves, table_path, sep = "\t", quote = FALSE, compress = "gzip")
fwrite(adjustment_qc, file.path(stage_dir,
  "cytokine_complexity_adjustment_qc.tsv"), sep = "\t", quote = FALSE)

median_raw_cor <- median(abs(adjustment_qc$raw_spearman_nFeature_RNA))
median_residual_cor <- median(abs(
  adjustment_qc$residual_spearman_nFeature_RNA))
writeLines(c("# Cytokine raw/residual-UCell/original-z curve comparison", "",
  "- All 38 cytokine signals are shown in the established cytokine UCell k-means order.",
  "- Each gene-set/method curve was independently min-max scaled to 0-1 for shape comparison only.",
  "- Residualization model: `UCell ~ ns(log1p(nFeature_RNA), df=4)`, fitted independently per cytokine without pseudotime.",
  "- Existing permutation calls were not recalculated or changed.",
  paste0("- Reconstructed raw-UCell curves matched the stored curves with maximum absolute difference ",
    format(max_raw_difference, scientific = TRUE, digits = 3), "."),
  paste0("- Median absolute Spearman correlation with nFeature_RNA changed from ",
    signif(median_raw_cor, 3), " to ", signif(median_residual_cor, 3), ".")),
  file.path(stage_dir, "REPORT.md"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
write_json(list(status = "complete", cytokine_sets = 38L,
  adjustment = "lm(UCell ~ ns(log1p(nFeature_RNA), df=4))",
  scaling = paste("independent min-max 0-1 within each",
    "cytokine/method curve; visualization only"),
  raw_curve_maximum_absolute_difference = max_raw_difference,
  median_absolute_raw_spearman_nFeature_RNA = median_raw_cor,
  median_absolute_residual_spearman_nFeature_RNA = median_residual_cor,
  permutation_calls_changed = FALSE, bundle_unchanged = TRUE),
  file.path(stage_dir, "completion.json"), pretty = TRUE,
  auto_unbox = TRUE, digits = NA)

wb_close_bundle(bundle)
stopifnot(file.rename(stage_dir, output_dir))
message("Completed cytokine three-method curve comparison at: ", output_dir)
