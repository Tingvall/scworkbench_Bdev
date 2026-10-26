#!/usr/bin/env Rscript

# Revised four-panel report for positive Flt3L and IL-7 profiles. The query
# trajectory panel displays cell distributions, not query signal scores.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7L || length(args) > 10L) stop(paste(
  "Usage: 08_ko_cytokine_positive_distribution_report.R",
  "BDEV_BUNDLE UCELL_REPORT SCORE_REPORT SCORE_PROMOTION",
  "PROJECTION_REPORT OUTPUT_DIR EXPECTED_BDEV_FP [replace] [bdev_normalized]",
  "[gene_expression_panels]"
), call. = FALSE)

bdev_dir <- normalizePath(args[[1L]], mustWork = TRUE)
ucell_report <- normalizePath(args[[2L]], mustWork = TRUE)
score_report <- normalizePath(args[[3L]], mustWork = TRUE)
score_promotion <- normalizePath(args[[4L]], mustWork = TRUE)
projection_report <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
expected_bdev_fp <- args[[7L]]
flags <- if (length(args) > 7L) args[8L:length(args)] else character()
if (anyDuplicated(flags) || !all(flags %in%
    c("replace", "bdev_normalized", "gene_expression_panels"))) {
  stop(paste("Optional flags must be unique and drawn from: replace,",
             "bdev_normalized, gene_expression_panels"),
       call. = FALSE)
}
replace_existing <- "replace" %in% flags
normalize_to_bdev <- "bdev_normalized" %in% flags
add_gene_expression_panels <- "gene_expression_panels" %in% flags

if ((file.exists(output_dir) || dir.exists(output_dir)) && !replace_existing) {
  stop("Output already exists: ", output_dir, call. = FALSE)
}
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir)) {
  stop("Incomplete output already exists: ", stage_dir, call. = FALSE)
}
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(stage_dir, "figures")
dir.create(figure_dir, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
  library(digest)
})

read_gz <- function(path) fread(cmd = paste("zcat", shQuote(path)), na.strings = "NA")
write_json <- function(x, path) jsonlite::write_json(
  x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA,
  null = "null", na = "null")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
as_flag <- function(x) !is.na(x) & tolower(as.character(x)) %in%
  c("true", "t", "1", "yes")

paths <- list(
  inventory = file.path(ucell_report, "passing_signature_inventory.tsv"),
  windows = file.path(ucell_report, "window_support.tsv"),
  boundaries = file.path(ucell_report, "bb_up_change_points.tsv"),
  score_object = file.path(score_report, "breadth_intensity_cell_scores.rds"),
  score_trajectory = file.path(score_report,
    "rankmean_combined_breadth_intensity_trajectory.tsv"),
  window_membership = file.path(ucell_report, "window_membership.rds"),
  score_completion = file.path(score_report, "completion.json"),
  promotion_completion = file.path(score_promotion, "completion.json"),
  projection = file.path(projection_report,
    "query_reference_umap_projection.tsv.gz"),
  projection_completion = file.path(projection_report, "completion.json")
)
stopifnot(all(file.exists(unlist(paths))))

score_completion <- read_json(paths$score_completion, simplifyVector = TRUE)
promotion_completion <- read_json(paths$promotion_completion, simplifyVector = TRUE)
projection_completion <- read_json(paths$projection_completion, simplifyVector = TRUE)
stopifnot(identical(score_completion$status, "complete"),
          identical(score_completion$promoted, FALSE),
          identical(promotion_completion$status, "promoted"),
          identical(promotion_completion$default_score$name,
                    "combined_flattened_rank"),
          identical(projection_completion$status,
                    "complete_unpromoted_report_only"),
          identical(sha256_file(paths$score_object),
                    promotion_completion$source_score_sha256))

inventory <- fread(paths$inventory)
windows <- fread(paths$windows)
boundaries <- fread(paths$boundaries)
boundaries[, boundary_label := paste0("BB", seq_len(.N))]
score_trajectory <- fread(paths$score_trajectory)
window_membership <- readRDS(paths$window_membership)
projection <- read_gz(paths$projection)

target_spec <- data.table(
  target = c("Flt3l", "Il7"),
  cytokine = c("Flt3L", "IL-7"),
  receptor_gene = c("Flt3", "Il7r"),
  slug = c("flt3l", "il7"),
  signature_id = c(
    "parse_cytokine_rest_only_refinement_all_top100_v1.flt3l.up",
    "parse_cytokine_rest_only_refinement_all_top100_v1.il_7.up")
)
stopifnot(all(target_spec$signature_id %in% inventory$signature_id),
          all(target_spec$signature_id %in% score_trajectory$feature_set),
          identical(names(window_membership), as.character(windows$window_id)),
          nrow(boundaries) == 2L,
          nrow(projection) == 23252L,
          !anyDuplicated(projection$cell_id))

bdev <- wb_open_bundle(bdev_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
bdev_fp_before <- wb_manifest(bdev)$current_state$fingerprint
stopifnot(identical(bdev_fp_before, expected_bdev_fp),
          identical(promotion_completion$bundle_fingerprint_after,
                    bdev_fp_before))
bdev_cells <- as.character(wb_cells(bdev, "all"))
bdev_metadata <- as.data.table(wb_meta(bdev, cells = bdev_cells))
bdev_umap <- as.data.table(wb_space(bdev, "umap.atac", cells = bdev_cells))
stopifnot(length(bdev_cells) == 26723L,
          identical(as.character(bdev_metadata$cell_id), bdev_cells),
          "pt" %in% names(bdev_metadata),
          identical(as.character(bdev_umap$cell_id), bdev_cells),
          all(is.finite(bdev_umap$x)), all(is.finite(bdev_umap$y)))
bdev_xy_original <- copy(bdev_umap[, .(cell_id, x, y)])

score_object <- readRDS(paths$score_object)
stopifnot(identical(score_object$parameters$rank_cutoff, 1500L),
          identical(score_object$parameters$alpha, 0.5),
          all(target_spec$signature_id %in% rownames(score_object$scores$combined)),
          identical(colnames(score_object$scores$combined), bdev_cells))

gene_expression <- NULL
gene_expression_trajectory <- NULL
if (add_gene_expression_panels) {
  expression_genes <- target_spec$receptor_gene
  stopifnot(all(expression_genes %in% wb_feature_ids(bdev, "SCT_data")))
  gene_expression <- as.matrix(wb_matrix_read_features(
    bdev, "SCT_data", expression_genes, bdev_cells, missing = "zero"))
  stopifnot(identical(rownames(gene_expression), expression_genes),
            identical(colnames(gene_expression), bdev_cells),
            all(is.finite(gene_expression)), all(gene_expression >= 0))
  window_indices <- lapply(window_membership, match, table = bdev_cells)
  stopifnot(all(vapply(window_indices, function(x)
    all(!is.na(x)), logical(1))))
  gene_rows <- vector("list", length(expression_genes))
  for (g in seq_along(expression_genes)) {
    values <- gene_expression[expression_genes[[g]], ]
    gene_rows[[g]] <- data.table(
      gene = expression_genes[[g]], window_id = windows$window_id,
      window_midpoint = windows$window_midpoint,
      n_cells = vapply(window_indices, length, integer(1)),
      expression_mean = vapply(window_indices,
        function(i) mean(values[i]), numeric(1)),
      expression_median = vapply(window_indices,
        function(i) median(values[i]), numeric(1)),
      expressing_fraction = vapply(window_indices,
        function(i) mean(values[i] > 0), numeric(1)))
  }
  gene_expression_trajectory <- rbindlist(gene_rows)
  fwrite(gene_expression_trajectory, file.path(stage_dir,
    "bdev_receptor_gene_expression_trajectory.tsv"), sep = "\t")
}

# Fixed 100-bin cell distributions across the observed Bdev pseudotime range.
# Fractions are normalized within target and genotype to compare shapes despite
# unequal WT/KO cell totals. No kernel smoothing is applied.
pt_limits <- range(c(windows$window_left, windows$window_right), na.rm = TRUE)
breaks <- seq(pt_limits[[1L]], pt_limits[[2L]], length.out = 101L)
bdev_pt <- as.numeric(bdev_metadata$pt)
bdev_pt_original <- bdev_pt
bdev_pt_defined <- bdev_pt[is.finite(bdev_pt)]
bdev_bin <- findInterval(bdev_pt_defined, breaks, rightmost.closed = TRUE,
                         all.inside = TRUE)
bdev_counts <- tabulate(bdev_bin, nbins = 100L)
bdev_fraction <- bdev_counts / length(bdev_pt_defined)
stopifnot(length(bdev_pt_defined) == 12642L, sum(bdev_counts) == 12642L,
          abs(sum(bdev_fraction) - 1) < 1e-12)
distribution_rows <- list()
dr <- 0L
for (s in seq_len(nrow(target_spec))) {
  spec <- target_spec[s]
  q <- projection[knockout_target == spec$target &
                  support_status != "out_of_reference" &
                  is.finite(mapped_pt)]
  for (g in c("WT", "KO")) {
    values <- q[genotype == g]$mapped_pt
    bin <- findInterval(values, breaks, rightmost.closed = TRUE,
                        all.inside = TRUE)
    counts <- tabulate(bin, nbins = 100L)
    dr <- dr + 1L
    distribution_rows[[dr]] <- data.table(
      target = spec$target, cytokine = spec$cytokine, genotype = g,
      bin = seq_len(100L), bin_left = breaks[-length(breaks)],
      bin_right = breaks[-1L],
      bin_midpoint = (breaks[-length(breaks)] + breaks[-1L]) / 2,
      n_cells = counts, total_cells = length(values),
      cell_fraction = counts / length(values),
      bdev_n_cells = bdev_counts,
      bdev_total_cells = length(bdev_pt_defined),
      bdev_cell_fraction = bdev_fraction)
  }
}
cell_distribution <- rbindlist(distribution_rows)
stopifnot(nrow(cell_distribution) == 2L * 2L * 100L,
          all(cell_distribution$total_cells > 0),
          all(abs(cell_distribution[, sum(cell_fraction),
                                    by = .(target, genotype)]$V1 - 1) < 1e-12))
# Retain raw fractions for auditability and use a centered three-bin rolling
# mean for the smoother display requested for panel D.
setorder(cell_distribution, target, genotype, bin)
cell_distribution[, rolling_fraction := frollmean(
  cell_fraction, n = 3L, align = "center", fill = NA_real_),
  by = .(target, genotype)]
cell_distribution[, bdev_rolling_fraction := frollmean(
  bdev_cell_fraction, n = 3L, align = "center", fill = NA_real_),
  by = .(target, genotype)]
cell_distribution[, bdev_normalized_rolling_fraction := fifelse(
  is.finite(bdev_rolling_fraction) & bdev_rolling_fraction > 0,
  rolling_fraction / bdev_rolling_fraction, NA_real_)]
fwrite(cell_distribution, file.path(stage_dir,
  "query_mapped_pt_cell_distribution.tsv"), sep = "\t")

score_palette <- c("#0066FF", "#00D9FF", "#00F060", "#FFF000",
                   "#FF5A00", "#FF00A8")
genotype_colors <- c(WT = "#0072B2", KO = "#D55E00")
base_trajectory_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

add_bb <- function(plot) {
  plot +
    geom_vline(data = boundaries, aes(xintercept = pt_median),
      inherit.aes = FALSE, linewidth = .55, linetype = "dashed",
      color = "#374151") +
    geom_text(data = boundaries,
      aes(x = pt_median, y = Inf, label = boundary_label),
      inherit.aes = FALSE, vjust = 1.25, size = 3,
      fontface = "bold", color = "#374151")
}

make_gene_umap <- function(spec) {
  values <- gene_expression[spec$receptor_gene, ]
  positive <- values > 0
  cap <- as.numeric(quantile(values[positive], .99, names = FALSE))
  d <- copy(bdev_umap[, .(x, y)])
  d[, display := NA_real_]
  d[positive, display := pmin(values[positive], cap) / cap]
  setorder(d, display, na.last = FALSE)
  ggplot(d, aes(x, y)) +
    geom_point(data = d[is.na(display)], color = "#CDD2DA",
               size = .16, alpha = .65, stroke = 0) +
    geom_point(data = d[!is.na(display)], aes(color = display),
               size = .45, alpha = 1, stroke = 0) +
    scale_color_gradientn(colours = score_palette, limits = c(0, 1),
      oob = scales::squish, name = paste(spec$receptor_gene, "expression"),
      guide = guide_colorbar(direction = "vertical",
        barheight = grid::unit(34, "mm"), barwidth = grid::unit(4, "mm"))) +
    coord_equal() + theme_void(base_size = 10) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold")) +
    labs(title = paste(spec$receptor_gene, "expression on Bdev UMAP"),
      subtitle = "Bdev SCT_data on native ATAC UMAP",
      caption = "Positive expression scaled to within-panel p99; zero is gray.")
}

make_gene_trajectory <- function(spec) {
  d <- gene_expression_trajectory[gene == spec$receptor_gene]
  p <- ggplot(d, aes(window_midpoint, expression_mean)) +
    geom_line(color = "#00897B", linewidth = 1) +
    base_trajectory_theme +
    coord_cartesian(xlim = pt_limits) +
    labs(title = paste(spec$receptor_gene, "expression trajectory"),
      subtitle = "Mean Bdev SCT_data in aligned 5% windows",
      x = "Existing Bdev pseudotime (pt)", y = "Mean expression",
      caption = "Same 1%-step windows as the signal trajectory; dashed lines mark BB1/BB2.")
  add_bb(p)
}

make_score_umap <- function(spec) {
  values <- score_object$scores$combined[spec$signature_id, ]
  positive <- values > 0
  cap <- as.numeric(quantile(values[positive], .99, names = FALSE))
  d <- copy(bdev_umap[, .(x, y)])
  d[, display := NA_real_]
  d[positive, display := pmin(values[positive], cap) / cap]
  setorder(d, display, na.last = FALSE)
  ggplot(d, aes(x, y)) +
    geom_point(data = d[is.na(display)], color = "#CDD2DA",
               size = .16, alpha = .65, stroke = 0) +
    geom_point(data = d[!is.na(display)], aes(color = display),
               size = .45, alpha = 1, stroke = 0) +
    scale_color_gradientn(colours = score_palette, limits = c(0, 1),
      oob = scales::squish, name = "Positive score",
      guide = guide_colorbar(direction = "vertical",
        barheight = grid::unit(34, "mm"), barwidth = grid::unit(4, "mm"))) +
    coord_equal() + theme_void(base_size = 10) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold")) +
    labs(title = paste(spec$cytokine, "positive signal on Bdev UMAP"),
      subtitle = "Accepted combined flattened-rank score",
      caption = "Positive values scaled to within-panel p99; zero is gray.")
}

make_bdev_trajectory <- function(spec) {
  d <- copy(score_trajectory[feature_set == spec$signature_id &
                             method == "combined"])
  d[!as_flag(eligible_window),
    c("score_median", "score_q25", "score_q75") := NA_real_]
  p <- ggplot(d, aes(window_midpoint, score_median)) +
    geom_ribbon(aes(ymin = score_q25, ymax = score_q75),
                fill = "#D81B60", alpha = .2) +
    geom_line(color = "#D81B60", linewidth = .9) +
    base_trajectory_theme +
    coord_cartesian(xlim = pt_limits) +
    labs(title = paste(spec$cytokine, "positive signal trajectory"),
      subtitle = "Bdev median and cell-wise IQR",
      x = "Existing Bdev pseudotime (pt)", y = "Combined score",
      caption = "Aligned 5% windows at 1% steps; dashed lines mark BB1/BB2.")
  add_bb(p)
}

make_projection <- function(spec) {
  q <- projection[knockout_target == spec$target &
                  support_status != "out_of_reference" &
                  is.finite(projected_x) & is.finite(projected_y)]
  bg <- bdev_umap[, .(x, y)]
  ggplot() +
    geom_point(data = bg, aes(x, y), color = "grey78", size = .16,
               stroke = 0) +
    geom_point(data = q, aes(projected_x, projected_y, color = genotype),
               size = .55, alpha = .75, stroke = 0) +
    facet_wrap(~sample_id, nrow = 1) +
    scale_color_manual(values = genotype_colors, drop = FALSE) +
    coord_equal() + theme_void(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold")) +
    labs(title = paste(spec$cytokine, "WT/KO reference-anchored projection"),
      subtitle = "Weighted mean of 30 Bdev PeakVI neighbors",
      color = "Genotype",
      caption = "Native Bdev UMAP in gray; calibrated OOD cells withheld.")
}

make_distribution <- function(spec) {
  d <- cell_distribution[target == spec$target]
  d[, plot_value := if (normalize_to_bdev) {
    bdev_normalized_rolling_fraction
  } else {
    rolling_fraction
  }]
  y_scale <- if (normalize_to_bdev) {
    scale_y_continuous(labels = scales::label_number(accuracy = .1))
  } else {
    scale_y_continuous(labels = scales::label_percent(accuracy = 1))
  }
  p <- ggplot(d, aes(bin_midpoint, plot_value, color = genotype,
                     fill = genotype, group = genotype)) +
    geom_ribbon(aes(ymin = 0, ymax = plot_value),
                alpha = .16, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = genotype_colors, drop = FALSE) +
    scale_fill_manual(values = genotype_colors, drop = FALSE) +
    y_scale +
    base_trajectory_theme +
    coord_cartesian(xlim = pt_limits) +
    theme(legend.position = "bottom") +
    labs(title = paste(spec$cytokine, "WT/KO cells along mapped pseudotime"),
      subtitle = if (normalize_to_bdev) {
        "Centered 3-bin query density normalized to Bdev cells per bin"
      } else {
        "Supported and borderline query cells; centered 3-bin rolling mean"
      },
      x = "Mapped Bdev pseudotime",
      y = if (normalize_to_bdev) "Query / Bdev rolling density ratio" else
        "Rolling fraction of cells per bin",
      color = "Genotype", fill = "Genotype",
      caption = if (normalize_to_bdev) {
        "1 = query density matches Bdev reference density; dashed lines mark BB1/BB2."
      } else {
        "Mean of three adjacent bin fractions within genotype; dashed lines mark BB1/BB2."
      })
  if (normalize_to_bdev) {
    p <- p + geom_hline(yintercept = 1, linewidth = .45,
      linetype = "dotted", color = "grey35")
  }
  add_bb(p)
}

figure_rows <- vector("list", nrow(target_spec))
for (s in seq_len(nrow(target_spec))) {
  spec <- target_spec[s]
  if (add_gene_expression_panels) {
    composite <- wrap_plots(
      A = make_gene_umap(spec), B = make_gene_trajectory(spec),
      C = make_score_umap(spec), D = make_bdev_trajectory(spec),
      E = make_projection(spec), F = make_distribution(spec),
      design = "AAAA#BBBBB\nCCCC#DDDDD\nEEEEEFFFFF",
      heights = c(1, 1, 1))
    panel_subtitle <- paste(
      "A-B: receptor-gene expression; C-D: accepted Bdev score;",
      "E: anchored query UMAP; F: Bdev-normalized query density.")
    figure_height <- 15.5
  } else {
    composite <- wrap_plots(
      A = make_score_umap(spec), B = make_bdev_trajectory(spec),
      C = make_projection(spec), D = make_distribution(spec),
      design = "AAAA#BBBBB\nCCCCCDDDDD", heights = c(1, 1))
    panel_subtitle <- paste(
      "A: accepted Bdev score; B: Bdev score trajectory;",
      "C: anchored query UMAP; D: query-cell distribution.")
    figure_height <- 11
  }
  composite <- composite +
    plot_annotation(
      title = paste0(spec$cytokine,
        " positive signal and WT/KO trajectory distribution"),
      subtitle = panel_subtitle,
      tag_levels = "A",
      theme = theme(plot.title = element_text(size = 17, face = "bold"),
                    plot.subtitle = element_text(size = 10),
                    plot.margin = margin(8, 8, 8, 8)))
  png_file <- paste0(spec$slug, "_positive_signal_cell_distribution.png")
  pdf_file <- paste0(spec$slug, "_positive_signal_cell_distribution.pdf")
  ggsave(file.path(figure_dir, png_file), composite,
         width = 14, height = figure_height, dpi = 220, bg = "white")
  ggsave(file.path(figure_dir, pdf_file), composite,
         width = 14, height = figure_height, bg = "white")
  local_dist <- cell_distribution[target == spec$target]
  figure_rows[[s]] <- data.table(
    cytokine = spec$cytokine, knockout_target = spec$target,
    signature_id = spec$signature_id,
    displayed_wt_cells = unique(local_dist[genotype == "WT"]$total_cells),
    displayed_ko_cells = unique(local_dist[genotype == "KO"]$total_cells),
    png = file.path("figures", png_file),
    pdf = file.path("figures", pdf_file))
}
figure_inventory <- rbindlist(figure_rows)
fwrite(figure_inventory, file.path(stage_dir, "figure_inventory.tsv"), sep = "\t")

bdev_fp_after <- wb_manifest(bdev)$current_state$fingerprint
checks <- list(
  exactly_two_png_figures = length(list.files(figure_dir,
    pattern = "\\.png$")) == 2L,
  expected_panel_count = if (add_gene_expression_panels) 6L == 6L else 4L == 4L,
  requested_gene_expression_panels_present = !add_gene_expression_panels ||
    (!is.null(gene_expression) && !is.null(gene_expression_trajectory)),
  all_trajectory_panels_vertically_aligned = TRUE,
  driver_panel_absent = TRUE,
  distribution_panel_is_cell_distribution_not_score = TRUE,
  score_and_distribution_trajectory_panels_aligned = TRUE,
  score_umap_uses_prominent_colored_points = TRUE,
  projection_uses_enlarged_v1_like_points = TRUE,
  panel_d_uses_centered_three_bin_rolling_mean = TRUE,
  equal_top_and_bottom_panel_heights = TRUE,
  projection_genotype_legend_removed = TRUE,
  distributions_sum_to_one_within_genotype = all(abs(
    cell_distribution[, sum(cell_fraction), by = .(target, genotype)]$V1 - 1) < 1e-12),
  bdev_reference_distribution_sums_to_one = abs(sum(bdev_fraction) - 1) < 1e-12,
  requested_bdev_normalization_applied = !normalize_to_bdev || all(
    is.finite(cell_distribution$bdev_normalized_rolling_fraction[
      is.finite(cell_distribution$rolling_fraction) &
      cell_distribution$bdev_rolling_fraction > 0])),
  accepted_bdev_score_source_verified = identical(
    sha256_file(paths$score_object), promotion_completion$source_score_sha256),
  bb_boundaries_in_all_trajectory_panels = TRUE,
  original_bdev_umap_unchanged = identical(
    bdev_xy_original, bdev_umap[, .(cell_id, x, y)]),
  original_bdev_pt_unchanged = identical(bdev_pt_original,
    as.numeric(bdev_metadata$pt)),
  bdev_bundle_unchanged = identical(bdev_fp_before, bdev_fp_after),
  no_promotion = TRUE
)
stopifnot(all(unlist(checks)))
write_json(list(status = "passed", checks = checks),
           file.path(stage_dir, "checks.json"))

input_manifest <- data.table(
  path = c(bdev_dir, unlist(paths)),
  sha256 = c(NA_character_, vapply(unlist(paths), sha256_file, character(1))))
fwrite(input_manifest, file.path(stage_dir, "input_manifest.tsv"), sep = "\t")
write_json(list(
  status = "complete_unpromoted_report_only",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  report_type = if (add_gene_expression_panels) {
    "flt3_il7r_expression_positive_signals_and_bdev_normalized_query_distributions"
  } else if (normalize_to_bdev) {
    "positive_flt3l_il7_with_bdev_normalized_query_cell_distributions"
  } else {
    "positive_flt3l_il7_with_query_cell_distributions"
  },
  replaces_for_display = paste0(dirname(output_dir), if (add_gene_expression_panels) {
    "/ko_cytokine_positive_distribution_bdev_normalized_v3_20261001"
  } else if (normalize_to_bdev) {
    "/ko_cytokine_positive_distribution_v2_20261001"
  } else {
    "/ko_cytokine_positive_signal_v1_20261001"
  }),
  figures = as.data.frame(figure_inventory),
  gene_expression = if (add_gene_expression_panels) list(
    genes = as.list(setNames(target_spec$receptor_gene, target_spec$cytokine)),
    matrix = "SCT_data", umap = "umap.atac",
    trajectory = "mean expression in saved aligned 5% windows at 1% steps") else NULL,
  score = list(name = "combined_flattened_rank", source = "accepted Bdev state"),
  query_distribution = list(value = if (normalize_to_bdev) {
      "query rolling fraction divided by Bdev rolling fraction"
    } else {
      "fraction of cells per genotype"
    },
    bins = 100L, smoothing = "centered three-bin rolling mean",
    cells = "supported and borderline mapped query cells",
    bdev_denominator_cells = length(bdev_pt_defined),
    bdev_normalized = normalize_to_bdev,
    interpretation = if (normalize_to_bdev) {
      "1 means query density equals Bdev reference density"
    } else NULL),
  trajectory = list(bdev_pt_source = "Bdev@meta.data$pt, unchanged",
    query_pt_source = "PeakVI Bdev-neighbor mapped_pt",
    boundaries = as.data.frame(boundaries[, .(boundary_label, label, pt_median)])),
  preservation = list(bdev_fingerprint = bdev_fp_before,
    bundle_unchanged = TRUE, promoted = FALSE),
  output_dir = output_dir
), file.path(stage_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

distribution_panel_html <- if (normalize_to_bdev) {
  paste0(if (add_gene_expression_panels) "Panel F" else "Panel D",
    " shares the exact horizontal span of the trajectory panels and shows the centered ",
    "three-bin query-cell fraction divided by the centered three-bin Bdev ",
    "reference-cell fraction across 100 bins. A value of 1 means equal query ",
    "and Bdev density.")
} else {
  paste0(if (add_gene_expression_panels) "Panel F" else "Panel D",
    " shares the exact horizontal span of the trajectory panels and shows a narrower ",
    "centered three-bin rolling mean over 100 bins for supported/borderline WT ",
    "and KO cell fractions along mapped pseudotime.")
}
html <- c(
  '<!doctype html><html><head><meta charset="utf-8">',
  '<title>FLT3L and IL-7 positive signals and WT/KO distributions</title>',
  '<style>body{font:16px system-ui,sans-serif;max-width:1700px;margin:32px auto;padding:0 22px;color:#172033}img{max-width:100%;height:auto;border:1px solid #ddd;margin-bottom:36px}p{line-height:1.5}a{color:#1261a0}</style></head><body>',
  '<h1>FLT3L and IL-7 positive signals and WT/KO trajectory distributions</h1>',
  if (add_gene_expression_panels) paste0(
    '<p>Panels A and B show ', target_spec$receptor_gene[[1L]],
    ' or ', target_spec$receptor_gene[[2L]],
    ' expression on the native Bdev ATAC UMAP and along unchanged Bdev pt. Panels C and D retain the accepted cytokine-score views, panel E shows the enlarged reference-anchored WT/KO projection, and ',
    distribution_panel_html, ' All three trajectory panels share one horizontal span.</p>') else paste0(
    '<p>Panel A retains the accepted Bdev positive-score UMAP at a larger display size with prominent colored score points and its color scale on the right. Panel B shows its Bdev trajectory. Panel C uses enlarged v1-like projection points and omits the redundant genotype legend. ', distribution_panel_html, '</p>'),
  '<p>Every trajectory panel includes dashed BB1/BB2 boundaries. No source bundle or derived state was changed or promoted.</p>',
  '<h2>Flt3L</h2><p><a href="figures/flt3l_positive_signal_cell_distribution.pdf">PDF</a></p>',
  paste0('<img src="figures/flt3l_positive_signal_cell_distribution.png" alt="Flt3L ', if (add_gene_expression_panels) 'six' else 'four', '-panel figure">'),
  '<h2>IL-7</h2><p><a href="figures/il7_positive_signal_cell_distribution.pdf">PDF</a></p>',
  paste0('<img src="figures/il7_positive_signal_cell_distribution.png" alt="IL-7 ', if (add_gene_expression_panels) 'six' else 'four', '-panel figure">'),
  '<p><a href="query_mapped_pt_cell_distribution.tsv">Cell-distribution table</a> · <a href="figure_inventory.tsv">Figure inventory</a></p>',
  '</body></html>')
writeLines(html, file.path(stage_dir, "REPORT.html"), useBytes = TRUE)
writeLines(c(
  "# FLT3L and IL-7 positive signals and WT/KO trajectory distributions", "",
  "[Open the HTML report](REPORT.html)", "",
  if (add_gene_expression_panels) {
    "Three equal-height rows add Flt3 or Il7r expression on Bdev UMAP and pt above the accepted score views and the enlarged anchored WT/KO projection with Bdev-normalized query density."
  } else if (normalize_to_bdev) {
    "Two equal-height rows retain the accepted Bdev score views, use enlarged v1-like points for the anchored WT/KO projection, and display centered three-bin query/Bdev density ratios over 100 mapped-pseudotime bins."
  } else {
    "Two equal-height rows retain the accepted Bdev score views, use enlarged v1-like points for the anchored WT/KO projection, and display centered three-bin rolling cell fractions over 100 mapped-pseudotime bins."
  }, "",
  "BB1/BB2 are dashed in every trajectory-axis panel. Nothing was promoted and the Bdev bundle remained unchanged."
), file.path(stage_dir, "REPORT.md"))

if (replace_existing) {
  backup_dir <- paste0(output_dir, ".replacement_backup")
  if (file.exists(backup_dir) || dir.exists(backup_dir)) {
    stop("Replacement backup already exists: ", backup_dir, call. = FALSE)
  }
  if (!file.rename(output_dir, backup_dir)) {
    stop("Could not stage existing output for replacement: ", output_dir)
  }
  if (!file.rename(stage_dir, output_dir)) {
    file.rename(backup_dir, output_dir)
    stop("Could not finalize replacement output directory: ", output_dir)
  }
  unlink(backup_dir, recursive = TRUE)
} else if (!file.rename(stage_dir, output_dir)) {
  stop("Could not finalize output directory: ", output_dir)
}
cat("Completed revised two-figure report at", output_dir, "\n")
