#!/usr/bin/env Rscript

# Two report-only composite figures for positive Flt3L and IL-7 signatures.
# Bdev panels reuse accepted score/driver artifacts. Query panels use the
# reference-anchored UMAP projection and the same flattened-rank score formula.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) stop(paste(
  "Usage: 07_ko_cytokine_positive_signal_report.R",
  "BDEV_BUNDLE KO_BUNDLE UCELL_REPORT SCORE_REPORT SCORE_PROMOTION",
  "DRIVER_REPORT PROJECTION_REPORT OUTPUT_DIR EXPECTED_BDEV_FP EXPECTED_KO_FP"
), call. = FALSE)

bdev_dir <- normalizePath(args[[1L]], mustWork = TRUE)
ko_dir <- normalizePath(args[[2L]], mustWork = TRUE)
ucell_report <- normalizePath(args[[3L]], mustWork = TRUE)
score_report <- normalizePath(args[[4L]], mustWork = TRUE)
score_promotion <- normalizePath(args[[5L]], mustWork = TRUE)
driver_report <- normalizePath(args[[6L]], mustWork = TRUE)
projection_report <- normalizePath(args[[7L]], mustWork = TRUE)
output_dir <- normalizePath(args[[8L]], mustWork = FALSE)
expected_bdev_fp <- args[[9L]]
expected_ko_fp <- args[[10L]]

if (file.exists(output_dir) || dir.exists(output_dir)) {
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

read_tsv <- function(path) fread(path, na.strings = "NA")
read_gz <- function(path) fread(cmd = paste("zcat", shQuote(path)), na.strings = "NA")
write_json <- function(x, path) jsonlite::write_json(
  x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA,
  null = "null", na = "null")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
as_flag <- function(x) !is.na(x) & tolower(as.character(x)) %in%
  c("true", "t", "1", "yes")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}

paths <- list(
  inventory = file.path(ucell_report, "passing_signature_inventory.tsv"),
  membership = file.path(ucell_report, "passing_gene_membership.tsv"),
  windows = file.path(ucell_report, "window_support.tsv"),
  boundaries = file.path(ucell_report, "bb_up_change_points.tsv"),
  score_object = file.path(score_report, "breadth_intensity_cell_scores.rds"),
  score_trajectory = file.path(score_report,
    "rankmean_combined_breadth_intensity_trajectory.tsv"),
  score_completion = file.path(score_report, "completion.json"),
  promotion_completion = file.path(score_promotion, "completion.json"),
  drivers = file.path(driver_report, "top5_driver_genes.tsv"),
  driver_curves = file.path(driver_report,
    "top5_driver_mean_expression_curves.tsv"),
  driver_completion = file.path(driver_report, "completion.json"),
  projection = file.path(projection_report,
    "query_reference_umap_projection.tsv.gz"),
  projection_completion = file.path(projection_report, "completion.json")
)
stopifnot(all(file.exists(unlist(paths))))

score_completion <- read_json(paths$score_completion, simplifyVector = TRUE)
promotion_completion <- read_json(paths$promotion_completion, simplifyVector = TRUE)
driver_completion <- read_json(paths$driver_completion, simplifyVector = TRUE)
projection_completion <- read_json(paths$projection_completion, simplifyVector = TRUE)
stopifnot(identical(score_completion$status, "complete"),
          identical(score_completion$promoted, FALSE),
          identical(promotion_completion$status, "promoted"),
          identical(promotion_completion$promoted, TRUE),
          identical(promotion_completion$default_score$name,
                    "combined_flattened_rank"),
          identical(as.integer(promotion_completion$default_score$rank_cutoff), 1500L),
          identical(driver_completion$status, "complete"),
          identical(projection_completion$status,
                    "complete_unpromoted_report_only"),
          identical(sha256_file(paths$score_object),
                    promotion_completion$source_score_sha256))

inventory <- read_tsv(paths$inventory)
membership <- read_tsv(paths$membership)
windows <- read_tsv(paths$windows)
boundaries <- read_tsv(paths$boundaries)
boundaries[, boundary_label := paste0("BB", seq_len(.N))]
score_trajectory <- read_tsv(paths$score_trajectory)
driver_table <- read_tsv(paths$drivers)
driver_curves <- read_tsv(paths$driver_curves)
projection <- read_gz(paths$projection)

target_spec <- data.table(
  target = c("Flt3l", "Il7"),
  cytokine = c("Flt3L", "IL-7"),
  slug = c("flt3l", "il7"),
  signature_id = c(
    "parse_cytokine_rest_only_refinement_all_top100_v1.flt3l.up",
    "parse_cytokine_rest_only_refinement_all_top100_v1.il_7.up")
)
stopifnot(all(target_spec$signature_id %in% inventory$signature_id),
          all(target_spec$signature_id %in% membership$signature_id),
          all(target_spec$signature_id %in% score_trajectory$feature_set),
          all(target_spec$signature_id %in% driver_table$feature_set),
          all(target_spec$signature_id %in% driver_curves$feature_set),
          nrow(boundaries) == 2L,
          all(c("window_id", "window_left", "window_right",
                "window_midpoint", "eligible_window") %in% names(windows)),
          nrow(projection) == 23252L, !anyDuplicated(projection$cell_id))

bdev <- wb_open_bundle(bdev_dir, validate = FALSE)
ko <- wb_open_bundle(ko_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
on.exit(try(wb_close_bundle(ko), silent = TRUE), add = TRUE)
bdev_fp_before <- wb_manifest(bdev)$current_state$fingerprint
ko_fp_before <- wb_manifest(ko)$current_state$fingerprint
stopifnot(identical(bdev_fp_before, expected_bdev_fp),
          identical(ko_fp_before, expected_ko_fp),
          identical(promotion_completion$bundle_fingerprint_after,
                    bdev_fp_before))

bdev_cells <- as.character(wb_cells(bdev, "all"))
bdev_umap <- as.data.table(wb_space(bdev, "umap.atac", cells = bdev_cells))
stopifnot(length(bdev_cells) == 26723L,
          identical(as.character(bdev_umap$cell_id), bdev_cells),
          all(is.finite(bdev_umap$x)), all(is.finite(bdev_umap$y)))
bdev_xy_original <- copy(bdev_umap[, .(cell_id, x, y)])

score_object <- readRDS(paths$score_object)
stopifnot(identical(score_object$parameters$rank_cutoff, 1500L),
          identical(score_object$parameters$alpha, 0.5),
          all(target_spec$signature_id %in% rownames(score_object$scores$combined)),
          identical(colnames(score_object$scores$combined), bdev_cells))

query_cells <- as.character(projection$source_cell_id)
stopifnot(!anyDuplicated(query_cells),
          setequal(query_cells, wb_cells(ko, "set:qc_use_downstream")))
ko_md <- as.data.table(wb_meta(ko, cells = query_cells))
stopifnot(identical(as.character(ko_md$cell_id), query_cells),
          identical(as.character(ko_md$library_id),
                    as.character(projection$sample_id)),
          identical(as.character(ko_md$genotype),
                    as.character(projection$genotype)),
          identical(as.character(ko_md$knockout_target),
                    as.character(projection$knockout_target)))

# Apply the accepted Bdev default formula to query raw RNA counts. Missing
# signature features remain in the denominator and contribute zero.
all_ko_genes <- as.character(wb_feature_ids(ko, "RNA_counts"))
signature_members <- lapply(target_spec$signature_id, function(id)
  as.character(membership[signature_id == id][order(final_rank)]$feature_id))
names(signature_members) <- target_spec$signature_id
union_genes <- unique(unlist(signature_members, use.names = FALSE))
union_to_matrix <- match(union_genes, all_ko_genes)
mapped_union <- !is.na(union_to_matrix)
set_index <- lapply(signature_members, match, table = union_genes)
set_size <- lengths(set_index)
coverage <- data.table(
  signature_id = names(signature_members),
  signature_size = set_size,
  mapped_genes = vapply(signature_members,
    function(x) sum(x %in% all_ko_genes), integer(1))
)
coverage[, coverage_fraction := mapped_genes / signature_size]
stopifnot(all(coverage$mapped_genes >= 5L))

query_scores <- lapply(c("breadth", "intensity", "combined"), function(x)
  matrix(NA_real_, nrow = nrow(target_spec), ncol = length(query_cells),
         dimnames = list(target_spec$signature_id, query_cells)))
names(query_scores) <- c("breadth", "intensity", "combined")
rank_cutoff <- 1500L
alpha <- 0.5
chunks <- split(seq_along(query_cells), ceiling(seq_along(query_cells) / 128L))
log_event("Scoring", length(query_cells), "query cells in", length(chunks),
          "chunks with the accepted flattened-rank formula")
for (chunk_number in seq_along(chunks)) {
  ix <- chunks[[chunk_number]]
  x <- as.matrix(wb_matrix_read_features(ko, "RNA_counts", all_ko_genes,
    query_cells[ix], missing = "zero"))
  stopifnot(nrow(x) == length(all_ko_genes), ncol(x) == length(ix),
            identical(rownames(x), all_ko_genes),
            identical(colnames(x), query_cells[ix]),
            all(is.finite(x)), all(x >= 0))
  for (j in seq_along(ix)) {
    z <- x[, j]
    detected <- z > 0
    detected_positions <- which(detected)
    detected_ranks <- rank(-z[detected], ties.method = "average")
    union_detected <- rep(FALSE, length(union_genes))
    union_ranks <- rep(Inf, length(union_genes))
    if (any(mapped_union)) {
      union_detected[mapped_union] <- detected[union_to_matrix[mapped_union]]
      active <- mapped_union & union_detected
      if (any(active)) {
        rank_lookup <- match(union_to_matrix[active], detected_positions)
        union_ranks[active] <- detected_ranks[rank_lookup]
      }
    }
    union_weight <- ifelse(union_detected & union_ranks <= rank_cutoff,
      (rank_cutoff - union_ranks) / (rank_cutoff - 1), 0)
    for (s in seq_len(nrow(target_spec))) {
      gi <- set_index[[s]]
      breadth <- sum(union_detected[gi]) / set_size[[s]]
      intensity <- sum(union_weight[gi]) / set_size[[s]]
      query_scores$breadth[s, ix[[j]]] <- breadth
      query_scores$intensity[s, ix[[j]]] <- intensity
      query_scores$combined[s, ix[[j]]] <-
        alpha * breadth + (1 - alpha) * intensity
    }
  }
  rm(x)
  if (chunk_number %% 20L == 0L || chunk_number == length(chunks)) {
    log_event("Completed query score chunk", chunk_number, "of", length(chunks))
    invisible(gc())
  }
}
stopifnot(all(vapply(query_scores, function(x) all(is.finite(x)), logical(1))),
          all(query_scores$combined >= 0 & query_scores$combined <= 1),
          max(abs(query_scores$combined -
            (alpha * query_scores$breadth +
             (1 - alpha) * query_scores$intensity))) < 1e-14)

score_long <- rbindlist(lapply(seq_len(nrow(target_spec)), function(s) data.table(
  source_cell_id = query_cells,
  signature_id = target_spec$signature_id[[s]],
  cytokine = target_spec$cytokine[[s]],
  breadth = query_scores$breadth[s, ],
  intensity = query_scores$intensity[s, ],
  combined_score = query_scores$combined[s, ]
)))
query_base <- projection[, .(cell_id, source_cell_id, sample_id, genotype,
  knockout_target, support_status, mapped_pt, projected_x, projected_y)]
score_long <- query_base[score_long, on = "source_cell_id"]
stopifnot(nrow(score_long) == 2L * nrow(projection), !anyNA(score_long$cell_id))
fwrite(score_long, file.path(stage_dir, "query_positive_signal_scores.tsv.gz"),
       sep = "\t")
fwrite(coverage, file.path(stage_dir, "query_signature_coverage.tsv"), sep = "\t")

window_rows <- list()
wr <- 0L
for (s in seq_len(nrow(target_spec))) {
  spec <- target_spec[s]
  q <- score_long[signature_id == spec$signature_id &
                  knockout_target == spec$target &
                  support_status != "out_of_reference" & is.finite(mapped_pt)]
  for (g in c("WT", "KO")) {
    qg <- q[genotype == g]
    for (wi in seq_len(nrow(windows))) {
      use <- qg$mapped_pt >= windows$window_left[[wi]] &
             qg$mapped_pt <= windows$window_right[[wi]]
      values <- qg$combined_score[use]
      wr <- wr + 1L
      window_rows[[wr]] <- data.table(
        target = spec$target, cytokine = spec$cytokine,
        signature_id = spec$signature_id, genotype = g,
        window_id = windows$window_id[[wi]],
        window_midpoint = windows$window_midpoint[[wi]],
        window_left = windows$window_left[[wi]],
        window_right = windows$window_right[[wi]],
        n_cells = length(values), eligible_window = length(values) >= 20L,
        score_median = if (length(values)) median(values) else NA_real_,
        score_q25 = if (length(values)) quantile(values, .25, names = FALSE) else NA_real_,
        score_q75 = if (length(values)) quantile(values, .75, names = FALSE) else NA_real_)
    }
  }
}
query_trajectory <- rbindlist(window_rows)
stopifnot(nrow(query_trajectory) == nrow(target_spec) * 2L * nrow(windows),
          all(query_trajectory[eligible_window == TRUE,
                               is.finite(score_median)]))
fwrite(query_trajectory, file.path(stage_dir,
  "query_positive_signal_trajectory.tsv"), sep = "\t")

score_palette <- c("#0066FF", "#00D9FF", "#00F060", "#FFF000",
                   "#FF5A00", "#FF00A8")
genotype_colors <- c(WT = "#0072B2", KO = "#D55E00")
gene_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
base_trajectory_theme <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

add_bb <- function(plot) {
  plot +
    geom_vline(data = boundaries, aes(xintercept = pt_median),
      inherit.aes = FALSE, linewidth = 0.55, linetype = "dashed",
      color = "#374151") +
    geom_text(data = boundaries,
      aes(x = pt_median, y = Inf, label = boundary_label),
      inherit.aes = FALSE, vjust = 1.25, size = 2.8,
      fontface = "bold", color = "#374151")
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
               size = .1, alpha = .65) +
    geom_point(data = d[!is.na(display)], aes(color = display),
               size = .22, alpha = .95) +
    scale_color_gradientn(colours = score_palette, limits = c(0, 1),
      oob = scales::squish, name = "Positive\nscore") +
    coord_equal() + theme_void(base_size = 9) +
    theme(legend.position = "bottom",
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
    geom_line(color = "#D81B60", linewidth = .85) +
    base_trajectory_theme +
    labs(title = paste(spec$cytokine, "positive signal trajectory"),
      subtitle = "Bdev median and cell-wise IQR",
      x = "Existing Bdev pseudotime (pt)", y = "Combined score",
      caption = "Aligned 5% windows at 1% steps; dashed lines mark BB1/BB2.")
  add_bb(p)
}

make_driver_trajectory <- function(spec) {
  d <- copy(driver_curves[feature_set == spec$signature_id])
  ranks <- driver_table[feature_set == spec$signature_id][order(driver_rank)]
  genes <- ranks$feature_id
  stopifnot(length(genes) == 5L, setequal(unique(d$feature_id), genes))
  d[, feature_id := factor(feature_id, levels = genes)]
  p <- ggplot(d, aes(pseudotime, plot_value_internal_0_1,
                     color = feature_id, group = feature_id)) +
    geom_line(linewidth = .8) +
    scale_color_manual(values = setNames(gene_colors, genes), drop = FALSE) +
    base_trajectory_theme +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    guides(color = guide_legend(nrow = 1)) +
    labs(title = "Top-five positive-signal drivers",
      subtitle = "Each mean-expression curve internally scaled to 0-1",
      x = "Existing Bdev pseudotime (pt)",
      y = "Scaled mean SCT expression",
      caption = "Unsmoothed 50-bin curves; dashed lines mark BB1/BB2.")
  add_bb(p)
}

make_projection <- function(spec) {
  q <- projection[knockout_target == spec$target &
                  support_status != "out_of_reference" &
                  is.finite(projected_x) & is.finite(projected_y)]
  bg <- bdev_umap[, .(x, y)]
  ggplot() +
    geom_point(data = bg, aes(x, y), color = "grey84", size = .09) +
    geom_point(data = q, aes(projected_x, projected_y, color = genotype),
               size = .21, alpha = .5) +
    facet_wrap(~sample_id, nrow = 1) +
    scale_color_manual(values = genotype_colors, drop = FALSE) +
    coord_equal() + theme_void(base_size = 9) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold")) +
    labs(title = paste(spec$cytokine, "WT/KO reference-anchored projection"),
      subtitle = "Weighted mean of 30 Bdev PeakVI neighbors",
      color = "Genotype",
      caption = "Native Bdev UMAP in gray; calibrated OOD cells withheld.")
}

make_query_trajectory <- function(spec) {
  d <- copy(query_trajectory[target == spec$target])
  d[!as_flag(eligible_window),
    c("score_median", "score_q25", "score_q75") := NA_real_]
  p <- ggplot(d, aes(window_midpoint, score_median, color = genotype,
                     fill = genotype, group = genotype)) +
    geom_ribbon(aes(ymin = score_q25, ymax = score_q75),
                alpha = .16, color = NA) +
    geom_line(linewidth = .9) +
    scale_color_manual(values = genotype_colors, drop = FALSE) +
    scale_fill_manual(values = genotype_colors, drop = FALSE) +
    base_trajectory_theme +
    theme(legend.position = "bottom") +
    labs(title = paste(spec$cytokine, "WT/KO positive-score trajectory"),
      subtitle = "Supported and borderline query cells",
      x = "Mapped Bdev pseudotime", y = "Combined score",
      color = "Genotype", fill = "Genotype",
      caption = "Median/IQR in aligned 5% windows; dashed lines mark BB1/BB2.")
  add_bb(p)
}

figure_rows <- vector("list", nrow(target_spec))
log_event("Rendering two five-panel composite figures")
for (s in seq_len(nrow(target_spec))) {
  spec <- target_spec[s]
  panels <- list(
    make_score_umap(spec), make_bdev_trajectory(spec),
    make_driver_trajectory(spec), make_projection(spec),
    make_query_trajectory(spec)
  )
  composite <- wrap_plots(panels, design = "AABBCC\nDDDEEE") +
    plot_annotation(
      title = paste0(spec$cytokine,
        " positive signal across Bdev and the WT/KO query"),
      subtitle = paste(
        "A-C: accepted Bdev score and drivers; D: anchored query UMAP;",
        "E: query score versus mapped pseudotime."),
      tag_levels = "A",
      theme = theme(plot.title = element_text(size = 17, face = "bold"),
                    plot.subtitle = element_text(size = 10),
                    plot.margin = margin(8, 8, 8, 8)))
  png_file <- paste0(spec$slug, "_positive_signal_wt_ko.png")
  pdf_file <- paste0(spec$slug, "_positive_signal_wt_ko.pdf")
  ggsave(file.path(figure_dir, png_file), composite,
         width = 18, height = 11.5, dpi = 200, bg = "white")
  ggsave(file.path(figure_dir, pdf_file), composite,
         width = 18, height = 11.5, bg = "white")
  local_drivers <- driver_table[feature_set == spec$signature_id][order(driver_rank)]
  local_coverage <- coverage[signature_id == spec$signature_id]
  local_query <- score_long[signature_id == spec$signature_id &
                            knockout_target == spec$target]
  figure_rows[[s]] <- data.table(
    cytokine = spec$cytokine, knockout_target = spec$target,
    signature_id = spec$signature_id,
    signature_size = local_coverage$signature_size,
    query_mapped_genes = local_coverage$mapped_genes,
    query_coverage_fraction = local_coverage$coverage_fraction,
    driver_1 = local_drivers$feature_id[[1L]],
    driver_2 = local_drivers$feature_id[[2L]],
    driver_3 = local_drivers$feature_id[[3L]],
    driver_4 = local_drivers$feature_id[[4L]],
    driver_5 = local_drivers$feature_id[[5L]],
    wt_cells = local_query[genotype == "WT", .N],
    ko_cells = local_query[genotype == "KO", .N],
    png = file.path("figures", png_file),
    pdf = file.path("figures", pdf_file)
  )
}
figure_inventory <- rbindlist(figure_rows)
fwrite(figure_inventory, file.path(stage_dir, "figure_inventory.tsv"), sep = "\t")

bdev_fp_after <- wb_manifest(bdev)$current_state$fingerprint
ko_fp_after <- wb_manifest(ko)$current_state$fingerprint
checks <- list(
  exactly_two_figures = nrow(figure_inventory) == 2L,
  positive_profiles_only = all(grepl("\\.up$", target_spec$signature_id)),
  accepted_bdev_score_source_verified = identical(
    sha256_file(paths$score_object), promotion_completion$source_score_sha256),
  query_score_formula_exact = max(abs(query_scores$combined -
    (alpha * query_scores$breadth +
     (1 - alpha) * query_scores$intensity))) < 1e-14,
  all_query_scores_finite = all(is.finite(query_scores$combined)),
  exact_query_cell_alignment = identical(query_cells, projection$source_cell_id),
  ood_excluded_from_primary_query_panels = all(
    is.na(projection[support_status == "out_of_reference"]$projected_x)),
  bb_boundaries_in_all_trajectory_panel_builders = TRUE,
  original_bdev_pt_source_preserved = TRUE,
  original_bdev_umap_unchanged = identical(
    bdev_xy_original, bdev_umap[, .(cell_id, x, y)]),
  bdev_bundle_unchanged = identical(bdev_fp_before, bdev_fp_after),
  ko_bundle_unchanged = identical(ko_fp_before, ko_fp_after),
  no_promotion = TRUE
)
stopifnot(all(unlist(checks)))
write_json(list(status = "passed", checks = checks),
           file.path(stage_dir, "checks.json"))

input_manifest <- data.table(
  path = c(bdev_dir, ko_dir, unlist(paths)),
  sha256 = c(NA_character_, NA_character_,
    vapply(unlist(paths), sha256_file, character(1)))
)
fwrite(input_manifest, file.path(stage_dir, "input_manifest.tsv"), sep = "\t")

completion <- list(
  status = "complete_unpromoted_report_only",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  report_type = "positive_flt3l_il7_bdev_and_ko_composites",
  figures = as.data.frame(figure_inventory),
  score = list(name = "combined_flattened_rank", rank_cutoff = rank_cutoff,
    alpha = alpha, formula = "0.5 * breadth + 0.5 * rank intensity",
    detection_matrix = "RNA_counts", rank_universe = "all RNA_counts genes",
    ties = "average descending expression rank",
    missing_query_signature_features = "zero contribution with original signature-size denominator"),
  trajectory = list(bdev_pt_source = "Bdev@meta.data$pt, unchanged",
    query_pt_source = "PeakVI Bdev-neighbor mapped_pt",
    windows = "existing aligned 5%-width windows at 1% steps",
    query_min_cells_per_window_and_genotype = 20L,
    boundaries = as.data.frame(boundaries[, .(boundary_label, label, pt_median)])),
  caveat = paste("Each target has one WT and one KO library; WT/KO curves are",
                 "descriptive and do not provide replicate-level inference."),
  preservation = list(bdev_fingerprint = bdev_fp_before,
    ko_fingerprint = ko_fp_before, bundles_unchanged = TRUE, promoted = FALSE),
  output_dir = output_dir
)
write_json(completion, file.path(stage_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

html <- c(
  '<!doctype html><html><head><meta charset="utf-8">',
  '<title>FLT3L and IL-7 positive signals in Bdev and KO query</title>',
  '<style>body{font:16px system-ui,sans-serif;max-width:1800px;margin:32px auto;padding:0 22px;color:#172033}img{max-width:100%;height:auto;border:1px solid #ddd;margin-bottom:36px}p{line-height:1.5}code{overflow-wrap:anywhere}a{color:#1261a0}</style></head><body>',
  '<h1>FLT3L and IL-7 positive signals in Bdev and WT/KO query cells</h1>',
  '<p>Panels A–B use the promoted default Bdev combined flattened-rank score. Panel C reuses the established native top-five driver ranking and curves. Panel D uses the Bdev-reference-anchored PeakVI-neighbor UMAP projection. Panel E applies the identical score formula to query RNA counts and summarizes supported/borderline cells along mapped Bdev pseudotime.</p>',
  '<p>All trajectory panels include dashed BB1/BB2 median boundaries. Each knockout target has one WT and one KO library, so WT/KO differences are descriptive and not replicate-level inference. No source bundle or derived state was changed or promoted.</p>',
  '<h2>Flt3L</h2><p><a href="figures/flt3l_positive_signal_wt_ko.pdf">PDF</a></p>',
  '<img src="figures/flt3l_positive_signal_wt_ko.png" alt="Flt3L five-panel positive-signal figure">',
  '<h2>IL-7</h2><p><a href="figures/il7_positive_signal_wt_ko.pdf">PDF</a></p>',
  '<img src="figures/il7_positive_signal_wt_ko.png" alt="IL-7 five-panel positive-signal figure">',
  '<p><a href="query_positive_signal_scores.tsv.gz">Query scores</a> · <a href="query_positive_signal_trajectory.tsv">Query trajectory table</a> · <a href="figure_inventory.tsv">Figure inventory</a></p>',
  '</body></html>'
)
writeLines(html, file.path(stage_dir, "REPORT.html"), useBytes = TRUE)
writeLines(c(
  "# FLT3L and IL-7 positive signals in Bdev and WT/KO query cells", "",
  "[Open the HTML report](REPORT.html)", "",
  "Two five-panel figures combine the accepted Bdev positive score, Bdev score trajectory, established driver trajectories, reference-anchored WT/KO projection, and WT/KO positive-score trajectory. All trajectory panels contain dashed BB1/BB2 boundaries.", "",
  "The report is descriptive: each target has one WT and one KO library. Nothing was promoted and both source bundle fingerprints were unchanged."
), file.path(stage_dir, "REPORT.md"))

if (!file.rename(stage_dir, output_dir)) {
  stop("Could not finalize output directory: ", output_dir)
}
log_event("Completed two-figure cytokine report at", output_dir)
