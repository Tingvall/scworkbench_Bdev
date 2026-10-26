#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 09_all_cytokine_refined_driver_report.R",
  "<bdev_bundle> <existing_38_profile_report> <output_dir>",
  "<expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
source_report <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
expected_fingerprint <- args[[4L]]
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(stage_dir, "profile_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(patchwork)
})

read_tsv <- function(path) read.delim(path, check.names = FALSE,
  na.strings = "NA", stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
slug <- function(x) gsub("(^_+|_+$)", "",
  tolower(gsub("[^A-Za-z0-9]+", "_", x)))
as_flag <- function(x) !is.na(x) & tolower(as.character(x)) %in%
  c("true", "t", "1", "yes")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}

required_inputs <- c("completion.json", "passing_signature_inventory.tsv",
  "passing_gene_membership.tsv", "ucell_scores.rds",
  "ucell_score_summary.tsv", "sliding_window_summary.tsv",
  "window_support.tsv", "trajectory_parameters.json",
  "bb_up_change_points.tsv")
input_paths <- file.path(source_report, required_inputs)
stopifnot(all(file.exists(input_paths)))
source_completion <- read_json(input_paths[[1L]], simplifyVector = TRUE)
stopifnot(identical(source_completion$status, "complete"),
  identical(source_completion$promoted, FALSE),
  identical(as.integer(source_completion$passing_signatures), 38L),
  identical(as.integer(source_completion$defined_pt_cells), 12642L),
  identical(as.integer(source_completion$missing_pt_cells), 14081L),
  identical(as.numeric(source_completion$width_fraction), 0.05),
  identical(as.numeric(source_completion$step_fraction), 0.01))

inventory <- read_tsv(file.path(source_report,
  "passing_signature_inventory.tsv"))
members <- read_tsv(file.path(source_report, "passing_gene_membership.tsv"))
score_family <- readRDS(file.path(source_report, "ucell_scores.rds"))
scores <- data.frame(key = as.character(score_family$key),
  cell_id = as.character(score_family$cell_id),
  value = as.numeric(score_family$value), stringsAsFactors = FALSE)
score_summary <- read_tsv(file.path(source_report, "ucell_score_summary.tsv"))
trajectory <- read_tsv(file.path(source_report,
  "sliding_window_summary.tsv"))
window_support <- read_tsv(file.path(source_report, "window_support.tsv"))
trajectory_parameters <- read_json(file.path(source_report,
  "trajectory_parameters.json"), simplifyVector = TRUE)
change_points <- read_tsv(file.path(source_report, "bb_up_change_points.tsv"))
change_points$boundary_label <- paste0("BB", seq_len(nrow(change_points)))

stopifnot(nrow(inventory) == 38L, all(as_flag(inventory$final_overall_pass)),
  all(as.integer(inventory$final_size) >= 20L),
  !anyDuplicated(inventory$signature_id),
  setequal(unique(members$signature_id), inventory$signature_id),
  !anyDuplicated(members[c("signature_id", "feature_id")]),
  setequal(unique(scores$key), inventory$signature_id),
  setequal(unique(trajectory$feature_set), inventory$signature_id),
  all(is.finite(scores$value)), all(scores$value >= 0 & scores$value <= 1),
  identical(trajectory_parameters$center, "median"),
  identical(trajectory_parameters$band,
    "25th to 75th percentiles (IQR)"))
inventory <- inventory[order(inventory$cytokine,
  match(inventory$direction, c("up", "down"))), , drop = FALSE]
members <- members[order(match(members$signature_id,
  inventory$signature_id), members$final_rank), , drop = FALSE]

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint),
  identical(source_completion$bdev_bundle_fingerprint, fingerprint_before))
cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "pt"))
stopifnot(identical(as.character(metadata$cell_id), cells),
  length(cells) == 26723L, sum(is.finite(metadata$pt)) == 12642L,
  sum(is.na(metadata$pt)) == 14081L)
sct_features <- as.character(wb_feature_ids(bundle, "SCT_data"))
stopifnot(all(unique(members$feature_id) %in% sct_features))

# Stage the exact reviewed memberships on SCT_data and let the native
# trajectory gene-set operation define and summarize the five drivers.
namespace <- "parse_cytokine_rest_only_refinement_all_top100_v1_drivers"
import_table <- data.frame(
  feature_set = members$signature_id,
  set_name = paste(members$cytokine, toupper(members$direction)),
  feature_id = members$feature_id,
  description = paste("Reviewed all-cytokine competitor-refined",
    members$cytokine, members$direction, "signature"),
  stringsAsFactors = FALSE)
import_csv <- file.path(stage_dir, "passing_gene_membership_import.csv")
write.csv(import_table, import_csv, row.names = FALSE)

session_root <- file.path(dirname(bundle_dir), "sessions",
  "all_cytokine_refined_top100_v1")
session_id <- "all_cytokine_refined_top100_drivers_v1"
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) {
  wb_session_resume(bundle, workspace)
} else {
  wb_session_create(bundle, workspace_root = session_root,
    session_id = session_id)
}
if (!is.null(session$durable$promotion_receipt)) {
  stop("Driver review session was unexpectedly promoted")
}

find_run <- function(session, label, op) {
  marker <- paste0("[", label, "] ")
  ids <- names(session$runs)[vapply(session$runs, function(x)
    identical(x$op_name, op) && startsWith(as.character(x$notes %||% ""),
      marker), logical(1))]
  stopifnot(length(ids) <= 1L)
  if (length(ids)) session$runs[[ids[[1L]]]] else NULL
}

import_run <- find_run(session, "driver_membership_import",
  "feature_set_import")
if (is.null(import_run)) {
  log_event("Staging the 38 reviewed memberships on SCT_data")
  imported <- wb_session_run_op(session, bundle, "feature_set_import",
    inputs = list(matrix = "SCT_data", cell_set = "all", min_cells = 10L,
      source_path = import_csv, source_format = "table",
      namespace = namespace, collection = "competitor_refined_drivers",
      species = "Mus musculus", organism = "Mus musculus",
      gene_id_type = "bundle_feature_id", set_id_col = "feature_set",
      set_name_col = "set_name", gene_col = "feature_id",
      description_col = "description", split_by_collection = FALSE,
      overwrite = FALSE, strict = TRUE,
      created_by = "scworkbench_Bdev driver review",
      source_id = "all_cytokine_competitor_refinement_top100_v1",
      source_version = "2026-09-25", origin = "derived",
      is_readonly = TRUE),
    notes = paste0("[driver_membership_import] Exact 38 passing memberships; ",
      "staged only for native driver analysis; no promotion."))
  session <- imported$session
  import_run <- imported$run
} else {
  log_event("Reusing completed driver membership import")
}
import_manifest <- import_run$outputs$import_manifest$table
stopifnot(nrow(import_manifest) == 38L,
  setequal(import_manifest$feature_set, inventory$signature_id),
  all(import_manifest$n_genes_raw == import_manifest$n_genes_mapped_to_bundle))

pt_use <- is.finite(metadata$pt)
pt <- wb_prim_cell_vector_family(rep("bdev_pt", sum(pt_use)),
  cells[pt_use], as.numeric(metadata$pt[pt_use]),
  meta = list(source_bundle = "bdev_reference", source_field = "pt",
    authoritative_values = "exact Bdev metadata values", direction = "as_is",
    fitted_model = FALSE))
driver_run <- find_run(session, "top5_mean_ranked_trajectory_drivers",
  "trajectory_geneset_dynamics")
if (is.null(driver_run)) {
  log_event("Calculating native mean-ranked top-five trajectory drivers for 38 profiles")
  driver_result <- wb_session_run_op(session, bundle,
    "trajectory_geneset_dynamics", inputs = list(matrix = "SCT_data",
      pseudotime = pt, pseudotime_key = "bdev_pt",
      membership = pt$cell_id, cell_set = "all", namespace = namespace,
      feature_sets = inventory$signature_id, score_method = "mean",
      pseudotime_direction = "as_is", n_bins = 50L,
      summary_fun = "mean", min_cells_per_bin = 20L,
      min_set_size = 20L, max_set_size = 5000L, top_n = 38L,
      top_driver_genes = 5L, smooth_method = "none", bootstrap_n = 0L,
      n_permutations = 0L, bootstrap_seed = 1729L,
      permutation_seed = 1729L),
    notes = paste0("[top5_mean_ranked_trajectory_drivers] Native driver ranking within ",
      "each reviewed signature on Bdev SCT_data. Drivers rank by dynamic ",
      "score then amplitude from unsmoothed per-bin means on exact pt."))
  session <- driver_result$session
  driver_run <- driver_result$run
} else {
  log_event("Reusing completed native driver analysis")
}

driver_table <- driver_run$outputs$driver_table$table
mean_driver_curves <- driver_run$outputs$driver_curve_table$table
driver_counts <- table(factor(driver_table$feature_set,
  levels = inventory$signature_id))
stopifnot(identical(as.integer(driver_counts), rep(5L, 38L)),
  !anyDuplicated(driver_table[c("feature_set", "feature_id")]),
  all(is.finite(driver_table$dynamic_score)),
  all(driver_table$dynamic_score > 0))

# Re-register just the mean-ranked drivers, then use a second native operation
# to calculate the requested median-expression lines without changing ranks.
median_namespace <- paste0(namespace, "_mean_ranked_top5")
driver_id_map <- setNames(paste0(inventory$signature_id,
  ".mean_ranked_top5"), inventory$signature_id)
ranked_import <- merge(driver_table,
  inventory[, c("signature_id", "cytokine", "direction")],
  by.x = "feature_set", by.y = "signature_id", sort = FALSE)
ranked_import <- ranked_import[order(match(ranked_import$feature_set,
  inventory$signature_id), ranked_import$driver_rank), ]
ranked_import$source_feature_set <- ranked_import$feature_set
ranked_import$feature_set <- unname(driver_id_map[ranked_import$feature_set])
ranked_import$set_name <- paste(ranked_import$cytokine,
  toupper(ranked_import$direction), "mean-ranked top five drivers")
ranked_import$description <- paste("Mean-ranked top five Bdev trajectory drivers for",
  ranked_import$cytokine, ranked_import$direction)
ranked_import_csv <- file.path(stage_dir,
  "mean_ranked_top5_driver_membership_import.csv")
write.csv(ranked_import[, c("feature_set", "set_name", "feature_id",
  "description")], ranked_import_csv, row.names = FALSE)

ranked_import_run <- find_run(session, "mean_ranked_top5_driver_import",
  "feature_set_import")
if (is.null(ranked_import_run)) {
  log_event("Staging mean-ranked top-five driver memberships for median curves")
  imported_drivers <- wb_session_run_op(session, bundle,
    "feature_set_import", inputs = list(matrix = "SCT_data",
      cell_set = "all", min_cells = 10L, source_path = ranked_import_csv,
      source_format = "table", namespace = median_namespace,
      collection = "mean_ranked_top5_drivers", species = "Mus musculus",
      organism = "Mus musculus", gene_id_type = "bundle_feature_id",
      set_id_col = "feature_set", set_name_col = "set_name",
      gene_col = "feature_id", description_col = "description",
      split_by_collection = FALSE, overwrite = FALSE, strict = TRUE,
      created_by = "scworkbench_Bdev driver review",
      source_id = "all_cytokine_mean_ranked_top5_drivers",
      source_version = "2026-09-28", origin = "derived",
      is_readonly = TRUE),
    notes = paste0("[mean_ranked_top5_driver_import] Exact top five from the ",
      "native mean-expression driver ranking; staged for median curves only."))
  session <- imported_drivers$session
  ranked_import_run <- imported_drivers$run
} else {
  log_event("Reusing mean-ranked top-five membership import")
}
stopifnot(nrow(ranked_import_run$outputs$import_manifest$table) == 38L)

median_run <- find_run(session, "mean_ranked_top5_median_curves",
  "trajectory_geneset_dynamics")
if (is.null(median_run)) {
  log_event("Calculating native median curves for the mean-ranked drivers")
  median_result <- wb_session_run_op(session, bundle,
    "trajectory_geneset_dynamics", inputs = list(matrix = "SCT_data",
      pseudotime = pt, pseudotime_key = "bdev_pt",
      membership = pt$cell_id, cell_set = "all",
      namespace = median_namespace,
      feature_sets = unname(driver_id_map[inventory$signature_id]),
      score_method = "mean", pseudotime_direction = "as_is",
      n_bins = 50L, summary_fun = "median", min_cells_per_bin = 20L,
      min_set_size = 5L, max_set_size = 5L, top_n = 38L,
      top_driver_genes = 5L, smooth_method = "none", bootstrap_n = 0L,
      n_permutations = 0L, bootstrap_seed = 1729L,
      permutation_seed = 1729L),
    notes = paste0("[mean_ranked_top5_median_curves] Native unsmoothed median ",
      "SCT_data curves for the independently mean-ranked driver genes."))
  session <- median_result$session
  median_run <- median_result$run
} else {
  log_event("Reusing native median driver curves")
}
median_driver_curves <- median_run$outputs$driver_curve_table$table
driver_fit <- median_run$outputs$fit_summary$table
reverse_driver_map <- setNames(names(driver_id_map), unname(driver_id_map))
median_driver_curves$feature_set <- unname(
  reverse_driver_map[median_driver_curves$feature_set])
median_driver_curves <- median_driver_curves[, setdiff(names(median_driver_curves),
  "driver_rank"), drop = FALSE]
median_driver_curves <- merge(median_driver_curves,
  driver_table[, c("feature_set", "feature_id", "driver_rank")],
  by = c("feature_set", "feature_id"), all.x = TRUE, sort = FALSE)
mean_driver_curves <- mean_driver_curves[, setdiff(names(mean_driver_curves),
  "driver_rank"), drop = FALSE]
mean_driver_curves <- merge(mean_driver_curves,
  driver_table[, c("feature_set", "feature_id", "driver_rank")],
  by = c("feature_set", "feature_id"), all = FALSE, sort = FALSE)
internal_minmax <- function(x) {
  limits <- range(x)
  if (diff(limits) == 0) return(rep(0, length(x)))
  (x - limits[[1L]]) / diff(limits)
}
curve_group <- interaction(mean_driver_curves$feature_set,
  mean_driver_curves$feature_id, drop = TRUE)
mean_driver_curves$plot_value_internal_0_1 <- ave(mean_driver_curves$value,
  curve_group, FUN = internal_minmax)
scaled_min <- tapply(mean_driver_curves$plot_value_internal_0_1, curve_group, min)
scaled_max <- tapply(mean_driver_curves$plot_value_internal_0_1, curve_group, max)
stopifnot(setequal(unique(median_driver_curves$feature_set), inventory$signature_id),
  setequal(unique(mean_driver_curves$feature_set), inventory$signature_id),
  !anyNA(median_driver_curves$driver_rank),
  !anyNA(mean_driver_curves$driver_rank),
  all(is.finite(median_driver_curves$value)),
  all(is.finite(mean_driver_curves$value)),
  all(is.finite(mean_driver_curves$plot_value_internal_0_1)),
  all(mean_driver_curves$plot_value_internal_0_1 >= 0 &
    mean_driver_curves$plot_value_internal_0_1 <= 1),
  all(abs(scaled_min) < 1e-12),
  all(abs(scaled_max) < 1e-12 | abs(scaled_max - 1) < 1e-12),
  identical(unique(driver_fit$summary_fun), "median"),
  identical(unique(driver_fit$smooth_method), "none"),
  identical(unique(as.integer(driver_fit$top_driver_genes)), 5L),
  identical(unique(as.integer(driver_fit$n_bins)), 50L))
write_tsv(import_manifest, file.path(stage_dir, "membership_import_manifest.tsv"))
write_tsv(driver_table, file.path(stage_dir, "top5_driver_genes.tsv"))
write_tsv(median_driver_curves, file.path(stage_dir,
  "top5_driver_median_expression_curves.tsv"))
write_tsv(mean_driver_curves, file.path(stage_dir,
  "top5_driver_mean_expression_curves.tsv"))
write_tsv(driver_fit, file.path(stage_dir, "driver_fit_summary.tsv"))
write_json(wb_session_dag(session, bundle = bundle, include_artifacts = TRUE,
  include_active = TRUE), file.path(stage_dir, "native_session_dag.json"))

space <- wb_space(bundle, "umap.atac", cells = cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells))
  coordinates <- data.frame(cell_id = cells, x = as.numeric(space$x),
    y = as.numeric(space$y), stringsAsFactors = FALSE)
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  coordinates <- data.frame(cell_id = cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]), stringsAsFactors = FALSE)
}
score_matrix <- lapply(inventory$signature_id, function(id) {
  local <- scores[scores$key == id, ]
  stopifnot(identical(as.character(local$cell_id), cells))
  local
})
names(score_matrix) <- inventory$signature_id

profile_palette <- c(up = "#D73027", down = "#2166AC")
gene_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
unsupported <- window_support[!window_support$eligible_window, , drop = FALSE]
unsupported$region_left <- unsupported$window_left
unsupported$region_right <- unsupported$window_right

make_umap <- function(id, title) {
  local <- score_matrix[[id]]
  values <- as.numeric(local$value)
  positive <- values > 0
  cap <- as.numeric(quantile(values[positive], 0.99, names = FALSE, type = 7))
  display <- rep(NA_real_, length(values))
  display[positive] <- pmin(values[positive], cap) / cap
  data <- coordinates
  data$display_score <- display
  data <- data[order(ifelse(is.na(data$display_score), -Inf,
    data$display_score)), ]
  ggplot(data, aes(x, y)) +
    geom_point(data = data[is.na(data$display_score), ], color = "#CDD2DA",
      size = 0.12, alpha = 0.7, stroke = 0) +
    geom_point(data = data[!is.na(data$display_score), ],
      aes(color = display_score), size = 0.28, alpha = 1, stroke = 0) +
    scale_color_gradientn(colours = c("#0066FF", "#00D9FF", "#00F060",
      "#FFF000", "#FF5A00", "#FF00A8"), limits = c(0, 1),
      oob = scales::squish, name = "UCell\n(panel scale)") +
    coord_equal() + theme_void(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
    labs(title = "UCell score on Bdev umap.atac",
      subtitle = title, caption = "Positive scores scaled to within-profile p99; zero is gray.")
}

make_trajectory <- function(id, direction) {
  data <- trajectory[trajectory$feature_set == id &
    trajectory$eligible_window, , drop = FALSE]
  ggplot(data, aes(window_midpoint, score_median)) +
    geom_rect(data = unsupported,
      aes(xmin = region_left, xmax = region_right, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "#64748B", alpha = 0.13) +
    geom_vline(data = change_points, aes(xintercept = pt_median),
      inherit.aes = FALSE, linewidth = 0.5, linetype = "dotted",
      color = "#111827") +
    geom_ribbon(aes(ymin = score_q25, ymax = score_q75),
      fill = profile_palette[[direction]], alpha = 0.22, colour = NA) +
    geom_line(color = profile_palette[[direction]], linewidth = 0.9) +
    labs(title = "UCell trajectory: 5% rolling windows",
      x = "Existing Bdev pseudotime (pt)", y = "UCell score",
      caption = "Median and IQR; 1% steps; ≥20 cells; no smoothing.") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"))
}

make_drivers <- function(id) {
  data <- mean_driver_curves[mean_driver_curves$feature_set == id, , drop = FALSE]
  ranks <- driver_table[driver_table$feature_set == id,
    c("feature_id", "driver_rank"), drop = FALSE]
  ranks <- ranks[order(ranks$driver_rank), ]
  genes <- ranks$feature_id
  colors <- setNames(gene_palette[seq_along(genes)], genes)
  data$feature_id <- factor(data$feature_id, levels = genes)
  ends <- do.call(rbind, lapply(genes, function(gene) {
    x <- data[as.character(data$feature_id) == gene, , drop = FALSE]
    x[which.max(x$pseudotime), , drop = FALSE]
  }))
  span <- diff(range(data$pseudotime))
  ggplot(data, aes(pseudotime, plot_value_internal_0_1, color = feature_id,
      group = feature_id)) +
    geom_vline(data = change_points, aes(xintercept = pt_median),
      inherit.aes = FALSE, linewidth = 0.7, linetype = "dashed",
      color = "#374151") +
    geom_text(data = change_points,
      aes(x = pt_median, y = Inf, label = boundary_label),
      inherit.aes = FALSE, vjust = 1.35, fontface = "bold", size = 3.2,
      color = "#374151") +
    geom_line(linewidth = 1) +
    geom_text_repel(data = ends, aes(label = feature_id), direction = "y",
      hjust = 0, nudge_x = 0.02 * span, seed = 1729,
      box.padding = 0.2, point.padding = 0.15, min.segment.length = 0,
      max.overlaps = Inf, show.legend = FALSE) +
    scale_color_manual(values = colors, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.20))) +
    coord_cartesian(clip = "off") +
    labs(title = "Top five member-gene drivers",
      subtitle = paste("Native mean-expression dynamic-score rank;",
        "each mean-expression curve scaled internally to 0-1"),
      x = "Existing Bdev pseudotime (pt)",
      y = "Within-gene scaled mean expression",
      caption = paste("Dashed lines mark BB1/BB2 medians.",
        "50 equal-width bins; ≥20 cells; unsmoothed means.",
        "Raw mean and median values remain in output tables.")) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
      legend.position = "none", plot.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 80, 5.5, 5.5))
}

log_event("Rendering 38 combined UMAP, rolling-IQR, and driver figures")
pdf_path <- file.path(figure_dir,
  "all_38_profiles_umap_iqr_top5_drivers.pdf")
grDevices::pdf(pdf_path, width = 15, height = 10.5, onefile = TRUE)
figure_rows <- vector("list", nrow(inventory))
for (i in seq_len(nrow(inventory))) {
  row <- inventory[i, ]
  id <- as.character(row$signature_id)
  profile_label <- paste(row$cytokine, toupper(row$direction), "—",
    row$final_size, "genes")
  combined <- (make_umap(id, profile_label) |
    make_trajectory(id, row$direction)) / make_drivers(id) +
    plot_layout(heights = c(1, 0.9)) +
    plot_annotation(title = paste("Competitor-refined profile:", profile_label),
      theme = theme(plot.title = element_text(size = 17, face = "bold")))
  print(combined)
  file <- paste0(sprintf("%02d_", i), slug(paste(row$cytokine,
    row$direction)), "_umap_iqr_top5_drivers.png")
  ggsave(file.path(figure_dir, file), combined, width = 15, height = 10.5,
    dpi = 220, units = "in", bg = "white")
  local_drivers <- driver_table[driver_table$feature_set == id, ]
  local_drivers <- local_drivers[order(local_drivers$driver_rank), ]
  figure_rows[[i]] <- data.frame(signature_id = id,
    cytokine = row$cytokine, direction = row$direction,
    signature_size = as.integer(row$final_size),
    driver_1 = local_drivers$feature_id[[1L]],
    driver_2 = local_drivers$feature_id[[2L]],
    driver_3 = local_drivers$feature_id[[3L]],
    driver_4 = local_drivers$feature_id[[4L]],
    driver_5 = local_drivers$feature_id[[5L]], file = file,
    stringsAsFactors = FALSE)
}
grDevices::dev.off()
figure_inventory <- do.call(rbind, figure_rows)
write_tsv(figure_inventory, file.path(stage_dir, "figure_inventory.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
figure_paths <- file.path(figure_dir, figure_inventory$file)
checks <- list(source_has_38_passing_profiles = nrow(inventory) == 38L,
  every_signature_has_at_least_20_genes = all(inventory$final_size >= 20L),
  exact_source_pt_preserved = identical(pt$value,
    as.numeric(metadata$pt[pt_use])),
  source_missing_pt_preserved = sum(is.na(metadata$pt)) == 14081L,
  five_unique_drivers_per_profile = all(driver_counts == 5L),
  driver_curves_are_native_medians =
    identical(unique(driver_fit$summary_fun), "median"),
  driver_plot_values_are_internal_0_1 =
    all(mean_driver_curves$plot_value_internal_0_1 >= 0 &
      mean_driver_curves$plot_value_internal_0_1 <= 1),
  no_driver_smoothing = identical(unique(driver_fit$smooth_method), "none"),
  every_profile_figure_exists = all(file.exists(figure_paths)),
  multipage_pdf_exists = file.exists(pdf_path),
  bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after), result_promotions = 0L)
stopifnot(all(unlist(checks[names(checks) != "result_promotions"])),
  checks$result_promotions == 0L)
write_json(list(status = "passed", checks = checks),
  file.path(stage_dir, "checks.json"))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
provenance_inputs <- c(input_paths, file.path(bundle_dir, "manifest.json"),
  script_path, "figure_scripts/bdev_reference/run_09_all_cytokine_refined_driver_report.sh",
  "config/bdev_reference.yml")
input_manifest <- data.frame(path = normalizePath(provenance_inputs,
  mustWork = TRUE), sha256 = vapply(provenance_inputs, wb_source_sha256,
  character(1)), stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

html_rows <- vapply(seq_len(nrow(figure_inventory)), function(i) {
  x <- figure_inventory[i, ]
  paste0("<h2>", x$cytokine, " ", toupper(x$direction), "</h2>",
    "<p>Drivers: ", paste(unlist(x[paste0("driver_", 1:5)]),
      collapse = ", "), "</p>",
    '<img src="profile_figures/', x$file, '" alt="', x$cytokine, " ",
    x$direction, ' combined profile figure">')
}, character(1))
html <- c('<!doctype html><html><head><meta charset="utf-8">',
  '<title>All-cytokine refined profile drivers</title>',
  '<style>body{font:16px system-ui,sans-serif;max-width:1400px;margin:32px auto;padding:0 22px;color:#172033}img{max-width:100%;height:auto;border:1px solid #ddd}p{line-height:1.5}code{overflow-wrap:anywhere}a{color:#1261a0}</style></head><body>',
  '<h1>Thirty-eight competitor-refined cytokine profiles in Bdev</h1>',
  '<p>Each figure combines the native RNA-count UCell score on the unchanged <code>umap.atac</code>, the existing 5% rolling-window median/IQR trajectory, and five member-gene drivers ranked by native scWorkbench trajectory dynamic score on <code>SCT_data</code>. Driver curves show unsmoothed mean expression in 50 equal-width bins, independently min-max scaled to 0-1 within each gene so sparse, low-expression trajectories remain visible. Raw mean and median curves remain in separate tables.</p>',
  '<p><a href="profile_figures/all_38_profiles_umap_iqr_top5_drivers.pdf">Download the 38-page PDF</a> · <a href="figure_inventory.tsv">Figure inventory</a> · <a href="top5_driver_genes.tsv">Driver table</a></p>',
  html_rows, '</body></html>')
writeLines(html, file.path(stage_dir, "REPORT.html"))
writeLines(c("# All-cytokine refined profile drivers", "",
  "[Open the report](REPORT.html)", "",
  "Thirty-eight passing profiles; five native scWorkbench-ranked member-gene drivers per profile.",
  "Each combined figure contains the UCell UMAP, 5% rolling-window median/IQR curve, and internally 0-1 scaled mean SCT_data driver-expression curves."),
  file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", promoted = FALSE,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  source_report = source_report,
  source_report_sha256 = wb_source_sha256(file.path(source_report,
    "completion.json")), bundle = bundle_dir,
  bundle_fingerprint = fingerprint_before, profiles = 38L,
  cytokines = length(unique(inventory$cytokine)), drivers_per_profile = 5L,
  driver_matrix = "SCT_data",
  driver_ranking = "scWorkbench mean-expression dynamic_score, amplitude tie-breaker",
  driver_summary = "mean", driver_median_values_retained = TRUE,
  driver_bins = 50L,
  driver_min_cells_per_bin = 20L, driver_smoothing = "none",
  driver_plot_scale = "per-gene min-max scaling of mean trajectory to [0,1]",
  driver_raw_mean_and_median_values_retained = TRUE,
  ucell_matrix = "RNA_counts", embedding = "umap.atac",
  rolling_width_fraction = 0.05, rolling_step_fraction = 0.01,
  rolling_min_cells = 20L, defined_pt_cells = 12642L,
  missing_pt_cells = 14081L, result_promotions = 0L,
  driver_session_workspace = workspace,
  multipage_pdf = file.path("profile_figures", basename(pdf_path)),
  figures = figure_inventory$file,
  report_sha256 = wb_source_sha256(file.path(stage_dir, "REPORT.html"))),
  file.path(stage_dir, "completion.json"))

if (!file.rename(stage_dir, output_dir)) {
  stop("Failed to atomically finalize output directory")
}
log_event("Completed 38-profile driver report at", output_dir)
