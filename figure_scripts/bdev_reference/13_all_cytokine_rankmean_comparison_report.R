#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 13_all_cytokine_rankmean_comparison_report.R",
  "<bdev_bundle> <ucell_report> <gene_change_report>",
  "<output_dir> <expected_bundle_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
ucell_report <- normalizePath(args[[2L]], mustWork = TRUE)
change_report <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_fingerprint <- args[[5L]]
if (dir.exists(output_dir)) stop("Output already exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(stage_dir, "profile_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
})

read_tsv <- function(path) read.delim(path, check.names = FALSE,
  na.strings = "NA", stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt"); on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
slug <- function(x) gsub("(^_+|_+$)", "", tolower(gsub("[^A-Za-z0-9]+", "_", x)))
as_flag <- function(x) !is.na(x) & tolower(as.character(x)) %in%
  c("true", "t", "1", "yes")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n"); flush.console()
}

ucell_inputs <- c("completion.json", "passing_signature_inventory.tsv",
  "passing_gene_membership.tsv", "ucell_scores.rds",
  "sliding_window_summary.tsv", "window_support.tsv",
  "window_membership.rds", "bb_up_change_points.tsv")
change_inputs <- c("completion.json", "neighborhood_rank_change_cell_scores.rds",
  "signature_neighborhood_rank_change_trajectory.tsv", "checks.json")
ucell_paths <- file.path(ucell_report, ucell_inputs)
change_paths <- file.path(change_report, change_inputs)
stopifnot(all(file.exists(ucell_paths)), all(file.exists(change_paths)))
ucell_completion <- read_json(ucell_paths[[1L]], simplifyVector = TRUE)
change_completion <- read_json(change_paths[[1L]], simplifyVector = TRUE)
stopifnot(identical(ucell_completion$status, "complete"),
  identical(ucell_completion$promoted, FALSE),
  identical(as.integer(ucell_completion$passing_signatures), 38L),
  identical(change_completion$status, "complete"),
  identical(change_completion$promoted, FALSE),
  identical(as.integer(change_completion$profiles), 38L),
  identical(change_completion$trajectory_change$aggregation,
    "mean across all signature genes"))

inventory <- read_tsv(file.path(ucell_report,
  "passing_signature_inventory.tsv"))
members <- read_tsv(file.path(ucell_report, "passing_gene_membership.tsv"))
ucell_family <- readRDS(file.path(ucell_report, "ucell_scores.rds"))
ucell_scores <- data.frame(key = as.character(ucell_family$key),
  cell_id = as.character(ucell_family$cell_id),
  value = as.numeric(ucell_family$value), stringsAsFactors = FALSE)
ucell_trajectory <- read_tsv(file.path(ucell_report,
  "sliding_window_summary.tsv"))
window_support <- read_tsv(file.path(ucell_report, "window_support.tsv"))
window_membership <- readRDS(file.path(ucell_report, "window_membership.rds"))
change_points <- read_tsv(file.path(ucell_report, "bb_up_change_points.tsv"))
change_points$boundary_label <- paste0("BB", seq_len(nrow(change_points)))
change_scores <- readRDS(file.path(change_report,
  "neighborhood_rank_change_cell_scores.rds"))
change_trajectory <- read_tsv(file.path(change_report,
  "signature_neighborhood_rank_change_trajectory.tsv"))

stopifnot(nrow(inventory) == 38L, all(as_flag(inventory$final_overall_pass)),
  all(inventory$final_size >= 20L), !anyDuplicated(inventory$signature_id),
  setequal(unique(members$signature_id), inventory$signature_id),
  setequal(unique(ucell_scores$key), inventory$signature_id),
  setequal(unique(change_scores$key), inventory$signature_id),
  setequal(unique(change_trajectory$feature_set), inventory$signature_id),
  nrow(window_support) == length(window_membership),
  identical(as.character(window_support$window_id), names(window_membership)))
inventory <- inventory[order(inventory$cytokine,
  match(inventory$direction, c("up", "down"))), , drop = FALSE]
members <- members[order(match(members$signature_id,
  inventory$signature_id), members$final_rank), , drop = FALSE]

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint),
  identical(ucell_completion$bdev_bundle_fingerprint, fingerprint_before),
  identical(change_completion$bundle_fingerprint, fingerprint_before))
cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "pt"))
stopifnot(identical(as.character(metadata$cell_id), cells),
  length(cells) == 26723L, sum(is.finite(metadata$pt)) == 12642L,
  sum(is.na(metadata$pt)) == 14081L)
for (id in inventory$signature_id) {
  stopifnot(identical(ucell_scores$cell_id[ucell_scores$key == id], cells),
    identical(change_scores$cell_id[change_scores$key == id], cells))
}

# Stage the exact memberships in a separate durable session and run the native
# rank_mean trajectory operation. Nothing is promoted to the accepted bundle.
namespace <- "parse_cytokine_rest_only_refinement_rankmean_comparison_v1"
import_table <- data.frame(feature_set = members$signature_id,
  set_name = paste(members$cytokine, toupper(members$direction)),
  feature_id = members$feature_id,
  description = paste("Reviewed competitor-aware cytokine signature",
    members$cytokine, members$direction), stringsAsFactors = FALSE)
import_csv <- file.path(stage_dir, "passing_gene_membership_import.csv")
write.csv(import_table, import_csv, row.names = FALSE)

session_root <- file.path(dirname(bundle_dir), "sessions",
  "all_cytokine_native_rankmean_comparison_v1")
session_id <- "all_cytokine_native_rankmean_comparison_v1"
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) wb_session_resume(bundle, workspace) else
  wb_session_create(bundle, workspace_root = session_root,
    session_id = session_id)
if (!is.null(session$durable$promotion_receipt)) {
  stop("Rank-mean comparison session was unexpectedly promoted")
}
find_run <- function(session, label, op) {
  marker <- paste0("[", label, "] ")
  ids <- names(session$runs)[vapply(session$runs, function(x) {
    note <- if (is.null(x$notes)) "" else as.character(x$notes)
    identical(x$op_name, op) && startsWith(note, marker)
  }, logical(1))]
  stopifnot(length(ids) <= 1L)
  if (length(ids)) session$runs[[ids[[1L]]]] else NULL
}
import_run <- find_run(session, "rankmean_membership_import",
  "feature_set_import")
if (is.null(import_run)) {
  log_event("Staging 38 signatures for native rank_mean")
  result <- wb_session_run_op(session, bundle, "feature_set_import",
    inputs = list(matrix = "SCT_data", cell_set = "all", min_cells = 10L,
      source_path = import_csv, source_format = "table",
      namespace = namespace, collection = "rankmean_comparison",
      species = "Mus musculus", organism = "Mus musculus",
      gene_id_type = "bundle_feature_id", set_id_col = "feature_set",
      set_name_col = "set_name", gene_col = "feature_id",
      description_col = "description", split_by_collection = FALSE,
      overwrite = FALSE, strict = TRUE,
      created_by = "scworkbench_Bdev rankmean comparison",
      source_id = "all_cytokine_competitor_refinement_top100_v1",
      source_version = "2026-09-25", origin = "derived",
      is_readonly = TRUE),
    notes = paste0("[rankmean_membership_import] Exact 38 reviewed sets; ",
      "staged for native rank_mean comparison; no promotion."))
  session <- result$session; import_run <- result$run
} else log_event("Reusing staged rank_mean memberships")
stopifnot(nrow(import_run$outputs$import_manifest$table) == 38L,
  setequal(import_run$outputs$import_manifest$table$feature_set,
    inventory$signature_id))

pt_use <- is.finite(metadata$pt)
pt_family <- wb_prim_cell_vector_family(rep("bdev_pt", sum(pt_use)),
  cells[pt_use], as.numeric(metadata$pt[pt_use]),
  meta = list(source_bundle = "bdev_reference", source_field = "pt",
    authoritative_values = "exact Bdev metadata values", direction = "as_is",
    fitted_model = FALSE))
rankmean_run <- find_run(session, "native_rankmean_50bin",
  "trajectory_geneset_dynamics")
if (is.null(rankmean_run)) {
  log_event("Running native trajectory_geneset_dynamics rank_mean")
  result <- wb_session_run_op(session, bundle,
    "trajectory_geneset_dynamics", inputs = list(matrix = "SCT_data",
      pseudotime = pt_family, pseudotime_key = "bdev_pt",
      membership = pt_family$cell_id, cell_set = "all",
      namespace = namespace, feature_sets = inventory$signature_id,
      score_method = "rank_mean", pseudotime_direction = "as_is",
      n_bins = 50L, summary_fun = "median", min_cells_per_bin = 20L,
      min_set_size = 20L, max_set_size = 5000L, top_n = 38L,
      top_driver_genes = 5L, smooth_method = "none", bootstrap_n = 0L,
      n_permutations = 0L, bootstrap_seed = 1729L,
      permutation_seed = 1729L),
    notes = paste0("[native_rankmean_50bin] Built-in rank_mean on SCT_data; ",
      "exact Bdev pt; median in 50 equal-width bins; no promotion."))
  session <- result$session; rankmean_run <- result$run
} else log_event("Reusing native rank_mean operation")
native_curve <- rankmean_run$outputs$signature_table$table
native_dynamics <- rankmean_run$outputs$dynamics_table$table
native_fit <- rankmean_run$outputs$fit_summary$table
write_tsv(native_curve, file.path(stage_dir,
  "native_rankmean_50bin_trajectory.tsv"))
write_tsv(native_dynamics, file.path(stage_dir,
  "native_rankmean_dynamics.tsv"))
write_tsv(native_fit, file.path(stage_dir, "native_rankmean_fit_summary.tsv"))
stopifnot(setequal(unique(native_curve$feature_set), inventory$signature_id),
  nrow(native_fit) == 1L, identical(native_fit$score_method[[1L]], "rank_mean"))

# trajectory_geneset_dynamics does not export its per-cell rank_mean values.
# Reproduce its implementation exactly: rank the union of selected genes within
# each cell, then average those ranks for every feature set.
background_genes <- unique(members$feature_id)
stopifnot(all(background_genes %in% wb_feature_ids(bundle, "SCT_data")))
log_event("Reconstructing exact per-cell native rank_mean values for UMAP")
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data",
  background_genes, cells, missing = "zero"))
rank_matrix <- apply(expression, 2L, function(z)
  rank(z, ties.method = "average", na.last = "keep"))
if (is.null(dim(rank_matrix))) rank_matrix <- matrix(rank_matrix, ncol = 1L)
rownames(rank_matrix) <- rownames(expression); colnames(rank_matrix) <- cells
rankmean_matrix <- matrix(NA_real_, nrow = nrow(inventory), ncol = length(cells),
  dimnames = list(inventory$signature_id, cells))
for (i in seq_len(nrow(inventory))) {
  id <- inventory$signature_id[[i]]
  set_genes <- members$feature_id[members$signature_id == id]
  rankmean_matrix[i, ] <- colMeans(rank_matrix[set_genes, , drop = FALSE])
}
rm(expression, rank_matrix); invisible(gc())
stopifnot(all(is.finite(rankmean_matrix)), all(rankmean_matrix >= 1),
  all(rankmean_matrix <= length(background_genes)))
rankmean_scores <- do.call(rbind, lapply(seq_len(nrow(inventory)), function(i)
  data.frame(key = inventory$signature_id[[i]], cell_id = cells,
    value = as.numeric(rankmean_matrix[i, ]), stringsAsFactors = FALSE)))
saveRDS(rankmean_scores, file.path(stage_dir,
  "rankmean_cell_scores.rds"), compress = "gzip")
write_tsv_gz(rankmean_scores, file.path(stage_dir,
  "rankmean_cell_scores.tsv.gz"))

# Verify the reconstructed cell scores against every native 50-bin median.
pt_values <- metadata$pt[pt_use]
native_bins <- as.integer(cut(pt_values, breaks = 50L,
  include.lowest = TRUE, labels = FALSE))
validation_rows <- vector("list", nrow(inventory))
for (i in seq_len(nrow(inventory))) {
  id <- inventory$signature_id[[i]]
  expected <- aggregate(rankmean_matrix[id, pt_use],
    by = list(bin = native_bins), FUN = median)
  names(expected)[[2L]] <- "reconstructed_value"
  observed <- native_curve[native_curve$feature_set == id,
    c("bin", "value_raw"), drop = FALSE]
  merged <- merge(observed, expected, by = "bin", all.x = TRUE, sort = TRUE)
  validation_rows[[i]] <- data.frame(feature_set = id,
    bins_compared = nrow(merged),
    missing_reconstructed_bins = sum(!is.finite(merged$reconstructed_value)),
    maximum_absolute_difference = max(abs(merged$value_raw -
      merged$reconstructed_value)), stringsAsFactors = FALSE)
}
native_validation <- do.call(rbind, validation_rows)
write_tsv(native_validation, file.path(stage_dir,
  "native_rankmean_per_cell_reconstruction_validation.tsv"))
stopifnot(all(native_validation$bins_compared >= 5L),
  all(native_validation$missing_reconstructed_bins == 0L),
  max(native_validation$maximum_absolute_difference) < 1e-10)

# Use the exact existing rolling memberships for an apples-to-apples trajectory
# display of UCell and rank_mean. Gene-wise change retains its gene-wise IQR.
window_index <- lapply(window_membership, match, table = cells)
stopifnot(all(vapply(window_index, function(x) !anyNA(x), logical(1))))
rankmean_rows <- vector("list", nrow(inventory))
for (i in seq_len(nrow(inventory))) {
  id <- inventory$signature_id[[i]]; values <- rankmean_matrix[id, ]
  rankmean_rows[[i]] <- data.frame(feature_set = id,
    window_id = window_support$window_id,
    window_midpoint = window_support$window_midpoint,
    eligible_window = as_flag(window_support$eligible_window),
    score_median = vapply(window_index, function(ix) median(values[ix]), numeric(1)),
    score_q25 = vapply(window_index, function(ix)
      quantile(values[ix], 0.25, names = FALSE, type = 7), numeric(1)),
    score_q75 = vapply(window_index, function(ix)
      quantile(values[ix], 0.75, names = FALSE, type = 7), numeric(1)),
    stringsAsFactors = FALSE)
}
rankmean_trajectory <- do.call(rbind, rankmean_rows)
write_tsv(rankmean_trajectory, file.path(stage_dir,
  "rankmean_rolling_trajectory.tsv"))

method_summary <- do.call(rbind, lapply(seq_len(nrow(inventory)), function(i) {
  id <- inventory$signature_id[[i]]
  u <- ucell_scores$value[ucell_scores$key == id]
  r <- rankmean_matrix[id, ]
  g <- change_scores$value[change_scores$key == id]
  keep <- is.finite(g)
  data.frame(feature_set = id,
    spearman_rankmean_vs_ucell_all_cells = cor(r, u, method = "spearman"),
    spearman_genechange_vs_ucell_scoreable_cells = cor(g[keep], u[keep],
      method = "spearman"),
    spearman_genechange_vs_rankmean_scoreable_cells = cor(g[keep], r[keep],
      method = "spearman"),
    rankmean_min = min(r), rankmean_max = max(r),
    genechange_scoreable_cells = sum(keep), stringsAsFactors = FALSE)
}))
write_tsv(method_summary, file.path(stage_dir, "method_comparison.tsv"))

space <- wb_space(bundle, "umap.atac", cells = cells)
coordinates <- if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells))
  data.frame(cell_id = cells, x = as.numeric(space$x), y = as.numeric(space$y))
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  data.frame(cell_id = cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]))
}

positive_map <- function(values, title, subtitle, palette, legend_title) {
  finite <- is.finite(values)
  limits <- quantile(values[finite], c(0.01, 0.99), names = FALSE, type = 7)
  if (diff(limits) <= 0) limits <- range(values[finite])
  data <- coordinates; data$value <- values
  data <- data[order(ifelse(is.finite(data$value), data$value, -Inf)), ]
  ggplot(data, aes(x, y)) +
    geom_point(color = "#CDD2DA", size = 0.1, alpha = 0.65, stroke = 0) +
    geom_point(data = data[is.finite(data$value), ], aes(color = value),
      size = 0.24, alpha = 1, stroke = 0) +
    scale_color_gradientn(colours = palette, limits = limits,
      oob = scales::squish, name = legend_title) +
    coord_equal() + theme_void(base_size = 9) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
    labs(title = title, subtitle = subtitle, caption = "Raw score; colors clipped at p1/p99.")
}
change_map <- function(values, title, subtitle) {
  finite <- is.finite(values)
  cap <- quantile(abs(values[finite]), 0.99, names = FALSE, type = 7)
  if (!is.finite(cap) || cap <= 0) cap <- 0.01
  data <- coordinates; data$value <- values
  data <- data[order(ifelse(is.finite(data$value), abs(data$value), -Inf)), ]
  ggplot(data, aes(x, y)) +
    geom_point(color = "#B8BEC8", size = 0.1, alpha = 0.65, stroke = 0) +
    geom_point(data = data[is.finite(data$value), ], aes(color = value),
      size = 0.24, alpha = 1, stroke = 0) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-cap, cap), oob = scales::squish,
      name = "Mean rank\nchange") +
    coord_equal() + theme_void(base_size = 9) +
    theme(panel.background = element_rect(fill = "#E5E7EB", colour = NA),
      legend.position = "bottom", plot.title = element_text(face = "bold")) +
    labs(title = title, subtitle = subtitle,
      caption = "Blue: decrease; white: zero; red: increase; gray: unavailable.")
}
trajectory_plot <- function(data, value, low, high, title, y_label, color,
                            ribbon_caption) {
  data <- data[is.finite(data[[value]]), , drop = FALSE]
  ggplot(data, aes(x = window_midpoint, y = .data[[value]])) +
    geom_hline(yintercept = if (grepl("change", value)) 0 else NA_real_,
      color = "#6B7280", linewidth = 0.4) +
    geom_vline(data = change_points, aes(xintercept = pt_median),
      inherit.aes = FALSE, linewidth = 0.5, linetype = "dashed",
      color = "#374151") +
    geom_ribbon(aes(ymin = .data[[low]], ymax = .data[[high]]),
      fill = color, alpha = 0.2) +
    geom_line(color = color, linewidth = 0.85) +
    labs(title = title, x = "Existing Bdev pseudotime (pt)", y = y_label,
      caption = ribbon_caption) +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"))
}

log_event("Rendering 38 UCell/rank_mean/gene-change comparisons")
pdf_path <- file.path(figure_dir,
  "all_38_profiles_ucell_rankmean_genechange.pdf")
grDevices::pdf(pdf_path, width = 18, height = 10, onefile = TRUE)
figure_rows <- vector("list", nrow(inventory))
for (i in seq_len(nrow(inventory))) {
  row <- inventory[i, ]; id <- row$signature_id
  label <- paste(row$cytokine, toupper(row$direction), "-", row$final_size,
    "genes")
  u_values <- ucell_scores$value[ucell_scores$key == id]
  r_values <- rankmean_matrix[id, ]
  g_values <- change_scores$value[change_scores$key == id]
  u_curve <- ucell_trajectory[ucell_trajectory$feature_set == id, ]
  r_curve <- rankmean_trajectory[rankmean_trajectory$feature_set == id, ]
  g_curve <- change_trajectory[change_trajectory$feature_set == id, ]
  top <- positive_map(u_values, "Native UCell", label,
      c("#E8EAF6", "#3F51B5", "#D81B60"), "UCell") |
    positive_map(r_values, "Built-in rank_mean", label,
      c("#F3E5F5", "#7B1FA2", "#240046"), "Mean gene\nrank") |
    change_map(g_values, "Gene-wise neighborhood rank change", label)
  bottom <- trajectory_plot(u_curve, "score_median", "score_q25", "score_q75",
      "UCell trajectory", "UCell score", "#D81B60",
      "Median and cell-wise IQR; existing 5% windows.") |
    trajectory_plot(r_curve, "score_median", "score_q25", "score_q75",
      "rank_mean trajectory", "Mean within-cell gene rank", "#7B1FA2",
      "Median and cell-wise IQR; same 5% windows.") |
    trajectory_plot(g_curve, "score_mean", "score_q25", "score_q75",
      "Gene-wise change trajectory", "Mean percentile-rank change", "#00796B",
      "Mean across genes; ribbon is gene-wise IQR; five-window lag.")
  combined <- top / bottom + plot_layout(heights = c(1, 0.9)) +
    plot_annotation(title = paste("Three gene-set scoring views:", label),
      subtitle = paste("UCell: within-cell enrichment; rank_mean: mean rank among",
        length(background_genes), "set-union genes; gene-wise: k=30 change versus five windows earlier"),
      theme = theme(plot.title = element_text(size = 16, face = "bold")))
  print(combined)
  file <- paste0(sprintf("%02d_", i), slug(paste(row$cytokine,
    row$direction)), "_ucell_rankmean_genechange.png")
  ggsave(file.path(figure_dir, file), combined, width = 18, height = 10,
    dpi = 220, units = "in", bg = "white")
  comparison <- method_summary[method_summary$feature_set == id, ]
  figure_rows[[i]] <- data.frame(signature_id = id, cytokine = row$cytokine,
    direction = row$direction, signature_size = as.integer(row$final_size),
    spearman_rankmean_vs_ucell =
      comparison$spearman_rankmean_vs_ucell_all_cells,
    spearman_genechange_vs_ucell =
      comparison$spearman_genechange_vs_ucell_scoreable_cells,
    spearman_genechange_vs_rankmean =
      comparison$spearman_genechange_vs_rankmean_scoreable_cells,
    file = file, stringsAsFactors = FALSE)
}
grDevices::dev.off()
figure_inventory <- do.call(rbind, figure_rows)
write_tsv(figure_inventory, file.path(stage_dir, "figure_inventory.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
checks <- list(source_has_38_passing_profiles = nrow(inventory) == 38L,
  all_signatures_have_at_least_20_genes = all(inventory$final_size >= 20L),
  native_rankmean_operation = identical(rankmean_run$op_name,
    "trajectory_geneset_dynamics"),
  native_score_method_rankmean = identical(native_fit$score_method[[1L]],
    "rank_mean"),
  native_bins_retained_for_every_profile =
    all(table(native_curve$feature_set) >= 5L),
  per_cell_reconstruction_exact =
    max(native_validation$maximum_absolute_difference) < 1e-10,
  rank_background_is_set_union = setequal(background_genes,
    unique(members$feature_id)),
  exact_window_membership_reused = nrow(rankmean_trajectory) ==
    nrow(inventory) * nrow(window_support),
  exact_source_pt_preserved = identical(as.numeric(metadata$pt),
    as.numeric(change_scores$pt[change_scores$key == inventory$signature_id[[1L]]])),
  source_missing_pt_preserved = sum(is.na(metadata$pt)) == 14081L,
  every_figure_exists = all(file.exists(file.path(figure_dir,
    figure_inventory$file))), multipage_pdf_exists = file.exists(pdf_path),
  bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after), result_promotions = 0L)
stopifnot(all(unlist(checks[names(checks) != "result_promotions"])),
  checks$result_promotions == 0L,
  is.null(session$durable$promotion_receipt))
write_json(list(status = "passed", checks = checks),
  file.path(stage_dir, "checks.json"))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
provenance_inputs <- c(ucell_paths, change_paths,
  file.path(bundle_dir, "manifest.json"), script_path,
  "figure_scripts/bdev_reference/run_13_all_cytokine_rankmean_comparison_report.sh",
  "config/bdev_reference.yml")
input_manifest <- data.frame(path = normalizePath(provenance_inputs,
  mustWork = TRUE), sha256 = vapply(provenance_inputs, wb_source_sha256,
  character(1)), stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

html_rows <- vapply(seq_len(nrow(figure_inventory)), function(i) {
  x <- figure_inventory[i, ]
  paste0("<h2>", x$cytokine, " ", toupper(x$direction), "</h2>",
    "<p>Spearman: rank_mean vs UCell ",
    sprintf("%.3f", x$spearman_rankmean_vs_ucell),
    "; gene-change vs UCell ", sprintf("%.3f", x$spearman_genechange_vs_ucell),
    "; gene-change vs rank_mean ",
    sprintf("%.3f", x$spearman_genechange_vs_rankmean), ".</p>",
    '<img src="profile_figures/', x$file, '" alt="', x$cytokine, " ",
    x$direction, ' score comparison">')
}, character(1))
html <- c('<!doctype html><html><head><meta charset="utf-8">',
  '<title>UCell, native rank_mean, and gene-wise rank-change comparison</title>',
  '<style>body{font:16px system-ui,sans-serif;max-width:1700px;margin:32px auto;padding:0 22px;color:#172033}img{max-width:100%;height:auto;border:1px solid #ddd}p,li{line-height:1.5}code{overflow-wrap:anywhere}a{color:#1261a0}</style></head><body>',
  '<h1>UCell, native rank_mean, and gene-wise rank-change comparison</h1>',
  '<p>Native scWorkbench <code>trajectory_geneset_dynamics(score_method="rank_mean")</code> ranks the union of the 753 signature genes within each cell and averages member-gene ranks. Its native 50-bin output is retained. Because the native operation does not export per-cell scores, the documented formula was reconstructed for UMAP display and validated exactly against every native bin.</p>',
  '<p>The main trajectory figures place all methods on the same saved 5%-wide, 1%-step windows. UCell and rank_mean ribbons are cell-wise IQRs. The gene-wise method first performs k=30 smoothing, ranks each gene across cells, differences it against five windows earlier, and then averages genes; its ribbon is the gene-wise IQR. These scores answer different questions and their numeric scales are not directly comparable.</p>',
  '<p><a href="profile_figures/all_38_profiles_ucell_rankmean_genechange.pdf">Download 38-page PDF</a> · <a href="method_comparison.tsv">Method correlations</a> · <a href="rankmean_cell_scores.tsv.gz">Per-cell rank_mean scores</a> · <a href="native_rankmean_50bin_trajectory.tsv">Native trajectory output</a> · <a href="native_rankmean_per_cell_reconstruction_validation.tsv">Exact validation</a></p>',
  html_rows, '</body></html>')
writeLines(html, file.path(stage_dir, "REPORT.html"))
writeLines(c("# UCell, rank_mean, and gene-wise rank-change comparison", "",
  "[Open report](REPORT.html)", "",
  "All 38 profiles contain three UMAPs and three matched-window trajectories."),
  file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", promoted = FALSE,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  analysis_scope = paste("native scWorkbench rank_mean comparison with saved",
    "UCell and user-authorized gene-wise rank-change scores; no promotion"),
  source_ucell_report = ucell_report, source_gene_change_report = change_report,
  bundle = bundle_dir, bundle_fingerprint = fingerprint_before,
  profiles = 38L, cytokines = length(unique(inventory$cytokine)),
  rankmean = list(native_operation = "trajectory_geneset_dynamics",
    matrix = "SCT_data", score_method = "rank_mean",
    background = "union of all selected signature genes",
    background_gene_count = length(background_genes),
    native_trajectory = "50 equal-width bins; median; no smoothing",
    per_cell_export = "exact implementation reconstruction",
    maximum_native_validation_difference =
      max(native_validation$maximum_absolute_difference)),
  comparison_trajectory = list(window_width_fraction = 0.05,
    window_step_fraction = 0.01, min_cells = 20L,
    window_membership = "exact saved UCell report memberships"),
  gene_change = list(neighborhood_k = 30L, lag_windows = 5L,
    aggregation = "mean gene-level percentile-rank change"),
  defined_pt_cells = 12642L, missing_pt_cells = 14081L,
  result_promotions = 0L, native_session_workspace = workspace,
  multipage_pdf = file.path("profile_figures", basename(pdf_path)),
  figures = figure_inventory$file,
  report_sha256 = wb_source_sha256(file.path(stage_dir, "REPORT.html"))),
  file.path(stage_dir, "completion.json"))

if (!file.rename(stage_dir, output_dir)) stop("Failed to finalize output")
log_event("Completed rank_mean comparison report at", output_dir)
