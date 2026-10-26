#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 200)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 14_requested_cluster_driver_genes.R",
  "<bundle_dir> <preflight_dir> <trajectory_report_root> <cluster_dir>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
preflight_dir <- normalizePath(args[[2L]], mustWork = TRUE)
report_root <- normalizePath(args[[3L]], mustWork = TRUE)
cluster_dir <- normalizePath(args[[4L]], mustWork = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(splines)
})
read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE, stringsAsFactors = FALSE))
}
targets <- data.table(analysis_group = c("bp", "bp", "mf", "mf", "cc", "cc", "cc"),
  cluster_id = c("C06", "C12", "C01", "C06", "C05", "C06", "C09"))

cluster_membership <- fread(file.path(cluster_dir,
  "gene_set_cluster_membership.tsv"))
centroids <- fread(file.path(cluster_dir, "cluster_centroid_profiles.tsv"))
members <- read_tsv_gz(file.path(preflight_dir,
  "all_group_eligible_mapped_memberships.tsv.gz"))
target_sets <- merge(cluster_membership, targets,
  by = c("analysis_group", "cluster_id"))
stopifnot(nrow(target_sets) > 0L,
  nrow(unique(target_sets[, .(analysis_group, cluster_id)])) == nrow(targets))
target_members <- merge(members, target_sets[, .(analysis_group, feature_set,
  cluster_id)], by = c("analysis_group", "feature_set"))
set_sizes <- target_members[, .(set_size = uniqueN(feature_id)),
  by = .(analysis_group, cluster_id, feature_set)]
target_members <- merge(target_members, set_sizes,
  by = c("analysis_group", "cluster_id", "feature_set"), sort = FALSE)
genes <- sort(unique(target_members$feature_id))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
cells <- as.character(wb_cells(bundle, "all"))
metadata <- as.data.table(wb_meta(bundle, cells = cells))
stopifnot(identical(as.character(metadata$cell_id), cells))
pt_all <- as.numeric(metadata$pt)
pt_defined <- is.finite(pt_all)
pt <- pt_all[pt_defined]
pt_cells <- cells[pt_defined]
depth <- log1p(as.numeric(metadata$nCount_RNA[pt_defined]))
depth_z <- as.numeric(scale(depth))
depth_design <- cbind(1, ns(pt, df = 8L), depth_z)
depth_qr <- qr(depth_design)

curve_reference <- read_tsv_gz(file.path(report_root, "bp", "window_curves.tsv.gz"))
windows <- unique(curve_reference[, .(window_index, grid_window_index,
  window_center, window_left, window_right, n_cells)])
setorder(windows, window_index)
pt_limits <- range(pt)
pt_span <- diff(pt_limits)
window_width <- 0.05 * pt_span
windows[, reconstructed_center := pt_limits[[1L]] +
  (grid_window_index - 1L) * 0.01 * pt_span]
windows[, reconstructed_left := pmax(pt_limits[[1L]],
  reconstructed_center - window_width / 2)]
windows[, reconstructed_right := pmin(pt_limits[[2L]],
  reconstructed_center + window_width / 2)]
membership_index <- lapply(seq_len(nrow(windows)), function(i)
  which(pt >= windows$reconstructed_left[[i]] &
    pt <= windows$reconstructed_right[[i]]))
observed_window_counts <- as.integer(lengths(membership_index))
expected_window_counts <- as.integer(windows$n_cells)
if (!identical(observed_window_counts, expected_window_counts)) {
  bad <- which(observed_window_counts != expected_window_counts)
  stop("Window membership mismatch at: ", paste0(windows$window_index[bad],
    " expected=", expected_window_counts[bad], " observed=",
    observed_window_counts[bad], collapse = "; "))
}
weight_matrix <- matrix(0, nrow(windows), length(pt))
for (i in seq_len(nrow(windows)))
  weight_matrix[i, membership_index[[i]]] <- 1 / length(membership_index[[i]])

expression <- as.matrix(wb_matrix_read_features(bundle, "RNA_log1p_cp10k",
  genes, pt_cells, missing = "NA"))
stopifnot(setequal(rownames(expression), genes), all(is.finite(expression)))
expression <- expression[genes, , drop = FALSE]
coefficients <- qr.coef(depth_qr, t(expression))
depth_beta <- as.numeric(coefficients[ncol(depth_design), ])
adjusted <- expression - depth_beta %o% depth_z
gene_means <- rowMeans(adjusted)
gene_sd <- apply(adjusted, 1L, sd)
gene_curves <- adjusted %*% t(weight_matrix)
deviation_curves <- sweep(gene_curves, 1L, gene_means, "-")
finite_population_factor <- sqrt((1 - windows$n_cells / ncol(adjusted)) /
  windows$n_cells)
gene_z <- deviation_curves / (gene_sd %o% finite_population_factor)
gene_z[!is.finite(gene_z)] <- NA_real_

ranking_rows <- list(); cursor <- 0L
for (i in seq_len(nrow(targets))) {
  g <- targets$analysis_group[[i]]
  cid <- targets$cluster_id[[i]]
  sets <- target_sets[analysis_group == g & cluster_id == cid, feature_set]
  cluster_members <- target_members[analysis_group == g & cluster_id == cid]
  cluster_genes <- sort(unique(cluster_members$feature_id))
  n_sets <- length(sets)
  gene_weights <- cluster_members[, .(
    gene_set_membership_count = uniqueN(feature_set),
    gene_set_membership_fraction = uniqueN(feature_set) / n_sets,
    mean_set_score_weight = sum(1 / set_size) / n_sets), by = feature_id]
  centroid <- centroids[analysis_group == g & cluster_id == cid]
  setorder(centroid, window_index)
  stopifnot(identical(centroid$window_index, windows$window_index))
  peak_window <- which.max(centroid$mean_complete_z)
  gi <- match(cluster_genes, genes)
  correlations <- vapply(gi, function(j) {
    x <- gene_z[j, ]
    if (all(!is.finite(x)) || sd(x, na.rm = TRUE) == 0) return(NA_real_)
    cor(x, centroid$mean_complete_z, use = "pairwise.complete.obs")
  }, numeric(1))
  driver <- gene_weights[match(cluster_genes, feature_id)]
  driver[, `:=`(
    analysis_group = g, cluster_id = cid, n_gene_sets_in_cluster = n_sets,
    centroid_peak_pt = windows$window_center[[peak_window]],
    gene_z_at_centroid_peak = gene_z[gi, peak_window],
    gene_profile_correlation = correlations,
    weighted_peak_contribution = mean_set_score_weight *
      deviation_curves[gi, peak_window])]
  driver[, aligned_driver_score := weighted_peak_contribution *
    pmax(gene_profile_correlation, 0)]
  positive_total <- sum(pmax(driver$aligned_driver_score, 0), na.rm = TRUE)
  driver[, positive_driver_fraction := if (positive_total > 0)
    pmax(aligned_driver_score, 0) / positive_total else NA_real_]
  setorder(driver, -aligned_driver_score, -weighted_peak_contribution,
    -gene_z_at_centroid_peak, feature_id, na.last = TRUE)
  driver[, driver_rank := seq_len(.N)]
  cursor <- cursor + 1L
  ranking_rows[[cursor]] <- driver
}
rankings <- rbindlist(ranking_rows, use.names = TRUE, fill = TRUE)
setcolorder(rankings, c("analysis_group", "cluster_id", "driver_rank",
  "feature_id", "n_gene_sets_in_cluster", "gene_set_membership_count",
  "gene_set_membership_fraction", "mean_set_score_weight",
  "centroid_peak_pt", "gene_z_at_centroid_peak",
  "gene_profile_correlation", "weighted_peak_contribution",
  "aligned_driver_score", "positive_driver_fraction"))
fwrite(rankings, file.path(cluster_dir,
  "requested_cluster_gene_driver_rankings.tsv"), sep = "\t", quote = FALSE,
  na = "NA")
top <- rankings[driver_rank <= 15L]
fwrite(top, file.path(cluster_dir,
  "requested_cluster_top15_driver_genes.tsv"), sep = "\t", quote = FALSE,
  na = "NA")

group_label <- c(bp = "GO BP", mf = "GO MF", cc = "GO CC")
report <- c("# Requested cluster driver genes", "",
  "Drivers are descriptive contributions, not independent significance tests. Genes are ranked by their depth-adjusted peak contribution to the mean cluster gene-set score, multiplied by positive correlation with the cluster's complete z-score centroid. Gene-set recurrence is weighted by inverse set size.", "")
for (i in seq_len(nrow(targets))) {
  g <- targets$analysis_group[[i]]; cid <- targets$cluster_id[[i]]
  x <- rankings[analysis_group == g & cluster_id == cid & driver_rank <= 10L]
  report <- c(report, paste0("## ", group_label[[g]], " ", cid), "",
    paste0("Cluster sets: ", unique(x$n_gene_sets_in_cluster),
      "; centroid peak pt: ", signif(unique(x$centroid_peak_pt), 6)), "",
    "| Rank | Gene | Set membership | Peak gene z | Profile correlation | Positive driver fraction |",
    "|---:|---|---:|---:|---:|---:|",
    vapply(seq_len(nrow(x)), function(j) paste0("| ", x$driver_rank[[j]],
      " | ", x$feature_id[[j]], " | ", x$gene_set_membership_count[[j]],
      "/", x$n_gene_sets_in_cluster[[j]], " | ",
      signif(x$gene_z_at_centroid_peak[[j]], 4), " | ",
      signif(x$gene_profile_correlation[[j]], 4), " | ",
      signif(x$positive_driver_fraction[[j]], 4), " |"), character(1)), "")
}
writeLines(report, file.path(cluster_dir,
  "REQUESTED_CLUSTER_DRIVER_GENES.md"))
fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
wb_close_bundle(bundle)
