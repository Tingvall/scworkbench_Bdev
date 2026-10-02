#!/usr/bin/env Rscript

# External, read-only mapping QC for the completed scWorkbench PeakVI stage.
# RANN is used because the registered scWorkbench neighbor operation is not
# suitable for reference-restricted query-to-reference searches at this scale.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: 11_peakvi_mapping_qc.R BUNDLE_DIR STAGE_BUNDLE MAPPING_JSON REPORT_DIR",
       call. = FALSE)
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
stage_dir <- normalizePath(args[[2L]], mustWork = TRUE)
mapping_json <- normalizePath(args[[3L]], mustWork = TRUE)
report_dir <- normalizePath(args[[4L]], mustWork = FALSE)
if (file.exists(report_dir) || dir.exists(report_dir)) {
  stop("QC report directory already exists; refusing to overwrite: ", report_dir,
       call. = FALSE)
}

options(width = 180, warn = 1)
suppressPackageStartupMessages({
  library(scworkbench)
  library(RANN)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

seed <- 1729L
ks <- c(15L, 30L, 50L)
primary_k <- 30L

write_json <- function(x, path) jsonlite::write_json(
  x, path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null",
  digits = NA
)

log_event <- function(...) {
  msg <- paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), paste(..., collapse = "\t"), sep = "\t")
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(report_dir, "progress.log"), append = TRUE)
}

weighted_categorical <- function(index, distance, labels, k) {
  index <- index[, seq_len(k), drop = FALSE]
  distance <- distance[, seq_len(k), drop = FALSE]
  labs <- matrix(labels[index], nrow = nrow(index))
  sigma <- distance[, k]
  sigma[sigma <= 0 | !is.finite(sigma)] <- 1
  weights <- exp(-(distance / sigma)^2)
  lev <- sort(unique(labels[!is.na(labels)]))
  scores <- vapply(lev, function(z) rowSums(weights * (labs == z), na.rm = TRUE), numeric(nrow(index)))
  if (is.null(dim(scores))) scores <- matrix(scores, ncol = 1L, dimnames = list(NULL, lev))
  colnames(scores) <- lev
  totals <- rowSums(scores)
  ord <- t(apply(scores, 1L, order, decreasing = TRUE))
  winner <- lev[ord[, 1L]]
  top1 <- scores[cbind(seq_len(nrow(scores)), ord[, 1L])]
  top2 <- if (ncol(scores) > 1L) scores[cbind(seq_len(nrow(scores)), ord[, 2L])] else rep(0, nrow(scores))
  p <- scores / pmax(totals, .Machine$double.eps)
  entropy <- -rowSums(ifelse(p > 0, p * log(p), 0)) / log(ncol(scores))
  if (ncol(scores) == 1L) entropy[] <- 0
  data.table(prediction = winner, confidence = top1 / totals,
             margin = (top1 - top2) / totals, entropy = entropy)
}

weighted_quantiles <- function(index, distance, values, k, probs = c(.25, .5, .75)) {
  index <- index[, seq_len(k), drop = FALSE]
  distance <- distance[, seq_len(k), drop = FALSE]
  vals <- matrix(values[index], nrow = nrow(index))
  sigma <- distance[, k]
  sigma[sigma <= 0 | !is.finite(sigma)] <- 1
  weights <- exp(-(distance / sigma)^2)
  out <- matrix(NA_real_, nrow(index), length(probs))
  for (i in seq_len(nrow(index))) {
    keep <- is.finite(vals[i, ]) & is.finite(weights[i, ]) & weights[i, ] > 0
    if (!any(keep)) next
    o <- order(vals[i, keep])
    v <- vals[i, keep][o]
    w <- weights[i, keep][o]
    cw <- cumsum(w) / sum(w)
    out[i, ] <- v[vapply(probs, function(p) which(cw >= p)[1L], integer(1))]
  }
  colnames(out) <- paste0("q", probs * 100)
  out
}

neighbor_table <- function(query_ids, reference_ids, result) {
  k <- ncol(result$nn.idx)
  data.table(
    query_cell_id = rep(query_ids, each = k), rank = rep(seq_len(k), times = length(query_ids)),
    reference_cell_id = reference_ids[as.vector(t(result$nn.idx))],
    distance = as.vector(t(result$nn.dists))
  )
}

categorical_metrics <- function(truth, pred) {
  truth <- as.character(truth); pred <- as.character(pred)
  lev <- sort(unique(truth))
  per <- data.table(label = lev, n = vapply(lev, function(z) sum(truth == z), integer(1)),
                    accuracy = vapply(lev, function(z) mean(pred[truth == z] == z), numeric(1)))
  list(accuracy = mean(pred == truth), balanced_accuracy = mean(per$accuracy), per_label = per)
}

plot_save <- function(p, name, width = 8, height = 6) {
  ggsave(file.path(report_dir, "figures", name), p, width = width, height = height,
         dpi = 160, limitsize = FALSE)
}

main <- function() {
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(report_dir, "figures"), showWarnings = FALSE)
  on.exit({
    if (!file.exists(file.path(report_dir, "completion.json"))) {
      writeLines("incomplete", file.path(report_dir, "INCOMPLETE"))
    }
  }, add = TRUE)

  log_event("open", "Opening immutable source and staged bundles")
  source_bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
  staged <- wb_open_bundle(stage_dir, validate = FALSE)
  on.exit(try(wb_close_bundle(source_bundle), silent = TRUE), add = TRUE)
  on.exit(try(wb_close_bundle(staged), silent = TRUE), add = TRUE)
  manifest_before <- wb_manifest(source_bundle)
  fingerprint <- manifest_before$current_state$fingerprint
  stopifnot(identical(manifest_before$bundle_id, "bdev_integrated_ko_mapping_v1"))

  ref_ids <- wb_cells(source_bundle, "set:bdev_reference")
  query_ids <- wb_cells(source_bundle, "set:ko_query")
  all_ids <- c(ref_ids, query_ids)
  md <- as.data.table(wb_meta(source_bundle, cells = all_ids))
  stopifnot(identical(as.character(md$cell_id), all_ids), length(ref_ids) == 26723L,
            length(query_ids) == 23252L)
  ref_md <- md[seq_along(ref_ids)]
  query_md <- md[length(ref_ids) + seq_along(query_ids)]
  reference_pt_original <- ref_md$reference_pt

  joint <- as.matrix(wb_space(staged, "peakvi.full.joint", cells = all_ids))
  umap <- wb_space(staged, "umap.peakvi.full.joint", cells = all_ids)
  stopifnot(identical(rownames(joint), all_ids), all(is.finite(joint)),
            identical(as.character(umap$cell_id), all_ids))
  ref_x <- joint[seq_along(ref_ids), , drop = FALSE]
  query_x <- joint[length(ref_ids) + seq_along(query_ids), , drop = FALSE]

  log_event("query_neighbors", "Exact query-to-Bdev k=50 RANN search")
  query_nn <- RANN::nn2(ref_x, query_x, k = max(ks), searchtype = "standard", eps = 0)
  fwrite(neighbor_table(query_ids, ref_ids, query_nn),
         file.path(report_dir, "query_bdev_neighbors_k50.tsv.gz"), sep = "\t")

  log_event("cross_validation", "Five-fold cluster-stratified Bdev calibration")
  set.seed(seed)
  folds <- integer(length(ref_ids))
  for (z in sort(unique(ref_md$reference_peak_cluster))) {
    ii <- which(ref_md$reference_peak_cluster == z)
    folds[ii[sample.int(length(ii))]] <- rep(seq_len(5L), length.out = length(ii))
  }
  cv_idx <- matrix(NA_integer_, length(ref_ids), max(ks))
  cv_dist <- matrix(NA_real_, length(ref_ids), max(ks))
  for (fold in seq_len(5L)) {
    test <- which(folds == fold); train <- which(folds != fold)
    ans <- RANN::nn2(ref_x[train, , drop = FALSE], ref_x[test, , drop = FALSE],
                     k = max(ks), searchtype = "standard", eps = 0)
    cv_idx[test, ] <- matrix(train[ans$nn.idx], nrow = length(test))
    cv_dist[test, ] <- ans$nn.dists
    log_event("cross_validation", paste0("fold_", fold, "_complete"))
  }

  cv <- data.table(cell_id = ref_ids, fold = folds,
                   truth_cluster = ref_md$reference_peak_cluster,
                   truth_phase = ref_md$reference_phase,
                   nearest_distance = cv_dist[, 1L], median30_distance = apply(cv_dist[, 1:30], 1L, median),
                   kth30_distance = cv_dist[, 30L])
  query <- data.table(cell_id = query_ids,
                      source_cell_id = query_md$source_cell_id,
                      sample_id = query_md$sample_id, genotype = query_md$genotype,
                      knockout_target = query_md$knockout_target, tissue = query_md$tissue,
                      nearest_distance = query_nn$nn.dists[, 1L],
                      median30_distance = apply(query_nn$nn.dists[, 1:30], 1L, median),
                      kth30_distance = query_nn$nn.dists[, 30L])

  log_event("label_transfer", "Cluster and phase transfer at k=15,30,50")
  for (k in ks) {
    cc <- weighted_categorical(cv_idx, cv_dist, ref_md$reference_peak_cluster, k)
    cp <- weighted_categorical(cv_idx, cv_dist, ref_md$reference_phase, k)
    qc <- weighted_categorical(query_nn$nn.idx, query_nn$nn.dists, ref_md$reference_peak_cluster, k)
    qp <- weighted_categorical(query_nn$nn.idx, query_nn$nn.dists, ref_md$reference_phase, k)
    setnames(cc, names(cc), paste0(c("predicted_cluster", "cluster_confidence", "cluster_margin", "cluster_entropy"), "_k", k))
    setnames(cp, names(cp), paste0(c("predicted_phase", "phase_confidence", "phase_margin", "phase_entropy"), "_k", k))
    setnames(qc, names(qc), paste0(c("predicted_cluster", "cluster_confidence", "cluster_margin", "cluster_entropy"), "_k", k))
    setnames(qp, names(qp), paste0(c("predicted_phase", "phase_confidence", "phase_margin", "phase_entropy"), "_k", k))
    cv <- cbind(cv, cc, cp); query <- cbind(query, qc, qp)
  }
  cv[, cluster_k_agreement := predicted_cluster_k15 == predicted_cluster_k30 & predicted_cluster_k30 == predicted_cluster_k50]
  cv[, phase_k_agreement := predicted_phase_k15 == predicted_phase_k30 & predicted_phase_k30 == predicted_phase_k50]
  query[, cluster_k_agreement := predicted_cluster_k15 == predicted_cluster_k30 & predicted_cluster_k30 == predicted_cluster_k50]
  query[, phase_k_agreement := predicted_phase_k15 == predicted_phase_k30 & predicted_phase_k30 == predicted_phase_k50]

  q95 <- unname(quantile(cv$median30_distance, .95, type = 8))
  q99 <- unname(quantile(cv$median30_distance, .99, type = 8))
  query[, support_status := fifelse(median30_distance <= q95, "supported",
                            fifelse(median30_distance <= q99, "borderline", "out_of_reference"))]

  log_event("pseudotime", "Exact query-to-defined-pseudotime Bdev search and calibration")
  pt_ref <- which(is.finite(ref_md$reference_pt))
  pt_nn <- RANN::nn2(ref_x[pt_ref, , drop = FALSE], query_x, k = max(ks),
                     searchtype = "standard", eps = 0)
  fwrite(neighbor_table(query_ids, ref_ids[pt_ref], pt_nn),
         file.path(report_dir, "query_defined_pt_neighbors_k50.tsv.gz"), sep = "\t")
  for (k in ks) {
    pq <- weighted_quantiles(pt_nn$nn.idx, pt_nn$nn.dists, ref_md$reference_pt[pt_ref], k)
    query[[paste0("mapped_pt_raw_k", k)]] <- pq[, "q50"]
    if (k == primary_k) {
      query[, `:=`(mapped_pt_q25 = pq[, "q25"], mapped_pt_q75 = pq[, "q75"])]
    }
  }
  query[, mapped_pt_sensitivity_range := do.call(pmax, c(.SD, na.rm = TRUE)) - do.call(pmin, c(.SD, na.rm = TRUE)),
        .SDcols = paste0("mapped_pt_raw_k", ks)]
  query[, mapped_pt_iqr := mapped_pt_q75 - mapped_pt_q25]
  query[, mapped_pt := fifelse(support_status == "out_of_reference", NA_real_, mapped_pt_raw_k30)]

  pt_cv <- data.table(cell_id = ref_ids[pt_ref], truth_pt = ref_md$reference_pt[pt_ref], fold = folds[pt_ref])
  pt_cv_idx <- matrix(NA_integer_, length(pt_ref), max(ks)); pt_cv_dist <- matrix(NA_real_, length(pt_ref), max(ks))
  for (fold in seq_len(5L)) {
    test <- which(folds[pt_ref] == fold); train <- which(folds[pt_ref] != fold)
    ans <- RANN::nn2(ref_x[pt_ref[train], , drop = FALSE], ref_x[pt_ref[test], , drop = FALSE],
                     k = max(ks), searchtype = "standard", eps = 0)
    pt_cv_idx[test, ] <- matrix(train[ans$nn.idx], nrow = length(test))
    pt_cv_dist[test, ] <- ans$nn.dists
  }
  for (k in ks) pt_cv[[paste0("predicted_pt_k", k)]] <- weighted_quantiles(
    pt_cv_idx, pt_cv_dist, ref_md$reference_pt[pt_ref], k)[, "q50"]
  pt_cv[, sensitivity_range := do.call(pmax, c(.SD, na.rm = TRUE)) - do.call(pmin, c(.SD, na.rm = TRUE)),
        .SDcols = paste0("predicted_pt_k", ks)]

  fwrite(query, file.path(report_dir, "query_mapping.tsv.gz"), sep = "\t")
  fwrite(cv, file.path(report_dir, "reference_cv_predictions.tsv.gz"), sep = "\t")
  fwrite(pt_cv, file.path(report_dir, "reference_pt_cv_predictions.tsv.gz"), sep = "\t")

  cluster_metrics <- categorical_metrics(cv$truth_cluster, cv$predicted_cluster_k30)
  phase_metrics <- categorical_metrics(cv$truth_phase, cv$predicted_phase_k30)
  fwrite(cluster_metrics$per_label, file.path(report_dir, "cluster_cv_per_label.tsv"), sep = "\t")
  cluster_confusion <- as.data.table(as.table(table(truth = cv$truth_cluster, prediction = cv$predicted_cluster_k30)))
  fwrite(cluster_confusion, file.path(report_dir, "cluster_cv_confusion.tsv"), sep = "\t")
  err <- pt_cv$predicted_pt_k30 - pt_cv$truth_pt
  pt_metrics <- list(n = nrow(pt_cv), spearman = cor(pt_cv$truth_pt, pt_cv$predicted_pt_k30, method = "spearman"),
                     pearson = cor(pt_cv$truth_pt, pt_cv$predicted_pt_k30), mae = mean(abs(err)),
                     rmse = sqrt(mean(err^2)), median_absolute_error = median(abs(err)))

  sample_summary <- query[, .(n_cells = .N, supported_fraction = mean(support_status == "supported"),
                              borderline_fraction = mean(support_status == "borderline"),
                              out_of_reference_fraction = mean(support_status == "out_of_reference"),
                              cluster_k_agreement = mean(cluster_k_agreement), phase_k_agreement = mean(phase_k_agreement),
                              mapped_pt_fraction = mean(is.finite(mapped_pt)), mapped_pt_median = median(mapped_pt, na.rm = TRUE),
                              mapped_pt_q25 = quantile(mapped_pt, .25, na.rm = TRUE), mapped_pt_q75 = quantile(mapped_pt, .75, na.rm = TRUE)),
                          by = .(sample_id, genotype, knockout_target, tissue)]
  fwrite(sample_summary, file.path(report_dir, "query_summary_by_sample.tsv"), sep = "\t")
  fwrite(query[, .N, by = .(sample_id, genotype, knockout_target, predicted_cluster_k30)][order(sample_id, genotype, predicted_cluster_k30)],
         file.path(report_dir, "query_cluster_counts.tsv"), sep = "\t")

  log_event("figures", "Rendering mapping QC figures")
  plot_dt <- data.table(cell_id = all_ids, x = umap$x, y = umap$y,
                        role = c(rep("Bdev reference", length(ref_ids)), rep("KO query", length(query_ids))),
                        reference_pt = c(ref_md$reference_pt, rep(NA_real_, length(query_ids))),
                        mapped_pt = c(rep(NA_real_, length(ref_ids)), query$mapped_pt),
                        cluster = c(ref_md$reference_peak_cluster, query$predicted_cluster_k30),
                        support = c(rep("reference", length(ref_ids)), query$support_status))
  base_theme <- theme_classic(base_size = 11)
  plot_save(ggplot(plot_dt, aes(x, y, color = role)) + geom_point(size = .12, alpha = .55) + base_theme,
            "joint_umap_role.png")
  pt_plot <- copy(plot_dt); pt_plot[, display_pt := fifelse(role == "Bdev reference", reference_pt, mapped_pt)]
  plot_save(ggplot(pt_plot[is.finite(display_pt)], aes(x, y, color = display_pt)) + geom_point(size = .16) +
              scale_color_viridis_c(name = "pt") + facet_wrap(~role) + base_theme,
            "joint_umap_pseudotime.png", 10, 5)
  plot_save(ggplot(plot_dt[role == "KO query"], aes(x, y, color = cluster)) + geom_point(size = .2) +
              guides(color = guide_legend(override.aes = list(size = 2))) + base_theme,
            "query_umap_predicted_cluster.png")
  plot_save(ggplot(plot_dt[role == "KO query"], aes(x, y, color = support)) + geom_point(size = .2) + base_theme,
            "query_umap_support.png")
  plot_save(ggplot(pt_cv, aes(truth_pt, predicted_pt_k30)) + geom_point(size = .25, alpha = .25) +
              geom_abline(slope = 1, intercept = 0, color = "red") + coord_equal() + base_theme,
            "reference_pt_cv.png")
  plot_save(ggplot(cluster_confusion, aes(prediction, truth, fill = N)) + geom_tile() +
              scale_fill_viridis_c() + theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 90)),
            "cluster_cv_confusion.png", 9, 8)
  plot_save(ggplot(query[is.finite(mapped_pt)], aes(sample_id, mapped_pt, fill = genotype)) +
              geom_violin(scale = "width", trim = TRUE) + coord_flip() + base_theme,
            "mapped_pt_by_sample.png", 8, 6)

  stopifnot(identical(reference_pt_original, ref_md$reference_pt))
  manifest_after <- wb_manifest(source_bundle)
  stopifnot(identical(fingerprint, manifest_after$current_state$fingerprint))
  mapping_receipt <- jsonlite::read_json(mapping_json, simplifyVector = TRUE)
  metrics <- list(
    cluster = c(list(n = nrow(cv)), cluster_metrics[c("accuracy", "balanced_accuracy")],
                    list(k_agreement = mean(cv$cluster_k_agreement))),
    phase = c(list(n = nrow(cv)), phase_metrics[c("accuracy", "balanced_accuracy")],
                  list(k_agreement = mean(cv$phase_k_agreement))),
    pseudotime = pt_metrics,
    distance_calibration = list(statistic = "median Euclidean distance to first 30 Bdev neighbors",
                                q95 = q95, q99 = q99,
                                note = "Five-fold held-out Bdev calibration uses 80% reference candidates per fold and is conservative relative to the full reference search."),
    query_support = as.list(table(query$support_status)),
    query_cluster_k_agreement = mean(query$cluster_k_agreement),
    query_phase_k_agreement = mean(query$phase_k_agreement)
  )
  write_json(metrics, file.path(report_dir, "metrics.json"))
  methods <- list(
    status = "complete_unpromoted", scope = "read-only PeakVI mapping QC and Bdev-restricted transfer",
    implementation = "approved external RANN calculation because registered scWorkbench neighbors are brute-force and do not support this reference-restricted workflow",
    nearest_neighbors = list(package = "RANN", version = as.character(packageVersion("RANN")), exact = TRUE,
                             eps = 0, distance = "Euclidean", k = ks, primary_k = primary_k),
    cross_validation = list(folds = 5L, seed = seed, stratified_by = "reference_peak_cluster"),
    weights = "exp(-(distance/kth_distance)^2)", pseudotime = "weighted median; weighted q25/q75 at k=30",
    ood = "median first-30 distance; supported <= CV q95, borderline <= CV q99, out_of_reference > CV q99",
    inputs = list(bundle = bundle_dir, stage_bundle = stage_dir, full_mapping_report = mapping_json,
                  bundle_fingerprint = fingerprint, full_mapping_runs = mapping_receipt$runs),
    preservation = list(reference_pt_unchanged = TRUE, source_bundle_unchanged = TRUE, promoted = FALSE)
  )
  write_json(methods, file.path(report_dir, "methods.json"))
  writeLines(c(
    "# Bdev-restricted PeakVI mapping QC", "",
    "Status: complete and unpromoted.", "",
    sprintf("The exact k-nearest-neighbor analysis covers %s KO query cells and %s Bdev reference cells in the 20-dimensional PeakVI space.", format(length(query_ids), big.mark=","), format(length(ref_ids), big.mark=",")),
    sprintf("Five-fold Bdev calibration gave primary k=30 cluster accuracy %.3f (balanced %.3f), phase accuracy %.3f, and pseudotime Spearman rho %.3f.", cluster_metrics$accuracy, cluster_metrics$balanced_accuracy, phase_metrics$accuracy, pt_metrics$spearman),
    sprintf("OOD support from held-out Bdev median-30 distances: q95 %.6g and q99 %.6g. Query calls: %s.", q95, q99, paste(names(table(query$support_status)), as.integer(table(query$support_status)), collapse=", ")),
    "", "Mapped pseudotime is a weighted median over Bdev cells with defined reference_pt. Raw estimates and k-sensitivity are retained; the primary mapped_pt is withheld for out-of-reference cells. Original Bdev reference_pt values and missing entries were not changed.",
    "", "All outputs are report artifacts only. Nothing was promoted into the source or staged bundle."
  ), file.path(report_dir, "REPORT.md"))
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
  write_json(list(status = "complete_unpromoted", completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                  cells = list(reference = length(ref_ids), query = length(query_ids), reference_pt_defined = length(pt_ref)),
                  bundle_fingerprint_unchanged = TRUE, reference_pt_preserved = TRUE, metrics = metrics),
             file.path(report_dir, "completion.json"))
  log_event("complete", "Mapping QC complete and unpromoted")
}

main()
