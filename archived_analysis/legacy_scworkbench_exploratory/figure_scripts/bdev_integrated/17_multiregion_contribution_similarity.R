#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 17_multiregion_contribution_similarity.R",
  "<bundle_dir> <memberships.tsv.gz> <source_root> <output_dir> <primary_threshold>"))

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
membership_path <- normalizePath(args[[2L]], mustWork = TRUE)
source_root <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
primary_threshold <- as.numeric(args[[5L]])
thresholds <- c(0.40, 0.50, 0.60, 0.70, 0.80)
stopifnot(primary_threshold %in% thresholds)
primary_label <- formatC(primary_threshold, format = "f", digits = 2)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(jsonlite)
})
read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA")
write_tsv_gz <- function(x, path) {
  if (nrow(x)) {
    fwrite(x, path, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
  } else {
    con <- gzfile(path, "wt"); on.exit(close(con), add = TRUE)
    writeLines(paste(names(x), collapse = "\t"), con)
  }
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ..., "\n")
  flush.console()
}
weighted_jaccard <- function(a_idx, a_val, b_idx, b_val) {
  common <- intersect(a_idx, b_idx)
  if (!length(common)) return(0)
  ai <- match(common, a_idx); bi <- match(common, b_idx)
  overlap <- sum(pmin(a_val[ai], b_val[bi]))
  overlap / (2 - overlap)
}
greedy_assign <- function(sets, edges, threshold) {
  n <- nrow(sets)
  rank_order <- order(-sets$best_peak_z, sets$best_q,
    -sets$best_prominence, sets$feature_set)
  assigned <- integer(n); similarity_to_rep <- rep(NA_real_, n)
  representatives <- integer(n); group_number <- 0L
  use <- edges[full_similarity >= threshold]
  adjacency <- vector("list", n)
  if (nrow(use)) for (k in seq_len(nrow(use))) {
    a <- use$i[[k]]; b <- use$j[[k]]; s <- use$full_similarity[[k]]
    adjacency[[a]] <- rbind(adjacency[[a]], c(b, s))
    adjacency[[b]] <- rbind(adjacency[[b]], c(a, s))
  }
  for (rep_i in rank_order) {
    if (assigned[[rep_i]] != 0L) next
    group_number <- group_number + 1L
    representatives[[group_number]] <- rep_i
    assigned[[rep_i]] <- group_number; similarity_to_rep[[rep_i]] <- 1
    neighbors <- adjacency[[rep_i]]
    if (!is.null(neighbors)) {
      neighbors <- as.matrix(neighbors)
      open <- neighbors[, 1L][assigned[as.integer(neighbors[, 1L])] == 0L]
      if (length(open)) {
        assigned[as.integer(open)] <- group_number
        similarity_to_rep[as.integer(open)] <-
          neighbors[match(open, neighbors[, 1L]), 2L]
      }
    }
  }
  rep_by_group <- setNames(representatives[seq_len(group_number)],
    seq_len(group_number))
  representative_i <- unname(rep_by_group[as.character(assigned)])
  data.table(set_row = seq_len(n), redundancy_group = assigned,
    representative_row = representative_i,
    similarity_to_representative = similarity_to_rep,
    retained_representative = seq_len(n) == representative_i)
}

groups <- c("cytokines", "bp", "mf", "cc")
members_all <- read_tsv_gz(membership_path)
bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
cells <- as.character(wb_cells(bundle, "all"))
metadata <- as.data.table(wb_meta(bundle, cells = cells))
stopifnot(identical(as.character(metadata$cell_id), cells),
  all(is.finite(metadata$nCount_RNA)), all(metadata$nCount_RNA > 0))
pt_all <- as.numeric(metadata$pt); pt_defined <- is.finite(pt_all)
pt <- pt_all[pt_defined]; pt_cells <- cells[pt_defined]
stopifnot(length(pt) == 12642L, sum(!pt_defined) == 14081L)
log_depth <- log1p(as.numeric(metadata$nCount_RNA)[pt_defined])
depth_z <- (log_depth - mean(log_depth)) / sd(log_depth)

# Exact shared 20-cell analysis grid.
pt_limits <- range(pt); pt_span <- diff(pt_limits)
centers <- pt_limits[[1L]] + seq(0, 1, by = 0.01) * pt_span
window_width <- 0.05 * pt_span
windows <- data.table(grid_window_index = seq_along(centers),
  window_center = centers,
  window_left = pmax(pt_limits[[1L]], centers - window_width / 2),
  window_right = pmin(pt_limits[[2L]], centers + window_width / 2))
window_members <- lapply(seq_len(nrow(windows)), function(i)
  which(pt >= windows$window_left[[i]] & pt <= windows$window_right[[i]]))
windows[, n_cells := lengths(window_members)]
windows <- windows[n_cells >= 20L]
window_members <- lapply(seq_len(nrow(windows)), function(i)
  which(pt >= windows$window_left[[i]] & pt <= windows$window_right[[i]]))
weight_matrix <- matrix(0, nrow(windows), length(pt))
for (i in seq_len(nrow(windows)))
  weight_matrix[i, window_members[[i]]] <- 1 / length(window_members[[i]])
stopifnot(nrow(windows) == 97L, min(windows$n_cells) == 23L)

all_summary <- list(); all_primary_mapping <- list()
all_retained_sets <- list(); all_retained_regions <- list()
all_signature_rows <- list(); group_completion <- list()
for (group in groups) {
  log_event("Building multi-region gene-by-pseudotime maps for", group)
  group_dir <- file.path(source_root, group)
  regions <- fread(file.path(group_dir, "significant_positive_regions.tsv"))
  transforms <- read_tsv_gz(file.path(group_dir,
    "gene_standardization_diagnostics.tsv.gz"))[standardizable == TRUE]
  setorder(transforms, feature_id); genes <- transforms$feature_id
  final_sets <- unique(regions$feature_set)
  members <- members_all[analysis_group == group &
    feature_set %chin% final_sets & feature_id %chin% genes]
  expression <- as.matrix(wb_matrix_read_features(bundle,
    "RNA_log1p_cp10k", genes, pt_cells, missing = "NA"))
  adjusted <- expression - tcrossprod(transforms$conditional_depth_beta, depth_z)
  standardized <- (adjusted - transforms$adjusted_center) /
    transforms$adjusted_scale
  gene_windows <- standardized %*% t(weight_matrix)
  rownames(gene_windows) <- genes
  rm(expression, adjusted, standardized); invisible(gc())

  ranked_regions <- copy(regions)
  setorder(ranked_regions, feature_set, -peak_z, q_region,
    -relative_peak_prominence, peak_pseudotime)
  strongest_region <- ranked_regions[, .SD[1L], by = feature_set]
  set_stats <- strongest_region[,
    .(analysis_group = group, feature_set, set_name,
      best_region_id = region_id, best_peak_z = peak_z, best_q = q_region,
      best_prominence = relative_peak_prominence,
      best_peak_pseudotime = peak_pseudotime)]
  region_counts <- regions[, .(n_significant_regions = .N), by = feature_set]
  set_stats[region_counts, on = "feature_set",
    n_significant_regions := i.n_significant_regions]
  setorder(set_stats, feature_set)
  set_integer <- setNames(seq_len(nrow(set_stats)), set_stats$feature_set)
  gene_integer <- setNames(seq_along(genes), genes)

  signatures <- vector("list", nrow(set_stats)); signature_rows <- list()
  for (i in seq_len(nrow(set_stats))) {
    fs <- set_stats$feature_set[[i]]
    fs_regions <- regions[feature_set == fs]
    significant_windows <- which(vapply(seq_len(nrow(windows)), function(w)
      any(windows$window_center[[w]] >= fs_regions$region_left &
        windows$window_center[[w]] <= fs_regions$region_right), logical(1)))
    fs_genes <- members[feature_set == fs, feature_id]
    contribution <- pmax(gene_windows[fs_genes, significant_windows,
      drop = FALSE], 0)
    nz <- which(contribution > 0, arr.ind = TRUE)
    stopifnot(nrow(nz) > 0L)
    raw_weight <- contribution[nz]
    weight <- raw_weight / sum(raw_weight)
    gene_i <- unname(gene_integer[fs_genes[nz[, 1L]]])
    window_i <- significant_windows[nz[, 2L]]
    full_i <- gene_i + (window_i - 1L) * length(genes)
    gene_marginal <- rowsum(weight, gene_i, reorder = TRUE)
    window_marginal <- rowsum(weight, window_i, reorder = TRUE)
    signatures[[i]] <- list(full_index = full_i, full_weight = weight,
      gene_index = as.integer(rownames(gene_marginal)),
      gene_weight = as.numeric(gene_marginal),
      window_index = as.integer(rownames(window_marginal)),
      window_weight = as.numeric(window_marginal),
      significant_windows = significant_windows)
    signature_rows[[i]] <- data.table(analysis_group = group,
      feature_set = fs, feature_id = fs_genes[nz[, 1L]],
      window_index = window_i,
      window_center = windows$window_center[window_i],
      positive_contribution = raw_weight,
      normalized_map_weight = weight)
  }
  all_signature_rows[[group]] <- rbindlist(signature_rows)

  # Candidate sets must share at least one gene. Full similarity additionally
  # requires that shared genes contribute in the same pseudotime windows.
  gene_to_sets <- split(unname(set_integer[members$feature_set]),
    unname(gene_integer[members$feature_id]))
  edge_chunks <- vector("list", nrow(set_stats))
  for (i in seq_len(nrow(set_stats) - 1L)) {
    candidate <- unique(unlist(gene_to_sets[
      as.character(signatures[[i]]$gene_index)], use.names = FALSE))
    candidate <- candidate[candidate > i]
    if (!length(candidate)) next
    full_sim <- vapply(candidate, function(j) weighted_jaccard(
      signatures[[i]]$full_index, signatures[[i]]$full_weight,
      signatures[[j]]$full_index, signatures[[j]]$full_weight), numeric(1))
    keep <- full_sim >= min(thresholds)
    if (any(keep)) {
      j_keep <- candidate[keep]; sim_keep <- full_sim[keep]
      gene_sim <- vapply(j_keep, function(j) weighted_jaccard(
        signatures[[i]]$gene_index, signatures[[i]]$gene_weight,
        signatures[[j]]$gene_index, signatures[[j]]$gene_weight), numeric(1))
      time_sim <- vapply(j_keep, function(j) weighted_jaccard(
        signatures[[i]]$window_index, signatures[[i]]$window_weight,
        signatures[[j]]$window_index, signatures[[j]]$window_weight), numeric(1))
      shared_windows <- vapply(j_keep, function(j) length(intersect(
        signatures[[i]]$significant_windows,
        signatures[[j]]$significant_windows)), integer(1))
      edge_chunks[[i]] <- data.table(i = i, j = j_keep,
        full_similarity = sim_keep,
        shared_contribution_mass = 2 * sim_keep / (1 + sim_keep),
        gene_marginal_similarity = gene_sim,
        temporal_marginal_similarity = time_sim,
        shared_significant_windows = shared_windows,
        significant_windows_a = length(signatures[[i]]$significant_windows),
        significant_windows_b = vapply(j_keep, function(j)
          length(signatures[[j]]$significant_windows), integer(1)))
    }
    if (i %% 500L == 0L) log_event(group, "multi-region similarity row", i,
      "of", nrow(set_stats))
  }
  edges <- rbindlist(Filter(Negate(is.null), edge_chunks))
  if (!nrow(edges)) edges <- data.table(i = integer(), j = integer(),
    full_similarity = numeric(), shared_contribution_mass = numeric(),
    gene_marginal_similarity = numeric(), temporal_marginal_similarity = numeric(),
    shared_significant_windows = integer(), significant_windows_a = integer(),
    significant_windows_b = integer())
  edges[, `:=`(feature_set_a = set_stats$feature_set[i],
    feature_set_b = set_stats$feature_set[j],
    best_z_a = set_stats$best_peak_z[i], best_z_b = set_stats$best_peak_z[j],
    n_regions_a = set_stats$n_significant_regions[i],
    n_regions_b = set_stats$n_significant_regions[j])]
  write_tsv_gz(edges, file.path(stage_dir,
    paste0(group, "_multiregion_edges_ge0.40.tsv.gz")))

  assignments <- lapply(thresholds, function(threshold) {
    a <- greedy_assign(set_stats, edges, threshold)
    a[, similarity_threshold := threshold]
    a
  })
  for (k in seq_along(thresholds)) {
    threshold <- thresholds[[k]]; a <- assignments[[k]]
    reps <- set_stats[a$retained_representative]
    retained_region_count <- regions[feature_set %chin% reps$feature_set, .N]
    all_summary[[paste(group, threshold)]] <- data.table(
      analysis_group = group, similarity_threshold = threshold,
      input_gene_sets = nrow(set_stats), input_regions = nrow(regions),
      similarity_edges = sum(edges$full_similarity >= threshold),
      retained_gene_sets = nrow(reps), retained_regions = retained_region_count,
      removed_gene_sets = nrow(set_stats) - nrow(reps),
      removed_regions = nrow(regions) - retained_region_count)
  }

  primary <- assignments[[which(thresholds == primary_threshold)]]
  mapping <- cbind(copy(set_stats), primary[, .(set_row, redundancy_group,
    representative_row, similarity_to_representative,
    retained_representative)])
  mapping[, `:=`(representative_feature_set =
      set_stats$feature_set[representative_row],
    representative_best_peak_z = set_stats$best_peak_z[representative_row],
    representative_best_peak_pseudotime =
      set_stats$best_peak_pseudotime[representative_row],
    edge_i = pmin(set_row, representative_row),
    edge_j = pmax(set_row, representative_row))]
  edge_metrics <- edges[, .(edge_i = i, edge_j = j,
    shared_contribution_mass, gene_marginal_similarity,
    temporal_marginal_similarity, shared_significant_windows,
    member_significant_windows = significant_windows_a,
    representative_significant_windows = significant_windows_b)]
  # Edge orientation can be opposite to member/representative orientation, but
  # the symmetric similarities and shared-window count remain valid. Window
  # totals are recomputed directly below for unambiguous reporting.
  mapping[edge_metrics, on = c("edge_i", "edge_j"), `:=`(
    shared_contribution_mass = i.shared_contribution_mass,
    gene_marginal_similarity = i.gene_marginal_similarity,
    temporal_marginal_similarity = i.temporal_marginal_similarity,
    shared_significant_windows = i.shared_significant_windows)]
  mapping[retained_representative == TRUE, `:=`(
    shared_contribution_mass = 1, gene_marginal_similarity = 1,
    temporal_marginal_similarity = 1,
    shared_significant_windows = lengths(lapply(set_row, function(x)
      signatures[[x]]$significant_windows)))]
  mapping[, `:=`(
    member_significant_windows = lengths(lapply(set_row, function(x)
      signatures[[x]]$significant_windows)),
    representative_significant_windows = lengths(lapply(representative_row,
      function(x) signatures[[x]]$significant_windows)))]
  mapping[, c("set_row", "representative_row", "edge_i", "edge_j") := NULL]
  setorder(mapping, redundancy_group, -retained_representative, -best_peak_z)
  retained_sets <- mapping[retained_representative == TRUE]
  retained_regions <- regions[feature_set %chin% retained_sets$feature_set]
  all_primary_mapping[[group]] <- mapping
  all_retained_sets[[group]] <- retained_sets
  all_retained_regions[[group]] <- retained_regions
  group_completion[[group]] <- list(input_gene_sets = nrow(set_stats),
    retained_gene_sets = nrow(retained_sets), input_regions = nrow(regions),
    retained_regions = nrow(retained_regions), edges_ge_0_40 = nrow(edges))
  rm(gene_windows); invisible(gc())
}

summary <- rbindlist(all_summary)
setorder(summary, analysis_group, similarity_threshold)
mapping <- rbindlist(all_primary_mapping, fill = TRUE)
retained_sets <- rbindlist(all_retained_sets, fill = TRUE)
retained_regions <- rbindlist(all_retained_regions, fill = TRUE)
signature_rows <- rbindlist(all_signature_rows)
write_tsv(summary, file.path(stage_dir,
  "multiregion_threshold_sensitivity.tsv"))
write_tsv_gz(mapping, file.path(stage_dir, paste0("primary_", primary_label,
  "_gene_set_to_representative.tsv.gz")))
write_tsv(retained_sets, file.path(stage_dir, paste0("primary_", primary_label,
  "_retained_gene_sets.tsv")))
write_tsv_gz(retained_regions, file.path(stage_dir, paste0("primary_",
  primary_label, "_retained_regions.tsv.gz")))
write_tsv_gz(signature_rows, file.path(stage_dir,
  "multiregion_gene_by_window_signatures.tsv.gz"))

primary_summary <- summary[similarity_threshold == primary_threshold]
table_lines <- c(
  "| Group | Input sets | Retained sets | Input regions | Retained regions |",
  "|---|---:|---:|---:|---:|",
  vapply(seq_len(nrow(primary_summary)), function(i) paste0("| ",
    primary_summary$analysis_group[[i]], " | ",
    primary_summary$input_gene_sets[[i]], " | ",
    primary_summary$retained_gene_sets[[i]], " | ",
    primary_summary$input_regions[[i]], " | ",
    primary_summary$retained_regions[[i]], " |"), character(1)))
writeLines(c("# Multi-region pseudotime-aware contribution similarity", "",
  paste0("Primary weighted-Jaccard threshold: ", primary_threshold), "",
  "Each gene set is represented by one normalized gene-by-pseudotime map. Positive standardized contributions are retained only in the union of that set's significant regions; all other gene-window entries are zero.", "",
  "Similarity therefore integrates every significant region while requiring shared genes to contribute at shared pseudotime windows. Missing weak regions can be tolerated, whereas missing regions carrying substantial contribution mass reduce similarity.", "",
  "Greedy non-transitive groups retain the gene set with the highest maximum region z-score. All significant regions belonging to retained representatives are preserved.", "",
  table_lines), file.path(stage_dir, "REPORT.md"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after),
  all(mapping$representative_best_peak_z >= mapping$best_peak_z),
  all(mapping$similarity_to_representative >= primary_threshold),
  all(mapping$shared_contribution_mass >=
    2 * primary_threshold / (1 + primary_threshold) - 1e-12))
completion <- list(status = "complete", source_analysis = source_root,
  primary_similarity_threshold = primary_threshold,
  similarity = paste0("weighted Jaccard of normalized positive gene-by-window ",
    "contributions within the union of each set's significant regions"),
  grouping = "greedy direct-to-highest-z gene-set representative; non-transitive",
  significant_window_definition = paste0("window center falls within any final ",
    "significant region for the gene set"),
  group_results = group_completion,
  pseudotime_source = "exact Bdev metadata pt",
  bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = fingerprint_after, bundle_unchanged = TRUE)
write_json(completion, file.path(stage_dir, "completion.json"))
wb_close_bundle(bundle)
stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed multi-region similarity at", output_dir)
