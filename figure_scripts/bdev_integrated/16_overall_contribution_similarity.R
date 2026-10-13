#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 16_overall_contribution_similarity.R",
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
  use <- edges[similarity >= threshold]
  adjacency <- vector("list", n)
  if (nrow(use)) for (k in seq_len(nrow(use))) {
    a <- use$i[[k]]; b <- use$j[[k]]; s <- use$similarity[[k]]
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
stopifnot(identical(as.character(metadata$cell_id), cells))
pt_all <- as.numeric(metadata$pt); pt_defined <- is.finite(pt_all)
pt <- pt_all[pt_defined]; pt_cells <- cells[pt_defined]
log_depth <- log1p(as.numeric(metadata$nCount_RNA)[pt_defined])
depth_z <- (log_depth - mean(log_depth)) / sd(log_depth)

pt_limits <- range(pt); pt_span <- diff(pt_limits)
centers <- pt_limits[[1L]] + seq(0, 1, by = 0.01) * pt_span
window_width <- 0.05 * pt_span
windows <- data.table(window_center = centers,
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
stopifnot(nrow(windows) == 97L)

all_summary <- list(); all_primary_mapping <- list()
all_retained_sets <- list(); all_retained_regions <- list()
all_signature_rows <- list(); group_completion <- list()
for (group in groups) {
  log_event("Building overall driver signatures for", group)
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
  # Each gene's overall weight is its total positive activity across all usable
  # windows. Thus this signature ignores where peaks occur while retaining the
  # relative strength of each gene's positive trajectory contribution.
  gene_overall <- rowSums(pmax(gene_windows, 0))
  names(gene_overall) <- genes
  rm(expression, adjusted, standardized, gene_windows); invisible(gc())

  ranked_regions <- copy(regions)
  setorder(ranked_regions, feature_set, -peak_z, q_region,
    -relative_peak_prominence, peak_pseudotime)
  set_stats <- ranked_regions[, .SD[1L], by = feature_set][,
    .(analysis_group = group, feature_set, set_name, best_region_id = region_id,
      best_peak_z = peak_z, best_q = q_region,
      best_prominence = relative_peak_prominence,
      best_peak_pseudotime = peak_pseudotime)]
  setorder(set_stats, feature_set)
  set_integer <- setNames(seq_len(nrow(set_stats)), set_stats$feature_set)
  gene_integer <- setNames(seq_along(genes), genes)
  signatures <- vector("list", nrow(set_stats)); signature_rows <- list()
  for (i in seq_len(nrow(set_stats))) {
    fs <- set_stats$feature_set[[i]]
    fs_genes <- members[feature_set == fs, feature_id]
    weight <- gene_overall[fs_genes]
    keep <- is.finite(weight) & weight > 0
    fs_genes <- fs_genes[keep]; weight <- weight[keep]
    weight <- weight / sum(weight)
    signatures[[i]] <- list(index = unname(gene_integer[fs_genes]), weight = weight)
    signature_rows[[i]] <- data.table(analysis_group = group,
      feature_set = fs, feature_id = fs_genes,
      total_positive_window_contribution = gene_overall[fs_genes],
      normalized_overall_weight = weight)
  }
  all_signature_rows[[group]] <- rbindlist(signature_rows)

  # Candidate pairs need at least one shared gene; unlike the region analysis,
  # there is deliberately no pseudotime-distance restriction.
  gene_to_sets <- split(unname(set_integer[members$feature_set]),
    unname(gene_integer[members$feature_id]))
  edge_chunks <- vector("list", nrow(set_stats))
  for (i in seq_len(nrow(set_stats) - 1L)) {
    candidate <- unique(unlist(gene_to_sets[as.character(signatures[[i]]$index)],
      use.names = FALSE))
    candidate <- candidate[candidate > i]
    if (!length(candidate)) next
    sim <- vapply(candidate, function(j) weighted_jaccard(
      signatures[[i]]$index, signatures[[i]]$weight,
      signatures[[j]]$index, signatures[[j]]$weight), numeric(1))
    keep <- sim >= min(thresholds)
    if (any(keep)) edge_chunks[[i]] <- data.table(i = i,
      j = candidate[keep], similarity = sim[keep])
    if (i %% 500L == 0L) log_event(group, "overall similarity row", i,
      "of", nrow(set_stats))
  }
  edges <- rbindlist(Filter(Negate(is.null), edge_chunks))
  if (!nrow(edges)) edges <- data.table(i = integer(), j = integer(),
    similarity = numeric())
  edges[, `:=`(feature_set_a = set_stats$feature_set[i],
    feature_set_b = set_stats$feature_set[j],
    best_z_a = set_stats$best_peak_z[i], best_z_b = set_stats$best_peak_z[j])]
  write_tsv_gz(edges, file.path(stage_dir,
    paste0(group, "_overall_similarity_edges_ge0.40.tsv.gz")))

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
      similarity_edges = sum(edges$similarity >= threshold),
      retained_gene_sets = nrow(reps), retained_regions = retained_region_count,
      removed_gene_sets = nrow(set_stats) - nrow(reps),
      removed_regions = nrow(regions) - retained_region_count)
  }

  primary <- assignments[[which(thresholds == primary_threshold)]]
  mapping <- cbind(copy(set_stats), primary[, .(redundancy_group,
    representative_row, similarity_to_representative,
    retained_representative)])
  mapping[, `:=`(
    representative_feature_set = set_stats$feature_set[representative_row],
    representative_best_peak_z = set_stats$best_peak_z[representative_row],
    representative_best_peak_pseudotime =
      set_stats$best_peak_pseudotime[representative_row])]
  mapping[, representative_row := NULL]
  setorder(mapping, redundancy_group, -retained_representative, -best_peak_z)
  retained_sets <- mapping[retained_representative == TRUE]
  retained_regions <- regions[feature_set %chin% retained_sets$feature_set]
  all_primary_mapping[[group]] <- mapping
  all_retained_sets[[group]] <- retained_sets
  all_retained_regions[[group]] <- retained_regions
  group_completion[[group]] <- list(input_gene_sets = nrow(set_stats),
    retained_gene_sets = nrow(retained_sets), input_regions = nrow(regions),
    retained_regions = nrow(retained_regions), edges_ge_0_40 = nrow(edges))
}

summary <- rbindlist(all_summary)
setorder(summary, analysis_group, similarity_threshold)
mapping <- rbindlist(all_primary_mapping, fill = TRUE)
retained_sets <- rbindlist(all_retained_sets, fill = TRUE)
retained_regions <- rbindlist(all_retained_regions, fill = TRUE)
signature_rows <- rbindlist(all_signature_rows)
write_tsv(summary, file.path(stage_dir, "overall_threshold_sensitivity.tsv"))
write_tsv_gz(mapping, file.path(stage_dir, paste0("primary_", primary_label,
  "_gene_set_to_representative.tsv.gz")))
write_tsv(retained_sets, file.path(stage_dir, paste0("primary_", primary_label,
  "_retained_gene_sets.tsv")))
write_tsv_gz(retained_regions, file.path(stage_dir, paste0("primary_",
  primary_label, "_retained_regions.tsv.gz")))
write_tsv_gz(signature_rows, file.path(stage_dir,
  "overall_positive_contribution_signatures.tsv.gz"))

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
writeLines(c("# Overall gene-set contribution similarity", "",
  paste0("Primary weighted-Jaccard threshold: ", primary_threshold), "",
  "For each gene, positive standardized contribution was summed across all 97 usable pseudotime windows. These gene weights were normalized within each gene set and compared without any region-overlap or peak-distance restriction.", "",
  "Greedy groups retain the gene set with the highest z-score among directly similar, still-unassigned sets. When a gene set is removed, all of its regions are removed from the reduced overall list.", "",
  table_lines, "",
  "This is deliberately complementary to the region-level result: it detects common overall drivers even when their enriched regions occur at different pseudotimes."),
  file.path(stage_dir, "REPORT.md"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after),
  all(mapping$representative_best_peak_z >= mapping$best_peak_z),
  all(mapping$similarity_to_representative >= primary_threshold))
completion <- list(status = "complete", source_analysis = source_root,
  primary_similarity_threshold = primary_threshold,
  similarity = paste0("weighted Jaccard of gene weights from total positive ",
    "contribution across all 97 pseudotime windows"),
  temporal_restriction = "none", grouping =
    "greedy direct-to-highest-z gene-set representative; non-transitive",
  group_results = group_completion,
  pseudotime_source = "exact Bdev metadata pt",
  bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = fingerprint_after, bundle_unchanged = TRUE)
write_json(completion, file.path(stage_dir, "completion.json"))
wb_close_bundle(bundle)
stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed overall contribution similarity at", output_dir)
