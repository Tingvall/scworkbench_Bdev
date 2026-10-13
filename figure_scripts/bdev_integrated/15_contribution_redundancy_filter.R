#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 15_contribution_redundancy_filter.R",
  "<bundle_dir> <memberships.tsv.gz> <source_root> <output_dir> <primary_threshold>"))

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
membership_path <- normalizePath(args[[2L]], mustWork = TRUE)
source_root <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
primary_threshold <- as.numeric(args[[5L]])
thresholds <- c(0.40, 0.50, 0.60, 0.70, 0.80)
minimum_threshold <- min(thresholds)
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

groups <- c("cytokines", "bp", "mf", "cc")
stopifnot(all(file.exists(file.path(source_root, groups, "completion.json"))))
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

# Recreate the exact 20-cell window grid used by the source analysis.
pt_limits <- range(pt); pt_span <- diff(pt_limits)
window_width_fraction <- 0.05; step_fraction <- 0.01
window_width <- window_width_fraction * pt_span
centers <- pt_limits[[1L]] + seq(0, 1, by = step_fraction) * pt_span
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

weighted_jaccard <- function(a_idx, a_val, b_idx, b_val) {
  common <- intersect(a_idx, b_idx)
  if (!length(common)) return(0)
  ai <- match(common, a_idx); bi <- match(common, b_idx)
  overlap <- sum(pmin(a_val[ai], b_val[bi]))
  overlap / (2 - overlap)
}

greedy_assign <- function(regions, edges, threshold) {
  n <- nrow(regions)
  rank_order <- order(-regions$peak_z, regions$q_region,
    -regions$relative_peak_prominence, regions$feature_set,
    regions$peak_pseudotime)
  assigned <- integer(n); representative_similarity <- rep(NA_real_, n)
  group_number <- 0L; representatives <- integer(n)
  edge_use <- edges[similarity >= threshold]
  adjacency <- vector("list", n)
  if (nrow(edge_use)) for (k in seq_len(nrow(edge_use))) {
    a <- edge_use$i[[k]]; b <- edge_use$j[[k]]; s <- edge_use$similarity[[k]]
    adjacency[[a]] <- rbind(adjacency[[a]], c(b, s))
    adjacency[[b]] <- rbind(adjacency[[b]], c(a, s))
  }
  for (rep_i in rank_order) {
    if (assigned[[rep_i]] != 0L) next
    group_number <- group_number + 1L
    representatives[[group_number]] <- rep_i
    assigned[[rep_i]] <- group_number
    representative_similarity[[rep_i]] <- 1
    neighbors <- adjacency[[rep_i]]
    if (!is.null(neighbors)) {
      neighbors <- as.matrix(neighbors)
      open <- neighbors[, 1L][assigned[as.integer(neighbors[, 1L])] == 0L]
      if (length(open)) {
        assigned[as.integer(open)] <- group_number
        representative_similarity[as.integer(open)] <-
          neighbors[match(open, neighbors[, 1L]), 2L]
      }
    }
  }
  rep_by_group <- setNames(representatives[seq_len(group_number)],
    seq_len(group_number))
  representative_i <- unname(rep_by_group[as.character(assigned)])
  data.table(region_row = seq_len(n), redundancy_group = assigned,
    representative_row = representative_i,
    similarity_to_representative = representative_similarity,
    retained_representative = seq_len(n) == representative_i)
}

all_threshold_summary <- list(); all_primary_mapping <- list()
all_primary_representatives <- list(); all_signature_rows <- list()
group_completion <- list()
for (group in groups) {
  log_event("Building contribution signatures for", group)
  group_dir <- file.path(source_root, group)
  regions <- fread(file.path(group_dir, "significant_positive_regions.tsv"))
  transforms <- read_tsv_gz(file.path(group_dir,
    "gene_standardization_diagnostics.tsv.gz"))[standardizable == TRUE]
  members <- members_all[analysis_group == group &
    feature_id %chin% transforms$feature_id]
  genes <- transforms$feature_id
  setorder(transforms, feature_id); genes <- transforms$feature_id
  expression <- as.matrix(wb_matrix_read_features(bundle,
    "RNA_log1p_cp10k", genes, pt_cells, missing = "NA"))
  stopifnot(identical(rownames(expression), genes), all(is.finite(expression)))
  adjusted <- expression - tcrossprod(transforms$conditional_depth_beta, depth_z)
  standardized <- (adjusted - transforms$adjusted_center) /
    transforms$adjusted_scale
  gene_windows <- standardized %*% t(weight_matrix)
  rownames(gene_windows) <- genes
  rm(expression, adjusted, standardized); invisible(gc())

  gene_integer <- setNames(seq_along(genes), genes)
  signatures <- vector("list", nrow(regions))
  signature_rows <- vector("list", nrow(regions))
  for (i in seq_len(nrow(regions))) {
    fs <- regions$feature_set[[i]]
    fs_genes <- members[feature_set == fs, feature_id]
    contribution <- gene_windows[fs_genes, regions$peak_window[[i]], drop = TRUE]
    positive <- is.finite(contribution) & contribution > 0
    fs_genes <- fs_genes[positive]; contribution <- contribution[positive]
    stopifnot(length(contribution) >= 1L, sum(contribution) > 0)
    weight <- contribution / sum(contribution)
    signatures[[i]] <- list(index = unname(gene_integer[fs_genes]), weight = weight)
    signature_rows[[i]] <- data.table(analysis_group = group,
      region_id = regions$region_id[[i]], feature_set = fs,
      peak_pseudotime = regions$peak_pseudotime[[i]], feature_id = fs_genes,
      positive_contribution = contribution, normalized_positive_weight = weight)
  }
  signature_table <- rbindlist(signature_rows)
  all_signature_rows[[group]] <- signature_table

  # Calculate each temporally eligible pair once. Save edges down to the lowest
  # all requested threshold sensitivities use the same similarity calculations.
  log_event("Calculating temporally local similarities for", nrow(regions),
    group, "regions")
  edge_chunks <- vector("list", nrow(regions))
  for (i in seq_len(nrow(regions) - 1L)) {
    j <- seq.int(i + 1L, nrow(regions))
    temporal <- (regions$region_left[j] <= regions$region_right[[i]] &
      regions$region_right[j] >= regions$region_left[[i]]) |
      abs(regions$peak_pseudotime[j] - regions$peak_pseudotime[[i]]) <=
        window_width
    j <- j[temporal]
    if (!length(j)) next
    sim <- vapply(j, function(k) weighted_jaccard(
      signatures[[i]]$index, signatures[[i]]$weight,
      signatures[[k]]$index, signatures[[k]]$weight), numeric(1))
    keep <- sim >= minimum_threshold
    if (any(keep)) edge_chunks[[i]] <- data.table(i = i, j = j[keep],
      similarity = sim[keep])
    if (i %% 500L == 0L) log_event(group, "similarity row", i, "of",
      nrow(regions))
  }
  edges <- rbindlist(Filter(Negate(is.null), edge_chunks))
  if (!nrow(edges)) edges <- data.table(i = integer(), j = integer(),
    similarity = numeric())
  edges[, `:=`(region_id_a = regions$region_id[i],
    region_id_b = regions$region_id[j],
    feature_set_a = regions$feature_set[i],
    feature_set_b = regions$feature_set[j],
    peak_pt_a = regions$peak_pseudotime[i],
    peak_pt_b = regions$peak_pseudotime[j])]
  write_tsv_gz(edges, file.path(stage_dir,
    paste0(group, "_similarity_edges_ge0.40.tsv.gz")))

  assignments <- lapply(thresholds, function(threshold) {
    a <- greedy_assign(regions, edges, threshold)
    a[, similarity_threshold := threshold]
    a
  })
  for (k in seq_along(thresholds)) {
    threshold <- thresholds[[k]]; a <- assignments[[k]]
    reps <- regions[a$retained_representative]
    all_threshold_summary[[paste(group, threshold)]] <- data.table(
      analysis_group = group, similarity_threshold = threshold,
      input_gene_sets = uniqueN(regions$feature_set),
      input_regions = nrow(regions), similarity_edges =
        sum(edges$similarity >= threshold),
      retained_gene_sets = uniqueN(reps$feature_set),
      retained_regions = nrow(reps), removed_regions = nrow(regions) - nrow(reps))
  }

  primary <- assignments[[which(thresholds == primary_threshold)]]
  mapping <- cbind(copy(regions), primary[, .(redundancy_group,
    representative_row, similarity_to_representative,
    retained_representative)])
  mapping[, `:=`(
    representative_region_id = regions$region_id[representative_row],
    representative_feature_set = regions$feature_set[representative_row],
    representative_peak_z = regions$peak_z[representative_row],
    representative_peak_pseudotime = regions$peak_pseudotime[representative_row])]
  mapping[, representative_row := NULL]
  setorder(mapping, redundancy_group, -retained_representative, -peak_z)
  representatives <- mapping[retained_representative == TRUE]
  all_primary_mapping[[group]] <- mapping
  all_primary_representatives[[group]] <- representatives
  group_completion[[group]] <- list(input_regions = nrow(regions),
    input_gene_sets = uniqueN(regions$feature_set),
    retained_regions = nrow(representatives),
    retained_gene_sets = uniqueN(representatives$feature_set),
    similarity_edges_ge_0_40 = nrow(edges))
  rm(gene_windows); invisible(gc())
}

threshold_summary <- rbindlist(all_threshold_summary)
setorder(threshold_summary, analysis_group, similarity_threshold)
primary_mapping <- rbindlist(all_primary_mapping, fill = TRUE)
primary_representatives <- rbindlist(all_primary_representatives, fill = TRUE)
signature_table <- rbindlist(all_signature_rows)
write_tsv(threshold_summary, file.path(stage_dir,
  "threshold_sensitivity_summary.tsv"))
write_tsv_gz(primary_mapping, file.path(stage_dir,
  paste0("primary_", primary_label, "_region_to_representative.tsv.gz")))
write_tsv_gz(primary_representatives, file.path(stage_dir,
  paste0("primary_", primary_label,
    "_retained_representative_regions.tsv.gz")))
write_tsv_gz(signature_table, file.path(stage_dir,
  "positive_gene_contribution_signatures.tsv.gz"))

primary_summary <- threshold_summary[similarity_threshold == primary_threshold]
table_lines <- c(
  "| Group | Input sets | Input regions | Retained sets | Retained regions | Removed regions |",
  "|---|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(primary_summary)), function(i) paste0("| ",
    primary_summary$analysis_group[[i]], " | ",
    primary_summary$input_gene_sets[[i]], " | ",
    primary_summary$input_regions[[i]], " | ",
    primary_summary$retained_gene_sets[[i]], " | ",
    primary_summary$retained_regions[[i]], " | ",
    primary_summary$removed_regions[[i]], " |"), character(1)))
writeLines(c("# Contribution-weighted gene-set redundancy filter", "",
  paste0("Primary weighted-Jaccard threshold: ", primary_threshold), "",
  "Positive gene contributions at each representative peak were normalized to sum to one. Regions were compared only within the same analysis group and when their intervals overlapped or their peaks were no more than one 5% pseudotime window apart.", "",
  "Regions were processed from highest to lowest peak z-score. Each retained region absorbed directly similar, still-unassigned regions; ties were resolved by q-value, prominence, gene-set identifier, and peak pseudotime. This non-transitive greedy rule prevents similarity chaining.", "",
  table_lines, "",
  "Threshold sensitivity for 0.40, 0.50, 0.60, 0.70, and 0.80 is recorded in `threshold_sensitivity_summary.tsv`. The complete representative mapping and contribution signatures are retained for audit."),
  file.path(stage_dir, "REPORT.md"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after),
  nrow(primary_mapping) == sum(primary_summary$input_regions),
  nrow(primary_representatives) == sum(primary_summary$retained_regions),
  all(primary_mapping$representative_peak_z >= primary_mapping$peak_z),
  all(primary_mapping$similarity_to_representative >= primary_threshold))
completion <- list(status = "complete", source_analysis = source_root,
  primary_similarity_threshold = primary_threshold,
  sensitivity_thresholds = as.list(thresholds),
  similarity = "weighted Jaccard of normalized positive gene contributions",
  temporal_eligibility = paste0("region overlap or peak distance <= one ",
    "5% pseudotime window"),
  grouping = "greedy direct-to-highest-z representative; non-transitive",
  group_results = group_completion,
  pseudotime_source = "exact Bdev metadata pt",
  bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = fingerprint_after,
  bundle_unchanged = TRUE)
write_json(completion, file.path(stage_dir, "completion.json"))
wb_close_bundle(bundle)
stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed contribution redundancy filter at", output_dir)
