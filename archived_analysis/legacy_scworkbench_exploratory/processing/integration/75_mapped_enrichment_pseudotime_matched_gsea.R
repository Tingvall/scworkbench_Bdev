#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "Usage: 75_mapped_enrichment_pseudotime_matched_gsea.R",
  "<query_bundle> <session_workspace> <projection.tsv.gz>",
  "<density.tsv> <regions.tsv> <all_gene_de.tsv>",
  "<cytokine_membership.tsv> <output_dir> <seed>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
session_workspace <- normalizePath(args[[2L]], mustWork = TRUE)
projection_path <- normalizePath(args[[3L]], mustWork = TRUE)
density_path <- normalizePath(args[[4L]], mustWork = TRUE)
region_path <- normalizePath(args[[5L]], mustWork = TRUE)
gene_path <- normalizePath(args[[6L]], mustWork = TRUE)
cytokine_path <- normalizePath(args[[7L]], mustWork = TRUE)
output_dir <- normalizePath(args[[8L]], mustWork = FALSE)
seed <- as.integer(args[[9L]])
stopifnot(is.finite(seed))
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
  library(jsonlite)
})

write_json <- function(x, path) jsonlite::write_json(x, path,
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
md5 <- function(path) unname(tools::md5sum(path))
pretty_set <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  tools::toTitleCase(tolower(gsub("_", " ", x)))
}

regions <- fread(region_path)
stopifnot(nrow(regions) == 2L,
  identical(sort(regions$region_id), sort(c("Flt3l_shared_peak1",
    "Il7_shared_peak2"))))
regions[, comparison_id := fifelse(target == "Flt3l",
  "Flt3L_KO_vs_WT", "IL7_KO_vs_WT")]
regions[, comparison := fifelse(target == "Flt3l",
  "FLT3LKO vs FLT3L WT", "IL7KO vs IL-7 WT")]

density <- fread(density_path)
required_density <- c("target", "genotype", "bin", "bin_left", "bin_right",
  "bin_midpoint", "n_cells", "bdev_n_cells")
stopifnot(all(required_density %in% names(density)))
bin_defs <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
  reg <- regions[i]
  x <- unique(density[target == reg$target & genotype == "WT" &
      bin >= reg$first_bin & bin <= reg$last_bin,
    .(target, bin, bin_left, bin_right, bin_midpoint, bdev_n_cells)])
  stopifnot(nrow(x) == reg$n_bins, all(x$bin_left >= reg$pt_left - 1e-10),
    all(x$bin_right <= reg$pt_right + 1e-10))
  x[, `:=`(region_id = reg$region_id,
    comparison_id = reg$comparison_id,
    bdev_weight = bdev_n_cells / sum(bdev_n_cells))]
  x
}))

projection <- fread(cmd = paste("gzip -cd", shQuote(projection_path)))[
  support_status %chin% c("supported", "borderline") & is.finite(mapped_pt),
  .(source_cell_id, sample_id, genotype, knockout_target, support_status,
    mapped_pt)]
stopifnot(!anyDuplicated(projection$source_cell_id),
  setequal(unique(projection$genotype), c("WT", "KO")))
assignments <- rbindlist(lapply(seq_len(nrow(bin_defs)), function(i) {
  b <- bin_defs[i]
  projection[knockout_target == b$target & mapped_pt >= b$bin_left &
      mapped_pt < b$bin_right,
    .(source_cell_id, sample_id, genotype, knockout_target, support_status,
      mapped_pt, region_id = b$region_id, comparison_id = b$comparison_id,
      bin = b$bin, bin_left = b$bin_left, bin_right = b$bin_right,
      bin_midpoint = b$bin_midpoint, bdev_n_cells = b$bdev_n_cells,
      bdev_weight = b$bdev_weight)]
}))
stopifnot(!anyDuplicated(assignments$source_cell_id),
  nrow(assignments[region_id == "Flt3l_shared_peak1" & genotype == "WT"]) == 4394L,
  nrow(assignments[region_id == "Flt3l_shared_peak1" & genotype == "KO"]) == 3342L,
  nrow(assignments[region_id == "Il7_shared_peak2" & genotype == "WT"]) == 1912L,
  nrow(assignments[region_id == "Il7_shared_peak2" & genotype == "KO"]) == 1355L)
bin_counts <- assignments[, .(n_cells = .N),
  by = .(region_id, comparison_id, bin, bin_left, bin_right, bin_midpoint,
    bdev_n_cells, bdev_weight, genotype)]
bin_counts_wide <- dcast(bin_counts,
  region_id + comparison_id + bin + bin_left + bin_right + bin_midpoint +
    bdev_n_cells + bdev_weight ~ genotype, value.var = "n_cells", fill = 0L)
setnames(bin_counts_wide, c("WT", "KO"), c("n_wt", "n_ko"))
bin_counts_wide[, reliable_20 := n_wt >= 20L & n_ko >= 20L]
stopifnot(all(bin_counts_wide$n_wt > 0L), all(bin_counts_wide$n_ko > 0L),
  sum(!bin_counts_wide$reliable_20) == 2L)
fwrite(bin_counts_wide, file.path(stage_dir, "pseudotime_bin_cell_counts.tsv"),
  sep = "\t")

all_genes <- fread(gene_path, select = "feature_id")$feature_id
stopifnot(length(all_genes) == 33696L, !anyDuplicated(all_genes))
bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
staged_bundle_dir <- file.path(session_workspace, "stage")
staged <- wb_open_bundle(staged_bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(staged), silent = TRUE), add = TRUE)
matrix_id <- "RNA_log1p_cp10k_wt_enriched_peak_de_v1"
stopifnot(matrix_id %in% wb_matrices(staged)$matrix_id)
cell_order <- assignments$source_cell_id

message("Calculating bin-matched expression means for ", length(all_genes),
  " genes and ", length(cell_order), " cells")
blocks <- split(all_genes, ceiling(seq_along(all_genes) / 128L))
bin_means <- rbindlist(lapply(seq_along(blocks), function(block_index) {
  if (block_index %% 25L == 0L)
    message("Expression block ", block_index, "/", length(blocks))
  genes <- blocks[[block_index]]
  m <- wb_matrix_read_features(staged, matrix_id, genes, cell_order,
    missing = "NA")
  stopifnot(!anyNA(m), identical(rownames(m), genes),
    identical(colnames(m), cell_order))
  linear <- expm1(as.matrix(m))
  rbindlist(lapply(seq_len(nrow(bin_counts_wide)), function(i) {
    b <- bin_counts_wide[i]
    wt <- which(assignments$region_id == b$region_id &
      assignments$bin == b$bin & assignments$genotype == "WT")
    ko <- which(assignments$region_id == b$region_id &
      assignments$bin == b$bin & assignments$genotype == "KO")
    stopifnot(length(wt) == b$n_wt, length(ko) == b$n_ko)
    data.table(region_id = b$region_id, comparison_id = b$comparison_id,
      bin = b$bin, bin_midpoint = b$bin_midpoint,
      bdev_n_cells = b$bdev_n_cells, bdev_weight = b$bdev_weight,
      n_wt = b$n_wt, n_ko = b$n_ko, reliable_20 = b$reliable_20,
      feature_id = genes,
      mean_cp10k_wt = rowMeans(linear[, wt, drop = FALSE]),
      mean_cp10k_ko = rowMeans(linear[, ko, drop = FALSE]))
  }))
}), use.names = TRUE)
bin_means[, bin_ko_vs_wt_log2fc :=
  log2((mean_cp10k_ko + 1) / (mean_cp10k_wt + 1))]
fwrite(bin_means, file.path(stage_dir,
  "pseudotime_bin_gene_expression_means.tsv.gz"), sep = "\t",
  compress = "gzip", na = "NA")

rankings <- bin_means[, {
  primary_weight <- bdev_weight / sum(bdev_weight)
  sensitivity_weight <- fifelse(reliable_20, bdev_weight, 0)
  sensitivity_weight <- sensitivity_weight / sum(sensitivity_weight)
  .(weighted_mean_cp10k_wt = sum(primary_weight * mean_cp10k_wt),
    weighted_mean_cp10k_ko = sum(primary_weight * mean_cp10k_ko),
    bdev_weighted_ko_vs_wt_log2fc = sum(primary_weight * bin_ko_vs_wt_log2fc),
    sensitivity_ko_vs_wt_log2fc = sum(sensitivity_weight * bin_ko_vs_wt_log2fc),
    bins_used = .N, sensitivity_bins_used = sum(reliable_20))
}, by = .(region_id, comparison_id, feature_id)]

break_ties <- function(d, primary, secondary) {
  setorderv(d, c(primary, secondary, "feature_id"), c(1L, 1L, 1L))
  values <- d[[primary]]
  unique_scores <- sort(unique(values))
  positive_gaps <- diff(unique_scores)
  minimum_gap <- if (length(positive_gaps)) min(positive_gaps) else 1
  d[, rank_score := get(primary)]
  d[, tie_size := .N, by = primary]
  for (value in d[tie_size > 1L, unique(get(primary))]) {
    idx <- which(d[[primary]] == value)
    local_gaps <- c(value - max(unique_scores[unique_scores < value],
      na.rm = TRUE), min(unique_scores[unique_scores > value],
      na.rm = TRUE) - value)
    local_gaps <- local_gaps[is.finite(local_gaps) & local_gaps > 0]
    margin <- if (length(local_gaps)) min(local_gaps) / 4 else
      minimum_gap / 4
    d$rank_score[idx] <- value + seq(-margin, margin,
      length.out = length(idx))
  }
  stopifnot(!anyDuplicated(d$rank_score))
  d
}

ranking_list <- lapply(split(rankings, rankings$comparison_id), function(d) {
  d <- copy(d[weighted_mean_cp10k_wt + weighted_mean_cp10k_ko > 0])
  d <- break_ties(d, "bdev_weighted_ko_vs_wt_log2fc",
    "sensitivity_ko_vs_wt_log2fc")
  setorder(d, -rank_score)
  d[, rank := seq_len(.N)]
  d
})
rankings <- rbindlist(ranking_list)
fwrite(rankings, file.path(stage_dir,
  "pseudotime_matched_full_detected_gene_rankings.tsv.gz"), sep = "\t",
  compress = "gzip", na = "NA")

message("Loading mouse-native MSigDB Hallmark and GO memberships")
hallmark <- as.data.table(msigdbr(db_species = "MM", species = "Mus musculus",
  collection = "MH"))
go <- as.data.table(msigdbr(db_species = "MM", species = "Mus musculus",
  collection = "M5"))[gs_subcollection %chin% c("GO:BP", "GO:MF", "GO:CC")]
cytokine <- fread(cytokine_path)
stopifnot(uniqueN(hallmark$gs_name) == 50L,
  setequal(unique(go$gs_subcollection), c("GO:BP", "GO:MF", "GO:CC")),
  uniqueN(cytokine$signature_id) == 38L)
memberships <- rbindlist(list(
  unique(hallmark[, .(collection = "Hallmark", set_id = gs_name,
    set_label = gs_name, gene = gene_symbol)]),
  unique(go[, .(collection = sub("GO:", "GO_", gs_subcollection),
    set_id = gs_name, set_label = gs_name, gene = gene_symbol)]),
  unique(cytokine[, .(collection = "Cytokine", set_id = signature_id,
    set_label = paste(cytokine, direction), gene = feature_id)])))
setorder(memberships, collection, set_id, gene)
fwrite(memberships, file.path(stage_dir, "gene_set_memberships.tsv.gz"),
  sep = "\t", compress = "gzip")

run_one <- function(ranking, collection_name, seed_offset,
    score_column = "rank_score") {
  member <- memberships[collection == collection_name]
  pathways <- split(member$gene, member$set_id)
  stats <- ranking[[score_column]]
  names(stats) <- ranking$feature_id
  stats <- sort(stats, decreasing = TRUE)
  set.seed(seed + seed_offset)
  ans <- as.data.table(fgseaMultilevel(pathways = pathways, stats = stats,
    minSize = 15L, maxSize = 500L, eps = 0, scoreType = "std", nproc = 1L))
  setnames(ans, "pathway", "set_id")
  ans[, `:=`(comparison_id = ranking$comparison_id[[1L]],
    region_id = ranking$region_id[[1L]], collection = collection_name,
    enrichment_direction = fifelse(NES > 0, "KO higher", "WT higher"),
    leading_edge_genes = vapply(leadingEdge, paste, collapse = ",",
      FUN.VALUE = character(1)))]
  ans[, leadingEdge := NULL]
  setcolorder(ans, c("comparison_id", "region_id", "collection", "set_id",
    "enrichment_direction", "size", "ES", "NES", "pval", "padj",
    "log2err", "leading_edge_genes"))
  ans[, absolute_NES_for_order := abs(NES)]
  setorder(ans, padj, -absolute_NES_for_order, set_id)
  ans[, absolute_NES_for_order := NULL]
  ans
}

collections <- c("Hallmark", "GO_BP", "GO_MF", "GO_CC", "Cytokine")
gsea_results <- list()
offset <- 0L
for (comparison_id in names(ranking_list)) {
  for (collection_name in collections) {
    offset <- offset + 1L
    message("GSEA: ", comparison_id, " / ", collection_name)
    key <- paste(comparison_id, collection_name, sep = "__")
    gsea_results[[key]] <- run_one(ranking_list[[comparison_id]],
      collection_name, offset)
  }
}
results <- rbindlist(gsea_results, fill = TRUE)
fwrite(results, file.path(stage_dir,
  "all_collections_pseudotime_matched_gsea.tsv.gz"), sep = "\t",
  compress = "gzip", na = "NA")

# Hallmark sensitivity analysis excludes bins with fewer than 20 cells in
# either genotype, while retaining Bdev-density weights over remaining bins.
sensitivity_results <- rbindlist(lapply(seq_along(ranking_list), function(i) {
  d <- copy(ranking_list[[i]])
  d[, sensitivity_rank_score := sensitivity_ko_vs_wt_log2fc]
  d <- break_ties(d, "sensitivity_rank_score",
    "bdev_weighted_ko_vs_wt_log2fc")
  run_one(d, "Hallmark", 100L + i, "rank_score")
}))
setnames(sensitivity_results, c("enrichment_direction", "ES", "NES", "pval",
  "padj", "log2err", "leading_edge_genes"),
  c("sensitivity_enrichment_direction", "sensitivity_ES",
    "sensitivity_NES", "sensitivity_pval", "sensitivity_padj",
    "sensitivity_log2err", "sensitivity_leading_edge_genes"))
hallmark_primary <- results[collection == "Hallmark"]
hallmark_sensitivity <- merge(hallmark_primary, sensitivity_results,
  by = c("comparison_id", "region_id", "collection", "set_id", "size"),
  all = TRUE)
fwrite(hallmark_sensitivity, file.path(stage_dir,
  "hallmark_sparse_bin_sensitivity.tsv"), sep = "\t", na = "NA")

summary <- results[, .(sets_tested = .N, fdr_below_0_05 = sum(padj < 0.05),
  ko_higher_fdr_below_0_05 = sum(padj < 0.05 & NES > 0),
  wt_higher_fdr_below_0_05 = sum(padj < 0.05 & NES < 0)),
  by = .(comparison_id, region_id, collection)]
fwrite(summary, file.path(stage_dir, "gsea_summary.tsv"), sep = "\t")

make_heatmap <- function(collection_name, max_each = Inf) {
  d <- copy(results[collection == collection_name])
  if (is.finite(max_each)) {
    selected <- unique(d[order(padj, -abs(NES)), head(set_id, max_each),
      by = comparison_id]$V1)
    d <- d[set_id %chin% selected]
  }
  d[, display := if (collection_name == "Cytokine")
    sub("^parse_cytokine_rest_only_refinement_all_top100_v1\\.", "", set_id) else
    pretty_set(set_id)]
  d[, max_abs := max(abs(NES)), by = set_id]
  d[, mean_nes := mean(NES), by = set_id]
  setorder(d, max_abs, mean_nes, set_id)
  d[, display := factor(display, levels = unique(display))]
  d[, comparison_label := factor(comparison_id,
    levels = c("Flt3L_KO_vs_WT", "IL7_KO_vs_WT"),
    labels = c("FLT3LKO / WT", "IL7KO / WT"))]
  d[, mark := fifelse(padj < 0.001, "***", fifelse(padj < 0.01, "**",
    fifelse(padj < 0.05, "*", "")))]
  lim <- max(abs(d$NES))
  p <- ggplot(d, aes(comparison_label, display, fill = NES)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(aes(label = mark), size = 2.4) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, limits = c(-lim, lim), name = "NES") +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6.6),
      axis.title = element_blank(), plot.title = element_text(face = "bold")) +
    labs(title = paste(collection_name, "pseudotime-matched GSEA"),
      subtitle = "Positive/blue: KO higher; negative/red: WT higher; stars mark within-comparison collection FDR")
  stub <- tolower(collection_name)
  height <- max(6, min(18, 2.5 + 0.19 * uniqueN(d$set_id)))
  ggsave(file.path(stage_dir, "figures", paste0(stub, "_gsea_heatmap.png")),
    p, width = 8.5, height = height, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(stub, "_gsea_heatmap.pdf")),
    p, width = 8.5, height = height, device = cairo_pdf, bg = "white")
}
make_heatmap("Hallmark")
make_heatmap("GO_BP", 20L)
make_heatmap("GO_MF", 20L)
make_heatmap("GO_CC", 20L)
make_heatmap("Cytokine")

rank_cor <- rankings[, .(spearman_primary_vs_sensitivity = cor(
    bdev_weighted_ko_vs_wt_log2fc, sensitivity_ko_vs_wt_log2fc,
    method = "spearman")), by = comparison_id]
hallmark_cor <- hallmark_sensitivity[, .(
  hallmark_NES_primary_vs_sensitivity = cor(NES, sensitivity_NES,
    method = "spearman")), by = comparison_id]
audit <- merge(rank_cor, hallmark_cor, by = "comparison_id")
fwrite(audit, file.path(stage_dir, "ranking_sensitivity_audit.tsv"), sep = "\t")

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after))
diagnostics <- list(status = "complete_unpromoted_report_only_outside_scworkbench",
  rationale = paste("scWorkbench pathway_enrich is hypergeometric",
    "over-representation analysis and does not implement preranked GSEA"),
  user_approved_outside_analysis = TRUE,
  region_policy = "mapped shared enrichment",
  ranking = paste("Bdev-density-weighted mean across within-bin",
    "log2((mean CP10K KO + 1)/(mean CP10K WT + 1))"),
  comparisons = as.list(regions[, .(comparison_id, region_id, pt_left,
    pt_right)]),
  bins = list(Flt3L = 8L, IL7 = 11L,
    sparse_sensitivity_exclusion = "n_WT <20 or n_KO <20"),
  gene_sets = list(database = "Mouse MSigDB 2026.1.Mm plus accepted cytokine signatures",
    collections = collections),
  gsea = list(method = "fgseaMultilevel", min_size = 15L,
    max_size = 500L, eps = 0, seed = seed,
    p_adjustment = "BH separately within comparison and collection"),
  summary = summary, sensitivity = audit,
  biological_replicates_per_genotype = 1L,
  inference = paste("descriptive; cells and pseudotime bins are not",
    "biological replicates"),
  normalized_matrix = matrix_id, bundle_fingerprint_before = fingerprint_before,
  bundle_fingerprint_after = fingerprint_after, bundle_unchanged = TRUE,
  inputs = list(projection = projection_path, density = density_path,
    regions = region_path, genes = gene_path, cytokines = cytokine_path,
    projection_md5 = md5(projection_path), density_md5 = md5(density_path),
    regions_md5 = md5(region_path), genes_md5 = md5(gene_path),
    cytokines_md5 = md5(cytokine_path)),
  package_versions = list(R = R.version.string,
    fgsea = as.character(packageVersion("fgsea")),
    msigdbr = as.character(packageVersion("msigdbr"))))
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# Pseudotime-matched WT-versus-KO preranked GSEA", "",
  "The primary analysis uses the mapped shared-enrichment regions selected before ranking: `Flt3l_shared_peak1` and `Il7_shared_peak2`.",
  "", "Within each region, normalized expression was summarized separately for WT and KO in the existing mapped-pseudotime bins. Each gene was ranked by the Bdev-density-weighted mean of bin-wise KO-versus-WT log2 fold changes. Positive NES denotes KO-higher enrichment; negative NES denotes WT-higher enrichment.",
  "", "Hallmark, GO BP/MF/CC, and the 38 accepted cytokine signatures were tested. BH correction is separate within each comparison and collection.",
  "", "A Hallmark sensitivity analysis excludes bins with fewer than 20 cells in either genotype and renormalizes the Bdev weights.",
  "", "There is one biological sample per genotype. Results describe these measured samples; neither cells nor pseudotime bins provide biological replication."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

wb_close_bundle(staged)
wb_close_bundle(bundle)
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
message("Completed pseudotime-matched GSEA at: ", output_dir)
