#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 76_mapped_enrichment_dynamic_change_gsea.R",
  "<pseudotime_matched_gsea_dir> <output_dir> <seed>"))
input_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
seed <- as.integer(args[[3L]])
stopifnot(is.finite(seed))
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(fgsea)
  library(jsonlite)
})

write_json <- function(x, path) jsonlite::write_json(x, path,
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
md5 <- function(path) unname(tools::md5sum(path))
pretty_set <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^GOBP_", "", x)
  x <- sub("^GOMF_", "", x)
  x <- sub("^GOCC_", "", x)
  x <- sub("^parse_cytokine_rest_only_refinement_all_top100_v1\\.", "", x)
  tools::toTitleCase(tolower(gsub("[_.]", " ", x)))
}

means_path <- file.path(input_dir, "pseudotime_bin_gene_expression_means.tsv.gz")
membership_path <- file.path(input_dir, "gene_set_memberships.tsv.gz")
counts_path <- file.path(input_dir, "pseudotime_bin_cell_counts.tsv")
stopifnot(file.exists(means_path), file.exists(membership_path),
  file.exists(counts_path))
message("Reading matched-bin expression means and frozen memberships")
means <- fread(cmd = paste("gzip -cd", shQuote(means_path)))
memberships <- fread(cmd = paste("gzip -cd", shQuote(membership_path)))
counts <- fread(counts_path)
required <- c("region_id", "comparison_id", "bin", "bin_midpoint",
  "bdev_weight", "n_wt", "n_ko", "reliable_20", "feature_id",
  "mean_cp10k_wt", "mean_cp10k_ko")
stopifnot(all(required %in% names(means)),
  all(c("collection", "set_id", "set_label", "gene") %in% names(memberships)),
  setequal(unique(memberships$collection),
    c("Hallmark", "GO_BP", "GO_MF", "GO_CC", "Cytokine")))

# The primary endpoints use the first two and last two bins having at least 20
# cells in each genotype. This avoids calling loss of terminal KO coverage an
# expression trajectory. The sensitivity endpoints use the literal first and
# last two mapped-region bins, including the sparse terminal bin.
bin_windows <- unique(means[, .(region_id, comparison_id, bin, bin_midpoint,
  bdev_weight, n_wt, n_ko, reliable_20)])
setorder(bin_windows, comparison_id, bin)
bin_windows[, primary_window := "middle"]
bin_windows[reliable_20 == TRUE, primary_window := {
  z <- rep("middle", .N)
  z[seq_len(min(2L, .N))] <- "early"
  z[seq.int(max(1L, .N - 1L), .N)] <- "late"
  z
}, by = comparison_id]
bin_windows[, sensitivity_window := {
  z <- rep("middle", .N)
  z[seq_len(min(2L, .N))] <- "early"
  z[seq.int(max(1L, .N - 1L), .N)] <- "late"
  z
}, by = comparison_id]
stopifnot(bin_windows[, all(c("early", "late") %chin% primary_window),
  by = comparison_id]$V1 |> all())
fwrite(bin_windows, file.path(stage_dir, "endpoint_bin_definitions.tsv"),
  sep = "\t")
means <- merge(means, bin_windows[, .(region_id, comparison_id, bin,
  primary_window, sensitivity_window)],
  by = c("region_id", "comparison_id", "bin"), all.x = TRUE)
means[, `:=`(log2_wt = log2(mean_cp10k_wt + 1),
  log2_ko = log2(mean_cp10k_ko + 1))]

weighted_window <- function(value, weight, window, label) {
  use <- window == label
  sum(value[use] * weight[use]) / sum(weight[use])
}
weighted_slope <- function(value, x, weight, use) {
  value <- value[use]
  x <- x[use]
  weight <- weight[use] / sum(weight[use])
  xbar <- sum(weight * x)
  sum(weight * (x - xbar) * value) / sum(weight * (x - xbar)^2)
}

message("Calculating gene-level genotype-by-pseudotime change statistics")
gene_stats <- means[, {
  x <- (bin_midpoint - min(bin_midpoint)) /
    (max(bin_midpoint) - min(bin_midpoint))
  wt_early <- weighted_window(log2_wt, bdev_weight, primary_window, "early")
  wt_late <- weighted_window(log2_wt, bdev_weight, primary_window, "late")
  ko_early <- weighted_window(log2_ko, bdev_weight, primary_window, "early")
  ko_late <- weighted_window(log2_ko, bdev_weight, primary_window, "late")
  swt_early <- weighted_window(log2_wt, bdev_weight, sensitivity_window, "early")
  swt_late <- weighted_window(log2_wt, bdev_weight, sensitivity_window, "late")
  sko_early <- weighted_window(log2_ko, bdev_weight, sensitivity_window, "early")
  sko_late <- weighted_window(log2_ko, bdev_weight, sensitivity_window, "late")
  reliable <- reliable_20
  .(wt_early_log2 = wt_early, wt_late_log2 = wt_late,
    ko_early_log2 = ko_early, ko_late_log2 = ko_late,
    wt_change = wt_late - wt_early, ko_change = ko_late - ko_early,
    differential_change = (wt_late - wt_early) - (ko_late - ko_early),
    sensitivity_wt_change = swt_late - swt_early,
    sensitivity_ko_change = sko_late - sko_early,
    sensitivity_differential_change =
      (swt_late - swt_early) - (sko_late - sko_early),
    wt_reliable_slope = weighted_slope(log2_wt, x, bdev_weight, reliable),
    ko_reliable_slope = weighted_slope(log2_ko, x, bdev_weight, reliable),
    slope_interaction = weighted_slope(log2_wt, x, bdev_weight, reliable) -
      weighted_slope(log2_ko, x, bdev_weight, reliable),
    maximum_mean_cp10k = max(mean_cp10k_wt, mean_cp10k_ko))
}, by = .(region_id, comparison_id, feature_id)]
gene_stats <- gene_stats[maximum_mean_cp10k > 0]

break_ties <- function(d, primary, secondary) {
  setorderv(d, c(primary, secondary, "feature_id"), c(1L, 1L, 1L))
  values <- d[[primary]]
  unique_scores <- sort(unique(values))
  gaps <- diff(unique_scores)
  minimum_gap <- if (length(gaps)) min(gaps) else 1
  d[, rank_score := get(primary)]
  d[, tie_size := .N, by = primary]
  for (value in d[tie_size > 1L, unique(get(primary))]) {
    idx <- which(d[[primary]] == value)
    local_gaps <- c(value - max(unique_scores[unique_scores < value],
      na.rm = TRUE), min(unique_scores[unique_scores > value],
      na.rm = TRUE) - value)
    local_gaps <- local_gaps[is.finite(local_gaps) & local_gaps > 0]
    margin <- if (length(local_gaps)) min(local_gaps) / 4 else minimum_gap / 4
    d$rank_score[idx] <- value + seq(-margin, margin, length.out = length(idx))
  }
  stopifnot(!anyDuplicated(d$rank_score))
  d
}

ranking_list <- lapply(split(gene_stats, gene_stats$comparison_id), function(d) {
  d <- break_ties(copy(d), "differential_change", "slope_interaction")
  setorder(d, -rank_score)
  d[, rank := seq_len(.N)]
  d
})
gene_stats <- rbindlist(ranking_list)
fwrite(gene_stats, file.path(stage_dir,
  "gene_dynamic_change_rankings.tsv.gz"), sep = "\t", compress = "gzip")

run_one <- function(ranking, collection_name, seed_offset,
    score_column = "rank_score") {
  pathways <- split(memberships[collection == collection_name]$gene,
    memberships[collection == collection_name]$set_id)
  stats <- ranking[[score_column]]
  names(stats) <- ranking$feature_id
  stats <- sort(stats, decreasing = TRUE)
  set.seed(seed + seed_offset)
  ans <- as.data.table(fgseaMultilevel(pathways = pathways, stats = stats,
    minSize = 15L, maxSize = 500L, eps = 0, scoreType = "std", nproc = 1L))
  setnames(ans, "pathway", "set_id")
  ans[, `:=`(comparison_id = ranking$comparison_id[[1L]],
    region_id = ranking$region_id[[1L]], collection = collection_name,
    enrichment_direction = fifelse(NES > 0,
      "WT stronger increase", "KO stronger increase"),
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
message("Running differential-change preranked GSEA")
offset <- 0L
results_list <- list()
for (comparison_id in names(ranking_list)) {
  for (collection_name in collections) {
    offset <- offset + 1L
    message("GSEA: ", comparison_id, " / ", collection_name)
    results_list[[paste(comparison_id, collection_name)]] <- run_one(
      ranking_list[[comparison_id]], collection_name, offset)
  }
}
results <- rbindlist(results_list)

# Average member-gene changes are used only to give the directional GSEA calls
# a readable dynamics class. They are not substituted for the GSEA statistic.
set_changes <- merge(memberships[, .(collection, set_id, gene)],
  gene_stats, by.x = "gene", by.y = "feature_id", allow.cartesian = TRUE)
set_changes <- set_changes[, .(member_genes_detected = .N,
  mean_wt_change = mean(wt_change), mean_ko_change = mean(ko_change),
  mean_differential_change = mean(differential_change),
  mean_slope_interaction = mean(slope_interaction)),
  by = .(comparison_id, region_id, collection, set_id)]
results <- merge(results, set_changes,
  by = c("comparison_id", "region_id", "collection", "set_id"), all.x = TRUE)
results[, dynamics_class := fcase(
  mean_wt_change > 0 & mean_ko_change <= 0,
    "WT induction; KO flat/declining",
  mean_wt_change > mean_ko_change & mean_ko_change > 0,
    "shared induction; weaker in KO",
  mean_wt_change > 0 & mean_ko_change < mean_wt_change,
    "WT induction; attenuated in KO",
  mean_ko_change > 0 & mean_wt_change <= 0,
    "KO induction; WT flat/declining",
  mean_ko_change > mean_wt_change & mean_wt_change > 0,
    "shared induction; stronger in KO",
  mean_wt_change <= 0 & mean_ko_change < mean_wt_change,
    "KO-selective decline (not WT induction)",
  mean_ko_change <= 0 & mean_wt_change < mean_ko_change,
    "WT-selective decline (not KO induction)",
  default = "mixed/weak endpoint direction")]
results[, absolute_NES_for_order := abs(NES)]
setorder(results, comparison_id, collection, padj, -absolute_NES_for_order,
  set_id)
results[, absolute_NES_for_order := NULL]
fwrite(results, file.path(stage_dir,
  "all_collections_dynamic_change_gsea.tsv.gz"), sep = "\t",
  compress = "gzip", na = "NA")

summary <- results[, .(sets_tested = .N,
  fdr_below_0_05 = sum(padj < 0.05),
  wt_stronger_fdr_below_0_05 = sum(padj < 0.05 & NES > 0),
  ko_stronger_fdr_below_0_05 = sum(padj < 0.05 & NES < 0),
  wt_induction_ko_flat_or_declining = sum(padj < 0.05 & NES > 0 &
    dynamics_class == "WT induction; KO flat/declining"),
  shared_induction_weaker_in_ko = sum(padj < 0.05 & NES > 0 &
    dynamics_class == "shared induction; weaker in KO")),
  by = .(comparison_id, region_id, collection)]
fwrite(summary, file.path(stage_dir, "dynamic_gsea_summary.tsv"), sep = "\t")

# Sensitivity: literal terminal bins (including the n_KO < 20 terminal bin).
sensitivity_results <- list()
for (comparison_id in names(ranking_list)) {
  d <- copy(ranking_list[[comparison_id]])
  d <- break_ties(d, "sensitivity_differential_change", "differential_change")
  for (collection_name in collections) {
    offset <- offset + 1L
    sensitivity_results[[paste(comparison_id, collection_name)]] <-
      run_one(d, collection_name, offset)
  }
}
sensitivity <- rbindlist(sensitivity_results)
audit <- merge(results[, .(comparison_id, collection, set_id, NES, padj)],
  sensitivity[, .(comparison_id, collection, set_id,
    sensitivity_NES = NES, sensitivity_padj = padj)],
  by = c("comparison_id", "collection", "set_id"))
audit_summary <- audit[, .(NES_spearman = cor(NES, sensitivity_NES,
  method = "spearman"), primary_significant = sum(padj < 0.05),
  sensitivity_significant = sum(sensitivity_padj < 0.05)),
  by = .(comparison_id, collection)]
fwrite(audit, file.path(stage_dir, "terminal_sparse_bin_sensitivity.tsv.gz"),
  sep = "\t", compress = "gzip")
fwrite(audit_summary,
  file.path(stage_dir, "terminal_sparse_bin_sensitivity_summary.tsv"),
  sep = "\t")

make_heatmap <- function(collection_name, max_each = Inf) {
  d <- copy(results[collection == collection_name])
  if (is.finite(max_each)) {
    keep <- unique(d[order(padj, -abs(NES)), head(set_id, max_each),
      by = comparison_id]$V1)
    d <- d[set_id %chin% keep]
  }
  d[, display := pretty_set(set_id)]
  d[, maximum := max(abs(NES)), by = set_id]
  d[, average := mean(NES), by = set_id]
  setorder(d, maximum, average, set_id)
  d[, display := factor(display, levels = unique(display))]
  d[, comparison_label := factor(comparison_id,
    levels = c("Flt3L_KO_vs_WT", "IL7_KO_vs_WT"),
    labels = c("FLT3L", "IL-7"))]
  d[, mark := fifelse(padj < 0.001, "***", fifelse(padj < 0.01, "**",
    fifelse(padj < 0.05, "*", "")))]
  lim <- max(abs(d$NES))
  p <- ggplot(d, aes(comparison_label, display, fill = NES)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(aes(label = mark), size = 2.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-lim, lim), name = "NES") +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6.6),
      axis.title = element_blank(), plot.title = element_text(face = "bold")) +
    labs(title = paste(collection_name, "differential-change GSEA"),
      subtitle = paste("Red/positive: stronger early-to-late increase in WT;",
        "blue/negative: stronger increase in KO; stars: FDR"))
  stub <- tolower(collection_name)
  height <- max(6, min(18, 2.5 + 0.19 * uniqueN(d$set_id)))
  ggsave(file.path(stage_dir, "figures", paste0(stub,
    "_dynamic_change_gsea_heatmap.png")), p, width = 8.5, height = height,
    dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(stub,
    "_dynamic_change_gsea_heatmap.pdf")), p, width = 8.5, height = height,
    device = cairo_pdf, bg = "white")
}
make_heatmap("Hallmark")
make_heatmap("GO_BP", 20L)
make_heatmap("GO_MF", 20L)
make_heatmap("GO_CC", 20L)
make_heatmap("Cytokine")

# Curves for the most significant Hallmark and cytokine sets. Each member gene
# is centered to its own genotype-specific primary early baseline before the
# genes are averaged; zero therefore means that genotype's early state.
curve_selection <- results[collection %chin% c("Hallmark", "Cytokine"),
  head(.SD[order(padj, -abs(NES))], 12L),
  by = .(comparison_id, collection)]
curve_keys <- unique(curve_selection[, .(comparison_id, collection, set_id,
  padj, NES, dynamics_class)])
curve_members <- merge(curve_keys,
  memberships[, .(collection, set_id, gene)],
  by = c("collection", "set_id"), allow.cartesian = TRUE)
curve_data <- merge(curve_members,
  means[, .(comparison_id, bin, bin_midpoint, bdev_weight, reliable_20,
    feature_id, log2_wt, log2_ko, primary_window)],
  by.x = c("comparison_id", "gene"),
  by.y = c("comparison_id", "feature_id"), allow.cartesian = TRUE)
curve_baselines <- curve_data[primary_window == "early", .(
  wt_baseline = weighted.mean(log2_wt, bdev_weight),
  ko_baseline = weighted.mean(log2_ko, bdev_weight)),
  by = .(comparison_id, collection, set_id, gene)]
curve_data <- merge(curve_data, curve_baselines,
  by = c("comparison_id", "collection", "set_id", "gene"))
curve_data[, `:=`(wt_centered = log2_wt - wt_baseline,
  ko_centered = log2_ko - ko_baseline)]
curves <- curve_data[, .(WT = mean(wt_centered), KO = mean(ko_centered),
  member_genes = uniqueN(gene)),
  by = .(comparison_id, collection, set_id, bin, bin_midpoint, reliable_20,
    padj, NES, dynamics_class)]
curves <- melt(curves,
  id.vars = c("comparison_id", "collection", "set_id", "bin",
    "bin_midpoint", "reliable_20", "member_genes", "padj", "NES",
    "dynamics_class"), variable.name = "genotype",
  value.name = "baseline_centered_mean_log2_expression")
fwrite(curves, file.path(stage_dir,
  "selected_hallmark_cytokine_dynamic_curves.tsv.gz"), sep = "\t",
  compress = "gzip")

for (current_comparison in unique(curves$comparison_id)) {
  for (collection_name in c("Hallmark", "Cytokine")) {
    d <- copy(curves[comparison_id == current_comparison &
      collection == collection_name])
    labels <- unique(d[, .(set_id, padj, NES, dynamics_class)])
    labels[, facet := sprintf("%s\nNES %.2f, q %.2g\n%s",
      pretty_set(set_id), NES, padj, dynamics_class)]
    d <- merge(d, labels[, .(set_id, facet)], by = "set_id")
    d[, facet := factor(facet, levels = unique(facet))]
    p <- ggplot(d, aes(bin_midpoint,
        baseline_centered_mean_log2_expression, colour = genotype,
        group = genotype)) +
      geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
      geom_line(linewidth = 0.65) +
      geom_point(aes(shape = reliable_20), size = 1.4) +
      facet_wrap(~facet, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c(WT = "#B2182B", KO = "#2166AC")) +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
        name = "≥20 cells/genotype") +
      theme_bw(base_size = 8) +
      theme(strip.text = element_text(size = 6.5),
        legend.position = "bottom") +
      labs(title = paste(current_comparison, collection_name,
        "gene-set trajectories"),
        subtitle = paste("Member-gene log2 expression centered to each",
          "genotype's early baseline; open point marks sparse terminal bin"),
        x = "Mapped Bdev pseudotime",
        y = "Mean change from early baseline", colour = "Genotype")
    stub <- paste0(tolower(sub("_KO_vs_WT", "", current_comparison)), "_",
      tolower(collection_name), "_dynamic_curves")
    ggsave(file.path(stage_dir, "figures", paste0(stub, ".png")), p,
      width = 10.5, height = 11, dpi = 300, bg = "white")
    ggsave(file.path(stage_dir, "figures", paste0(stub, ".pdf")), p,
      width = 10.5, height = 11, device = cairo_pdf, bg = "white")
  }
}

top_failed <- results[padj < 0.05 & NES > 0 & mean_wt_change > 0 &
  mean_ko_change < mean_wt_change][order(comparison_id, collection, padj,
    -NES), head(.SD, 10L), by = .(comparison_id, collection)]
fwrite(top_failed, file.path(stage_dir,
  "significant_wt_induction_attenuated_in_ko.tsv"), sep = "\t", na = "NA")

report_lines <- c("# Differential trajectory-change WT-versus-KO GSEA", "",
  "This analysis asks which gene sets increase across the mapped enrichment region in WT but increase less, remain flat, or decline in the matched KO cells.", "",
  "## Design", "",
  "- Primary gene statistic: `(WT late - WT early) - (KO late - KO early)` on log2(mean linear CP10K + 1). Positive values mean stronger induction in WT.",
  "- Early and late are the first two and last two bins with at least 20 cells in each genotype. Bdev-density weights are renormalized within each endpoint.",
  "- A sensitivity analysis uses the literal first and last two mapped-region bins, including the terminal bin with fewer than 20 KO cells.",
  "- Hallmark, GO BP/MF/CC, and the 38 accepted cytokine signatures are tested using fgseaMultilevel; BH correction is separate within comparison and collection.",
  "- Gene-set curves are descriptive mean member-gene log-expression changes after centering each genotype to its own early baseline.", "",
  "## Significant results", "")
for (current_comparison in unique(summary$comparison_id)) {
  report_lines <- c(report_lines, paste0("### ", current_comparison), "",
    "| Collection | FDR < 0.05 | WT stronger | KO stronger | WT induction, KO flat/declining | Shared induction, weaker KO |",
    "|---|---:|---:|---:|---:|---:|")
  ss <- summary[comparison_id == current_comparison]
  report_lines <- c(report_lines, apply(ss, 1L, function(z) sprintf(
    "| %s | %s | %s | %s | %s | %s |", z[["collection"]],
    z[["fdr_below_0_05"]], z[["wt_stronger_fdr_below_0_05"]],
    z[["ko_stronger_fdr_below_0_05"]],
    z[["wt_induction_ko_flat_or_declining"]],
    z[["shared_induction_weaker_in_ko"]])), "")
  tf <- top_failed[comparison_id == current_comparison]
  if (nrow(tf)) {
    report_lines <- c(report_lines,
      "Top significant sets consistent with WT induction attenuated in KO:", "",
      vapply(seq_len(min(15L, nrow(tf))), function(i) sprintf(
        "- %s / %s: NES %.2f, q %.3g; WT change %.3f, KO change %.3f (%s).",
        tf$collection[[i]], pretty_set(tf$set_id[[i]]), tf$NES[[i]],
        tf$padj[[i]], tf$mean_wt_change[[i]], tf$mean_ko_change[[i]],
        tf$dynamics_class[[i]]), character(1)), "")
  }
}
report_lines <- c(report_lines, "## Interpretation limits", "",
  "There is one biological sample per genotype. These are descriptive within-sample trajectory contrasts; cells and pseudotime bins are not biological replicates. A positive GSEA result alone is not labeled failed induction unless the member-gene summary also increases in WT and increases less in KO.")
writeLines(report_lines, file.path(stage_dir, "REPORT.md"))

diagnostics <- list(
  status = "complete_unpromoted_report_only_outside_scworkbench",
  rationale = paste("scWorkbench pathway_enrich does not implement",
    "preranked genotype-by-pseudotime differential-change GSEA"),
  user_approved_outside_analysis = TRUE,
  source_run = input_dir, region_policy = "mapped shared enrichment",
  primary_statistic = "(WT_late-WT_early)-(KO_late-KO_early)",
  expression_scale = "log2(mean linear CP10K + 1) per bin and genotype",
  primary_endpoints = "first/last two bins with n_WT>=20 and n_KO>=20",
  sensitivity_endpoints = "literal first/last two region bins",
  curve_score = paste("mean member-gene log2 expression after each gene and",
    "genotype is centered to its primary early baseline"),
  gsea = list(method = "fgseaMultilevel", min_size = 15L,
    max_size = 500L, eps = 0, seed = seed,
    p_adjustment = "BH separately within comparison and collection"),
  summary = summary, sensitivity = audit_summary,
  biological_replicates_per_genotype = 1L,
  inference = "descriptive; cells and bins are not biological replicates",
  inputs = list(means = means_path, memberships = membership_path,
    counts = counts_path, means_md5 = md5(means_path),
    memberships_md5 = md5(membership_path), counts_md5 = md5(counts_path)),
  package_versions = list(R = R.version.string,
    fgsea = as.character(packageVersion("fgsea"))))
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
message("Completed differential trajectory-change GSEA at: ", output_dir)
