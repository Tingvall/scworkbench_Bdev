#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 40_fdr_supported_go_dynamic_curves.R",
  "<dynamic_gsea_dir> <top_go_dir> <output_dir>"))
dynamic_dir <- normalizePath(args[[1L]], mustWork = TRUE)
top_go_dir <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

top_path <- file.path(top_go_dir,
  "top_go_wt_induction_attenuated_in_ko.tsv")
endpoints_path <- file.path(dynamic_dir, "endpoint_bin_definitions.tsv")
diagnostics_path <- file.path(dynamic_dir, "diagnostics.json")
stopifnot(file.exists(top_path), file.exists(endpoints_path),
  file.exists(diagnostics_path))
top <- fread(top_path)
selected <- unique(top[terminal_sensitivity_supported == TRUE,
  .(comparison_id, region_id, collection, set_id, NES, padj,
    mean_wt_change, mean_ko_change, terminal_sensitivity_NES,
    terminal_sensitivity_padj)])
stopifnot(nrow(selected) > 0L,
  all(selected$padj < 0.05),
  all(selected$terminal_sensitivity_padj < 0.05),
  all(selected$NES > 0), all(selected$terminal_sensitivity_NES > 0))

dynamic_diagnostics <- jsonlite::read_json(diagnostics_path,
  simplifyVector = TRUE)
source_dir <- normalizePath(dynamic_diagnostics$source_run, mustWork = TRUE)
means_path <- file.path(source_dir,
  "pseudotime_bin_gene_expression_means.tsv.gz")
memberships_path <- file.path(source_dir, "gene_set_memberships.tsv.gz")
stopifnot(file.exists(means_path), file.exists(memberships_path))
means <- fread(cmd = paste("gzip -cd", shQuote(means_path)))
memberships <- fread(cmd = paste("gzip -cd", shQuote(memberships_path)))
endpoints <- fread(endpoints_path)
means <- merge(means,
  endpoints[, .(region_id, comparison_id, bin, primary_window)],
  by = c("region_id", "comparison_id", "bin"), all.x = TRUE)
means[, `:=`(log2_wt = log2(mean_cp10k_wt + 1),
  log2_ko = log2(mean_cp10k_ko + 1))]

selected_members <- merge(selected,
  memberships[, .(collection, set_id, gene)],
  by = c("collection", "set_id"), allow.cartesian = TRUE)
curve_data <- merge(selected_members,
  means[, .(comparison_id, region_id, bin, bin_midpoint, bdev_weight,
    reliable_20, primary_window, feature_id, log2_wt, log2_ko)],
  by.x = c("comparison_id", "region_id", "gene"),
  by.y = c("comparison_id", "region_id", "feature_id"),
  allow.cartesian = TRUE)
stopifnot(uniqueN(curve_data$set_id) == nrow(selected))

# Center every member gene to its genotype-specific, Bdev-density-weighted
# primary early baseline before averaging genes. This is identical to the
# curve construction used in the parent dynamic GSEA report.
baselines <- curve_data[primary_window == "early", .(
  wt_baseline = weighted.mean(log2_wt, bdev_weight),
  ko_baseline = weighted.mean(log2_ko, bdev_weight)),
  by = .(comparison_id, region_id, collection, set_id, gene)]
curve_data <- merge(curve_data, baselines,
  by = c("comparison_id", "region_id", "collection", "set_id", "gene"))
curve_data[, `:=`(wt_centered = log2_wt - wt_baseline,
  ko_centered = log2_ko - ko_baseline)]
curves <- curve_data[, .(
  WT = mean(wt_centered), KO = mean(ko_centered),
  member_genes = uniqueN(gene)),
  by = .(comparison_id, region_id, collection, set_id, bin, bin_midpoint,
    reliable_20, NES, padj, mean_wt_change, mean_ko_change,
    terminal_sensitivity_NES, terminal_sensitivity_padj)]
curves[, differential_curve := WT - KO]
fwrite(curves, file.path(stage_dir,
  "fdr_supported_go_dynamic_curves.tsv.gz"), sep = "\t", compress = "gzip")

pretty_set <- function(x) {
  x <- sub("^GOBP_", "", x)
  x <- sub("^GOMF_", "", x)
  x <- sub("^GOCC_", "", x)
  tools::toTitleCase(tolower(gsub("_", " ", x)))
}
collection_label <- c(GO_BP = "BP", GO_MF = "MF", GO_CC = "CC")
labels <- unique(curves[, .(collection, set_id, NES, padj,
  terminal_sensitivity_padj)])
labels[, facet := sprintf("%s: %s\nNES %.2f, q %.3g; terminal q %.3g",
  collection_label[collection], pretty_set(set_id), NES, padj,
  terminal_sensitivity_padj)]
curves <- merge(curves, labels[, .(collection, set_id, facet)],
  by = c("collection", "set_id"))
facet_order <- labels[order(match(collection, c("GO_BP", "GO_MF", "GO_CC")),
  padj)]$facet
curves[, facet := factor(facet, levels = facet_order)]

long <- melt(curves,
  id.vars = c("comparison_id", "region_id", "collection", "set_id", "bin",
    "bin_midpoint", "reliable_20", "member_genes", "NES", "padj",
    "mean_wt_change", "mean_ko_change", "terminal_sensitivity_NES",
    "terminal_sensitivity_padj", "differential_curve", "facet"),
  measure.vars = c("WT", "KO"), variable.name = "genotype",
  value.name = "baseline_centered_mean_log2_expression")

p <- ggplot(long, aes(bin_midpoint,
    baseline_centered_mean_log2_expression, colour = genotype,
    group = genotype)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_line(linewidth = 0.75) +
  geom_point(aes(shape = reliable_20), size = 1.8) +
  facet_wrap(~facet, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(WT = "#B2182B", KO = "#2166AC")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "≥20 cells/genotype", `FALSE` = "Sparse terminal"),
    name = "Bin coverage") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7.5), legend.position = "bottom",
    panel.grid.minor = element_blank()) +
  labs(title = "FDR-supported GO trajectories: WT induction attenuated in KO",
    subtitle = paste("Member-gene log2 expression centered to each genotype's",
      "early baseline; every term is supported with and without the sparse terminal bin"),
    x = "Mapped Bdev pseudotime",
    y = "Mean change from early baseline", colour = "Genotype")
height <- max(8, 2.6 + 2.5 * ceiling(uniqueN(curves$set_id) / 2))
ggsave(file.path(stage_dir, "fdr_supported_go_dynamic_curves.png"), p,
  width = 11, height = height, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "fdr_supported_go_dynamic_curves.pdf"), p,
  width = 11, height = height, device = cairo_pdf, bg = "white")

d <- unique(curves[, .(comparison_id, collection, set_id, member_genes,
  NES, padj, mean_wt_change, mean_ko_change, terminal_sensitivity_NES,
  terminal_sensitivity_padj)])
d[, group_order := match(collection, c("GO_BP", "GO_MF", "GO_CC"))]
setorder(d, comparison_id, group_order, padj)
d[, group_order := NULL]
report <- c("# FDR-supported GO trajectory curves", "",
  paste0("Curves are shown for ", nrow(d),
    " strict WT-induction/KO-attenuation GO terms that are FDR-significant in both the primary and sparse-terminal-bin sensitivity analyses."), "",
  "Each gene is centered to its genotype-specific early baseline before member genes are averaged. Open points mark terminal bins with fewer than 20 KO cells; those points are displayed but do not define the primary late endpoint.", "",
  "| GO group | Term | Genes | NES | FDR | WT change | KO change | Terminal FDR |",
  "|---|---|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(d)), function(i) sprintf(
    "| %s | %s | %d | %.2f | %.3g | %.3f | %.3f | %.3g |",
    d$collection[[i]], pretty_set(d$set_id[[i]]), d$member_genes[[i]],
    d$NES[[i]], d$padj[[i]], d$mean_wt_change[[i]],
    d$mean_ko_change[[i]], d$terminal_sensitivity_padj[[i]]), character(1)),
  "", "There is one biological sample per genotype; curves and enrichment results are descriptive within-sample trajectory contrasts.")
writeLines(report, file.path(stage_dir, "REPORT.md"))

diagnostics <- list(status = "complete_report_only",
  source_dynamic_gsea = dynamic_dir, source_top_go = top_go_dir,
  selection = paste("strict WT induction attenuated in KO with primary and",
    "terminal-sensitivity FDR < 0.05 and positive NES"),
  selected_terms = nrow(selected),
  score = paste("mean member-gene log2(mean CP10K + 1), each gene centered",
    "to genotype-specific early baseline"),
  sparse_terminal_bins_displayed = TRUE,
  sparse_terminal_bins_used_for_primary_endpoints = FALSE,
  biological_replicates_per_genotype = 1L)
jsonlite::write_json(diagnostics, file.path(stage_dir, "diagnostics.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
message("Completed FDR-supported GO curves at: ", output_dir)
