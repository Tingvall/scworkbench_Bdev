#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 78_summarize_go_failure_to_silence.R",
  "<dynamic_gsea_dir> <output_dir> <top_n>"))
input_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
top_n <- as.integer(args[[3L]])
stopifnot(is.finite(top_n), top_n >= 1L)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

results_path <- file.path(input_dir,
  "all_collections_dynamic_change_gsea.tsv.gz")
sensitivity_path <- file.path(input_dir,
  "terminal_sparse_bin_sensitivity.tsv.gz")
endpoints_path <- file.path(input_dir, "endpoint_bin_definitions.tsv")
diagnostics_path <- file.path(input_dir, "diagnostics.json")
stopifnot(file.exists(results_path), file.exists(sensitivity_path),
  file.exists(endpoints_path), file.exists(diagnostics_path))
results <- fread(cmd = paste("gzip -cd", shQuote(results_path)))
sensitivity <- fread(cmd = paste("gzip -cd", shQuote(sensitivity_path)))
go_groups <- c("GO_BP", "GO_MF", "GO_CC")
go <- results[collection %chin% go_groups]
stopifnot(setequal(unique(go$collection), go_groups))

# Failure to silence: the set decreases in WT, but KO changes upward relative
# to WT (less decrease, flat, or induction). Negative NES is required because
# the parent rank is (WT change - KO change).
inventory <- go[, .(terms_tested = .N,
  fdr_below_0_05 = sum(padj < 0.05),
  ko_stronger_fdr_below_0_05 = sum(padj < 0.05 & NES < 0),
  strict_failure_to_silence = sum(padj < 0.05 & NES < 0 &
    mean_wt_change < 0 & mean_ko_change > mean_wt_change)),
  by = .(comparison_id, region_id, collection)]
inventory[, group_order := match(collection, go_groups)]
setorder(inventory, comparison_id, group_order)
inventory[, group_order := NULL]
fwrite(inventory, file.path(stage_dir, "go_all_term_test_inventory.tsv"),
  sep = "\t")

eligible <- go[padj < 0.05 & NES < 0 & mean_wt_change < 0 &
  mean_ko_change > mean_wt_change]
setorder(eligible, comparison_id, collection, padj, NES, set_id)
top <- eligible[, head(.SD, top_n), by = .(comparison_id, collection)]
top <- merge(top,
  sensitivity[, .(comparison_id, collection, set_id,
    terminal_sensitivity_NES = sensitivity_NES,
    terminal_sensitivity_padj = sensitivity_padj)],
  by = c("comparison_id", "collection", "set_id"), all.x = TRUE)
top[, terminal_sensitivity_supported := terminal_sensitivity_padj < 0.05 &
  terminal_sensitivity_NES < 0]
top[, group_order := match(collection, go_groups)]
setorder(top, comparison_id, group_order, padj, NES, set_id)
top[, group_order := NULL]
fwrite(top, file.path(stage_dir, "top_go_failure_to_silence.tsv"),
  sep = "\t", na = "NA")

pretty_set <- function(x) {
  x <- sub("^GOBP_", "", x)
  x <- sub("^GOMF_", "", x)
  x <- sub("^GOCC_", "", x)
  tools::toTitleCase(tolower(gsub("_", " ", x)))
}

if (nrow(top)) {
  plot_data <- copy(top)
  plot_data[, label := pretty_set(set_id)]
  plot_data[, facet := factor(collection, levels = go_groups,
    labels = c("GO Biological Process", "GO Molecular Function",
      "GO Cellular Component"))]
  plot_data[, comparison := factor(comparison_id,
    levels = c("Flt3L_KO_vs_WT", "IL7_KO_vs_WT"),
    labels = c("FLT3L", "IL-7"))]
  plot_data[, panel := paste(comparison, facet, sep = " — ")]
  setorder(plot_data, comparison_id, collection, -NES, -padj)
  plot_data[, label := factor(label, levels = unique(label))]
  p <- ggplot(plot_data, aes(-NES, label, size = -log10(padj),
      colour = terminal_sensitivity_supported)) +
    geom_point(alpha = 0.9) +
    facet_wrap(~panel, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B"),
      labels = c(`TRUE` = "Supported", `FALSE` = "Not FDR-supported"),
      name = "Sparse-terminal\nsensitivity") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), axis.title.y = element_blank(),
      strip.text = element_text(face = "bold"), legend.position = "bottom") +
    labs(title = "Top GO programs that fail to silence in KO",
      subtitle = paste("At most", top_n,
        "terms per GO group; all GO terms were tested"),
      x = "Failure-to-silence enrichment magnitude (-NES)",
      size = expression(-log[10](FDR)))
  height <- max(9, 3 + 0.35 * nrow(top))
  ggsave(file.path(stage_dir, "figures",
    "top_go_failure_to_silence.png"), p, width = 10, height = height,
    dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures",
    "top_go_failure_to_silence.pdf"), p, width = 10, height = height,
    device = cairo_pdf, bg = "white")
}

# Draw trajectories for the reported terms that remain FDR-supported when the
# sparse terminal bin is included.
supported <- unique(top[terminal_sensitivity_supported == TRUE,
  .(comparison_id, region_id, collection, set_id, NES, padj,
    mean_wt_change, mean_ko_change, terminal_sensitivity_NES,
    terminal_sensitivity_padj)])
if (nrow(supported)) {
  parent_diagnostics <- jsonlite::read_json(diagnostics_path,
    simplifyVector = TRUE)
  source_dir <- normalizePath(parent_diagnostics$source_run, mustWork = TRUE)
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
  supported_members <- merge(supported,
    memberships[, .(collection, set_id, gene)],
    by = c("collection", "set_id"), allow.cartesian = TRUE)
  curve_data <- merge(supported_members,
    means[, .(comparison_id, region_id, bin, bin_midpoint, bdev_weight,
      reliable_20, primary_window, feature_id, log2_wt, log2_ko)],
    by.x = c("comparison_id", "region_id", "gene"),
    by.y = c("comparison_id", "region_id", "feature_id"),
    allow.cartesian = TRUE)
  stopifnot(uniqueN(curve_data$set_id) == nrow(supported))
  baselines <- curve_data[primary_window == "early", .(
    wt_baseline = weighted.mean(log2_wt, bdev_weight),
    ko_baseline = weighted.mean(log2_ko, bdev_weight)),
    by = .(comparison_id, region_id, collection, set_id, gene)]
  curve_data <- merge(curve_data, baselines,
    by = c("comparison_id", "region_id", "collection", "set_id", "gene"))
  curve_data[, `:=`(wt_centered = log2_wt - wt_baseline,
    ko_centered = log2_ko - ko_baseline)]
  curves <- curve_data[, .(WT = mean(wt_centered), KO = mean(ko_centered),
    member_genes = uniqueN(gene)),
    by = .(comparison_id, region_id, collection, set_id, bin, bin_midpoint,
      reliable_20, NES, padj, mean_wt_change, mean_ko_change,
      terminal_sensitivity_NES, terminal_sensitivity_padj)]
  fwrite(curves, file.path(stage_dir,
    "fdr_supported_failure_to_silence_curves.tsv.gz"), sep = "\t",
    compress = "gzip")
  collection_label <- c(GO_BP = "BP", GO_MF = "MF", GO_CC = "CC")
  curve_labels <- unique(curves[, .(collection, set_id, NES, padj,
    terminal_sensitivity_padj)])
  curve_labels[, facet := sprintf(
    "%s: %s\nNES %.2f, q %.3g; terminal q %.3g",
    collection_label[collection], pretty_set(set_id), NES, padj,
    terminal_sensitivity_padj)]
  curves <- merge(curves,
    curve_labels[, .(collection, set_id, facet)],
    by = c("collection", "set_id"))
  facet_order <- curve_labels[order(match(collection, go_groups), padj)]$facet
  curves[, facet := factor(facet, levels = facet_order)]
  long <- melt(curves,
    id.vars = c("comparison_id", "region_id", "collection", "set_id", "bin",
      "bin_midpoint", "reliable_20", "member_genes", "NES", "padj",
      "mean_wt_change", "mean_ko_change", "terminal_sensitivity_NES",
      "terminal_sensitivity_padj", "facet"),
    measure.vars = c("WT", "KO"), variable.name = "genotype",
    value.name = "baseline_centered_mean_log2_expression")
  curve_plot <- ggplot(long, aes(bin_midpoint,
      baseline_centered_mean_log2_expression, colour = genotype,
      group = genotype)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_line(linewidth = 0.72) +
    geom_point(aes(shape = reliable_20), size = 1.6) +
    facet_wrap(~facet, scales = "free_y", ncol = 2) +
    scale_colour_manual(values = c(WT = "#B2182B", KO = "#2166AC")) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
      labels = c(`TRUE` = "≥20 cells/genotype", `FALSE` = "Sparse terminal"),
      name = "Bin coverage") +
    theme_bw(base_size = 8.5) +
    theme(strip.text = element_text(size = 7), legend.position = "bottom",
      panel.grid.minor = element_blank()) +
    labs(title = "FDR-supported GO trajectories: failure to silence in KO",
      subtitle = paste("WT and KO member-gene expression centered to their",
        "own early baselines; open points mark the sparse terminal bin"),
      x = "Mapped Bdev pseudotime", y = "Mean change from early baseline",
      colour = "Genotype")
  curve_height <- max(10, 2.5 + 2.35 * ceiling(nrow(supported) / 2))
  ggsave(file.path(stage_dir, "figures",
    "fdr_supported_failure_to_silence_curves.png"), curve_plot,
    width = 11, height = curve_height, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures",
    "fdr_supported_failure_to_silence_curves.pdf"), curve_plot,
    width = 11, height = curve_height, device = cairo_pdf, bg = "white")
}

lines <- c("# GO trajectory programs that fail to silence in KO", "",
  paste0("All ", format(sum(inventory$terms_tested), big.mark = ","),
    " comparison-specific GO tests were retained in the parent analysis. This report shows at most ",
    top_n, " significant terms per GO group and comparison."), "",
  "A reported term must decrease from early to late in WT, change upward relative to WT in KO, have negative differential-change NES, and pass primary FDR < 0.05. Positive/negative here refers to trajectory change, not absolute expression.", "")
for (current_comparison in unique(inventory$comparison_id)) {
  lines <- c(lines, paste0("## ", current_comparison), "")
  for (current_group in go_groups) {
    inv <- inventory[comparison_id == current_comparison &
      collection == current_group]
    d <- top[comparison_id == current_comparison & collection == current_group]
    lines <- c(lines, paste0("### ", current_group), "",
      sprintf("Tested %d terms; %d were FDR-significant overall; %d changed more strongly in KO; %d met the strict failure-to-silence criterion.",
        inv$terms_tested, inv$fdr_below_0_05,
        inv$ko_stronger_fdr_below_0_05, inv$strict_failure_to_silence), "")
    if (!nrow(d)) {
      lines <- c(lines,
        "No term met the strict criterion at FDR < 0.05.", "")
    } else {
      lines <- c(lines,
        "| Term | NES | FDR | WT change | KO change | Terminal sensitivity FDR |",
        "|---|---:|---:|---:|---:|---:|",
        vapply(seq_len(nrow(d)), function(i) sprintf(
          "| %s | %.2f | %.3g | %.3f | %.3f | %.3g%s |",
          pretty_set(d$set_id[[i]]), d$NES[[i]], d$padj[[i]],
          d$mean_wt_change[[i]], d$mean_ko_change[[i]],
          d$terminal_sensitivity_padj[[i]],
          if (isTRUE(d$terminal_sensitivity_supported[[i]])) " ✓" else ""),
          character(1)), "")
    }
  }
}
lines <- c(lines, "## Limits", "",
  "There is one biological sample per genotype. These are descriptive within-sample trajectory contrasts; cells and pseudotime bins are not biological replicates. GO terms are redundant and should be interpreted as programs rather than independent discoveries.")
writeLines(lines, file.path(stage_dir, "REPORT.md"))

diagnostics <- list(status = "complete_report_only",
  source_run = input_dir, top_n_per_comparison_and_go_group = top_n,
  direction = "failure_to_silence_in_KO",
  filter = paste("padj < 0.05, NES < 0, mean WT change < 0,",
    "mean KO change > mean WT change"),
  all_go_terms_tested_in_parent = TRUE, inventory = inventory,
  reported_terms = nrow(top),
  reported_terms_terminal_sensitivity_supported = nrow(supported),
  biological_replicates_per_genotype = 1L)
jsonlite::write_json(diagnostics, file.path(stage_dir, "diagnostics.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
message("Completed GO failure-to-silence report at: ", output_dir)
