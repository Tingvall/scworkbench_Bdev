#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 77_summarize_dynamic_go_top_terms.R",
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
stopifnot(file.exists(results_path), file.exists(sensitivity_path))
results <- fread(cmd = paste("gzip -cd", shQuote(results_path)))
sensitivity <- fread(cmd = paste("gzip -cd", shQuote(sensitivity_path)))
go_groups <- c("GO_BP", "GO_MF", "GO_CC")
go <- results[collection %chin% go_groups]
stopifnot(setequal(unique(go$collection), go_groups))

inventory <- go[, .(terms_tested = .N,
  fdr_below_0_05 = sum(padj < 0.05),
  wt_stronger_fdr_below_0_05 = sum(padj < 0.05 & NES > 0),
  strict_wt_induction_attenuated_in_ko = sum(padj < 0.05 & NES > 0 &
  mean_wt_change > 0 & mean_ko_change < mean_wt_change)),
  by = .(comparison_id, region_id, collection)]
inventory[, group_order := match(collection, go_groups)]
setorder(inventory, comparison_id, group_order)
inventory[, group_order := NULL]
fwrite(inventory, file.path(stage_dir, "go_all_term_test_inventory.tsv"),
  sep = "\t")

# Restrict the report to significant terms that answer the stated biological
# question. All tested terms remain available in the parent run.
eligible <- go[padj < 0.05 & NES > 0 & mean_wt_change > 0 &
  mean_ko_change < mean_wt_change]
setorder(eligible, comparison_id, collection, padj, -NES, set_id)
top <- eligible[, head(.SD, top_n), by = .(comparison_id, collection)]
top <- merge(top,
  sensitivity[, .(comparison_id, collection, set_id,
    terminal_sensitivity_NES = sensitivity_NES,
    terminal_sensitivity_padj = sensitivity_padj)],
  by = c("comparison_id", "collection", "set_id"), all.x = TRUE)
top[, terminal_sensitivity_supported := terminal_sensitivity_padj < 0.05 &
  terminal_sensitivity_NES > 0]
top[, group_order := match(collection, go_groups)]
setorder(top, comparison_id, group_order, padj, -NES, set_id)
top[, group_order := NULL]
fwrite(top, file.path(stage_dir,
  "top_go_wt_induction_attenuated_in_ko.tsv"), sep = "\t", na = "NA")

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
  setorder(plot_data, comparison_id, collection, NES, padj)
  plot_data[, label := factor(label, levels = unique(label))]
  p <- ggplot(plot_data, aes(NES, label, size = -log10(padj),
      colour = terminal_sensitivity_supported)) +
    geom_point(alpha = 0.9) +
    facet_wrap(~panel, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B"),
      labels = c(`TRUE` = "Supported", `FALSE` = "Not FDR-supported"),
      name = "Sparse-terminal\nsensitivity") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), axis.title.y = element_blank(),
      strip.text = element_text(face = "bold"), legend.position = "bottom") +
    labs(title = "Top GO programs induced in WT and attenuated in KO",
      subtitle = paste("At most", top_n,
        "terms per GO group; all GO terms were tested"),
      x = "Differential-change GSEA NES",
      size = expression(-log[10](FDR)))
  height <- max(8, 3 + 0.35 * nrow(top))
  ggsave(file.path(stage_dir, "figures", "top_go_dynamic_change.png"), p,
    width = 10, height = height, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", "top_go_dynamic_change.pdf"), p,
    width = 10, height = height, device = cairo_pdf, bg = "white")
}

lines <- c("# Top GO differential-trajectory programs", "",
  paste0("All ", format(sum(inventory$terms_tested), big.mark = ","),
    " comparison-specific GO tests from BP, MF, and CC were retained in the parent result table. This focused report shows at most ",
    top_n, " significant terms per GO group and comparison."), "",
  "A term is reported only when the differential-change GSEA is significant and WT member-gene expression rises from early to late while KO rises less or declines. Positive enrichment caused only by a stronger KO decline is not reported here.", "")

for (current_comparison in unique(inventory$comparison_id)) {
  lines <- c(lines, paste0("## ", current_comparison), "")
  for (current_group in go_groups) {
    inv <- inventory[comparison_id == current_comparison &
      collection == current_group]
    d <- top[comparison_id == current_comparison & collection == current_group]
    lines <- c(lines, paste0("### ", current_group), "",
      sprintf("Tested %d terms; %d were FDR-significant overall; %d had stronger WT change; %d met the strict WT-induction/KO-attenuation criterion.",
        inv$terms_tested, inv$fdr_below_0_05,
        inv$wt_stronger_fdr_below_0_05,
        inv$strict_wt_induction_attenuated_in_ko), "")
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
  "There is one biological sample per genotype. These are descriptive within-sample trajectory contrasts; cells and pseudotime bins are not biological replicates. GO terms are redundant, so related terms should be interpreted as programs rather than independent discoveries.")
writeLines(lines, file.path(stage_dir, "REPORT.md"))

diagnostics <- list(status = "complete_report_only",
  source_run = input_dir, top_n_per_comparison_and_go_group = top_n,
  reporting_filter = paste("padj < 0.05, NES > 0, mean WT change > 0,",
    "mean KO change < mean WT change"),
  all_go_terms_tested_in_parent = TRUE, inventory = inventory,
  reported_terms = nrow(top),
  biological_replicates_per_genotype = 1L)
jsonlite::write_json(diagnostics, file.path(stage_dir, "diagnostics.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
message("Completed focused top-GO report at: ", output_dir)
