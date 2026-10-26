#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 33_plot_flt3l_peak_linked_gene_volcano.R LOCAL_DEG_REPORT",
  "OUTPUT_DIR EXPECTED_LOCAL_COMPLETION_SHA256 SCRIPT_VERSION LABEL_GENES"
), call. = FALSE)

local_report <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
expected_local_sha <- args[[3L]]
script_version <- args[[4L]]
label_genes <- as.integer(args[[5L]])
stopifnot(is.finite(label_genes), label_genes > 0L)

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(digest)
})
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

completion_path <- file.path(local_report, "completion.json")
stopifnot(identical(sha256_file(completion_path), expected_local_sha))
gene_path <- file.path(local_report,
  "all_unique_peak_linked_significant_genes.tsv")
genes <- fread(gene_path)
stopifnot(nrow(genes) > 0L, !anyDuplicated(genes$gene))
genes[, abs_log2fc := abs(log2fc_flt3lko_over_wt)]
setorder(genes, p_adj_bh, p_value, -abs_log2fc, gene)
genes[, significance_rank := seq_len(.N)]
genes[, neg_log10_padj := -log10(pmax(p_adj_bh, 1e-300))]
label_table <- copy(genes)
setorder(label_table, -abs_log2fc, p_adj_bh, p_value, gene)
label_table <- head(label_table, label_genes)
label_table[, fold_change_label_rank := seq_len(.N)]
genes[, fold_change_label_rank := NA_integer_]
genes[label_table, on = "gene",
  fold_change_label_rank := i.fold_change_label_rank]
genes[, label := fifelse(!is.na(fold_change_label_rank), gene, NA_character_)]
stopifnot(all(genes$p_adj_bh < 0.05),
  all(genes$direction == fifelse(genes$log2fc_flt3lko_over_wt > 0,
    "FLT3LKO_higher", "WT_higher")))
write_tsv(genes, file.path(stage_dir,
  "all_unique_peak_linked_significant_genes_volcano_data.tsv"))
write_tsv(label_table, file.path(stage_dir,
  paste0("top", label_genes,
    "_labeled_genes_by_absolute_log2fc.tsv")))

genes[, direction := factor(direction,
  levels = c("WT_higher", "FLT3LKO_higher"))]
fold_line <- log2(1.5)
p <- ggplot(genes, aes(log2fc_flt3lko_over_wt, neg_log10_padj,
  colour = direction)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2,
    colour = "grey55", linewidth = 0.45) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.45) +
  geom_vline(xintercept = c(-fold_line, fold_line), linetype = 3,
    colour = "grey55", linewidth = 0.45) +
  geom_point(alpha = 0.78, size = 2) +
  geom_text_repel(data = genes[!is.na(label)], aes(label = label),
    seed = 1729, size = 3.1, box.padding = 0.35, point.padding = 0.18,
    min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = c(WT_higher = "#D55E00",
    FLT3LKO_higher = "#0072B2"), labels = c(WT_higher = "WT higher",
      FLT3LKO_higher = "FLT3LKO higher")) +
  labs(x = "log2FC (FLT3LKO / WT; mean CP10K + 1)",
    y = "-log10(BH-adjusted P)", colour = "Direction",
    title = "FLT3L-responsive peak-linked genes",
    subtitle = paste0(nrow(genes),
      " unique significant genes within ±3 Mb of 229 validated peaks; top ",
      label_genes,
      " absolute fold changes labeled"),
    caption = paste0("All displayed genes already pass genome-wide BH P < 0.05. ",
      "Dotted vertical lines: |log2FC| = log2(1.5).")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank())
ggsave(file.path(stage_dir, "figures",
  "all_peak_linked_gene_volcano.png"), p, width = 11, height = 8,
  dpi = 240)
ggsave(file.path(stage_dir, "figures",
  "all_peak_linked_gene_volcano.pdf"), p, width = 11, height = 8,
  useDingbats = FALSE)

completion <- list(schema_version = "flt3l_peak_linked_gene_volcano.v1",
  script_version = script_version,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  inputs = list(local_deg_report = local_report,
    local_deg_completion_sha256 = expected_local_sha,
    unique_linked_genes = gene_path,
    unique_linked_genes_sha256 = sha256_file(gene_path)),
  plot = list(unique_genes = nrow(genes), labels = label_genes,
    label_rule = "largest absolute FLT3LKO/WT log2FC; ties by BH-adjusted P, raw P, gene",
    x = "FLT3LKO/WT standard log2FC", y = "-log10(BH-adjusted P)",
    colours = list(WT_higher = "#D55E00", FLT3LKO_higher = "#0072B2"),
    fold_reference = fold_line,
    significance_context = "input set prefiltered to genome-wide BH-adjusted P < 0.05",
    wt_higher = genes[direction == "WT_higher", .N],
    flt3lko_higher = genes[direction == "FLT3LKO_higher", .N]),
  inference = "descriptive cell-level; one biological sample per genotype")
write_json(completion, file.path(stage_dir, "completion.json"))
writeLines(c("# Volcano plot of FLT3L-responsive peak-linked genes", "",
  paste0("The volcano contains all ", nrow(genes),
    " unique significant genes linked to the 229 peak-centered ±3 Mb windows."),
  paste0("The ", label_genes,
    " unique genes with the largest absolute FLT3LKO/WT log2FC are labeled."), "",
  "The x-axis is FLT3LKO/WT standard log2FC: positive is FLT3LKO-higher and negative is WT-higher.",
  "WT-higher genes are orange; FLT3LKO-higher genes are blue.",
  "All displayed genes were prefiltered to genome-wide BH-adjusted P < 0.05.",
  "", "Outputs:", "",
  "- `figures/all_peak_linked_gene_volcano.png`",
  "- `figures/all_peak_linked_gene_volcano.pdf`",
  "- `all_unique_peak_linked_significant_genes_volcano_data.tsv`",
  paste0("- `top", label_genes,
    "_labeled_genes_by_absolute_log2fc.tsv`")),
  file.path(stage_dir, "REPORT.md"))
file.rename(stage_dir, output_dir)
cat("Completed linked-gene volcano:", output_dir, "\n")
