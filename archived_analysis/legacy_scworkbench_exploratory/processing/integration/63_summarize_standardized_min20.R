#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 63_summarize_standardized_min20.R <output_root>")
root <- normalizePath(args[[1L]], mustWork = TRUE)
suppressPackageStartupMessages(library(data.table))
groups <- c("cytokines", "bp", "mf", "cc")
stopifnot(all(file.exists(file.path(root, groups, "completion.json"))))
counts <- rbindlist(lapply(groups, function(g)
  fread(file.path(root, g, "filter_counts.tsv"))))
counts[, group_order := match(analysis_group, groups)]
setorder(counts, group_order)
counts[, group_order := NULL]
fwrite(counts, file.path(root, "group_filter_counts.tsv"), sep = "\t",
  quote = FALSE)
regions <- rbindlist(lapply(groups, function(g)
  fread(file.path(root, g, "significant_positive_regions.tsv"))), fill = TRUE)
fwrite(regions, file.path(root, "all_significant_positive_regions.tsv.gz"),
  sep = "\t", quote = FALSE, compress = "gzip")
lines <- c("# Standardized-gene min20 trajectory enrichment", "",
  "Per-gene depth adjustment and standardization before gene-set averaging; 1,000 permutations; BH separately within each group.", "",
  "Positive peaks were discovered before downstream filtering. Retained peaks require q<0.01, >=10% relative prominence, at least 3 positively contributing genes and at least 20% of all mapped genes, and no gene contributing >50% of total absolute contribution. Retained peaks separated by an uncovered gap strictly less than one 5% window were merged.", "",
  "| Group | Tested sets | Significant peaks | After prominence | After multi-gene support | After dominance | Final sets | Final regions |",
  "|---|---:|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(counts)), function(i) paste0("| ",
    counts$analysis_group[[i]], " | ", counts$tested_gene_sets[[i]], " | ",
    counts$significant_peaks[[i]], " | ",
    counts$prominence_retained_peaks[[i]], " | ",
    counts$multi_gene_support_retained_peaks[[i]], " | ",
    counts$dominance_retained_constituent_peaks[[i]], " | ",
    counts$final_gene_sets[[i]], " | ", counts$final_regions[[i]], " |"),
    character(1)), "",
  paste0("Total final regions: ", sum(counts$final_regions),
    " across ", sum(counts$final_gene_sets),
    " group-specific gene-set entries."))
writeLines(lines, file.path(root, "REPORT.md"))
