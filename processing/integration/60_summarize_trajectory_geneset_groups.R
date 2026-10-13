#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L)
  stop("Usage: 60_summarize_trajectory_geneset_groups.R <report_root>")
report_root <- normalizePath(args[[1L]], mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

groups <- c("cytokines", "bp", "mf", "cc")
group_labels <- c(cytokines = "Cytokines", bp = "GO BP",
  mf = "GO MF", cc = "GO CC")
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = TRUE)
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE,
  na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt"); on.exit(close(con), add = TRUE)
  write.table(as.data.frame(x), con, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
}

completion <- setNames(lapply(groups, function(g) {
  x <- read_json(file.path(report_root, g, "completion.json"))
  stopifnot(identical(x$status, "complete"),
    identical(x$analysis_group, g), isTRUE(x$bundle_unchanged),
    grepl(paste0("within ", g, " candidate peaks"), x$multiple_testing,
      fixed = TRUE))
  x
}), groups)

region_list <- vector("list", length(groups)); names(region_list) <- groups
set_list <- vector("list", length(groups)); names(set_list) <- groups
top_list <- vector("list", length(groups)); names(top_list) <- groups
for (g in groups) {
  regions <- fread(file.path(report_root, g, "significant_positive_regions.tsv"))
  figures <- fread(file.path(report_root, g, "umap_figure_index.tsv"))
  stopifnot(nrow(regions) == completion[[g]]$significant_positive_regions,
    uniqueN(regions$feature_set) == completion[[g]]$sets_with_significant_positive_regions,
    nrow(figures) == min(20L, uniqueN(regions$feature_set)),
    all(file.exists(file.path(report_root, g, figures$figure))))
  setorder(regions, q_region, p_region, -peak_z, feature_set, region_left)
  regions[, region_rank_within_set := seq_len(.N), by = feature_set]
  best <- regions[region_rank_within_set == 1L,
    .(analysis_group, feature_set, set_name, description,
      best_region_id = region_id, best_q = q_region, best_p = p_region,
      best_peak_z = peak_z, best_peak_score = peak_score,
      best_peak_pseudotime = peak_pseudotime,
      best_region_left = region_left, best_region_right = region_right)]
  aggregate <- regions[, .(
    n_significant_regions = .N,
    region_ids = paste(region_id, collapse = ";"),
    peak_pseudotimes = paste(signif(peak_pseudotime, 7), collapse = ";"),
    region_ranges = paste0(signif(region_left, 7), "-",
      signif(region_right, 7), collapse = ";")), by = feature_set]
  sets <- merge(best, aggregate, by = "feature_set", sort = FALSE)
  setorder(sets, best_q, best_p, -best_peak_z, feature_set)
  sets[, significance_rank_within_group := seq_len(.N)]
  setcolorder(sets, c("analysis_group", "significance_rank_within_group",
    "feature_set", "set_name", "description", "n_significant_regions",
    "best_region_id", "best_q", "best_p", "best_peak_z",
    "best_peak_score", "best_peak_pseudotime", "best_region_left",
    "best_region_right", "region_ids", "peak_pseudotimes", "region_ranges"))
  top <- sets[significance_rank_within_group <= 20L]
  top <- merge(top, figures[, .(feature_set, figure)], by = "feature_set",
    all.x = TRUE, sort = FALSE)
  setorder(top, significance_rank_within_group)
  stopifnot(!anyNA(top$figure), identical(top$feature_set, figures$feature_set))
  regions[, region_rank_within_set := NULL]
  region_list[[g]] <- regions
  set_list[[g]] <- sets
  top_list[[g]] <- top
}

all_regions <- rbindlist(region_list, use.names = TRUE, fill = TRUE)
all_sets <- rbindlist(set_list, use.names = TRUE, fill = TRUE)
top_sets <- rbindlist(top_list, use.names = TRUE, fill = TRUE)
write_tsv_gz(all_regions,
  file.path(report_root, "all_significant_positive_regions.tsv.gz"))
write_tsv(all_sets,
  file.path(report_root, "significant_gene_sets_by_group.tsv"))
write_tsv(top_sets,
  file.path(report_root, "top20_significant_gene_sets_by_group.tsv"))

summary_table <- rbindlist(lapply(groups, function(g) data.table(
  analysis_group = g, label = unname(group_labels[[g]]),
  tested_sets = completion[[g]]$tested_sets,
  zero_variance_untestable_sets = completion[[g]]$zero_variance_untestable_sets,
  significant_sets = completion[[g]]$sets_with_significant_positive_regions,
  significant_regions = completion[[g]]$significant_positive_regions,
  figures = completion[[g]]$figures_generated_for_top_significant_sets)))
write_tsv(summary_table, file.path(report_root, "group_summary.tsv"))

table_lines <- c("| Group | Tested | Zero-variance | Significant sets | Final regions | Figures |",
  "|---|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(summary_table)), function(i) paste0("| ",
    summary_table$label[[i]], " | ", summary_table$tested_sets[[i]], " | ",
    summary_table$zero_variance_untestable_sets[[i]], " | ",
    summary_table$significant_sets[[i]], " | ",
    summary_table$significant_regions[[i]], " | ",
    summary_table$figures[[i]], " |"), character(1)))

top_lines <- unlist(lapply(groups, function(g) {
  x <- top_list[[g]]
  c(paste0("## ", group_labels[[g]], ": top 20 significant gene sets"), "",
    "| Rank | Gene set | Regions | Best peak pt | Best z | Best q |",
    "|---:|---|---:|---:|---:|---:|",
    vapply(seq_len(nrow(x)), function(i) paste0("| ",
      x$significance_rank_within_group[[i]], " | ", x$set_name[[i]], " | ",
      x$n_significant_regions[[i]], " | ",
      signif(x$best_peak_pseudotime[[i]], 5), " | ",
      signif(x$best_peak_z[[i]], 5), " | ", signif(x$best_q[[i]], 5),
      " |"), character(1)), "")
}), use.names = FALSE)

writeLines(c("# Positive trajectory gene-set enrichment", "",
  "## Summary", "", table_lines, "",
  "Multiple testing used Benjamini-Hochberg correction separately within cytokines, GO BP, GO MF, and GO CC. Significance required q < 0.01; reported peaks additionally required at least 10% relative prominence. Retained peaks were merged only when the uncovered gap was strictly less than one 5% pseudotime window.",
  "",
  "The complete set-level results are in `significant_gene_sets_by_group.tsv`; all final regions are in `all_significant_positive_regions.tsv.gz`. Figures are limited to the 20 most significant sets per group, ranked by best q, then p, then peak z.",
  "", top_lines), file.path(report_root, "REPORT.md"))

fingerprints <- unique(vapply(completion, function(x)
  x$bundle_fingerprint_after, character(1)))
stopifnot(length(fingerprints) == 1L)
jsonlite::write_json(list(status = "complete", groups = as.list(groups),
  total_tested_sets = sum(summary_table$tested_sets),
  total_significant_sets = sum(summary_table$significant_sets),
  total_significant_regions = sum(summary_table$significant_regions),
  multiple_testing = "BH separately within each analysis group",
  q_cutoff = 0.01, relative_prominence_cutoff = 0.10,
  figures_per_group_maximum = 20,
  bundle_fingerprint = fingerprints[[1L]], bundle_unchanged = TRUE),
  file.path(report_root, "completion.json"), pretty = TRUE, auto_unbox = TRUE)
