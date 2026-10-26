#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 68_select_top_go_ucell_dynamic_sets.R",
  "<go_dynamic_report> <output_dir> <top_n>"))
report_root <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
top_n <- as.integer(args[[3L]])
stopifnot(is.finite(top_n), top_n > 0L)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

completion <- read_json(file.path(report_root, "completion.json"),
  simplifyVector = TRUE)
stopifnot(identical(completion$status, "complete"),
  identical(completion$inference_mode, "dynamic_range"),
  completion$permutations == 1000L,
  isTRUE(completion$groupwise_multiple_testing))

groups <- c("bp", "mf", "cc")
selected <- lapply(groups, function(g) {
  x <- fread(file.path(report_root, g, "significant_gene_sets.tsv"))
  stopifnot(nrow(x) >= top_n, all(x$q_dynamic < 0.01))
  x[, analysis_group := g]
  setorder(x, q_dynamic, p_dynamic, -observed_dynamic_range, feature_set)
  x <- x[seq_len(top_n)]
  x[, rank_within_group := seq_len(.N)]
  setcolorder(x, c("analysis_group", "rank_within_group", "feature_set",
    "set_name", "q_dynamic", "p_dynamic", "observed_dynamic_range",
    "n_regions", "best_region_fwer_p", "strongest_ucell_effect",
    "peak_pseudotimes"))
  fwrite(x, file.path(stage_dir, paste0(g, "_top", top_n,
    "_dynamic_gene_sets.tsv")), sep = "\t", quote = FALSE, na = "NA")
  x
})
selected_all <- rbindlist(selected, use.names = TRUE, fill = TRUE)
fwrite(selected_all, file.path(stage_dir, paste0("all_groups_top", top_n,
  "_dynamic_gene_sets.tsv")), sep = "\t", quote = FALSE, na = "NA")

fmt <- function(x) vapply(seq_len(nrow(x)), function(i) paste0("| ",
  x$rank_within_group[[i]], " | `", x$set_name[[i]], "` | ",
  signif(x$q_dynamic[[i]], 4), " | ", signif(x$p_dynamic[[i]], 4),
  " | ", signif(x$observed_dynamic_range[[i]], 5), " | ",
  x$n_regions[[i]], " |"), character(1))
report <- c("# Top dynamic GO gene sets", "",
  paste0("Top ", top_n,
    " dynamically significant gene sets were selected independently within BP, MF, and CC."),
  "Ranking: groupwise dynamic q, then permutation p, then descending native-UCell peak-to-trough range, then stable feature-set identifier.",
  "The UCell range is used only to order terms tied by permutation significance; it was not an inference cutoff.",
  "Positive-region count is reported but is not part of the ranking.", "")
for (g in groups) {
  x <- selected_all[analysis_group == g]
  report <- c(report, paste0("## ", toupper(g)), "",
    "| Rank | Gene set | Dynamic q | Dynamic p | UCell range | Positive regions |",
    "|---:|---|---:|---:|---:|---:|", fmt(x), "")
}
writeLines(report, file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", top_n_per_group = top_n,
  groups = groups, total_selected = nrow(selected_all),
  ranking = paste("q_dynamic, p_dynamic, descending observed native-UCell",
    "peak-to-trough range, feature_set"),
  region_metrics_used_for_ranking = FALSE,
  source_report = report_root), file.path(stage_dir, "completion.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)

stopifnot(file.rename(stage_dir, output_dir))
message("Completed top dynamic GO selection at: ", output_dir)
