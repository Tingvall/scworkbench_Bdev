#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 24_ebf1_il7_responsive_enhancers.R GLOBAL_REPORT OUTPUT_DIR",
  "SCRIPT_VERSION EXPECTED_GLOBAL_COMPLETION_SHA256"
), call. = FALSE)

global_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
script_version <- args[[3L]]
expected_completion_sha <- args[[4L]]
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
  library(jsonlite)
  library(digest)
})
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

completion_path <- file.path(global_dir, "completion.json")
stopifnot(identical(sha256_file(completion_path), expected_completion_sha))
all_da_path <- file.path(global_dir,
  "candidate_peaks_il7ko_vs_wt_differential_accessibility.tsv")
all_da <- fread(all_da_path)
validated <- all_da[globally_validated == TRUE]

ebf1 <- data.table(gene = "Ebf1", gene_id = "ENSMUSG00000057098",
  genome = "GRCm38/mm10", chromosome = "chr11", strand = "+",
  tss_1based = 44617317L, flank_bp = 3000000L)
ebf1[, `:=`(window_start_1based = tss_1based - flank_bp,
  window_end_1based = tss_1based + flank_bp)]
write_tsv(ebf1, file.path(stage_dir, "ebf1_tss_window_definition.tsv"))

local_all <- all_da[chromosome == ebf1$chromosome &
  peak_end_1based >= ebf1$window_start_1based &
  peak_start_1based <= ebf1$window_end_1based]
local_all[, `:=`(
  peak_midpoint_1based = (peak_start_1based + peak_end_1based) / 2,
  midpoint_relative_to_tss_bp = (peak_start_1based + peak_end_1based) / 2 -
    ebf1$tss_1based,
  interval_distance_to_tss_bp = fifelse(peak_start_1based <= ebf1$tss_1based &
    peak_end_1based >= ebf1$tss_1based, 0,
    fifelse(peak_end_1based < ebf1$tss_1based,
      peak_end_1based - ebf1$tss_1based,
      peak_start_1based - ebf1$tss_1based)))]
candidates <- local_all[globally_validated == TRUE]
candidates[, abs_interval_distance_to_tss_bp :=
  abs(interval_distance_to_tss_bp)]
setorder(candidates, p_val_adj, -standard_log2fc,
  -adjusted_spearman_rho, abs_interval_distance_to_tss_bp)
candidates[, abs_interval_distance_to_tss_bp := NULL]
candidates[, priority_rank := seq_len(.N)]
write_tsv(candidates, file.path(stage_dir,
  "ebf1_pm3mb_globally_validated_il7_responsive_enhancers.tsv"))
local_all[, abs_interval_distance_to_tss_bp :=
  abs(interval_distance_to_tss_bp)]
setorder(local_all, globally_validated, p_val_adj, -standard_log2fc,
  -adjusted_spearman_rho, abs_interval_distance_to_tss_bp)
local_all[, abs_interval_distance_to_tss_bp := NULL]
write_tsv(local_all, file.path(stage_dir,
  "ebf1_pm3mb_all_bdev_il7_responsive_peaks_with_query_da.tsv"))

if (nrow(local_all)) {
  pd <- copy(local_all)
  p <- ggplot(pd, aes(midpoint_relative_to_tss_bp / 1e6,
    adjusted_spearman_rho, size = abs(standard_log2fc),
    colour = -log10(pmax(p_val_adj, 1e-300)),
    shape = globally_validated)) + geom_point(alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey35") +
    scale_colour_viridis_c() +
    labs(x = "Peak midpoint relative to Ebf1 TSS (Mb)",
      y = "Bdev partial Spearman rho",
      size = "|WT/IL7KO log2FC|", colour = "-log10(query BH P)",
      shape = "Globally validated",
      title = "Ebf1 ±3 Mb IL7-responsive enhancer candidates") +
    theme_bw(base_size = 10)
  ggsave(file.path(stage_dir, "figures", "ebf1_candidate_locus_summary.png"),
    p, width = 9, height = 5.5, dpi = 180)
}

completion <- list(schema_version = "ebf1_il7_responsive_enhancers.v1",
  script_version = script_version,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  global_report = global_dir,
  global_completion_sha256 = expected_completion_sha,
  selection_rule = paste0("globally validated Bdev IL7-responsive and ",
    "WT-higher IL7KO-vs-WT BH<0.05 peak intersecting Ebf1 TSS ±3 Mb"),
  motif_used_for_selection = FALSE,
  ebf1 = as.list(ebf1[1]),
  counts = list(global_validated = nrow(validated),
    ebf1_pm3mb_bdev_responsive = nrow(local_all),
    ebf1_pm3mb_candidates = nrow(candidates)))
write_json(completion, file.path(stage_dir, "completion.json"))
writeLines(c("# Ebf1 IL7-responsive enhancer prioritization", "",
  paste0("This report intersects the frozen global validated set with Ebf1 ",
    "TSS ±3 Mb (chr11:", format(ebf1$window_start_1based, big.mark = ","),
    "–", format(ebf1$window_end_1based, big.mark = ","), "; GRCm38/mm10)."),
  "", paste0("Bdev-responsive peaks in the window: ", nrow(local_all)),
  paste0("Fully validated candidates retained: ", nrow(candidates)), "",
  "No EBF1 motif or any other motif result was used for selection or ranking.",
  "Global motif analysis is available only in the upstream global report.",
  "", "Primary table: `ebf1_pm3mb_globally_validated_il7_responsive_enhancers.tsv`",
  "Near-miss audit: `ebf1_pm3mb_all_bdev_il7_responsive_peaks_with_query_da.tsv`"),
  file.path(stage_dir, "REPORT.md"))
file.rename(stage_dir, output_dir)
cat("Completed Ebf1 follow-up:", output_dir, "\n")
