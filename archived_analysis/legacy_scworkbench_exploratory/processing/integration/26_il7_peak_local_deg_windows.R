#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "Usage: 26_il7_peak_local_deg_windows.R PEAK_ANNOTATION_REPORT",
  "DEG_REPORT BDEV_RDATA OUTPUT_DIR EXPECTED_PEAK_COMPLETION_SHA256",
  "EXPECTED_DEG_COMPLETION_SHA256 SCRIPT_VERSION WINDOW_BP PADJ_CUTOFF"
), call. = FALSE)

peak_report <- normalizePath(args[[1L]], mustWork = TRUE)
deg_report <- normalizePath(args[[2L]], mustWork = TRUE)
bdev_rdata <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_peak_sha <- args[[5L]]
expected_deg_sha <- args[[6L]]
script_version <- args[[7L]]
window_bp <- as.integer(args[[8L]])
padj_cutoff <- as.numeric(args[[9L]])
stopifnot(is.finite(window_bp), window_bp > 0L,
  is.finite(padj_cutoff), padj_cutoff > 0, padj_cutoff < 1)

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "per_peak"), recursive = TRUE,
  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
  library(SeuratObject)
  library(Signac)
})
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

peak_completion <- file.path(peak_report, "completion.json")
deg_completion <- file.path(deg_report, "completion.json")
stopifnot(identical(sha256_file(peak_completion), expected_peak_sha),
  identical(sha256_file(deg_completion), expected_deg_sha))

peaks <- fread(file.path(peak_report,
  "top51_globally_validated_peaks_closest_gene_annotation.tsv"))
stopifnot(nrow(peaks) == 51L,
  identical(peaks$significance_rank, seq_len(51L)),
  !anyDuplicated(peaks$feature_id))
peaks[, peak_center_1based :=
  (peak_start_1based + peak_end_1based) / 2]
peaks[, `:=`(
  window_start_1based = pmax(1, floor(peak_center_1based - window_bp)),
  window_end_1based = ceiling(peak_center_1based + window_bp))]

deg_path <- file.path(deg_report,
  "il7_shared_peak2_all_genes_normalized_counts_and_de.tsv")
deg_all <- fread(deg_path)
stopifnot(nrow(deg_all) == 33696L, !anyDuplicated(deg_all$feature_id),
  all(c("standard_log2fc", "p_val", "p_val_adj", "significant_bh05") %in%
    names(deg_all)))
deg <- deg_all[is.finite(p_val_adj) & p_val_adj < padj_cutoff &
  is.finite(standard_log2fc)]
stopifnot(nrow(deg) == sum(deg_all$significant_bh05 == TRUE, na.rm = TRUE))
deg <- deg[, .(gene = feature_id, mean_cp10k_wt, mean_cp10k_il7ko =
  mean_cp10k_ko, standard_log2fc_wt_over_il7ko = standard_log2fc,
  log2fc_il7ko_over_wt = -standard_log2fc, p_value = p_val,
  p_adj_bh = p_val_adj)]
deg[, direction := fifelse(log2fc_il7ko_over_wt > 0, "IL7KO_higher",
  fifelse(log2fc_il7ko_over_wt < 0, "WT_higher", "equal"))]

source_env <- new.env(parent = emptyenv())
con <- gzfile(bdev_rdata, "rb")
loaded <- load(con, envir = source_env)
close(con)
stopifnot("Bdev" %in% loaded)
annotation <- as.data.table(as.data.frame(
  Signac::Annotation(source_env$Bdev[["peaks"]])))
annotation <- annotation[!is.na(tx_id) & tx_id != "" &
  !is.na(gene_id) & gene_id != "" & !is.na(gene_name) & gene_name != "" &
  strand %in% c("+", "-")]
transcripts <- annotation[, .(
  transcript_start_1based = min(start),
  transcript_end_1based = max(end)),
  by = .(chromosome = as.character(seqnames), strand = as.character(strand),
    transcript_id = tx_id, gene_id, gene = gene_name,
    gene_biotype)]
transcripts[, gene_tss_1based := fifelse(strand == "+",
  transcript_start_1based, transcript_end_1based)]
rm(annotation, source_env); invisible(gc())

deg_annotated <- merge(deg, transcripts, by = "gene", allow.cartesian = TRUE)
mapped_genes <- intersect(deg$gene, transcripts$gene)
unmapped <- deg[!gene %in% mapped_genes]
write_tsv(unmapped, file.path(stage_dir,
  "significant_de_genes_without_bdev_tss_annotation.tsv"))

pair_tables <- vector("list", nrow(peaks))
for (i in seq_len(nrow(peaks))) {
  p <- peaks[i]
  local <- deg_annotated[chromosome == p$chromosome &
    gene_tss_1based >= p$window_start_1based &
    gene_tss_1based <= p$window_end_1based]
  if (nrow(local)) {
    local[, signed_tss_distance_from_peak_center_bp :=
      gene_tss_1based - p$peak_center_1based]
    local[, absolute_tss_distance_from_peak_center_bp :=
      abs(signed_tss_distance_from_peak_center_bp)]
    setorder(local, gene, absolute_tss_distance_from_peak_center_bp,
      transcript_id)
    local <- local[, .SD[1L], by = gene]
    local[, `:=`(peak_rank = p$significance_rank,
      peak_id = p$feature_id, peak_chromosome = p$chromosome,
      peak_start_1based = p$peak_start_1based,
      peak_end_1based = p$peak_end_1based,
      peak_center_1based = p$peak_center_1based,
      peak_closest_gene = p$closest_gene,
      peak_regulatory_class = p$regulatory_class,
      window_bp_each_side = window_bp,
      window_start_1based = p$window_start_1based,
      window_end_1based = p$window_end_1based)]
    local[, abs_log2fc_il7ko_over_wt := abs(log2fc_il7ko_over_wt)]
    setorder(local, p_adj_bh, -abs_log2fc_il7ko_over_wt,
      absolute_tss_distance_from_peak_center_bp, gene)
    local[, within_peak_rank := seq_len(.N)]
    local[, abs_log2fc_il7ko_over_wt := NULL]
  } else {
    local <- data.table()
  }
  pair_tables[[i]] <- local
}
all_pairs <- rbindlist(pair_tables, fill = TRUE)
stopifnot(nrow(all_pairs) > 0L,
  all(all_pairs$p_adj_bh < padj_cutoff),
  all(abs(all_pairs$signed_tss_distance_from_peak_center_bp) <= window_bp),
  all(unique(all_pairs$peak_id) %in% peaks$feature_id))
setcolorder(all_pairs, c("peak_rank", "peak_id", "peak_chromosome",
  "peak_start_1based", "peak_end_1based", "peak_center_1based",
  "peak_closest_gene", "peak_regulatory_class", "window_bp_each_side",
  "window_start_1based", "window_end_1based", "within_peak_rank", "gene",
  "gene_id", "gene_biotype", "transcript_id", "strand",
  "gene_tss_1based", "signed_tss_distance_from_peak_center_bp",
  "absolute_tss_distance_from_peak_center_bp", "mean_cp10k_wt",
  "mean_cp10k_il7ko", "standard_log2fc_wt_over_il7ko",
  "log2fc_il7ko_over_wt", "direction", "p_value", "p_adj_bh"))
fwrite(all_pairs, file.path(stage_dir,
  "all_51_peaks_pm3mb_significant_il7_peak2_de_genes.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip")

summary_stats <- all_pairs[, .(
  significant_genes = uniqueN(gene),
  il7ko_higher = uniqueN(gene[direction == "IL7KO_higher"]),
  wt_higher = uniqueN(gene[direction == "WT_higher"]),
  minimum_p_adj_bh = min(p_adj_bh),
  maximum_absolute_log2fc_il7ko_over_wt = max(abs(log2fc_il7ko_over_wt))),
  by = .(peak_rank, peak_id)]
summary <- merge(peaks[, .(peak_rank = significance_rank,
  peak_id = feature_id, peak_chromosome = chromosome,
  peak_start_1based, peak_end_1based, peak_center_1based,
  peak_closest_gene = closest_gene,
  peak_regulatory_class = regulatory_class,
  window_start_1based, window_end_1based)], summary_stats,
  by = c("peak_rank", "peak_id"), all.x = TRUE)
for (nm in c("significant_genes", "il7ko_higher", "wt_higher"))
  set(summary, which(is.na(summary[[nm]])), nm, 0L)
setorder(summary, peak_rank)
write_tsv(summary, file.path(stage_dir, "peak_local_deg_summary.tsv"))

manifest <- vector("list", nrow(peaks))
for (i in seq_len(nrow(peaks))) {
  p <- peaks[i]
  x <- all_pairs[peak_id == p$feature_id]
  safe_id <- gsub("[^A-Za-z0-9._-]", "_", p$feature_id)
  rel <- file.path("per_peak", sprintf("%02d_%s_pm3mb_significant_genes.tsv",
    p$significance_rank, safe_id))
  write_tsv(x, file.path(stage_dir, rel))
  manifest[[i]] <- data.table(peak_rank = p$significance_rank,
    peak_id = p$feature_id, significant_genes = nrow(x), file = rel)
}
manifest <- rbindlist(manifest)
write_tsv(manifest, file.path(stage_dir, "per_peak_file_manifest.tsv"))

completion <- list(schema_version = "il7_peak_local_deg_windows.v1",
  script_version = script_version,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  inputs = list(peak_annotation_report = peak_report,
    peak_annotation_completion_sha256 = expected_peak_sha,
    deg_report = deg_report, deg_completion_sha256 = expected_deg_sha,
    deg_table = deg_path, deg_table_sha256 = sha256_file(deg_path),
    bdev_rdata = bdev_rdata, bdev_rdata_sha256 = sha256_file(bdev_rdata)),
  rules = list(window = paste0("transcript TSS within ±", window_bp,
    " bp of peak center"),
    significance = paste0("genome-wide DEG BH-adjusted P < ", padj_cutoff),
    comparison = "IL7KO versus WT in exact shared IL7 peak2 cells",
    reported_log2fc = "log2((mean linear CP10K IL7KO + 1)/(mean linear CP10K WT + 1))",
    source_log2fc = "log2((mean linear CP10K WT + 1)/(mean linear CP10K IL7KO + 1))",
    duplicate_policy = "gene retained for every peak window it falls within; nearest transcript TSS retained within each peak-gene pair"),
  counts = list(peaks = nrow(peaks), tested_genes = nrow(deg_all),
    significant_genes = nrow(deg),
    significant_genes_with_bdev_tss = length(mapped_genes),
    significant_genes_without_bdev_tss = nrow(unmapped),
    peak_gene_pairs = nrow(all_pairs),
    unique_local_significant_genes = uniqueN(all_pairs$gene),
    peaks_with_at_least_one_gene = uniqueN(all_pairs$peak_id),
    il7ko_higher_pairs = all_pairs[direction == "IL7KO_higher", .N],
    wt_higher_pairs = all_pairs[direction == "WT_higher", .N]),
  inference = "descriptive cell-level; one biological sample per genotype")
write_json(completion, file.path(stage_dir, "completion.json"))

report <- c("# Significant IL7 peak2 DE genes within ±3 Mb of 51 peaks", "",
  paste0("Each window is centered on the peak midpoint and includes genes whose ",
    "annotated transcript TSS lies within ±", format(window_bp, big.mark = ","),
    " bp. DEG significance is genome-wide BH-adjusted P < ", padj_cutoff, "."),
  "", paste0("Peak–gene pairs: ", format(nrow(all_pairs), big.mark = ",")),
  paste0("Unique local significant genes: ", uniqueN(all_pairs$gene)),
  paste0("IL7KO-higher pairs: ", all_pairs[direction == "IL7KO_higher", .N]),
  paste0("WT-higher pairs: ", all_pairs[direction == "WT_higher", .N]), "",
  "The primary log2FC is IL7KO/WT; positive values are IL7KO-higher. The source WT/IL7KO log2FC is retained in a separate column.",
  "Genes are repeated when they fall within more than one overlapping peak window.",
  "", "## Outputs", "",
  "- `all_51_peaks_pm3mb_significant_il7_peak2_de_genes.tsv.gz`: all peak–gene pairs",
  "- `peak_local_deg_summary.tsv`: counts by peak and direction",
  "- `per_peak/`: one complete table for each peak",
  "- `per_peak_file_manifest.tsv`: per-peak file index",
  "", "Query DE inference is descriptive at the cell level because there is one biological sample per genotype.")
writeLines(report, file.path(stage_dir, "REPORT.md"))
file.rename(stage_dir, output_dir)
cat("Completed local DEG windows:", output_dir, "\n")
