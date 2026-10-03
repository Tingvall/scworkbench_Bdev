#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 25_annotate_il7_validated_peaks.R GLOBAL_REPORT BDEV_RDATA",
  "OUTPUT_DIR EXPECTED_GLOBAL_COMPLETION_SHA256 SCRIPT_VERSION PROMOTER_BP"
), call. = FALSE)

global_dir <- normalizePath(args[[1L]], mustWork = TRUE)
bdev_rdata <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
expected_global_sha <- args[[4L]]
script_version <- args[[5L]]
promoter_bp <- as.integer(args[[6L]])
stopifnot(is.finite(promoter_bp), promoter_bp > 0L)

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

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

global_completion <- file.path(global_dir, "completion.json")
stopifnot(identical(sha256_file(global_completion), expected_global_sha))
peaks <- fread(file.path(global_dir,
  "globally_validated_il7_responsive_peaks.tsv"))
stopifnot(nrow(peaks) == 51L, !anyDuplicated(peaks$feature_id),
  all(peaks$globally_validated == TRUE))

source_env <- new.env(parent = emptyenv())
con <- gzfile(bdev_rdata, "rb")
loaded <- load(con, envir = source_env)
close(con)
stopifnot("Bdev" %in% loaded)
annotation <- as.data.table(as.data.frame(
  Signac::Annotation(source_env$Bdev[["peaks"]])))
required <- c("seqnames", "start", "end", "strand", "tx_id",
  "gene_name", "gene_id", "gene_biotype")
stopifnot(all(required %in% names(annotation)))
annotation <- annotation[!is.na(tx_id) & tx_id != "" &
  !is.na(gene_id) & gene_id != "" & !is.na(gene_name) & gene_name != "" &
  strand %in% c("+", "-")]
transcripts <- annotation[, .(
  transcript_start_1based = min(start),
  transcript_end_1based = max(end)),
  by = .(chromosome = as.character(seqnames), strand = as.character(strand),
    transcript_id = tx_id, gene_id, gene_name, gene_biotype)]
transcripts[, transcript_tss_1based := fifelse(strand == "+",
  transcript_start_1based, transcript_end_1based)]
stopifnot(nrow(transcripts) > 10000L,
  !anyDuplicated(transcripts$transcript_id))

annotation_rows <- vector("list", nrow(peaks))
for (i in seq_len(nrow(peaks))) {
  p <- peaks[i]
  tx <- transcripts[chromosome == p$chromosome]
  stopifnot(nrow(tx) > 0L)
  tx[, interval_distance_to_peak_bp := fifelse(
    transcript_tss_1based < p$peak_start_1based,
    p$peak_start_1based - transcript_tss_1based,
    fifelse(transcript_tss_1based > p$peak_end_1based,
      transcript_tss_1based - p$peak_end_1based, 0))]
  tx[, midpoint_distance_to_tss_bp := abs(
    (p$peak_start_1based + p$peak_end_1based) / 2 -
      transcript_tss_1based)]
  setorder(tx, interval_distance_to_peak_bp, midpoint_distance_to_tss_bp,
    gene_name, transcript_id)
  nearest <- tx[1L]
  promoter_tx <- tx[interval_distance_to_peak_bp <= promoter_bp]
  promoter_genes <- paste(sort(unique(promoter_tx$gene_name)), collapse = ";")
  signed_distance <- ((p$peak_start_1based + p$peak_end_1based) / 2 -
    nearest$transcript_tss_1based) * ifelse(nearest$strand == "+", 1, -1)
  annotation_rows[[i]] <- data.table(feature_id = p$feature_id,
    closest_gene = nearest$gene_name,
    closest_gene_id = nearest$gene_id,
    closest_gene_biotype = nearest$gene_biotype,
    closest_transcript_id = nearest$transcript_id,
    closest_gene_strand = nearest$strand,
    closest_tss_1based = nearest$transcript_tss_1based,
    distance_peak_to_closest_tss_bp = nearest$interval_distance_to_peak_bp,
    signed_midpoint_distance_from_tss_bp = signed_distance,
    promoter_overlap_genes = ifelse(nrow(promoter_tx), promoter_genes, NA_character_),
    regulatory_class = ifelse(nrow(promoter_tx), "promoter",
      "putative_enhancer"))
}
peak_annotation <- rbindlist(annotation_rows)
annotated <- merge(peaks, peak_annotation, by = "feature_id", all.x = TRUE)
setorder(annotated, p_val_adj, p_val, -standard_log2fc, feature_id)
annotated[, significance_rank := seq_len(.N)]
setcolorder(annotated, c("significance_rank", "feature_id", "chromosome",
  "peak_start_1based", "peak_end_1based", "closest_gene",
  "closest_gene_id", "closest_gene_biotype", "closest_transcript_id",
  "closest_gene_strand", "closest_tss_1based",
  "distance_peak_to_closest_tss_bp", "signed_midpoint_distance_from_tss_bp",
  "regulatory_class", "promoter_overlap_genes", "standard_log2fc",
  "p_val", "p_val_adj", setdiff(names(annotated), c("significance_rank",
    "feature_id", "chromosome", "peak_start_1based", "peak_end_1based",
    "closest_gene", "closest_gene_id", "closest_gene_biotype",
    "closest_transcript_id", "closest_gene_strand", "closest_tss_1based",
    "distance_peak_to_closest_tss_bp", "signed_midpoint_distance_from_tss_bp",
    "regulatory_class", "promoter_overlap_genes", "standard_log2fc",
    "p_val", "p_val_adj"))))
stopifnot(nrow(annotated) == 51L, identical(annotated$significance_rank,
  seq_len(51L)), !anyNA(annotated$closest_gene),
  all(annotated$regulatory_class %in% c("promoter", "putative_enhancer")))
write_tsv(annotated, file.path(stage_dir,
  "top51_globally_validated_peaks_closest_gene_annotation.tsv"))
write_tsv(transcripts, file.path(stage_dir,
  "bdev_transcript_tss_annotation_used.tsv.gz"))

compact <- annotated[, .(significance_rank, feature_id, closest_gene,
  distance_peak_to_closest_tss_bp, regulatory_class, standard_log2fc,
  p_val_adj)]
write_tsv(compact, file.path(stage_dir,
  "top51_peak_gene_classification_compact.tsv"))
class_counts <- annotated[, .N, by = regulatory_class][order(regulatory_class)]
write_tsv(class_counts, file.path(stage_dir,
  "regulatory_classification_summary.tsv"))

completion <- list(schema_version = "il7_validated_peak_annotation.v1",
  script_version = script_version,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  global_report = global_dir,
  global_completion_sha256 = expected_global_sha,
  bdev_rdata = bdev_rdata,
  bdev_rdata_sha256 = sha256_file(bdev_rdata),
  annotation_source = "Signac::Annotation(Bdev[['peaks']]); GRCm38/mm10",
  closest_gene_rule = paste0("minimum interval distance from peak to any ",
    "annotated transcript TSS; ties by midpoint distance, gene name, transcript ID"),
  promoter_rule = paste0("peak interval within ", promoter_bp,
    " bp of any annotated transcript TSS"),
  enhancer_rule = paste0("no annotated transcript TSS within ", promoter_bp,
    " bp; labeled putative_enhancer"),
  ordering = "query IL7 WT-versus-IL7KO p_val_adj, then p_val",
  counts = list(peaks = nrow(annotated),
    promoters = annotated[regulatory_class == "promoter", .N],
    putative_enhancers = annotated[regulatory_class == "putative_enhancer", .N]))
write_json(completion, file.path(stage_dir, "completion.json"))
writeLines(c("# Annotation of the 51 globally validated IL7-responsive peaks", "",
  paste0("Peaks are ordered by query IL7 WT-versus-IL7KO BH-adjusted P value."),
  "Closest gene means closest annotated transcript TSS in the exact Bdev GRCm38 annotation.",
  paste0("Promoter: peak within ", promoter_bp, " bp of any transcript TSS. ",
    "All other peaks are labeled putative enhancer; this is a positional, not functional, classification."),
  "", paste0("Promoters: ", completion$counts$promoters),
  paste0("Putative enhancers: ", completion$counts$putative_enhancers), "",
  "Primary compact table: `top51_peak_gene_classification_compact.tsv`",
  "Full audit table: `top51_globally_validated_peaks_closest_gene_annotation.tsv`"),
  file.path(stage_dir, "REPORT.md"))
file.rename(stage_dir, output_dir)
cat("Completed peak annotation:", output_dir, "\n")
