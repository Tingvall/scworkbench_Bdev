#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) stop(paste(
  "Usage: 34_ebf1_putative_enhancers_three_way.R BDEV_BUNDLE BDEV_CELLS",
  "BDEV_PEAK_UNIVERSE TCF3_DA PLACSEQ SESSION_ROOT OUTPUT_DIR",
  "EXPECTED_BDEV_FP SCRIPT_VERSION EBF1_TSS_1BASED SEED"
), call. = FALSE)

bdev_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
bdev_cells_path <- normalizePath(args[[2L]], mustWork = TRUE)
bdev_peaks_path <- normalizePath(args[[3L]], mustWork = TRUE)
tcf3_da_path <- normalizePath(args[[4L]], mustWork = TRUE)
placseq_path <- normalizePath(args[[5L]], mustWork = TRUE)
session_root <- normalizePath(args[[6L]], mustWork = FALSE)
output_dir <- normalizePath(args[[7L]], mustWork = FALSE)
expected_bdev_fp <- args[[8L]]
script_version <- args[[9L]]
ebf1_tss <- as.integer(args[[10L]])
seed <- as.integer(args[[11L]])
stopifnot(ebf1_tss == 44617317L, is.finite(seed))

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(session_root, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(digest)
  library(splines)
})

log_event <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(..., collapse = " "), "\n")
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
collapse_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x)]))
  if (length(x)) paste(x, collapse = ";") else NA_character_
}

flank_bp <- 3000000L
promoter_flank_bp <- 2000L
anchor_slop_bp <- 20000L
window_start <- ebf1_tss - flank_bp
window_end <- ebf1_tss + flank_bp
promoter_start <- ebf1_tss - promoter_flank_bp
promoter_end <- ebf1_tss + promoter_flank_bp

peaks <- fread(bdev_peaks_path)
stopifnot(nrow(peaks) == 259L, !anyDuplicated(peaks$bdev_feature_id),
  all(peaks$mm10_seqname == "chr11"),
  all(peaks$mm10_end_1based >= window_start &
    peaks$mm10_start_1based <= window_end))
setnames(peaks, c("mm10_seqname", "mm10_start_1based", "mm10_end_1based"),
  c("chromosome", "peak_start_1based", "peak_end_1based"))
peaks[, `:=`(
  overlaps_ebf1_promoter = peak_end_1based >= promoter_start &
    peak_start_1based <= promoter_end,
  distal_putative_enhancer = !(peak_end_1based >= promoter_start &
    peak_start_1based <= promoter_end)
)]

# Route 1: retain every local GSE259355-to-Bdev overlap, then apply the frozen
# source-peak test. Coordinates in the DESeq2 table are recorded as 1-based.
log_event("Route 1: mapping GSE259355 differential peaks to Bdev peaks")
tcf3 <- fread(tcf3_da_path)
stopifnot(all(c("feature_id", "chromosome", "start_1based", "end_1based",
  "log2FoldChange", "pvalue", "padj_bh") %in% names(tcf3)))
tcf3[, tcf3_source_pass := is.finite(padj_bh) & padj_bh < 0.05 &
  is.finite(log2FoldChange) & abs(log2FoldChange) > 1]
tcf3[, tcf3_source_direction := fifelse(log2FoldChange < 0,
  "lower_after_TCF3_depletion", "higher_after_TCF3_depletion")]
tcf3_local <- tcf3[chromosome == "chr11" & end_1based >= window_start &
  start_1based <= window_end]

peak_intervals <- peaks[, .(chromosome, start = peak_start_1based,
  end = peak_end_1based, bdev_feature_id)]
source_intervals <- tcf3_local[, .(chromosome, start = start_1based,
  end = end_1based, tcf3_feature_id = feature_id, log2FoldChange, lfcSE,
  stat, pvalue, padj_bh, tcf3_source_pass, tcf3_source_direction)]
setkey(peak_intervals, chromosome, start, end)
tcf3_map <- foverlaps(source_intervals, peak_intervals,
  by.x = c("chromosome", "start", "end"),
  by.y = c("chromosome", "start", "end"), type = "any", nomatch = 0L)
setnames(tcf3_map, c("i.start", "i.end"),
  c("tcf3_start_1based", "tcf3_end_1based"))
tcf3_map[, `:=`(
  bdev_start_1based = start,
  bdev_end_1based = end,
  overlap_bp = pmax(0L, pmin(end, tcf3_end_1based) -
    pmax(start, tcf3_start_1based) + 1L)
)]
tcf3_map[, `:=`(
  bdev_overlap_fraction = overlap_bp / (bdev_end_1based - bdev_start_1based + 1),
  tcf3_overlap_fraction = overlap_bp /
    (tcf3_end_1based - tcf3_start_1based + 1)
)]
tcf3_map[, abs_log2fc_sort := abs(log2FoldChange)]
setorder(tcf3_map, bdev_feature_id, -tcf3_source_pass, padj_bh,
  -abs_log2fc_sort, tcf3_feature_id)
tcf3_map[, abs_log2fc_sort := NULL]
write_tsv(tcf3_map, file.path(stage_dir,
  "route1_all_tcf3_to_bdev_peak_overlaps.tsv"))

tcf3_summary <- tcf3_map[, {
  best <- order(!tcf3_source_pass, padj_bh, -abs(log2FoldChange),
    na.last = TRUE)[1L]
  list(tcf3_overlapping_source_peaks = .N,
    tcf3_passing_source_peaks = sum(tcf3_source_pass),
    tcf3_route_pass = any(tcf3_source_pass),
    tcf3_best_source_peak = tcf3_feature_id[best],
    tcf3_best_log2fc = log2FoldChange[best],
    tcf3_best_pvalue = pvalue[best],
    tcf3_best_padj = padj_bh[best],
    tcf3_best_direction = tcf3_source_direction[best],
    tcf3_max_overlap_bp = max(overlap_bp),
    tcf3_all_passing_source_peaks = collapse_unique(
      tcf3_feature_id[tcf3_source_pass]))
}, by = bdev_feature_id]

# Route 2: scWorkbench supplies native log1p-CP10K RNA and peak matrices in a
# durable report-only session. The engine has no feature-by-cell-vector screen,
# so the 259 local Spearman tests are computed reproducibly from typed reads.
log_event("Route 2: native scWorkbench RNA/ATAC normalization")
bdev <- wb_open_bundle(bdev_bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
bdev_fp <- wb_manifest(bdev)$current_state$fingerprint
stopifnot(identical(bdev_fp, expected_bdev_fp),
  "Ebf1" %in% wb_feature_ids(bdev, "RNA_counts"),
  all(peaks$bdev_feature_id %in% wb_feature_ids(bdev, "peaks_counts")))

session_id <- "ebf1_putative_enhancers_three_way_v1"
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) wb_session_resume(bdev, workspace) else
  wb_session_create(bdev, workspace_root = session_root, session_id = session_id)
stopifnot(is.null(session$durable$promotion_receipt))
find_run <- function(x, marker) {
  ids <- names(x$runs)[vapply(x$runs, function(y) {
    note <- if (is.null(y$notes)) "" else as.character(y$notes)
    identical(y$op_name, "matrix_log_normalize") &&
      startsWith(note, paste0("[", marker, "]"))
  }, logical(1L))]
  stopifnot(length(ids) <= 1L)
  if (length(ids)) x$runs[[ids[[1L]]]] else NULL
}
normalize_once <- function(session, marker, matrix_id, output_id) {
  run <- find_run(session, marker)
  if (!is.null(run)) return(session)
  result <- wb_session_run_op(session, bdev, "matrix_log_normalize",
    inputs = list(matrix = matrix_id, cell_set = "all",
      output_matrix_id = output_id, scale_factor = 10000,
      overwrite = FALSE),
    notes = paste0("[", marker, "] Ebf1 three-way enhancer report; ",
      "full native log1p-CP10K matrix; staged report-only; no promotion."))
  result$session
}
peak_norm_id <- "peaks_log1p_cp10k_ebf1_three_way_v1"
rna_norm_id <- "RNA_log1p_cp10k_ebf1_three_way_v1"
session <- normalize_once(session, "normalize_peaks", "peaks_counts", peak_norm_id)
session <- normalize_once(session, "normalize_rna", "RNA_counts", rna_norm_id)
stopifnot(is.null(session$durable$promotion_receipt))
view <- wb_session_view(session)
on.exit(try(wb_close_bundle(view), silent = TRUE), add = TRUE)

cell_input <- fread(bdev_cells_path)
stopifnot(nrow(cell_input) == 971L, !anyDuplicated(cell_input$cell_id),
  all(is.finite(cell_input$pt)))
cells <- as.character(cell_input$cell_id)
md <- as.data.table(wb_meta(bdev, cells = cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt", "nCount_RNA",
    "nCount_peaks")))
setnames(md, names(md)[1L], "cell_id")
md <- md[match(cells, cell_id)]
stopifnot(identical(md$cell_id, cells),
  max(abs(as.numeric(md$pt) - as.numeric(cell_input$pt))) < 1e-10,
  all(is.finite(md$nCount_RNA) & md$nCount_RNA > 0),
  all(is.finite(md$nCount_peaks) & md$nCount_peaks > 0))

ebf1_expression <- as.numeric(as.matrix(wb_matrix_read_features(view,
  rna_norm_id, "Ebf1", cells, missing = "NA"))[1L, ])
accessibility <- as.matrix(wb_matrix_read_features(view, peak_norm_id,
  peaks$bdev_feature_id, cells, missing = "NA"))
stopifnot(!anyNA(ebf1_expression), !anyNA(accessibility),
  identical(rownames(accessibility), peaks$bdev_feature_id))

expr_rank <- rank(ebf1_expression, ties.method = "average")
pt_rank <- rank(md$pt, ties.method = "average")
atac_depth_rank <- rank(log1p(md$nCount_peaks), ties.method = "average")
rna_depth_rank <- rank(log1p(md$nCount_RNA), ties.method = "average")
design <- cbind(1, ns(pt_rank, df = 3), atac_depth_rank, rna_depth_rank)
q <- qr.Q(qr(design))
expr_resid <- as.numeric(expr_rank - q %*% crossprod(q, expr_rank))
expr_resid_norm <- sqrt(sum(expr_resid^2))
stopifnot(expr_resid_norm > 0)
expr_resid <- expr_resid / expr_resid_norm
partial_df <- length(cells) - ncol(q) - 1L
p_from_rho <- function(r, df) {
  if (!is.finite(r)) return(NA_real_)
  r <- pmin(pmax(r, -1 + 1e-15), 1 - 1e-15)
  2 * pt(-abs(r) * sqrt(df / pmax(1 - r^2, 1e-30)), df = df)
}

correlation <- rbindlist(lapply(seq_len(nrow(accessibility)), function(i) {
  value <- accessibility[i, ]
  detected <- sum(value > 0)
  testable <- detected >= 20L && length(unique(value)) > 1L
  if (!testable) return(data.table(
    bdev_feature_id = rownames(accessibility)[i], detected_cells = detected,
    detection_rate = detected / length(value), mean_log1p_cp10k = mean(value),
    raw_spearman_rho = NA_real_, raw_pvalue = NA_real_,
    partial_spearman_rho = NA_real_, partial_pvalue = NA_real_, tested = FALSE))
  peak_rank <- rank(value, ties.method = "average")
  raw_rho <- suppressWarnings(cor(peak_rank, expr_rank))
  peak_resid <- as.numeric(peak_rank - q %*% crossprod(q, peak_rank))
  denom <- sqrt(sum(peak_resid^2))
  partial_rho <- if (denom > 0) sum((peak_resid / denom) * expr_resid) else NA_real_
  data.table(bdev_feature_id = rownames(accessibility)[i],
    detected_cells = detected, detection_rate = detected / length(value),
    mean_log1p_cp10k = mean(value), raw_spearman_rho = raw_rho,
    raw_pvalue = p_from_rho(raw_rho, length(value) - 2L),
    partial_spearman_rho = partial_rho,
    partial_pvalue = p_from_rho(partial_rho, partial_df), tested = TRUE)
}))
correlation[, raw_padj_bh := p.adjust(raw_pvalue, method = "BH")]
correlation[, partial_padj_bh := p.adjust(partial_pvalue, method = "BH")]
correlation[, correlation_direction := fifelse(partial_spearman_rho > 0,
  "positive", fifelse(partial_spearman_rho < 0, "negative", "zero_or_untested"))]
correlation[, correlation_significant := is.finite(partial_padj_bh) &
  partial_padj_bh < 0.05]
correlation[, correlation_route_pass := correlation_significant &
  partial_spearman_rho > 0]
setorder(correlation, partial_padj_bh, -partial_spearman_rho,
  bdev_feature_id, na.last = TRUE)
write_tsv(correlation, file.path(stage_dir,
  "route2_ebf1_expression_peak_accessibility_correlations.tsv"))
fwrite(data.table(cell_id = cells, pt = md$pt, nCount_RNA = md$nCount_RNA,
  nCount_peaks = md$nCount_peaks, Ebf1_log1p_cp10k = ebf1_expression),
  file.path(stage_dir, "route2_bdev_il7_peak2_ebf1_expression.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip")

# Route 3: exact promoter-anchor containment and +/-20 kb expansion of the
# opposite anchor before overlap with the frozen local Bdev peak universe.
log_event("Route 3: mapping Ebf1-promoter PLAC-seq contacts")
loops <- fread(cmd = paste("gzip -cd", shQuote(placseq_path)))
expected_loop_cols <- c("chr1", "s1", "e1", "chr2", "s2", "e2", "cc",
  "P-Value_Bias", "Q-Value_Bias", "bin1_low", "bin1_high", "bin2_low",
  "bin2_high", "sumCC", "StrongConn")
stopifnot(identical(names(loops), expected_loop_cols), nrow(loops) == 69898L)
loops[, interaction_id := sprintf("GSE162858_%06d", .I)]
anchor1_promoter <- loops$chr1 == "chr11" & loops$s1 <= ebf1_tss &
  loops$e1 >= ebf1_tss
anchor2_promoter <- loops$chr2 == "chr11" & loops$s2 <= ebf1_tss &
  loops$e2 >= ebf1_tss
promoter_loops <- loops[xor(anchor1_promoter, anchor2_promoter)]
promoter_loops[, promoter_anchor := ifelse(
  chr1 == "chr11" & s1 <= ebf1_tss & e1 >= ebf1_tss, "anchor1", "anchor2")]
promoter_loops[, `:=`(
  promoter_chr = "chr11",
  promoter_start = ifelse(promoter_anchor == "anchor1", s1, s2),
  promoter_end = ifelse(promoter_anchor == "anchor1", e1, e2),
  distal_chr = ifelse(promoter_anchor == "anchor1", chr2, chr1),
  distal_start = ifelse(promoter_anchor == "anchor1", s2, s1),
  distal_end = ifelse(promoter_anchor == "anchor1", e2, e1)
)]
promoter_loops[, `:=`(
  distal_expanded_start = pmax(1L, distal_start - anchor_slop_bp),
  distal_expanded_end = distal_end + anchor_slop_bp
)]
stopifnot(nrow(promoter_loops) > 0L,
  all(promoter_loops$`Q-Value_Bias` <= 0.05))
write_tsv(promoter_loops, file.path(stage_dir,
  "route3_all_ebf1_promoter_placseq_interactions.tsv"))

loop_intervals <- promoter_loops[, .(chromosome = distal_chr,
  start = distal_expanded_start, end = distal_expanded_end, interaction_id,
  promoter_anchor, promoter_start, promoter_end, distal_start, distal_end,
  distal_expanded_start, distal_expanded_end, cc,
  plac_pvalue = `P-Value_Bias`, plac_qvalue = `Q-Value_Bias`,
  bin1_low, bin1_high, bin2_low, bin2_high, sumCC, StrongConn)]
setkey(peak_intervals, chromosome, start, end)
plac_map <- foverlaps(loop_intervals, peak_intervals,
  by.x = c("chromosome", "start", "end"),
  by.y = c("chromosome", "start", "end"), type = "any", nomatch = 0L)
setnames(plac_map, c("i.start", "i.end"),
  c("expanded_anchor_start", "expanded_anchor_end"))
plac_map[, `:=`(bdev_start_1based = start, bdev_end_1based = end,
  overlap_bp = pmax(0L, pmin(end, expanded_anchor_end) -
    pmax(start, expanded_anchor_start) + 1L))]
setorder(plac_map, bdev_feature_id, plac_qvalue, -cc, interaction_id)
write_tsv(plac_map, file.path(stage_dir,
  "route3_all_placseq_to_bdev_peak_overlaps.tsv"))
plac_summary <- plac_map[, .(
  placseq_route_pass = TRUE,
  placseq_interaction_count = uniqueN(interaction_id),
  placseq_best_qvalue = min(plac_qvalue),
  placseq_max_contact_count = max(cc),
  placseq_max_sum_contact_count = max(sumCC),
  placseq_interaction_ids = collapse_unique(interaction_id),
  placseq_distal_anchors = collapse_unique(paste0(chromosome, ":",
    distal_start, "-", distal_end))
), by = bdev_feature_id]

# Master table: all 259 Bdev peaks, including zero-evidence and promoter peaks.
master <- merge(peaks, tcf3_summary, by = "bdev_feature_id", all.x = TRUE,
  sort = FALSE)
master <- merge(master, correlation, by = "bdev_feature_id", all.x = TRUE,
  sort = FALSE)
master <- merge(master, plac_summary, by = "bdev_feature_id", all.x = TRUE,
  sort = FALSE)
master[is.na(tcf3_overlapping_source_peaks), `:=`(
  tcf3_overlapping_source_peaks = 0L, tcf3_passing_source_peaks = 0L,
  tcf3_route_pass = FALSE)]
master[is.na(placseq_interaction_count), `:=`(
  placseq_route_pass = FALSE, placseq_interaction_count = 0L)]
stopifnot(nrow(master) == 259L, !anyNA(master$correlation_route_pass),
  !anyDuplicated(master$bdev_feature_id))
master[, evidence_count := as.integer(tcf3_route_pass) +
  as.integer(correlation_route_pass) + as.integer(placseq_route_pass)]
master[, evidence_routes := vapply(seq_len(.N), function(i) {
  x <- c(if (tcf3_route_pass[i]) "TCF3_DA",
    if (correlation_route_pass[i]) "Ebf1_correlation",
    if (placseq_route_pass[i]) "PLACseq")
  if (length(x)) paste(x, collapse = "+") else "none"
}, character(1L))]
master[, putative_ebf1_enhancer := distal_putative_enhancer & evidence_count > 0L]
master[, abs_distance_to_tss_sort := abs(interval_distance_to_tss_bp)]
setorder(master, -putative_ebf1_enhancer, -evidence_count,
  partial_padj_bh, tcf3_best_padj, placseq_best_qvalue,
  abs_distance_to_tss_sort, na.last = TRUE)
master[, abs_distance_to_tss_sort := NULL]
master[, priority_rank := fifelse(putative_ebf1_enhancer,
  seq_len(.N), NA_integer_)]
master[putative_ebf1_enhancer == FALSE, priority_rank := NA_integer_]
write_tsv(master, file.path(stage_dir,
  "all_259_bdev_ebf1_window_peaks_three_way_tests.tsv"))
write_tsv(master[putative_ebf1_enhancer == TRUE], file.path(stage_dir,
  "putative_ebf1_enhancer_candidates.tsv"))
write_tsv(master[evidence_count >= 2L & distal_putative_enhancer == TRUE],
  file.path(stage_dir, "putative_ebf1_enhancers_two_or_more_routes.tsv"))

combination <- master[, .(peaks = .N,
  distal_peaks = sum(distal_putative_enhancer),
  putative_enhancers = sum(putative_ebf1_enhancer)), by = evidence_routes]
setorder(combination, -putative_enhancers, -peaks, evidence_routes)
write_tsv(combination, file.path(stage_dir, "evidence_combination_counts.tsv"))

p_locus <- ggplot(master, aes(midpoint_relative_to_tss_bp / 1e6,
  evidence_count, colour = evidence_routes,
  shape = overlaps_ebf1_promoter)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey35") +
  scale_y_continuous(breaks = 0:3, limits = c(-0.15, 3.15)) +
  labs(title = "Three-way evidence for putative Ebf1 enhancers",
    x = "Bdev peak midpoint relative to Ebf1 TSS (Mb)",
    y = "Independent evidence routes", colour = "Evidence", shape = "Promoter overlap") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(stage_dir, "figures", "three_way_evidence_by_locus.png"),
  p_locus, width = 11, height = 5.8, dpi = 220)
ggsave(file.path(stage_dir, "figures", "three_way_evidence_by_locus.pdf"),
  p_locus, width = 11, height = 5.8)

plot_cor <- correlation[is.finite(partial_spearman_rho) &
  is.finite(partial_padj_bh)]
p_cor <- ggplot(plot_cor, aes(partial_spearman_rho,
  -log10(pmax(partial_padj_bh, 1e-300)),
  colour = correlation_route_pass)) + geom_point(alpha = 0.75) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#0072B2")) +
  labs(title = "Ebf1 expression–peak accessibility association",
    subtitle = "971 Bdev cells in Il7_shared_peak2; adjusted for pt and RNA/ATAC depth",
    x = "Partial Spearman rho", y = expression(-log[10](BH-adjusted~P)),
    colour = "Positive BH < 0.05") + theme_bw(base_size = 10)
ggsave(file.path(stage_dir, "figures", "ebf1_peak_correlation.png"),
  p_cor, width = 8, height = 5.5, dpi = 220)

counts <- list(
  bdev_window_peaks = nrow(master),
  promoter_overlapping_bdev_peaks = sum(master$overlaps_ebf1_promoter),
  distal_bdev_peaks = sum(master$distal_putative_enhancer),
  tcf3_local_source_peaks = nrow(tcf3_local),
  tcf3_local_source_peaks_passing = sum(tcf3_local$tcf3_source_pass),
  tcf3_route_bdev_peaks = sum(master$tcf3_route_pass),
  correlation_tested_bdev_peaks = sum(master$tested),
  correlation_significant_any_direction = sum(master$correlation_significant),
  correlation_positive_route_bdev_peaks = sum(master$correlation_route_pass),
  placseq_total_interactions = nrow(loops),
  placseq_ebf1_promoter_interactions = nrow(promoter_loops),
  placseq_route_bdev_peaks = sum(master$placseq_route_pass),
  putative_enhancer_union = sum(master$putative_ebf1_enhancer),
  putative_enhancers_two_or_more_routes = sum(master$distal_putative_enhancer &
    master$evidence_count >= 2L),
  putative_enhancers_all_three_routes = sum(master$distal_putative_enhancer &
    master$evidence_count == 3L)
)

completion <- list(
  schema_version = "ebf1_putative_enhancers_three_way.v2",
  status = "complete", completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script_version = script_version, seed = seed,
  bdev_bundle_fingerprint = bdev_fp,
  ebf1 = list(gene = "Ebf1", gene_id = "ENSMUSG00000057098",
    genome = "GRCm38/mm10", chromosome = "chr11", tss_1based = ebf1_tss,
    window_start_1based = window_start, window_end_1based = window_end,
    promoter_start_1based = promoter_start, promoter_end_1based = promoter_end),
  route1 = list(rule = "GSE259355 abs(log2FC)>1 and DESeq2 BH padj<0.05; direct mm10 interval overlap to Bdev peaks",
    log2fc_definition = "TCF3_depleted/control"),
  route2 = list(cells = "exact 971 Bdev cells in Il7_shared_peak2",
    matrices = c(rna_norm_id, peak_norm_id),
    method = "Spearman plus partial rank correlation controlling ns(pt rank, df=3), RNA depth rank, and ATAC depth rank",
    candidate_rule = "positive partial rho and BH padj<0.05 across 259 local peaks"),
  route3 = list(source = "GSE162858 230238 H3K4me3 PLAC-seq FitHiChIP",
    promoter_rule = "Ebf1 TSS contained in exactly one 5-kb interaction anchor",
    distal_peak_rule = "Bdev peak overlaps opposite anchor expanded 20 kb on each border"),
  counts = counts,
  input_sha256 = list(bdev_cells = sha256_file(bdev_cells_path),
    bdev_peak_universe = sha256_file(bdev_peaks_path),
    tcf3_da = sha256_file(tcf3_da_path), placseq = sha256_file(placseq_path))
)
write_json(completion, file.path(stage_dir, "completion.json"))

report <- c("# Three-way identification of putative Ebf1 enhancers", "",
  "## Scope", "",
  paste0("All ", counts$bdev_window_peaks, " Bdev peaks overlapping mm10 Ebf1 TSS ±3 Mb were evaluated and retained in the master table."),
  paste0("Promoter-overlapping peaks are reported but excluded from the putative-enhancer union; ", counts$distal_bdev_peaks, " peaks are distal by the ±2 kb promoter rule."), "",
  "## Evidence routes", "",
  paste0("1. **TCF3 differential accessibility:** ", counts$tcf3_route_bdev_peaks, " Bdev peaks overlap a GSE259355 peak with `abs(log2FC) > 1` and BH-adjusted P < 0.05."),
  paste0("2. **Ebf1 expression/accessibility correlation:** ", counts$correlation_positive_route_bdev_peaks, " peaks have positive partial Spearman correlation with BH-adjusted P < 0.05 in the exact 971-cell Bdev `Il7_shared_peak2` interval. ", counts$correlation_significant_any_direction, " peaks are significant in either direction."),
  paste0("3. **Ebf1 promoter PLAC-seq interaction:** ", counts$placseq_ebf1_promoter_interactions, " significant source interactions contain the Ebf1 TSS; their ±20 kb-expanded distal anchors overlap ", counts$placseq_route_bdev_peaks, " Bdev peaks."), "",
  "## Combined candidates", "",
  paste0("- Distal union supported by at least one route: ", counts$putative_enhancer_union),
  paste0("- Supported by at least two routes: ", counts$putative_enhancers_two_or_more_routes),
  paste0("- Supported by all three routes: ", counts$putative_enhancers_all_three_routes), "",
  "![Three-way evidence](figures/three_way_evidence_by_locus.png)", "",
  "![Ebf1 correlations](figures/ebf1_peak_correlation.png)", "",
  "## Primary outputs", "",
  "- `all_259_bdev_ebf1_window_peaks_three_way_tests.tsv`: complete one-row-per-Bdev-peak audit.",
  "- `putative_ebf1_enhancer_candidates.tsv`: distal union supported by one or more routes.",
  "- `putative_ebf1_enhancers_two_or_more_routes.tsv`: higher-confidence multi-route subset.",
  "- `route1_all_tcf3_to_bdev_peak_overlaps.tsv`: complete source-to-Bdev overlap mapping.",
  "- `route2_ebf1_expression_peak_accessibility_correlations.tsv`: complete 259-peak correlation test.",
  "- `route3_all_ebf1_promoter_placseq_interactions.tsv` and `route3_all_placseq_to_bdev_peak_overlaps.tsv`: complete promoter-loop and Bdev mapping audits.", "",
  "## Interpretation limits", "",
  "GSE259355 combines two 6-hour and one 8-hour depleted libraries. Correlations are cell-level associations within a restricted developmental interval and do not establish causality. GSE162858 PLAC-seq is orthogonal evidence from the 230-238 B-cell progenitor line, not the Bdev cells themselves. No result was promoted into a parent bundle.")
writeLines(report, file.path(stage_dir, "REPORT.md"))
html <- c("<!doctype html><html><head><meta charset='utf-8'><title>Ebf1 putative enhancers</title>",
  "<style>body{font-family:sans-serif;max-width:1100px;margin:2rem auto;line-height:1.45}img{max-width:100%}code{background:#eee;padding:.1rem .25rem}</style></head><body>",
  "<h1>Three-way identification of putative Ebf1 enhancers</h1>",
  paste0("<p>All <b>", counts$bdev_window_peaks, "</b> Bdev peaks in the Ebf1 ±3 Mb window were tested. The distal candidate union contains <b>", counts$putative_enhancer_union, "</b> peaks; <b>", counts$putative_enhancers_two_or_more_routes, "</b> have at least two evidence routes and <b>", counts$putative_enhancers_all_three_routes, "</b> have all three.</p>"),
  "<img src='figures/three_way_evidence_by_locus.png'><img src='figures/ebf1_peak_correlation.png'>",
  "<h2>Downloads</h2><ul>",
  "<li><a href='all_259_bdev_ebf1_window_peaks_three_way_tests.tsv'>All 259 Bdev peaks and tests</a></li>",
  "<li><a href='putative_ebf1_enhancer_candidates.tsv'>Candidate union</a></li>",
  "<li><a href='putative_ebf1_enhancers_two_or_more_routes.tsv'>Two-or-more-route candidates</a></li>",
  "<li><a href='route1_all_tcf3_to_bdev_peak_overlaps.tsv'>TCF3 overlap mapping</a></li>",
  "<li><a href='route2_ebf1_expression_peak_accessibility_correlations.tsv'>Ebf1 correlations</a></li>",
  "<li><a href='route3_all_placseq_to_bdev_peak_overlaps.tsv'>PLAC-seq mapping</a></li></ul>",
  "</body></html>")
writeLines(html, file.path(stage_dir, "REPORT.html"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

if (!file.rename(stage_dir, output_dir)) stop("Could not publish output directory")
cat("Completed Ebf1 three-way enhancer report: ", output_dir, "\n", sep = "")
