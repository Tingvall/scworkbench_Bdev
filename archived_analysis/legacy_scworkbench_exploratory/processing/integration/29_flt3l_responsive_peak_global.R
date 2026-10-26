#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L) stop(paste(
  "Usage: 29_flt3l_responsive_peak_global.R BDEV_BUNDLE QUERY_BUNDLE",
  "PEAK_MAPPING SHARED_REGION_REPORT SIGNATURE_MEMBERSHIP BDEV_RDATA",
  "OUTPUT_DIR BDEV_SESSION_ROOT QUERY_SESSION_ROOT EXPECTED_BDEV_FP",
  "EXPECTED_QUERY_FP SCRIPT_VERSION SEED"
), call. = FALSE)

bdev_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
query_bundle_dir <- normalizePath(args[[2L]], mustWork = TRUE)
peak_mapping_path <- normalizePath(args[[3L]], mustWork = TRUE)
shared_report_dir <- normalizePath(args[[4L]], mustWork = TRUE)
signature_path <- normalizePath(args[[5L]], mustWork = TRUE)
bdev_rdata_path <- normalizePath(args[[6L]], mustWork = TRUE)
output_dir <- normalizePath(args[[7L]], mustWork = FALSE)
bdev_session_root <- normalizePath(args[[8L]], mustWork = FALSE)
query_session_root <- normalizePath(args[[9L]], mustWork = FALSE)
expected_bdev_fp <- args[[10L]]
expected_query_fp <- args[[11L]]
script_version <- args[[12L]]
seed <- as.integer(args[[13L]])

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)
dir.create(bdev_session_root, recursive = TRUE, showWarnings = FALSE)
dir.create(query_session_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(jsonlite)
  library(digest)
  library(splines)
  library(SeuratObject)
  library(Signac)
})

write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
log_event <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(..., collapse = " "), "\n")

region <- fread(file.path(shared_report_dir, "shared_enriched_pt_regions.tsv"))[
  region_id == "Flt3l_shared_peak1"]
assignments <- fread(cmd = paste("zcat", shQuote(file.path(shared_report_dir,
  "shared_region_cell_assignments.tsv.gz"))))[region_id == "Flt3l_shared_peak1"]
stopifnot(nrow(region) == 1L, region$source_peak == "Flt3l_peak1",
  nrow(assignments) > 0L, all(assignments$mapped_pt >= region$pt_left),
  all(assignments$mapped_pt < region$pt_right))
write_tsv(region, file.path(stage_dir, "flt3l_peak1_region_definition.tsv"))
fwrite(assignments, file.path(stage_dir, "flt3l_peak1_query_cell_assignments.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip")

membership <- fread(signature_path)
flt3l_genes <- membership[cytokine == "Flt3L" & tolower(direction) == "up"]
setorder(flt3l_genes, final_rank)
stopifnot(nrow(flt3l_genes) >= 5L, !anyDuplicated(flt3l_genes$feature_id))
write_tsv(flt3l_genes, file.path(stage_dir, "flt3l_up_signature_genes.tsv"))

peak_mapping <- fread(cmd = paste("zcat", shQuote(peak_mapping_path)))[
  status == "accepted"]
stopifnot(!anyDuplicated(peak_mapping$bdev_feature_id))

bdev <- wb_open_bundle(bdev_bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
query <- wb_open_bundle(query_bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(query), silent = TRUE), add = TRUE)
bdev_fp <- wb_manifest(bdev)$current_state$fingerprint
query_fp <- wb_manifest(query)$current_state$fingerprint
stopifnot(identical(bdev_fp, expected_bdev_fp),
  identical(query_fp, expected_query_fp))

bdev_peak_ids <- wb_feature_ids(bdev, "peaks_counts")
query_peak_ids <- wb_feature_ids(query, "ATAC_bdev_mm10_peaks")
common_peaks <- intersect(bdev_peak_ids, query_peak_ids)
stopifnot(length(bdev_peak_ids) == 97846L, length(query_peak_ids) == 97788L,
  length(common_peaks) == 97788L,
  all(common_peaks %in% peak_mapping$bdev_feature_id))

bdev_cells_all <- wb_cells(bdev, "set:bdev_reference")
bdev_md <- as.data.table(wb_meta(bdev, cells = bdev_cells_all,
  cols = list(wb_ref_metadata_field("core:cell_id"), "pt", "nCount_RNA",
    "nCount_peaks")))
setnames(bdev_md, names(bdev_md)[1L], "cell_id")
bdev_md <- bdev_md[is.finite(pt) & pt >= region$pt_left & pt < region$pt_right]
stopifnot(nrow(bdev_md) >= 100L, !anyDuplicated(bdev_md$cell_id),
  all(is.finite(bdev_md$nCount_RNA) & bdev_md$nCount_RNA > 0),
  all(is.finite(bdev_md$nCount_peaks) & bdev_md$nCount_peaks > 0))
bdev_cells <- as.character(bdev_md$cell_id)
write_tsv(bdev_md, file.path(stage_dir, "bdev_flt3l_peak1_cells.tsv"))

mapped_genes <- intersect(flt3l_genes$feature_id,
  wb_feature_ids(bdev, "RNA_counts"))
stopifnot(length(mapped_genes) >= 5L)
gene_counts <- as.matrix(wb_matrix_read_features(bdev, "RNA_counts",
  mapped_genes, bdev_cells, missing = "NA"))
stopifnot(!anyNA(gene_counts))
rna_depth <- bdev_md$nCount_RNA[match(colnames(gene_counts), bdev_md$cell_id)]
flt3l_signal <- colMeans(log1p(t(t(gene_counts) / rna_depth) * 10000))
signal_table <- data.table(cell_id = names(flt3l_signal),
  pt = bdev_md$pt[match(names(flt3l_signal), bdev_md$cell_id)],
  nCount_RNA = rna_depth,
  nCount_peaks = bdev_md$nCount_peaks[match(names(flt3l_signal), bdev_md$cell_id)],
  flt3l_up_mean_log1p_cp10k = as.numeric(flt3l_signal))
fwrite(signal_table, file.path(stage_dir, "bdev_flt3l_peak1_flt3l_up_signal.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip")

write_json(wb_op_interface("matrix_log_normalize", bundle = bdev),
  file.path(stage_dir, "bdev_matrix_log_normalize_interface.json"))
bdev_session_id <- "flt3l_responsive_peak_global_v1_bdev"
bdev_workspace <- file.path(bdev_session_root, bdev_session_id)
bdev_session <- if (dir.exists(bdev_workspace))
  wb_session_resume(bdev, bdev_workspace) else
  wb_session_create(bdev, workspace_root = bdev_session_root,
    session_id = bdev_session_id)
stopifnot(is.null(bdev_session$durable$promotion_receipt))
find_run <- function(session, marker, operation) {
  ids <- names(session$runs)[vapply(session$runs, function(x) {
    note <- if (is.null(x$notes)) "" else as.character(x$notes)
    identical(x$op_name, operation) && startsWith(note, paste0("[", marker, "]"))
  }, logical(1))]
  stopifnot(length(ids) <= 1L)
  if (length(ids)) session$runs[[ids[[1L]]]] else NULL
}
bdev_norm_id <- "peaks_log1p_cp10k_flt3l_global_v1"
bdev_norm_run <- find_run(bdev_session, "normalize_bdev_peaks", "matrix_log_normalize")
if (is.null(bdev_norm_run)) {
  log_event("Native log1p-CP10K normalization of Bdev peaks_counts")
  result <- wb_session_run_op(bdev_session, bdev, "matrix_log_normalize",
    inputs = list(matrix = "peaks_counts", cell_set = "set:bdev_reference",
      min_cells = 10L, output_matrix_id = bdev_norm_id,
      scale_factor = 10000, overwrite = FALSE),
    notes = paste0("[normalize_bdev_peaks] Full Bdev peaks_counts log1p-CP10K; ",
      "staged report-only matrix; no promotion."))
  bdev_session <- result$session
  bdev_norm_run <- result$run
}
bdev_view <- wb_session_view(bdev_session)
on.exit(try(wb_close_bundle(bdev_view), silent = TRUE), add = TRUE)

# Vectorized report-only screen. scWorkbench supplies the typed matrix and
# normalization; it has no matrix-by-arbitrary-cell-vector association op.
score_rank <- rank(flt3l_signal, ties.method = "average")
pt_rank <- rank(signal_table$pt, ties.method = "average")
depth_rank <- rank(log1p(signal_table$nCount_peaks), ties.method = "average")
design <- cbind(1, ns(pt_rank, df = 3), depth_rank)
q <- qr.Q(qr(design))
score_resid <- as.numeric(score_rank - q %*% crossprod(q, score_rank))
score_resid <- score_resid / sqrt(sum(score_resid^2))
n_cells <- length(bdev_cells)
residual_df <- n_cells - ncol(q) - 1L
min_detected <- max(20L, ceiling(0.01 * n_cells))

screen_blocks <- split(common_peaks,
  ceiling(seq_along(common_peaks) / 128L))
screen_parts <- vector("list", length(screen_blocks))
for (bi in seq_along(screen_blocks)) {
  feats <- screen_blocks[[bi]]
  if (bi %% 50L == 1L) log_event("Correlation block", bi, "of",
    length(screen_blocks))
  m <- as.matrix(wb_matrix_read_features(bdev_view, bdev_norm_id, feats,
    bdev_cells, missing = "NA"))
  detected <- rowSums(m > 0)
  means <- rowMeans(m)
  raw_rho <- rep(NA_real_, nrow(m))
  adjusted_rho <- rep(NA_real_, nrow(m))
  keep <- which(detected >= min_detected & apply(m, 1L, function(z)
    length(unique(z)) > 1L))
  if (length(keep)) for (j in keep) {
    peak_rank <- rank(m[j, ], ties.method = "average")
    raw_rho[j] <- suppressWarnings(cor(peak_rank, score_rank))
    peak_resid <- peak_rank - q %*% crossprod(q, peak_rank)
    denom <- sqrt(sum(peak_resid^2))
    if (denom > 0) adjusted_rho[j] <- sum((peak_resid / denom) * score_resid)
  }
  p_from_rho <- function(r, df) {
    r <- pmin(pmax(r, -1 + 1e-15), 1 - 1e-15)
    2 * pt(-abs(r) * sqrt(df / pmax(1 - r^2, 1e-30)), df = df)
  }
  screen_parts[[bi]] <- data.table(feature_id = feats,
    detected_cells = detected, detection_rate = detected / n_cells,
    mean_log1p_cp10k = means, raw_spearman_rho = raw_rho,
    raw_p_value = p_from_rho(raw_rho, n_cells - 2L),
    adjusted_spearman_rho = adjusted_rho,
    adjusted_p_value = p_from_rho(adjusted_rho, residual_df))
}
screen <- rbindlist(screen_parts)
screen[, raw_p_adj_bh := p.adjust(raw_p_value, method = "BH")]
screen[, adjusted_p_adj_bh := p.adjust(adjusted_p_value, method = "BH")]
screen[, `:=`(tested = is.finite(adjusted_spearman_rho),
  flt3l_responsive = is.finite(adjusted_p_adj_bh) & adjusted_p_adj_bh < 0.05 &
    adjusted_spearman_rho > 0 & raw_spearman_rho > 0,
  high_confidence = is.finite(adjusted_p_adj_bh) & adjusted_p_adj_bh < 0.05 &
    adjusted_spearman_rho >= 0.1 & raw_spearman_rho > 0)]
screen <- merge(screen, peak_mapping[, .(feature_id = bdev_feature_id,
  chromosome = mm10_seqname, peak_start_1based = mm10_start_1based,
  peak_end_1based = mm10_end_1based)], by = "feature_id", all.x = TRUE)
setorder(screen, adjusted_p_adj_bh, -adjusted_spearman_rho, feature_id)
fwrite(screen, file.path(stage_dir, "all_bdev_peak_flt3l_associations.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip")
candidates <- screen[flt3l_responsive == TRUE]
write_tsv(candidates, file.path(stage_dir, "flt3l_responsive_peak_candidates.tsv"))
write_tsv(candidates[high_confidence == TRUE], file.path(stage_dir,
  "flt3l_responsive_peak_candidates_high_confidence_rho0p1.tsv"))
if (!nrow(candidates)) stop("No FLT3L-responsive peaks passed the frozen rule",
  call. = FALSE)

write_json(wb_op_interface("matrix_log_normalize", bundle = query),
  file.path(stage_dir, "query_matrix_log_normalize_interface.json"))
write_json(wb_op_interface("de_wilcox", bundle = query),
  file.path(stage_dir, "query_de_wilcox_interface.json"))
query_session_id <- "flt3l_responsive_peak_global_v1_query"
query_workspace <- file.path(query_session_root, query_session_id)
query_session <- if (dir.exists(query_workspace))
  wb_session_resume(query, query_workspace) else
  wb_session_create(query, workspace_root = query_session_root,
    session_id = query_session_id)
stopifnot(is.null(query_session$durable$promotion_receipt))
query_norm_id <- "ATAC_bdev_mm10_peaks_log1p_cp10k_flt3l_global_v1"
query_norm_run <- find_run(query_session, "normalize_query_peaks",
  "matrix_log_normalize")
if (is.null(query_norm_run)) {
  log_event("Native log1p-CP10K normalization of query Bdev-aligned peaks")
  result <- wb_session_run_op(query_session, query, "matrix_log_normalize",
    inputs = list(matrix = "ATAC_bdev_mm10_peaks", cell_set = "all",
      min_cells = 10L, output_matrix_id = query_norm_id,
      scale_factor = 10000, overwrite = FALSE),
    notes = paste0("[normalize_query_peaks] Full query/reference Bdev-aligned ",
      "peak matrix log1p-CP10K; staged report-only; no promotion."))
  query_session <- result$session
  query_norm_run <- result$run
}
query_bundle_cells <- wb_cells(query, "set:ko_query")
query_id_map <- as.data.table(wb_meta(query, cells = query_bundle_cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "source_cell_id",
    "sample_id", "genotype", "knockout_target")))
setnames(query_id_map, names(query_id_map)[1L], "cell_id")
stopifnot(!anyDuplicated(query_id_map$source_cell_id))
assignment_map <- merge(assignments, query_id_map,
  by = c("source_cell_id", "sample_id", "genotype", "knockout_target"),
  all.x = TRUE, sort = FALSE)
stopifnot(nrow(assignment_map) == nrow(assignments),
  !anyNA(assignment_map$cell_id), !anyDuplicated(assignment_map$cell_id))
fwrite(assignment_map, file.path(stage_dir,
  "flt3l_peak1_query_cell_id_mapping.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")
wt_cells <- as.character(assignment_map[genotype == "WT", cell_id])
ko_cells <- as.character(assignment_map[genotype == "KO", cell_id])
stopifnot(all(c(wt_cells, ko_cells) %in% query_bundle_cells))
de_run <- find_run(query_session, "flt3l_peak1_candidate_da", "de_wilcox")
if (is.null(de_run)) {
  log_event("Native differential accessibility for", nrow(candidates),
    "frozen candidates")
  result <- wb_session_run_op(query_session, query, "de_wilcox",
    inputs = list(matrix = query_norm_id, cell_set = "set:ko_query",
      min_cells = 10L, group_a = wt_cells, group_b = ko_cells,
      features = as.character(candidates$feature_id), min_pct = 0,
      p_adjust_method = "BH"),
    notes = paste0("[flt3l_peak1_candidate_da] Group A=FLT3L WT, group B=FLT3LKO; ",
      "exact saved shared peak1 cells; frozen Bdev FLT3L-responsive peaks; ",
      "descriptive cell-level test; no biological replicates."))
  query_session <- result$session
  de_run <- result$run
}
query_view <- wb_session_view(query_session)
on.exit(try(wb_close_bundle(query_view), silent = TRUE), add = TRUE)
da <- as.data.table(de_run$outputs$de_table$table)
stopifnot(setequal(da$feature_id, candidates$feature_id))
mean_blocks <- split(candidates$feature_id,
  ceiling(seq_len(nrow(candidates)) / 128L))
means <- rbindlist(lapply(mean_blocks, function(feats) {
  m <- expm1(as.matrix(wb_matrix_read_features(query_view, query_norm_id,
    feats, c(wt_cells, ko_cells), missing = "NA")))
  data.table(feature_id = feats,
    mean_cp10k_wt = rowMeans(m[, seq_along(wt_cells), drop = FALSE]),
    mean_cp10k_flt3lko = rowMeans(m[, length(wt_cells) + seq_along(ko_cells),
      drop = FALSE]))
}))
da <- merge(da, means, by = "feature_id", all.x = TRUE)
da[, standard_log2fc := log2((mean_cp10k_wt + 1) /
  (mean_cp10k_flt3lko + 1))]
da <- merge(da, candidates[, .(feature_id, chromosome, peak_start_1based,
  peak_end_1based, detected_cells, detection_rate, mean_log1p_cp10k,
  raw_spearman_rho, raw_p_value, raw_p_adj_bh, adjusted_spearman_rho,
  adjusted_p_value, adjusted_p_adj_bh, high_confidence)],
  by = "feature_id", all.x = TRUE)
da[, `:=`(direction = fifelse(standard_log2fc > 0, "WT_higher",
  fifelse(standard_log2fc < 0, "FLT3LKO_higher", "equal")),
  query_significant_bh05 = is.finite(p_val_adj) & p_val_adj < 0.05,
  globally_validated = is.finite(p_val_adj) & p_val_adj < 0.05 &
    is.finite(standard_log2fc) & standard_log2fc > 0)]
setorder(da, p_val_adj, -standard_log2fc, feature_id)
write_tsv(da, file.path(stage_dir,
  "candidate_peaks_flt3lko_vs_wt_differential_accessibility.tsv"))
validated <- da[globally_validated == TRUE]
write_tsv(validated, file.path(stage_dir,
  "globally_validated_flt3l_responsive_peaks.tsv"))

# Read-only recovery of the exact motif matches omitted from the bundle motif
# registry at import. This does not alter the Seurat object or any bundle.
log_event("Loading read-only Bdev motif-to-peak matrix")
source_env <- new.env(parent = emptyenv())
con <- gzfile(bdev_rdata_path, "rb")
loaded <- load(con, envir = source_env)
close(con)
stopifnot("Bdev" %in% loaded)
source_bdev <- source_env$Bdev
motif_matrix <- Signac::GetMotifData(source_bdev[["peaks"]], slot = "data")
peak_covariates <- as.data.table(source_bdev[["peaks"]][[]], keep.rownames =
  "feature_id")
stopifnot(nrow(motif_matrix) == length(bdev_peak_ids),
  all(common_peaks %in% rownames(motif_matrix)),
  all(c("GC.percent", "sequence.length") %in% names(peak_covariates)))
motif_ids <- colnames(motif_matrix)
motif_names <- tryCatch(Signac::ConvertMotifID(source_bdev[["peaks"]],
  id = motif_ids), error = function(e) setNames(motif_ids, motif_ids))
if (is.null(names(motif_names))) names(motif_names) <- motif_ids
rm(source_bdev, source_env); invisible(gc())

tested_cov <- merge(screen[tested == TRUE, .(feature_id, detection_rate,
  mean_log1p_cp10k)], peak_covariates[, .(feature_id, GC.percent,
  sequence.length)], by = "feature_id", all.x = TRUE)
stopifnot(!anyNA(tested_cov[, .(GC.percent, sequence.length)]))
make_background <- function(foreground_ids, label) {
  x <- copy(tested_cov)
  x[, foreground := feature_id %in% foreground_ids]
  bin <- function(z) as.integer(cut(z, unique(quantile(z,
    probs = seq(0, 1, 0.1), na.rm = TRUE)), include.lowest = TRUE,
    labels = FALSE))
  x[, `:=`(access_bin = bin(mean_log1p_cp10k), gc_bin = bin(GC.percent),
    length_bin = bin(log1p(sequence.length)))]
  x[, stratum := paste(access_bin, gc_bin, length_bin, sep = ":")]
  set.seed(seed + match(label, c("responsive", "validated")))
  fg <- x[foreground == TRUE]
  bg <- x[foreground == FALSE][fg[, .(n_fg = .N), by = stratum],
    on = "stratum", allow.cartesian = TRUE][, .SD[sample.int(.N,
      min(.N, 5L * n_fg[[1L]]))], by = stratum]
  unique(bg$feature_id)
}
motif_enrichment <- function(foreground_ids, label) {
  if (!length(foreground_ids)) return(data.table())
  background_ids <- make_background(foreground_ids, label)
  universe_ids <- c(foreground_ids, background_ids)
  mm <- motif_matrix[universe_ids, , drop = FALSE]
  fg_n <- length(foreground_ids); bg_n <- length(background_ids)
  fg_hits <- Matrix::colSums(mm[seq_len(fg_n), , drop = FALSE] > 0)
  bg_hits <- Matrix::colSums(mm[fg_n + seq_len(bg_n), , drop = FALSE] > 0)
  p <- vapply(seq_along(motif_ids), function(i) fisher.test(matrix(c(
    fg_hits[i], fg_n - fg_hits[i], bg_hits[i], bg_n - bg_hits[i]),
    nrow = 2L), alternative = "greater")$p.value, numeric(1))
  odds <- ((fg_hits + 0.5) / (fg_n - fg_hits + 0.5)) /
    ((bg_hits + 0.5) / (bg_n - bg_hits + 0.5))
  out <- data.table(set = label, motif_id = motif_ids,
    motif_name = as.character(motif_names[motif_ids]),
    foreground_peaks = fg_n, matched_background_peaks = bg_n,
    foreground_hits = as.numeric(fg_hits), background_hits = as.numeric(bg_hits),
    odds_ratio_haldane = odds, p_value = p, p_adj_bh = p.adjust(p, "BH"))
  write_tsv(data.table(set = label, feature_id = background_ids),
    file.path(stage_dir, paste0(label, "_motif_matched_background.tsv")))
  setorder(out, p_adj_bh, -odds_ratio_haldane, motif_id)
  out
}
motifs <- rbindlist(list(
  motif_enrichment(candidates$feature_id, "responsive"),
  motif_enrichment(validated$feature_id, "validated")), fill = TRUE)
write_tsv(motifs, file.path(stage_dir, "global_motif_enrichment.tsv"))

p1 <- ggplot(screen[tested == TRUE], aes(adjusted_spearman_rho,
  -log10(pmax(adjusted_p_adj_bh, 1e-300)), colour = flt3l_responsive)) +
  geom_point(size = 0.55, alpha = 0.55) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#009E73")) +
  labs(x = "Partial Spearman rho (adjusted for pt and ATAC depth)",
    y = "-log10(BH adjusted P)", colour = "FLT3L responsive",
    title = "Global Bdev FLT3L-responsive peak screen") + theme_bw(base_size = 10)
ggsave(file.path(stage_dir, "figures", "bdev_flt3l_peak_association.png"), p1,
  width = 8, height = 5.5, dpi = 180)

p2 <- ggplot(da, aes(standard_log2fc,
  -log10(pmax(p_val_adj, 1e-300)), colour = globally_validated)) +
  geom_point(size = 0.8, alpha = 0.65) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#009E73")) +
  labs(x = "Standard log2FC (WT / FLT3LKO; mean CP10K + 1)",
    y = "-log10(BH adjusted P)", colour = "Validated",
    title = "FLT3LKO versus WT accessibility in shared FLT3L peak1") +
  theme_bw(base_size = 10)
ggsave(file.path(stage_dir, "figures", "candidate_peak_da_volcano.png"), p2,
  width = 8, height = 5.5, dpi = 180)

if (nrow(motifs)) {
  mp <- motifs[p_adj_bh < 0.05 & is.finite(odds_ratio_haldane),
    head(.SD, 20L), by = set]
  if (nrow(mp)) {
    mp[, label := fifelse(is.na(motif_name) | motif_name == "", motif_id,
      paste0(motif_name, " (", motif_id, ")"))]
    p3 <- ggplot(mp, aes(log2(odds_ratio_haldane), reorder(label,
      log2(odds_ratio_haldane)), size = foreground_hits,
      colour = -log10(pmax(p_adj_bh, 1e-300)))) + geom_point() +
      facet_wrap(~set, scales = "free_y") +
      scale_colour_viridis_c() + labs(x = "log2 matched-background odds ratio",
        y = NULL, size = "Foreground hits", colour = "-log10(BH P)",
        title = "Motif enrichment in FLT3L-responsive elements") +
      theme_bw(base_size = 9)
    ggsave(file.path(stage_dir, "figures", "global_motif_enrichment.png"),
      p3, width = 10, height = 7, dpi = 180)
  }
}

summary <- list(schema_version = "flt3l_responsive_peak_global.v1",
  script_version = script_version, completed_at = format(Sys.time(),
    "%Y-%m-%dT%H:%M:%S%z"), seed = seed,
  inputs = list(bdev_bundle = bdev_bundle_dir, bdev_fingerprint = bdev_fp,
    query_bundle = query_bundle_dir, query_fingerprint = query_fp,
    peak_mapping = peak_mapping_path,
    peak_mapping_sha256 = sha256_file(peak_mapping_path),
    shared_region_report = shared_report_dir,
    signature_membership = signature_path,
    signature_membership_sha256 = sha256_file(signature_path),
    bdev_rdata = bdev_rdata_path,
    bdev_rdata_sha256 = sha256_file(bdev_rdata_path)),
  region = list(id = region$region_id, pt_left = region$pt_left,
    pt_right = region$pt_right),
  methods = list(flt3l_signal = paste0("arithmetic mean log1p CP10K over ",
    length(mapped_genes), " mapped Flt3L UP genes"),
    bdev_peak_normalization = "native scWorkbench matrix_log_normalize log1p CP10K",
    association = paste0("Spearman plus partial Spearman after residualizing ranks ",
      "against natural-spline(pt, df=3) and rank(log1p nCount_peaks)"),
    candidate_rule = "adjusted BH<0.05, adjusted rho>0, raw rho>0",
    high_confidence_rule = "candidate rule plus adjusted rho>=0.1",
    query_da = "native scWorkbench de_wilcox; group A WT, group B FLT3LKO",
    validated_rule = "Bdev candidate plus query BH<0.05 and standard log2FC>0",
    motif = paste0("exact source JASPAR motif matches; one-sided Fisher test ",
      "against deterministic up-to-5x controls matched within deciles of Bdev ",
      "mean accessibility, GC percent, and sequence length; BH correction"),
    inference = "descriptive cell-level; one biological sample per genotype"),
  counts = list(bdev_peak1_cells = n_cells,
    flt3l_up_genes = length(mapped_genes), peaks_in_common = length(common_peaks),
    peaks_tested = sum(screen$tested), candidates = nrow(candidates),
    high_confidence_candidates = sum(candidates$high_confidence),
    query_wt_cells = length(wt_cells), query_flt3lko_cells = length(ko_cells),
    globally_validated = nrow(validated)))
write_json(summary, file.path(stage_dir, "completion.json"))

report_lines <- c("# Global FLT3L-responsive peak analysis", "",
  paste0("Completed: ", summary$completed_at), "",
  "## Scope", "",
  paste0("Bdev discovery was restricted to `", region$region_id,
    "` (pt [", signif(region$pt_left, 6), ", ", signif(region$pt_right, 6),
    ")). The frozen candidates were then tested in the exact WT and FLT3LKO ",
    "query cells assigned to the same region."), "", "## Results", "",
  paste0("- Bdev cells: ", n_cells),
  paste0("- Peaks tested: ", sum(screen$tested)),
  paste0("- FLT3L-responsive candidates: ", nrow(candidates)),
  paste0("- High-confidence candidates (partial rho >= 0.1): ",
    sum(candidates$high_confidence)),
  paste0("- Query cells: ", length(wt_cells), " WT and ", length(ko_cells),
    " FLT3LKO"),
  paste0("- Globally validated WT-higher peaks: ", nrow(validated)), "",
  "## Primary outputs", "",
  "- `all_bdev_peak_flt3l_associations.tsv.gz`: complete global screen",
  "- `flt3l_responsive_peak_candidates.tsv`: frozen Bdev candidate list",
  "- `candidate_peaks_flt3lko_vs_wt_differential_accessibility.tsv`: full candidate DA",
  "- `globally_validated_flt3l_responsive_peaks.tsv`: final global intersection",
  "- `global_motif_enrichment.tsv`: motif analysis of responsive and validated sets",
  "", "## Interpretation limits", "",
  "The query DA is descriptive cell-level inference because each genotype has one biological sample. Motif enrichment characterizes the global responsive and validated peak sets and is not used to select peaks.")
writeLines(report_lines, file.path(stage_dir, "REPORT.md"))

file.rename(stage_dir, output_dir)
cat("Completed global FLT3L-responsive peak analysis:", output_dir, "\n")
