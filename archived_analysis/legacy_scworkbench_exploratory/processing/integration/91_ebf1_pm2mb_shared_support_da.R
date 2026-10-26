#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "Usage: 91_ebf1_pm2mb_shared_support_da.R <mapping_bundle>",
  "<peak_mapping> <support_report> <query_support_tsv_gz> <output_dir>",
  "<session_root> <expected_mapping_fingerprint> <seed> <script_version>"),
  call. = FALSE)
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
peak_mapping_path <- normalizePath(args[[2L]], mustWork = TRUE)
support_dir <- normalizePath(args[[3L]], mustWork = TRUE)
query_path <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
session_root <- normalizePath(args[[6L]], mustWork = FALSE)
expected_fp <- args[[7L]]
seed <- as.integer(args[[8L]])
script_version <- args[[9L]]
stopifnot(is.finite(seed))
if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE)
dir.create(session_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
log_event <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(..., collapse = " "), "\n")

ebf1 <- data.table(gene = "Ebf1", gene_id = "ENSMUSG00000057098",
  genome = "GRCm38/mm10", chromosome = "chr11", strand = "+",
  tss_1based = 44617317L, flank_bp = 2000000L)
ebf1[, `:=`(window_start_1based = tss_1based - flank_bp,
  window_end_1based = tss_1based + flank_bp)]
write_tsv(ebf1, file.path(stage_dir, "ebf1_tss_pm2mb_definition.tsv"))

peak_mapping <- fread(cmd = paste("zcat", shQuote(peak_mapping_path)))
stopifnot(all(c("bdev_feature_id", "mm10_seqname", "mm10_start_1based",
  "mm10_end_1based", "status") %in% names(peak_mapping)),
  !anyDuplicated(peak_mapping$bdev_feature_id))
window_peaks <- peak_mapping[status == "accepted" &
  mm10_seqname == ebf1$chromosome &
  mm10_end_1based >= ebf1$window_start_1based &
  mm10_start_1based <= ebf1$window_end_1based]
setorder(window_peaks, mm10_start_1based, mm10_end_1based)
window_peaks[, `:=`(
  peak_midpoint_1based = (mm10_start_1based + mm10_end_1based) / 2,
  midpoint_relative_to_tss_bp =
    (mm10_start_1based + mm10_end_1based) / 2 - ebf1$tss_1based,
  interval_distance_to_tss_bp = fifelse(
    mm10_start_1based <= ebf1$tss_1based &
      mm10_end_1based >= ebf1$tss_1based, 0,
    fifelse(mm10_end_1based < ebf1$tss_1based,
      mm10_end_1based - ebf1$tss_1based,
      mm10_start_1based - ebf1$tss_1based)))]
stopifnot(nrow(window_peaks) == 251L)
write_tsv(window_peaks, file.path(stage_dir,
  "ebf1_pm2mb_bdev_peak_universe.tsv"))

source_bins <- fread(file.path(support_dir, "source_bin_cell_counts.tsv"))
shared_bins <- source_bins[retained_run == TRUE & wt_cells >= 10L &
  ko_cells >= 10L & region_id %in%
    c("flt3l_wt_peak1_run1", "il7_wt_peak1_run1")]
setorder(shared_bins, target, bin)
shared_regions <- shared_bins[, .(first_bin = min(bin), last_bin = max(bin),
  pt_left = min(bin_left), pt_right = max(bin_right), n_bins = .N,
  minimum_wt_cells_per_bin = min(wt_cells),
  minimum_ko_cells_per_bin = min(ko_cells),
  expected_wt_cells = sum(wt_cells), expected_ko_cells = sum(ko_cells)),
  by = .(target, region_id)]
stopifnot(nrow(shared_regions) == 2L,
  identical(shared_regions[target == "Flt3l", c(first_bin, last_bin)],
    c(14L, 25L)),
  identical(shared_regions[target == "Il7", c(first_bin, last_bin)],
    c(22L, 30L)),
  shared_bins[, all(diff(bin) == 1L), by = target][, all(V1)])
write_tsv(shared_bins, file.path(stage_dir,
  "shared_wt_ko_supported_source_bins.tsv"))
write_tsv(shared_regions, file.path(stage_dir,
  "shared_wt_ko_supported_regions.tsv"))

query <- fread(cmd = paste("zcat", shQuote(query_path)))[
  current_accepted == TRUE & trajectory_support >= 0.25 & is.finite(mapped_pt)]
assignments <- rbindlist(lapply(seq_len(nrow(shared_bins)), function(i) {
  b <- shared_bins[i]
  query[knockout_target == b$target & mapped_pt >= b$bin_left &
    mapped_pt < b$bin_right,
    .(cell_id, source_cell_id, sample_id, genotype, knockout_target,
      mapped_pt, trajectory_support, bin = b$bin, bin_left = b$bin_left,
      bin_right = b$bin_right, region_id = b$region_id)]
}))
setorder(assignments, knockout_target, genotype, mapped_pt, cell_id)
assignment_counts <- assignments[, .(n_cells = .N),
  by = .(knockout_target, region_id, bin, bin_left, bin_right, genotype)]
count_wide <- dcast(assignment_counts,
  knockout_target + region_id + bin + bin_left + bin_right ~ genotype,
  value.var = "n_cells", fill = 0L)
expected <- shared_bins[, .(knockout_target = target, region_id, bin,
  expected_WT = wt_cells, expected_KO = ko_cells)]
count_wide <- merge(count_wide, expected,
  by = c("knockout_target", "region_id", "bin"), sort = FALSE)
stopifnot(all(count_wide$WT == count_wide$expected_WT),
  all(count_wide$KO == count_wide$expected_KO),
  all(count_wide$WT >= 10L), all(count_wide$KO >= 10L),
  !anyDuplicated(assignments$cell_id))
write_tsv(count_wide, file.path(stage_dir,
  "shared_region_bin_cell_counts.tsv"))
fwrite(assignments, file.path(stage_dir,
  "shared_region_cell_assignments.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_integrated_ko_mapping_v1"),
  identical(manifest$current_state$fingerprint, expected_fp),
  all(assignments$cell_id %in% wb_cells(bundle, "set:ko_query")),
  all(window_peaks$bdev_feature_id %in%
    wb_feature_ids(bundle, "ATAC_bdev_mm10_peaks")))
write_json(wb_op_interface("matrix_log_normalize", bundle = bundle),
  file.path(stage_dir, "matrix_log_normalize_interface.json"))
write_json(wb_op_interface("de_wilcox", bundle = bundle),
  file.path(stage_dir, "de_wilcox_interface.json"))

session_id <- "ebf1_pm2mb_shared_support_da_v1"
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) wb_session_resume(bundle, workspace) else
  wb_session_create(bundle, workspace_root = session_root,
    session_id = session_id)
stopifnot(is.null(session$durable$promotion_receipt))
find_run <- function(session, marker, operation) {
  ids <- names(session$runs)[vapply(session$runs, function(x) {
    note <- if (is.null(x$notes)) "" else as.character(x$notes)
    identical(x$op_name, operation) && startsWith(note, paste0("[", marker, "]"))
  }, logical(1L))]
  stopifnot(length(ids) <= 1L)
  if (length(ids)) session$runs[[ids[[1L]]]] else NULL
}

normalized_id <- "ATAC_bdev_mm10_peaks_log1p_cp10k_ebf1_pm2mb_shared_v1"
normal_run <- find_run(session, "normalize_full_atac", "matrix_log_normalize")
if (is.null(normal_run)) {
  log_event("Staging full-matrix log1p-CP10K normalization")
  result <- wb_session_run_op(session, bundle, "matrix_log_normalize",
    inputs = list(matrix = "ATAC_bdev_mm10_peaks", cell_set = "all",
      min_cells = 10L, output_matrix_id = normalized_id,
      scale_factor = 10000, overwrite = FALSE),
    notes = "[normalize_full_atac] Full 97,788-peak log1p-CP10K; staged only; no promotion.")
  session <- result$session
  normal_run <- result$run
} else log_event("Reusing staged full-matrix normalization")
stopifnot("matrix" %in% names(normal_run$outputs))

de_list <- list()
for (target in c("Flt3l", "Il7")) {
  a <- assignments[knockout_target == target]
  wt_cells <- as.character(a[genotype == "WT", cell_id])
  ko_cells <- as.character(a[genotype == "KO", cell_id])
  marker <- paste0("de_", tolower(target), "_shared_min10")
  de_run <- find_run(session, marker, "de_wilcox")
  if (is.null(de_run)) {
    log_event("Testing", nrow(window_peaks), "Ebf1-window peaks for", target,
      "WT", length(wt_cells), "KO", length(ko_cells))
    result <- wb_session_run_op(session, bundle, "de_wilcox",
      inputs = list(matrix = normalized_id, cell_set = "set:ko_query",
        min_cells = 10L, group_a = wt_cells, group_b = ko_cells,
        features = as.character(window_peaks$bdev_feature_id),
        min_pct = 0, p_adjust_method = "BH"),
      notes = paste0("[", marker, "] Group A=WT; group B=KO; accepted mapped cells with trajectory_support>=0.25; pooled only across consecutive original bins having >=10 cells in both genotypes; Ebf1 TSS +/-2 Mb; log1p-CP10K full-peak normalization; descriptive cell-level Wilcoxon; staged only; no promotion."))
    session <- result$session
    de_run <- result$run
  } else log_event("Reusing staged DA for", target)
  de <- as.data.table(de_run$outputs$de_table$table)
  stopifnot(nrow(de) == nrow(window_peaks),
    setequal(de$feature_id, window_peaks$bdev_feature_id))
  de <- merge(de, window_peaks[, .(feature_id = bdev_feature_id,
    chromosome = mm10_seqname, peak_start_1based = mm10_start_1based,
    peak_end_1based = mm10_end_1based, peak_midpoint_1based,
    midpoint_relative_to_tss_bp, interval_distance_to_tss_bp)],
    by = "feature_id", all.x = TRUE)
  target_value <- target
  reg <- shared_regions[target == target_value]
  de[, `:=`(target = target, region_id = reg$region_id,
    first_bin = reg$first_bin, last_bin = reg$last_bin,
    pt_left = reg$pt_left, pt_right = reg$pt_right,
    n_wt_cells = length(wt_cells), n_ko_cells = length(ko_cells),
    detection_fraction_difference_wt_minus_ko = pct_a - pct_b,
    direction = fifelse(avg_log1p_fc > 0, "WT_higher",
      fifelse(avg_log1p_fc < 0, "KO_higher", "equal")),
    significant_bh05 = is.finite(p_val_adj) & p_val_adj < 0.05)]
  de <- de[order(p_val_adj, -abs(avg_log1p_fc), feature_id)]
  write_tsv(de, file.path(stage_dir, paste0(tolower(target),
    "_ebf1_pm2mb_shared_min10_wt_vs_ko_da.tsv")))
  de_list[[target]] <- de
}

combined <- rbindlist(de_list, use.names = TRUE)
write_tsv(combined, file.path(stage_dir,
  "all_targets_ebf1_pm2mb_shared_min10_wt_vs_ko_da.tsv"))
significant <- combined[significant_bh05 == TRUE]
write_tsv(significant, file.path(stage_dir,
  "significant_ebf1_pm2mb_shared_min10_peaks_bh05.tsv"))
summary <- combined[, .(tested_peaks = .N,
  finite_tests = sum(is.finite(p_val)), significant_bh05 = sum(significant_bh05),
  significant_wt_higher = sum(significant_bh05 & direction == "WT_higher"),
  significant_ko_higher = sum(significant_bh05 & direction == "KO_higher"),
  minimum_p = min(p_val, na.rm = TRUE),
  minimum_bh_adjusted_p = min(p_val_adj, na.rm = TRUE)),
  by = .(target, region_id, first_bin, last_bin, pt_left, pt_right,
    n_wt_cells, n_ko_cells)]
write_tsv(summary, file.path(stage_dir, "differential_accessibility_summary.tsv"))
top <- combined[direction != "equal",
  head(.SD[order(p_val_adj, -abs(avg_log1p_fc))], 20L),
  by = .(target, direction)]
write_tsv(top, file.path(stage_dir, "top20_peaks_per_target_direction.tsv"))

combined[, plot_status := fifelse(significant_bh05 & direction == "WT_higher",
  "WT higher (BH < 0.05)", fifelse(significant_bh05 & direction == "KO_higher",
    "KO higher (BH < 0.05)", "Not significant"))]
combined[, target_label := factor(target, levels = c("Flt3l", "Il7"),
  labels = c("Flt3L", "IL-7"))]
cols <- c("WT higher (BH < 0.05)" = "#C84545",
  "KO higher (BH < 0.05)" = "#3B77A8", "Not significant" = "grey72")
p_location <- ggplot(combined,
    aes(midpoint_relative_to_tss_bp / 1e6,
      detection_fraction_difference_wt_minus_ko, colour = plot_status)) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 1.8, alpha = 0.85) +
  facet_wrap(~target_label, ncol = 1) + scale_colour_manual(values = cols) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank()) +
  labs(x = "Peak midpoint relative to Ebf1 TSS (Mb)",
    y = "Accessible-cell fraction: WT - KO", colour = NULL,
    title = "WT versus KO accessibility within Ebf1 TSS ±2 Mb",
    subtitle = "Only consecutive pseudotime bins with at least 10 WT and 10 KO cells")
ggsave(file.path(stage_dir, "figures", "ebf1_pm2mb_da_by_position.png"),
  p_location, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "figures", "ebf1_pm2mb_da_by_position.pdf"),
  p_location, width = 10, height = 8, device = cairo_pdf)

combined[, neglog10_padj := -log10(pmax(p_val_adj, 1e-300))]
p_volcano <- ggplot(combined[is.finite(neglog10_padj)],
    aes(avg_log1p_fc, neglog10_padj, colour = plot_status)) +
  geom_point(size = 1.8, alpha = 0.85) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  facet_wrap(~target_label) + scale_colour_manual(values = cols) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank()) +
  labs(x = "scWorkbench avg_log1p_fc (positive = WT higher)",
    y = expression(-log[10](BH~adjusted~italic(P))), colour = NULL,
    title = "Ebf1-window differential accessibility")
ggsave(file.path(stage_dir, "figures", "ebf1_pm2mb_da_volcano.png"),
  p_volcano, width = 11, height = 6.5, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "figures", "ebf1_pm2mb_da_volcano.pdf"),
  p_volcano, width = 11, height = 6.5, device = cairo_pdf)

write_json(list(status = "complete_unpromoted_report_only",
  script_version = script_version, mapping_bundle_id = manifest$bundle_id,
  mapping_bundle_fingerprint = expected_fp,
  peak_mapping = peak_mapping_path,
  peak_mapping_sha256 = sha256_file(peak_mapping_path),
  mapped_cell_selection = "current_accepted and trajectory_support >= 0.25",
  shared_region_rule = "consecutive original bins with WT >=10 and KO >=10 cells",
  normalization = "native scWorkbench full-peak log1p-CP10K",
  test = "native scWorkbench de_wilcox; WT group A versus KO group B; BH within target",
  ebf1 = as.list(ebf1), shared_regions = split(shared_regions,
    shared_regions$target), summary = split(summary, summary$target),
  biological_replicates_per_genotype = 1L, inference = "descriptive",
  bundle_mutated = FALSE, promotion = FALSE),
  file.path(stage_dir, "diagnostics.json"))

region_lines <- apply(shared_regions, 1L, function(x) paste0("- `", x[["target"]],
  "`: bins ", x[["first_bin"]], "–", x[["last_bin"]], ", pt [",
  sprintf("%.3f", as.numeric(x[["pt_left"]])), ", ",
  sprintf("%.3f", as.numeric(x[["pt_right"]])), "); ",
  x[["expected_wt_cells"]], " WT and ", x[["expected_ko_cells"]], " KO cells."))
summary_lines <- apply(summary, 1L, function(x) paste0("- `", x[["target"]],
  "`: ", x[["significant_bh05"]], " of ", x[["tested_peaks"]],
  " peaks BH-significant (", x[["significant_wt_higher"]], " WT-higher; ",
  x[["significant_ko_higher"]], " KO-higher)."))
writeLines(c("# Ebf1 ±2 Mb differential accessibility in shared WT/KO-supported trajectory regions", "",
  "The Ebf1 TSS is chr11:44,617,317 in GRCm38/mm10. The tested universe contains 251 Bdev peaks overlapping ±2 Mb.", "",
  "## Shared pseudotime domains", "", region_lines, "",
  "Every included original bin contains at least 10 accepted mapped WT cells and 10 accepted mapped KO cells. Cells require trajectory support ≥0.25.", "",
  "## Differential accessibility", "", summary_lines, "",
  "Accessibility was normalized using native scWorkbench full-matrix log1p-CP10K, then tested with native cell-level Wilcoxon tests separately for Flt3L and IL-7. BH correction is within each 251-peak target test.", "",
  "Because each genotype/target has one biological sample, these are descriptive cell-level differences rather than replicate-level inference.", "",
  "- [Complete results](all_targets_ebf1_pm2mb_shared_min10_wt_vs_ko_da.tsv)",
  "- [BH-significant peaks](significant_ebf1_pm2mb_shared_min10_peaks_bh05.tsv)",
  "- [Position figure](figures/ebf1_pm2mb_da_by_position.png)",
  "- [Volcano figure](figures/ebf1_pm2mb_da_volcano.png)"),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
wb_close_bundle(bundle); bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed Ebf1 ±2 Mb shared-support DA at", output_dir, "\n")
