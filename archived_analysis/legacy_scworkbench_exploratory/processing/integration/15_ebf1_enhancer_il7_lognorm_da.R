#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(paste(
  "Usage: 15_ebf1_enhancer_il7_lognorm_da.R MAPPING_BUNDLE PEAK_MAPPING",
  "ENRICHED_PT_REPORT RAW_DA_REPORT OUTPUT_DIR SESSION_ROOT",
  "EXPECTED_BUNDLE_FP SCRIPT_VERSION"
), call. = FALSE)

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
peak_mapping_path <- normalizePath(args[[2L]], mustWork = TRUE)
pt_report_dir <- normalizePath(args[[3L]], mustWork = TRUE)
raw_report_dir <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
session_root <- normalizePath(args[[6L]], mustWork = FALSE)
expected_fingerprint <- args[[7L]]
script_version <- args[[8L]]

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)
dir.create(session_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
log_event <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(..., collapse = " "), "\n")

ebf1 <- data.table(gene = "Ebf1", gene_id = "ENSMUSG00000057098",
  genome = "GRCm38/mm10", chromosome = "chr11", strand = "+",
  tss_1based = 44617317L, flank_bp = 3000000L)
ebf1[, `:=`(window_start_1based = tss_1based - flank_bp,
  window_end_1based = tss_1based + flank_bp)]
write_tsv(ebf1, file.path(stage_dir, "ebf1_tss_window_definition.tsv"))

peak_mapping <- fread(cmd = paste("zcat", shQuote(peak_mapping_path)))
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
stopifnot(nrow(window_peaks) == 259L,
  !anyDuplicated(window_peaks$bdev_feature_id))
write_tsv(window_peaks, file.path(stage_dir,
  "ebf1_pm3mb_bdev_peak_universe.tsv"))

region_path <- file.path(pt_report_dir, "wt_enriched_pt_regions.tsv")
assignment_path <- file.path(pt_report_dir, "cell_region_assignments.tsv.gz")
raw_de_path <- file.path(raw_report_dir,
  "all_il7_regions_ebf1_pm3mb_wt_vs_ko_differential_accessibility.tsv")
stopifnot(file.exists(region_path), file.exists(assignment_path),
  file.exists(raw_de_path))
regions <- fread(region_path)[target == "Il7"]
assignments <- fread(cmd = paste("zcat", shQuote(assignment_path)))[
  knockout_target == "Il7"]
setorder(regions, peak_number)
stopifnot(identical(regions$region_id, c("Il7_peak1", "Il7_peak2")),
  nrow(assignments) == 9406L)
write_tsv(regions, file.path(stage_dir, "il7_enriched_pt_regions.tsv"))
fwrite(assignments, file.path(stage_dir,
  "il7_enriched_pt_cell_assignments.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")
region_counts <- assignments[, .(n_cells = .N, n_samples = uniqueN(sample_id)),
  by = .(region_id, peak_number, genotype, pt_left, pt_right)]
setorder(region_counts, peak_number, genotype)
write_tsv(region_counts, file.path(stage_dir, "il7_region_cell_counts.tsv"))

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint),
  all(window_peaks$bdev_feature_id %in%
    wb_feature_ids(bundle, "ATAC_bdev_mm10_peaks")),
  all(assignments$cell_id %in% wb_cells(bundle, "set:ko_query")))

write_json(wb_op_interface("matrix_log_normalize", bundle = bundle),
  file.path(stage_dir, "matrix_log_normalize_interface.json"))
write_json(wb_op_interface("de_wilcox", bundle = bundle),
  file.path(stage_dir, "de_wilcox_interface.json"))
session_id <- "ebf1_pm3mb_il7_enriched_pt_lognorm_da_v1"
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) wb_session_resume(bundle, workspace) else
  wb_session_create(bundle, workspace_root = session_root,
    session_id = session_id)
stopifnot(is.null(session$durable$promotion_receipt))

find_run <- function(session, marker, operation) {
  ids <- names(session$runs)[vapply(session$runs, function(x) {
    note <- if (is.null(x$notes)) "" else as.character(x$notes)
    identical(x$op_name, operation) && startsWith(note, paste0("[", marker, "]"))
  }, logical(1))]
  stopifnot(length(ids) <= 1L)
  if (length(ids)) session$runs[[ids[[1L]]]] else NULL
}

normalized_matrix_id <- "ATAC_bdev_mm10_peaks_log1p_cp10k_ebf1_da_v1"
normal_run <- find_run(session, "normalize_atac_bdev_peaks", "matrix_log_normalize")
if (is.null(normal_run)) {
  log_event("Staging log1p-CP10K normalization of the full Bdev peak matrix")
  result <- wb_session_run_op(session, bundle, "matrix_log_normalize",
    inputs = list(matrix = "ATAC_bdev_mm10_peaks", cell_set = "all",
      min_cells = 10L, output_matrix_id = normalized_matrix_id,
      scale_factor = 10000, overwrite = FALSE),
    notes = paste0("[normalize_atac_bdev_peaks] Full-bundle log1p(count / ",
      "cell total across 97,788 Bdev peaks * 10,000); staged only; no promotion."))
  session <- result$session
  normal_run <- result$run
} else log_event("Reusing staged normalized ATAC matrix")
stopifnot("matrix" %in% names(normal_run$outputs))

de_tables <- list()
for (i in seq_len(nrow(regions))) {
  reg <- regions[i]
  region_cells <- assignments[region_id == reg$region_id]
  wt_cells <- as.character(region_cells[genotype == "WT", cell_id])
  ko_cells <- as.character(region_cells[genotype == "KO", cell_id])
  marker <- paste0("de_lognorm_", tolower(reg$region_id))
  de_run <- find_run(session, marker, "de_wilcox")
  if (is.null(de_run)) {
    log_event("Testing normalized accessibility for", reg$region_id,
      "WT", length(wt_cells), "KO", length(ko_cells))
    result <- wb_session_run_op(session, bundle, "de_wilcox",
      inputs = list(matrix = normalized_matrix_id,
        cell_set = "set:ko_query", min_cells = 10L,
        group_a = wt_cells, group_b = ko_cells,
        features = as.character(window_peaks$bdev_feature_id),
        min_pct = 0, p_adjust_method = "BH"),
      notes = paste0("[", marker, "] Group A=IL7 WT; group B=IL7 KO; ",
        "exact saved ", reg$region_id, " cells; log1p-CP10K normalized ",
        "across full Bdev peak matrix; all 259 Ebf1-window peaks supplied; ",
        "descriptive cell-level test; staged only; no promotion."))
    session <- result$session
    de_run <- result$run
  } else log_event("Reusing normalized DA run for", reg$region_id)
  de <- as.data.table(de_run$outputs$de_table$table)
  stopifnot(nrow(de) == 259L,
    setequal(de$feature_id, window_peaks$bdev_feature_id))
  de <- merge(de, window_peaks[, .(feature_id = bdev_feature_id,
    chromosome = mm10_seqname, peak_start_1based = mm10_start_1based,
    peak_end_1based = mm10_end_1based, peak_midpoint_1based,
    midpoint_relative_to_tss_bp, interval_distance_to_tss_bp)],
    by = "feature_id", all.x = TRUE)
  de[, `:=`(region_id = reg$region_id, peak_number = reg$peak_number,
    pt_left = reg$pt_left, pt_right = reg$pt_right,
    n_wt_cells = length(wt_cells), n_ko_cells = length(ko_cells),
    direction = fifelse(avg_log1p_fc > 0, "WT_higher",
      fifelse(avg_log1p_fc < 0, "KO_higher", "equal")),
    significant_bh05 = is.finite(p_val_adj) & p_val_adj < 0.05)]
  de <- de[order(p_val_adj, -abs(avg_log1p_fc), feature_id)]
  write_tsv(de, file.path(stage_dir, paste0(tolower(reg$region_id),
    "_ebf1_pm3mb_lognorm_wt_vs_ko_da.tsv")))
  de_tables[[reg$region_id]] <- de
}

normalized_de <- rbindlist(de_tables, use.names = TRUE)
write_tsv(normalized_de, file.path(stage_dir,
  "all_il7_regions_ebf1_pm3mb_lognorm_wt_vs_ko_da.tsv"))
significant <- normalized_de[significant_bh05 == TRUE]
write_tsv(significant, file.path(stage_dir,
  "all_il7_regions_ebf1_pm3mb_lognorm_significant_peaks_bh05.tsv"))
top <- normalized_de[direction != "equal",
  head(.SD[order(p_val_adj, -abs(avg_log1p_fc))], 20L),
  by = .(region_id, direction)]
write_tsv(top, file.path(stage_dir,
  "top20_lognorm_peaks_per_region_direction.tsv"))

raw <- fread(raw_de_path)
stopifnot(nrow(raw) == 518L)
comparison <- merge(
  normalized_de[, .(region_id, feature_id,
    norm_avg_log1p_fc = avg_log1p_fc, norm_p_val = p_val,
    norm_p_val_adj = p_val_adj, norm_significant_bh05 = significant_bh05,
    norm_direction = direction)],
  raw[, .(region_id, feature_id, raw_avg_log1p_fc = avg_log1p_fc,
    raw_p_val = p_val, raw_p_val_adj = p_val_adj,
    raw_significant_bh05 = significant_bh05, raw_direction = direction)],
  by = c("region_id", "feature_id"), all = TRUE)
comparison[, significance_change := fifelse(
  norm_significant_bh05 & raw_significant_bh05, "significant_both",
  fifelse(norm_significant_bh05, "normalized_only",
    fifelse(raw_significant_bh05, "raw_only", "neither")))]
write_tsv(comparison, file.path(stage_dir,
  "raw_vs_lognorm_differential_accessibility.tsv"))

summary_table <- normalized_de[, .(
  peaks_supplied = .N, peaks_with_finite_test = sum(is.finite(p_val)),
  significant_bh05 = sum(significant_bh05),
  significant_wt_higher = sum(significant_bh05 & direction == "WT_higher"),
  significant_ko_higher = sum(significant_bh05 & direction == "KO_higher")),
  by = .(region_id, peak_number, pt_left, pt_right, n_wt_cells, n_ko_cells)]
raw_summary <- comparison[, .(
  raw_significant = sum(raw_significant_bh05, na.rm = TRUE),
  significant_both = sum(significance_change == "significant_both"),
  raw_only = sum(significance_change == "raw_only"),
  normalized_only = sum(significance_change == "normalized_only")),
  by = region_id]
summary_table <- merge(summary_table, raw_summary, by = "region_id")
setorder(summary_table, peak_number)
write_tsv(summary_table, file.path(stage_dir,
  "lognorm_differential_accessibility_summary.tsv"))

normalized_de[, `:=`(
  neglog10_padj = -log10(pmax(p_val_adj, 1e-300)),
  plot_status = fifelse(significant_bh05 & avg_log1p_fc > 0,
    "WT higher (BH < 0.05)", fifelse(significant_bh05 & avg_log1p_fc < 0,
      "KO higher (BH < 0.05)", "Not significant")),
  region_label = factor(region_id, levels = c("Il7_peak1", "Il7_peak2"),
    labels = c("IL-7 peak 1", "IL-7 peak 2")))]
plot_colours <- c("WT higher (BH < 0.05)" = "#C84545",
  "KO higher (BH < 0.05)" = "#3B77A8", "Not significant" = "grey72")
volcano <- ggplot(normalized_de[is.finite(neglog10_padj)],
    aes(avg_log1p_fc, neglog10_padj, colour = plot_status)) +
  geom_point(size = 2, alpha = 0.85) + geom_vline(xintercept = 0,
    colour = "grey35", linewidth = 0.4) +
  facet_wrap(~region_label, ncol = 2) +
  scale_colour_manual(values = plot_colours) +
  labs(x = "Normalized avg_log1p_fc (positive = WT higher)",
    y = expression(-log[10](BH~adjusted~italic(P))), colour = NULL,
    title = "Depth-normalized Ebf1-window differential accessibility",
    subtitle = "Per-cell log1p-CP10K normalization across the full Bdev peak matrix") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank())
ggsave(file.path(stage_dir, "figures", "lognorm_da_volcano.pdf"), volcano,
  width = 10, height = 5.5, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "lognorm_da_volcano.png"), volcano,
  width = 10, height = 5.5, dpi = 300)

comparison[, region_label := factor(region_id,
  levels = c("Il7_peak1", "Il7_peak2"),
  labels = c("IL-7 peak 1", "IL-7 peak 2"))]
effect_compare <- ggplot(comparison,
    aes(raw_avg_log1p_fc, norm_avg_log1p_fc, colour = significance_change)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_point(size = 1.8, alpha = 0.8) + facet_wrap(~region_label) +
  scale_colour_manual(values = c(significant_both = "#7B3294",
    normalized_only = "#008837", raw_only = "#D95F02", neither = "grey72")) +
  labs(x = "Raw-count avg_log1p_fc", y = "Normalized avg_log1p_fc",
    colour = NULL, title = "Effect estimates before and after depth normalization") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank())
ggsave(file.path(stage_dir, "figures", "raw_vs_lognorm_effects.pdf"),
  effect_compare, width = 9, height = 5.5, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "raw_vs_lognorm_effects.png"),
  effect_compare, width = 9, height = 5.5, dpi = 300)

write_json(wb_session_dag(session, bundle = bundle, include_artifacts = TRUE,
  include_active = TRUE), file.path(stage_dir, "native_session_dag.json"))
fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
checks <- list(exact_two_regions = nrow(summary_table) == 2L,
  all_259_peaks_supplied_per_region = all(summary_table$peaks_supplied == 259L),
  normalized_matrix_staged = "matrix" %in% names(normal_run$outputs),
  complete_normalized_results = nrow(normalized_de) == 518L,
  complete_raw_comparison = nrow(comparison) == 518L,
  significant_table_matches = nrow(significant) ==
    sum(normalized_de$p_val_adj < 0.05, na.rm = TRUE),
  figures_exist = all(file.exists(file.path(stage_dir, "figures", c(
    "lognorm_da_volcano.png", "raw_vs_lognorm_effects.png")))),
  parent_bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after), result_promotions = 0L)
stopifnot(all(unlist(checks[names(checks) != "result_promotions"])),
  checks$result_promotions == 0L,
  is.null(session$durable$promotion_receipt))
write_json(list(status = "passed", checks = checks),
  file.path(stage_dir, "checks.json"))

answer_lines <- summary_table[, sprintf(
  "- %s: %s normalized-significant peaks (%s WT-higher; %s KO-higher), versus %s with raw counts; %s overlap, %s raw-only, %s normalized-only.",
  region_id, significant_bh05, significant_wt_higher,
  significant_ko_higher, raw_significant, significant_both, raw_only,
  normalized_only)]
report <- c("# Depth-normalized Ebf1-window differential accessibility", "",
  "## Results", "", answer_lines, "",
  "## Method", "",
  "- Same exact IL-7 peak-1 and peak-2 cell assignments as the raw-count analysis.",
  "- Native scWorkbench `matrix_log_normalize`: log1p(count / per-cell total across 97,788 Bdev peaks × 10,000).",
  "- Native scWorkbench `de_wilcox`: WT versus KO separately in each trajectory region, all 259 Ebf1-window peaks supplied, `min_pct=0`, BH adjustment within region.",
  "- The normalized matrix and test results remain staged in an unpromoted durable session.", "",
  "## Inference limitation", "",
  "Each genotype has one biological sample. Normalization addresses cell-depth differences but does not create biological replication; P values remain descriptive cell-level evidence.", "",
  "![Normalized volcano](figures/lognorm_da_volcano.png)", "",
  "![Raw versus normalized effects](figures/raw_vs_lognorm_effects.png)", "",
  "## Outputs", "",
  "- `all_il7_regions_ebf1_pm3mb_lognorm_wt_vs_ko_da.tsv`: complete normalized results.",
  "- `all_il7_regions_ebf1_pm3mb_lognorm_significant_peaks_bh05.tsv`: normalized significant peaks.",
  "- `raw_vs_lognorm_differential_accessibility.tsv`: direct peak-wise comparison.",
  "- `top20_lognorm_peaks_per_region_direction.tsv`: normalized ranking.")
writeLines(report, file.path(stage_dir, "REPORT.md"))
html_answers <- paste0("<li>", answer_lines, "</li>", collapse = "")
html <- c("<!doctype html><html><head><meta charset='utf-8'>",
  "<title>Normalized Ebf1 differential accessibility</title>",
  "<style>body{font-family:Arial,sans-serif;max-width:1100px;margin:30px auto;line-height:1.45;color:#222}img{max-width:100%}.warning{background:#fff3cd;border-left:5px solid #d79b00;padding:12px}</style></head><body>",
  "<h1>Depth-normalized Ebf1-window differential accessibility</h1>",
  paste0("<h2>Results</h2><ul>", html_answers, "</ul>"),
  "<div class='warning'>One sample per genotype: normalization addresses cell depth, but results remain descriptive cell-level evidence rather than replicate-level inference.</div>",
  "<h2>Method</h2><p>Native scWorkbench log1p-CP10K normalization across all 97,788 Bdev peaks per cell, followed by region-specific native Wilcoxon tests over all 259 Ebf1-window peaks.</p>",
  "<img src='figures/lognorm_da_volcano.png'>",
  "<img src='figures/raw_vs_lognorm_effects.png'>",
  "<h2>Downloads</h2><ul><li><a href='all_il7_regions_ebf1_pm3mb_lognorm_wt_vs_ko_da.tsv'>Complete normalized results</a></li><li><a href='all_il7_regions_ebf1_pm3mb_lognorm_significant_peaks_bh05.tsv'>Significant peaks</a></li><li><a href='raw_vs_lognorm_differential_accessibility.tsv'>Raw versus normalized comparison</a></li></ul>",
  "</body></html>")
writeLines(html, file.path(stage_dir, "REPORT.html"))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
input_manifest <- data.table(role = c("mapping_bundle_manifest",
  "accepted_peak_mapping", "enriched_pt_regions", "cell_assignments",
  "raw_da_results", "analysis_script"), path = normalizePath(c(
    file.path(bundle_dir, "manifest.json"), peak_mapping_path, region_path,
    assignment_path, raw_de_path, script_path), mustWork = TRUE))
input_manifest[, sha256 := vapply(path, sha256_file, character(1))]
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
write_json(list(status = "complete_unpromoted_report_only",
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  report = file.path(output_dir, "REPORT.html"), bundle = bundle_dir,
  bundle_fingerprint = fingerprint_before,
  bundle_fingerprint_unchanged = TRUE,
  source_matrix = "ATAC_bdev_mm10_peaks",
  normalized_matrix = normalized_matrix_id,
  normalization = "log1p_cp10k_across_full_97788_peak_matrix_per_cell",
  scale_factor = 10000L, tested_peaks_per_region = 259L,
  regions = split(summary_table, seq_len(nrow(summary_table))),
  session_workspace = workspace, promotions = 0L, checks = checks),
  file.path(stage_dir, "completion.json"))

if (!file.rename(stage_dir, output_dir))
  stop("Failed to finalize output directory")
log_event("Completed depth-normalized Ebf1 DA report", output_dir)
