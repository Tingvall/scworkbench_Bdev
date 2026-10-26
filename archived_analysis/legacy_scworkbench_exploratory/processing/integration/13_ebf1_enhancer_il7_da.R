#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 13_ebf1_enhancer_il7_da.R MAPPING_BUNDLE PEAK_MAPPING",
  "ENRICHED_PT_REPORT OUTPUT_DIR SESSION_ROOT EXPECTED_BUNDLE_FP SCRIPT_VERSION"
), call. = FALSE)

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
peak_mapping_path <- normalizePath(args[[2L]], mustWork = TRUE)
pt_report_dir <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
session_root <- normalizePath(args[[5L]], mustWork = FALSE)
expected_fingerprint <- args[[6L]]
script_version <- args[[7L]]

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

# The exact Bdev Signac annotation is GRCm38. It contains flattened Ebf1
# exon/CDS/UTR/gap records but no transcript/gene row. Because Ebf1 is on the
# positive strand, the minimum annotated coordinate is its TSS.
ebf1 <- data.table(gene = "Ebf1", gene_id = "ENSMUSG00000057098",
  genome = "GRCm38/mm10", chromosome = "chr11", strand = "+",
  tss_1based = 44617317L, flank_bp = 3000000L)
ebf1[, `:=`(window_start_1based = tss_1based - flank_bp,
  window_end_1based = tss_1based + flank_bp)]
write_tsv(ebf1, file.path(stage_dir, "ebf1_tss_window_definition.tsv"))

peak_mapping <- fread(cmd = paste("zcat", shQuote(peak_mapping_path)))
stopifnot(all(c("bdev_feature_id", "mm10_seqname", "mm10_start_1based",
  "mm10_end_1based", "status") %in% names(peak_mapping)),
  all(peak_mapping$status == "accepted"),
  !anyDuplicated(peak_mapping$bdev_feature_id))
window_peaks <- peak_mapping[mm10_seqname == ebf1$chromosome &
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
  all(window_peaks$mm10_seqname == "chr11"))
write_tsv(window_peaks, file.path(stage_dir,
  "ebf1_pm3mb_bdev_peak_universe.tsv"))

# Reuse the exact IL7 pseudotime-region definitions and cell assignments from
# the accepted enrichment-region analysis. Keep its two peaks separate.
region_path <- file.path(pt_report_dir, "wt_enriched_pt_regions.tsv")
assignment_path <- file.path(pt_report_dir, "cell_region_assignments.tsv.gz")
stopifnot(file.exists(region_path), file.exists(assignment_path))
regions <- fread(region_path)[target == "Il7"]
assignments <- fread(cmd = paste("zcat", shQuote(assignment_path)))[
  knockout_target == "Il7"]
setorder(regions, peak_number)
stopifnot(nrow(regions) == 2L,
  identical(regions$region_id, c("Il7_peak1", "Il7_peak2")),
  setequal(unique(assignments$region_id), regions$region_id),
  !anyDuplicated(assignments[, .(source_cell_id, region_id)]),
  all(assignments$mapped_pt >= assignments$pt_left &
    assignments$mapped_pt < assignments$pt_right))
write_tsv(regions, file.path(stage_dir, "il7_enriched_pt_regions.tsv"))
fwrite(assignments, file.path(stage_dir, "il7_enriched_pt_cell_assignments.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip")
region_counts <- assignments[, .(n_cells = .N, n_samples = uniqueN(sample_id),
  samples = paste(sort(unique(sample_id)), collapse = ";")),
  by = .(region_id, peak_number, genotype, pt_left, pt_right)]
setorder(region_counts, peak_number, genotype)
stopifnot(nrow(region_counts) == 4L, all(region_counts$n_samples == 1L),
  identical(region_counts[genotype == "WT"]$n_cells, c(4016L, 2021L)),
  identical(region_counts[genotype == "KO"]$n_cells, c(2009L, 1360L)))
write_tsv(region_counts, file.path(stage_dir, "il7_region_cell_counts.tsv"))

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint))
matrix_info <- wb_matrix_info(bundle, "ATAC_bdev_mm10_peaks")
stopifnot(nrow(matrix_info) == 1L, matrix_info$feature_type == "peak",
  matrix_info$layer == "counts",
  all(window_peaks$bdev_feature_id %in%
    wb_feature_ids(bundle, "ATAC_bdev_mm10_peaks")))
query_cells <- as.character(wb_cells(bundle, "set:ko_query"))
stopifnot(all(assignments$cell_id %in% query_cells))

write_json(wb_op_interface("de_wilcox", bundle = bundle),
  file.path(stage_dir, "de_wilcox_interface.json"))
session_id <- "ebf1_pm3mb_il7_enriched_pt_da_v1"
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) wb_session_resume(bundle, workspace) else
  wb_session_create(bundle, workspace_root = session_root,
    session_id = session_id)
if (!is.null(session$durable$promotion_receipt))
  stop("The Ebf1 report-only DA session was unexpectedly promoted")

de_tables <- list()
for (i in seq_len(nrow(regions))) {
  reg <- regions[i]
  region_cells <- assignments[region_id == reg$region_id]
  wt_cells <- as.character(region_cells[genotype == "WT", cell_id])
  ko_cells <- as.character(region_cells[genotype == "KO", cell_id])
  marker <- paste0("[ebf1_pm3mb_", tolower(reg$region_id),
    "_il7_wt_vs_ko] ")
  run_ids <- names(session$runs)[vapply(session$runs, function(x) {
    note <- if (is.null(x$notes)) "" else as.character(x$notes)
    identical(x$op_name, "de_wilcox") && startsWith(note, marker)
  }, logical(1))]
  stopifnot(length(run_ids) <= 1L)
  if (length(run_ids)) {
    log_event("Reusing completed run for", reg$region_id)
    de_run <- session$runs[[run_ids[[1L]]]]
  } else {
    log_event("Testing", nrow(window_peaks), "Ebf1-window peaks in",
      reg$region_id, ":", length(wt_cells), "WT and", length(ko_cells), "KO cells")
    result <- wb_session_run_op(session, bundle, "de_wilcox",
      inputs = list(matrix = "ATAC_bdev_mm10_peaks",
        cell_set = "set:ko_query", min_cells = 10L,
        group_a = wt_cells, group_b = ko_cells,
        features = as.character(window_peaks$bdev_feature_id),
        min_pct = 0, p_adjust_method = "BH"),
      notes = paste0(marker, "Group A=IL7 WT; group B=IL7 KO; exact ",
        reg$region_id, " saved cell assignments; pt [", reg$pt_left, ", ",
        reg$pt_right, "); all 259 Bdev mm10 peaks overlapping Ebf1 TSS ",
        "+/-3 Mb; raw fragment-count Wilcoxon; no feature prefilter; ",
        "descriptive cell-level test; staged only and not promoted."))
    session <- result$session
    de_run <- result$run
  }
  de <- as.data.table(de_run$outputs$de_table$table)
  stopifnot(nrow(de) == nrow(window_peaks),
    setequal(de$feature_id, window_peaks$bdev_feature_id))
  de <- merge(de, window_peaks[, .(feature_id = bdev_feature_id,
    chromosome = mm10_seqname, peak_start_1based = mm10_start_1based,
    peak_end_1based = mm10_end_1based, peak_midpoint_1based,
    midpoint_relative_to_tss_bp, interval_distance_to_tss_bp)],
    by = "feature_id", all.x = TRUE)
  de[, `:=`(region_id = reg$region_id, peak_number = reg$peak_number,
    pt_left = reg$pt_left, pt_right = reg$pt_right,
    n_wt_cells = length(wt_cells), n_ko_cells = length(ko_cells),
    detection_fraction_difference_wt_minus_ko = pct_a - pct_b,
    direction = fifelse(avg_log1p_fc > 0, "WT_higher",
      fifelse(avg_log1p_fc < 0, "KO_higher", "equal")),
    significant_bh05 = is.finite(p_val_adj) & p_val_adj < 0.05)]
  de <- de[order(p_val_adj, -abs(avg_log1p_fc), feature_id)]
  write_tsv(de, file.path(stage_dir, paste0(tolower(reg$region_id),
    "_ebf1_pm3mb_wt_vs_ko_differential_accessibility.tsv")))
  de_tables[[reg$region_id]] <- de
}

combined_de <- rbindlist(de_tables, use.names = TRUE)
write_tsv(combined_de, file.path(stage_dir,
  "all_il7_regions_ebf1_pm3mb_wt_vs_ko_differential_accessibility.tsv"))
significant <- combined_de[significant_bh05 == TRUE]
write_tsv(significant, file.path(stage_dir,
  "all_il7_regions_ebf1_pm3mb_significant_peaks_bh05.tsv"))
top_peaks <- combined_de[direction != "equal",
  head(.SD[order(p_val_adj, -abs(avg_log1p_fc))], 20L),
  by = .(region_id, direction)]
write_tsv(top_peaks, file.path(stage_dir,
  "top20_peaks_per_region_direction.tsv"))
summary_table <- combined_de[, .(tested_peaks = .N,
  peaks_with_finite_test = sum(is.finite(p_val)),
  peaks_without_test_statistic = sum(!is.finite(p_val)),
  significant_bh05 = sum(significant_bh05),
  significant_wt_higher = sum(significant_bh05 & direction == "WT_higher"),
  significant_ko_higher = sum(significant_bh05 & direction == "KO_higher"),
  minimum_p_value = min(p_val, na.rm = TRUE),
  minimum_bh_adjusted_p = min(p_val_adj, na.rm = TRUE)),
  by = .(region_id, peak_number, pt_left, pt_right, n_wt_cells, n_ko_cells)]
write_tsv(summary_table, file.path(stage_dir,
  "differential_accessibility_summary.tsv"))

combined_de[, neglog10_padj := -log10(pmax(p_val_adj, 1e-300))]
combined_de[, plot_status := fifelse(significant_bh05 & avg_log1p_fc > 0,
  "WT higher (BH < 0.05)", fifelse(significant_bh05 & avg_log1p_fc < 0,
    "KO higher (BH < 0.05)", "Not significant"))]
combined_de[, region_label := factor(region_id,
  levels = c("Il7_peak1", "Il7_peak2"), labels = c("IL-7 peak 1", "IL-7 peak 2"))]
plot_colours <- c("WT higher (BH < 0.05)" = "#C84545",
  "KO higher (BH < 0.05)" = "#3B77A8", "Not significant" = "grey72")

location_plot <- ggplot(combined_de,
    aes(midpoint_relative_to_tss_bp / 1e6,
      detection_fraction_difference_wt_minus_ko, colour = plot_status)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey35") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 2, alpha = 0.85) + facet_wrap(~region_label, ncol = 1) +
  scale_colour_manual(values = plot_colours) +
  labs(x = "Peak midpoint relative to Ebf1 TSS (Mb)",
    y = "Accessible-cell fraction: WT - KO", colour = NULL,
    title = "IL-7 WT versus KO accessibility around Ebf1",
    subtitle = "Exact WT-enriched pseudotime intervals are tested separately; dashed line marks the TSS") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank())
ggsave(file.path(stage_dir, "figures",
  "ebf1_pm3mb_accessibility_difference_by_position.pdf"), location_plot,
  width = 10, height = 8, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures",
  "ebf1_pm3mb_accessibility_difference_by_position.png"), location_plot,
  width = 10, height = 8, dpi = 300)

volcano <- ggplot(combined_de[is.finite(neglog10_padj)], aes(avg_log1p_fc, neglog10_padj,
    colour = plot_status)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey35") +
  facet_wrap(~region_label, ncol = 2) + scale_colour_manual(values = plot_colours) +
  labs(x = "scWorkbench avg_log1p_fc (positive = WT higher)",
    y = expression(-log[10](BH~adjusted~italic(P))), colour = NULL,
    title = "Ebf1-window differential accessibility",
    subtitle = "Native cell-level Wilcoxon test on raw Bdev peak fragment counts") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank())
ggsave(file.path(stage_dir, "figures", "ebf1_pm3mb_da_volcano.pdf"),
  volcano, width = 10, height = 5.5, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "ebf1_pm3mb_da_volcano.png"),
  volcano, width = 10, height = 5.5, dpi = 300)

write_json(wb_session_dag(session, bundle = bundle, include_artifacts = TRUE,
  include_active = TRUE), file.path(stage_dir, "native_session_dag.json"))
fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
checks <- list(
  exact_ebf1_tss_recorded = ebf1$tss_1based == 44617317L,
  window_is_plus_minus_3mb = ebf1$flank_bp == 3000000L,
  exact_two_il7_enriched_regions = identical(regions$region_id,
    c("Il7_peak1", "Il7_peak2")),
  exact_prior_cell_assignments_reused = nrow(assignments) == 9406L,
  all_overlapping_bdev_peaks_tested_per_region = nrow(combined_de) == 518L,
  complete_de_tables = all(combined_de[, .N, by = region_id]$N == 259L),
  exact_unestimable_peak_counts = identical(
    summary_table$peaks_without_test_statistic, c(4L, 1L)),
  significant_table_matches = nrow(significant) ==
    sum(combined_de$p_val_adj < 0.05, na.rm = TRUE),
  figures_exist = all(file.exists(file.path(stage_dir, "figures", c(
    "ebf1_pm3mb_accessibility_difference_by_position.png",
    "ebf1_pm3mb_da_volcano.png")))),
  parent_bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after), result_promotions = 0L)
stopifnot(all(unlist(checks[names(checks) != "result_promotions"])),
  checks$result_promotions == 0L, is.null(session$durable$promotion_receipt))
write_json(list(status = "passed", checks = checks),
  file.path(stage_dir, "checks.json"))

answer_lines <- summary_table[, sprintf(
  "- %s (PT %.3f to %.3f; %s WT and %s KO cells): %s significant among %s testable peaks (259 supplied; %s had no test statistic; WT higher: %s; KO higher: %s).",
  region_id, pt_left, pt_right, n_wt_cells, n_ko_cells, significant_bh05,
  peaks_with_finite_test, peaks_without_test_statistic,
  significant_wt_higher, significant_ko_higher)]
report_md <- c("# Ebf1 ±3 Mb differential accessibility in IL-7-enriched trajectory regions", "",
  "## Answer", "", answer_lines, "",
  "Significance is BH-adjusted P < 0.05 within each regional contrast; peaks with no test statistic are retained and reported as untestable.", "",
  "## Important inference limitation", "",
  paste("There is one IL-7 WT sample and one IL-7 KO sample. The Wilcoxon",
    "test treats cells as observations, so adjusted P values are descriptive",
    "cell-level evidence and are not biological-replicate-level genotype",
    "inference. Raw fragment-count ranks can also reflect per-cell sequencing depth."), "",
  "## Window and cell selection", "",
  "- Ebf1 TSS: GRCm38/mm10 chr11:44,617,317 (+ strand).",
  "- Window: chr11:41,617,317–47,617,317 (±3 Mb).",
  "- Tested universe: all 259 overlapping quantified Bdev peaks in each contrast.",
  "- Cells: exact `Il7_peak1` and `Il7_peak2` assignments from the prior WT-enriched pseudotime analysis; regions were kept separate.",
  "- Test: native scWorkbench `de_wilcox`, group A=WT, group B=KO, `min_pct=0`; all 259 peaks were supplied, with BH correction separately among peaks returning finite P values in each region.", "",
  "![Accessibility difference by position](figures/ebf1_pm3mb_accessibility_difference_by_position.png)", "",
  "![Differential-accessibility volcano](figures/ebf1_pm3mb_da_volcano.png)", "",
  "## Outputs", "",
  "- `ebf1_pm3mb_bdev_peak_universe.tsv`: complete tested peak universe.",
  "- `il7_enriched_pt_regions.tsv` and `il7_enriched_pt_cell_assignments.tsv.gz`: exact cell-selection provenance.",
  "- `all_il7_regions_ebf1_pm3mb_wt_vs_ko_differential_accessibility.tsv`: complete combined results.",
  "- `all_il7_regions_ebf1_pm3mb_significant_peaks_bh05.tsv`: BH-significant peaks.",
  "- `top20_peaks_per_region_direction.tsv`: ranked convenient summary.", "",
  "## Provenance", "", paste0("- Mapping bundle: `", bundle_dir, "`"),
  paste0("- Bundle fingerprint: `", fingerprint_before, "`"),
  paste0("- Enriched-PT source report: `", pt_report_dir, "`"),
  paste0("- Durable session: `", workspace, "`"),
  paste0("- Script version: `", script_version, "`"), "- Promotions: none.")
writeLines(report_md, file.path(stage_dir, "REPORT.md"))

html_answers <- paste0("<li>", answer_lines, "</li>", collapse = "")
html <- c("<!doctype html><html><head><meta charset='utf-8'>",
  "<title>Ebf1-window differential accessibility</title>",
  "<style>body{font-family:Arial,sans-serif;max-width:1100px;margin:30px auto;line-height:1.45;color:#222}code{background:#f1f1f1;padding:2px 4px}.warning{background:#fff3cd;border-left:5px solid #d79b00;padding:12px}img{max-width:100%}</style></head><body>",
  "<h1>Ebf1 ±3 Mb differential accessibility in IL-7-enriched trajectory regions</h1>",
  paste0("<h2>Answer</h2><ul>", html_answers, "</ul><p>Significance is BH-adjusted P &lt; 0.05 within each regional contrast; peaks with no test statistic are retained as untestable.</p>"),
  "<div class='warning'><b>Inference limitation:</b> each genotype is represented by one sample. These are descriptive cell-level Wilcoxon results, not biological-replicate-level genotype inference. Raw count ranks may also reflect cell sequencing depth.</div>",
  "<h2>Window and method</h2><p>Ebf1 TSS: GRCm38/mm10 chr11:44,617,317 (+ strand). All 259 overlapping Bdev peaks were tested separately in the exact two IL-7 WT-enriched pseudotime regions from the prior analysis.</p>",
  "<img src='figures/ebf1_pm3mb_accessibility_difference_by_position.png'>",
  "<img src='figures/ebf1_pm3mb_da_volcano.png'>",
  "<h2>Downloads</h2><ul><li><a href='ebf1_pm3mb_bdev_peak_universe.tsv'>Peak universe</a></li><li><a href='all_il7_regions_ebf1_pm3mb_wt_vs_ko_differential_accessibility.tsv'>Complete differential-accessibility table</a></li><li><a href='all_il7_regions_ebf1_pm3mb_significant_peaks_bh05.tsv'>BH-significant peaks</a></li></ul>",
  "</body></html>")
writeLines(html, file.path(stage_dir, "REPORT.html"))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
input_manifest <- data.table(role = c("mapping_bundle_manifest",
  "accepted_peak_mapping", "enriched_pt_regions", "enriched_pt_assignments",
  "analysis_script"), path = normalizePath(c(file.path(bundle_dir, "manifest.json"),
    peak_mapping_path, region_path, assignment_path, script_path), mustWork = TRUE))
input_manifest[, sha256 := vapply(path, sha256_file, character(1))]
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
completion <- list(status = "complete_unpromoted_report_only",
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  report = file.path(output_dir, "REPORT.html"), bundle = bundle_dir,
  bundle_fingerprint = fingerprint_before,
  bundle_fingerprint_unchanged = TRUE,
  matrix = "ATAC_bdev_mm10_peaks", matrix_layer = "counts",
  ebf1 = as.list(ebf1), tested_peaks_per_region = 259L,
  regions = split(summary_table, seq_len(nrow(summary_table))),
  operation = "de_wilcox", group_a = "IL7 WT", group_b = "IL7 KO",
  min_pct = 0, p_adjust_method = "BH_per_region",
  inference_scope = "descriptive_cell_level_only_one_sample_per_genotype",
  enriched_pt_source_report = pt_report_dir,
  session_workspace = workspace, promotions = 0L, checks = checks)
write_json(completion, file.path(stage_dir, "completion.json"))

if (!file.rename(stage_dir, output_dir))
  stop("Failed to finalize output directory")
log_event("Completed Ebf1-window enriched-PT DA report", output_dir)
