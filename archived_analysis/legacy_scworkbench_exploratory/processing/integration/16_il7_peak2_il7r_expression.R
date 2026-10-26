#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 16_il7_peak2_il7r_expression.R KO_BUNDLE PEAK2_TRACK_REPORT",
  "SESSION_ROOT OUTPUT_DIR EXPECTED_KO_FP SCRIPT_VERSION GENE"
), call. = FALSE)

ko_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
peak2_report <- normalizePath(args[[2L]], mustWork = TRUE)
session_root <- normalizePath(args[[3L]], mustWork = FALSE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_ko_fp <- args[[5L]]
script_version <- args[[6L]]
gene <- args[[7L]]
stopifnot(identical(gene, "Il7r"))

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output already exists: ", stage_dir, call. = FALSE)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage_dir, "figures"), showWarnings = FALSE)
dir.create(session_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

write_json <- function(x, path) jsonlite::write_json(
  x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA,
  null = "null", na = "null")
write_tsv <- function(x, path) data.table::fwrite(x, path, sep = "\t", na = "NA")
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 1e-300, "<1e-300",
    formatC(x, format = "e", digits = 2)))
}
find_run <- function(session, label, op) {
  marker <- paste0("[", label, "] ")
  ids <- names(session$runs)[vapply(session$runs, function(x) {
    note <- if (is.null(x$notes)) "" else as.character(x$notes)
    identical(x$op_name, op) && startsWith(note, marker)
  }, logical(1))]
  stopifnot(length(ids) <= 1L)
  if (length(ids)) session$runs[[ids[[1L]]]] else NULL
}

assignment_path <- file.path(peak2_report,
  "all_four_samples_il7_peak2_pt_cell_assignments.tsv.gz")
stopifnot(file.exists(assignment_path))
cells <- fread(cmd = paste("zcat", shQuote(assignment_path)))
needed <- c("source_cell_id", "sample_id", "genotype", "knockout_target",
  "support_status", "mapped_pt")
stopifnot(all(needed %in% names(cells)), !anyDuplicated(cells$source_cell_id),
  uniqueN(cells$sample_id) == 4L,
  setequal(cells$genotype, c("WT", "KO")),
  setequal(cells$knockout_target, c("Flt3l", "Il7")),
  all(cells$support_status %in% c("supported", "borderline")),
  all(is.finite(cells$mapped_pt)))

bundle <- wb_open_bundle(ko_bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_ko_fp))
downstream <- as.character(wb_cells(bundle, "set:qc_use_downstream"))
stopifnot(all(cells$source_cell_id %in% downstream))
meta <- as.data.table(wb_meta(bundle, cells = downstream))
alignment <- merge(cells, meta[, .(source_cell_id = cell_id,
  bundle_sample_id = library_id, bundle_genotype = genotype,
  bundle_target = knockout_target)], by = "source_cell_id", all.x = TRUE,
  sort = FALSE)
stopifnot(nrow(alignment) == nrow(cells), !anyNA(alignment$bundle_sample_id),
  identical(alignment$sample_id, alignment$bundle_sample_id),
  identical(alignment$genotype, alignment$bundle_genotype),
  identical(alignment$knockout_target, alignment$bundle_target))

write_json(wb_op_interface("matrix_log_normalize", bundle = bundle),
  file.path(stage_dir, "matrix_log_normalize_interface.json"))
write_json(wb_op_interface("de_wilcox", bundle = bundle),
  file.path(stage_dir, "de_wilcox_interface.json"))

session_id <- "il7_peak2_il7r_expression_v1"
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) wb_session_resume(bundle, workspace) else
  wb_session_create(bundle, workspace_root = session_root,
    session_id = session_id)
if (!is.null(session$durable$promotion_receipt))
  stop("The report-only session was unexpectedly promoted")

matrix_id <- "RNA_log1p_cp10k_il7_peak2_il7r_v1"
normal_run <- find_run(session, "normalize_rna", "matrix_log_normalize")
if (is.null(normal_run)) {
  result <- wb_session_run_op(session, bundle, "matrix_log_normalize",
    inputs = list(matrix = "RNA_counts", cell_set = "all", min_cells = 10L,
      output_matrix_id = matrix_id, scale_factor = 10000,
      overwrite = FALSE),
    notes = paste0("[normalize_rna] Full-bundle log1p CP10K matrix; ",
      "report-only Il7r analysis in the exact IL-7 peak-2 interval."))
  session <- result$session
  normal_run <- result$run
}
stopifnot("matrix" %in% names(normal_run$outputs))

contrast_specs <- data.table(
  contrast = c("FLT3L context: WT vs KO", "IL-7 context: WT vs KO"),
  target = c("Flt3l", "Il7"), label = c("flt3l", "il7"))
stats <- rbindlist(lapply(seq_len(nrow(contrast_specs)), function(i) {
  spec <- contrast_specs[i]
  wt <- cells[knockout_target == spec$target & genotype == "WT",
    source_cell_id]
  ko <- cells[knockout_target == spec$target & genotype == "KO",
    source_cell_id]
  stopifnot(length(wt) >= 100L, length(ko) >= 100L,
    !length(intersect(wt, ko)))
  run <- find_run(session, paste0("de_", spec$label), "de_wilcox")
  if (is.null(run)) {
    result <- wb_session_run_op(session, bundle, "de_wilcox",
      inputs = list(matrix = matrix_id, cell_set = "set:qc_use_downstream",
        min_cells = 10L, group_a = as.character(wt),
        group_b = as.character(ko), min_pct = 0,
        p_adjust_method = "BH"),
      notes = paste0("[de_", spec$label, "] Group A=WT, group B=KO; ",
        "cells restricted to the exact IL-7 peak-2 mapped-PT interval; ",
        "cell-level descriptive inference because n=1 sample per group."))
    session <<- result$session
    run <- result$run
  }
  de <- as.data.table(run$outputs$de_table$table)
  row <- de[feature_id == gene]
  stopifnot(nrow(row) == 1L)
  data.table(contrast = spec$contrast, knockout_target = spec$target,
    group_a = "WT", group_b = "KO", n_wt_cells = length(wt),
    n_ko_cells = length(ko), feature_id = row$feature_id,
    pct_wt = row$pct_a, pct_ko = row$pct_b,
    avg_expr_wt = row$avg_expr_a, avg_expr_ko = row$avg_expr_b,
    avg_log1p_fc = row$avg_log1p_fc, p_val = row$p_val,
    p_val_adj_genomewide_bh = row$p_val_adj)
}))

staged <- wb_session_view(session)
on.exit(try(wb_close_bundle(staged), silent = TRUE), add = TRUE)
cell_order <- as.character(cells$source_cell_id)
expr <- as.numeric(wb_matrix_read_features(staged, matrix_id, gene,
  cell_order, missing = "NA"))
stopifnot(length(expr) == nrow(cells), !anyNA(expr))
plot_data <- copy(cells)
plot_data[, log1p_cp10k := expr]
plot_data[, sample_label := factor(sample_id,
  levels = c("FLT3LKO_BM_WT", "FLT3LKO_BM_KO", "IL7KO_BM_WT", "IL7KO_BM_KO"),
  labels = c("FLT3L WT", "FLT3L KO", "IL7 WT", "IL7 KO"))]
plot_data[, sample_group := factor(paste(knockout_target, genotype),
  levels = c("Flt3l WT", "Flt3l KO", "Il7 WT", "Il7 KO"))]
write_tsv(plot_data[, .(source_cell_id, sample_id, genotype,
  knockout_target, support_status, mapped_pt, log1p_cp10k)],
  file.path(stage_dir, "il7r_cell_expression.tsv.gz"))
write_tsv(stats, file.path(stage_dir, "il7r_wt_vs_ko_statistics.tsv"))
write_tsv(cells[, .(n_cells = .N, median_mapped_pt = median(mapped_pt),
  minimum_mapped_pt = min(mapped_pt), maximum_mapped_pt = max(mapped_pt)),
  by = .(sample_id, genotype, knockout_target)],
  file.path(stage_dir, "sample_cell_counts.tsv"))

ann <- copy(stats)
ann[, x := c(1.5, 3.5)]
ann[, y := max(plot_data$log1p_cp10k) * c(1.08, 1.08)]
ann[, label := paste0("cell-level Wilcoxon\nBH P = ",
  fmt_p(p_val_adj_genomewide_bh))]
y_top <- max(ann$y) * 1.08
sample_colours <- c("FLT3L WT" = "#D55E00", "FLT3L KO" = "#E69F00",
  "IL7 WT" = "#0072B2", "IL7 KO" = "#CC79A7")
p <- ggplot(plot_data, aes(sample_label, log1p_cp10k,
    fill = sample_label)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.35,
    colour = "grey25") +
  geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.3,
    fill = "white", alpha = 0.85) +
  geom_segment(data = ann, aes(x = x - 0.42, xend = x + 0.42,
    y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = ann, aes(x = x, y = y, label = label),
    inherit.aes = FALSE, vjust = -0.25, size = 3.1) +
  scale_fill_manual(values = sample_colours) +
  scale_y_continuous(limits = c(0, y_top), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Il7r expression (log1p CP10K)",
    title = "Il7r expression within the IL-7 peak-2 pseudotime interval",
    subtitle = paste0("Exact peak-2 mapped-PT interval; descriptive ",
      "cell-level tests (one biological sample per group).")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9.5),
    plot.margin = margin(10, 12, 8, 10))
ggsave(file.path(stage_dir, "figures", "il7r_all_four_samples_violin.pdf"),
  p, width = 8.2, height = 5.5, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "il7r_all_four_samples_violin.png"),
  p, width = 8.2, height = 5.5, dpi = 300)

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
checks <- list(
  exact_four_samples = uniqueN(cells$sample_id) == 4L,
  exact_peak2_assignments = nrow(cells) == 11263L,
  supported_or_borderline_only = all(cells$support_status %in%
    c("supported", "borderline")),
  finite_mapped_pt_only = all(is.finite(cells$mapped_pt)),
  source_ids_align_to_bundle = all(cells$source_cell_id %in% downstream),
  metadata_alignment = nrow(alignment) == nrow(cells) &&
    !anyNA(alignment$bundle_sample_id),
  il7r_complete_for_all_cells = nrow(plot_data) == nrow(cells) &&
    !anyNA(plot_data$log1p_cp10k),
  native_pairwise_tests_complete = nrow(stats) == 2L &&
    all(is.finite(stats$p_val)),
  figure_exists = file.exists(file.path(stage_dir, "figures",
    "il7r_all_four_samples_violin.png")),
  parent_bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after),
  result_promotions_zero = is.null(session$durable$promotion_receipt))
stopifnot(all(unlist(checks)))

stats_lines <- vapply(seq_len(nrow(stats)), function(i) {
  x <- stats[i]
  paste0("<tr><td>", x$contrast, "</td><td>", x$n_wt_cells,
    "</td><td>", x$n_ko_cells, "</td><td>",
    sprintf("%+.4f", x$avg_log1p_fc), "</td><td>", fmt_p(x$p_val),
    "</td><td>", fmt_p(x$p_val_adj_genomewide_bh), "</td></tr>")
}, character(1))
html <- c("<!doctype html><html><head><meta charset='utf-8'>",
  "<title>IL-7 peak 2 Il7r expression</title>",
  "<style>body{font-family:sans-serif;max-width:1100px;margin:30px auto;line-height:1.45}img{max-width:100%}table{border-collapse:collapse}th,td{border:1px solid #bbb;padding:6px 9px;text-align:right}th:first-child,td:first-child{text-align:left}.note{background:#fff3cd;padding:12px}</style></head><body>",
  "<h1>Il7r expression in the IL-7 peak-2 pseudotime interval</h1>",
  "<img src='figures/il7r_all_four_samples_violin.png'>",
  "<h2>WT versus KO statistics</h2>",
  "<table><tr><th>Contrast</th><th>WT cells</th><th>KO cells</th><th>avg log1p effect</th><th>P</th><th>genome-wide BH P</th></tr>",
  stats_lines, "</table>",
  "<p class='note'><strong>Inference limitation:</strong> each condition contains one biological sample. The native scWorkbench Wilcoxon scores quantify cell-level separation and are useful descriptively, but cells are not biological replicates; these values do not establish replicate-level significance.</p>",
  "<h2>Method</h2><p>Cells are the exact supported/borderline query cells within the saved IL-7 peak-2 mapped-pseudotime bounds used for the corrected four-sample ATAC tracks. RNA counts were normalized with native scWorkbench <code>matrix_log_normalize</code> to log1p CP10K. Native <code>de_wilcox</code> compared WT with KO separately within the FLT3L and IL-7 sample pairs; BH adjustment is across all tested genes in each contrast. No result was promoted.</p>",
  "<h2>Downloads</h2><ul>",
  "<li><a href='figures/il7r_all_four_samples_violin.pdf'>Vector PDF</a></li>",
  "<li><a href='il7r_wt_vs_ko_statistics.tsv'>Statistics</a></li>",
  "<li><a href='il7r_cell_expression.tsv.gz'>Per-cell plotted values</a></li>",
  "<li><a href='sample_cell_counts.tsv'>Sample cell counts</a></li></ul>",
  "</body></html>")
writeLines(html, file.path(stage_dir, "REPORT.html"))

completion <- list(status = "complete_unpromoted_report_only",
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script_version = script_version, gene = gene,
  report = file.path(output_dir, "REPORT.html"),
  source_bundle = ko_bundle_dir,
  source_bundle_fingerprint = fingerprint_before,
  source_bundle_fingerprint_unchanged = TRUE,
  source_peak2_report = peak2_report,
  cell_definition = "all four samples; exact supported/borderline cells in the saved Il7_peak2 mapped-PT interval",
  expression_matrix = matrix_id, normalization = "log1p_cp10k",
  test = "native scWorkbench de_wilcox; WT versus KO within each cytokine context",
  inference_scope = "descriptive_cell_level_only_one_biological_sample_per_group",
  statistics = lapply(seq_len(nrow(stats)), function(i) as.list(stats[i])),
  session_workspace = workspace, promotions = 0L, checks = checks)
write_json(completion, file.path(stage_dir, "completion.json"))

wb_close_bundle(staged)
wb_close_bundle(bundle)
on.exit(NULL, add = FALSE)
if (!file.rename(stage_dir, output_dir))
  stop("Could not atomically finalize output directory")
cat("Completed Il7r peak-2 expression report:", output_dir, "\n")
