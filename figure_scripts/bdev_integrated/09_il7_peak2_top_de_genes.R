#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 09_il7_peak2_top_de_genes.R KO_BUNDLE PEAK2_TRACK_REPORT",
  "SESSION_WORKSPACE OUTPUT_DIR EXPECTED_KO_FP SCRIPT_VERSION"
), call. = FALSE)

ko_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
peak2_report <- normalizePath(args[[2L]], mustWork = TRUE)
session_workspace <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_ko_fp <- args[[5L]]
script_version <- args[[6L]]

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output already exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

write_json <- function(x, path) jsonlite::write_json(
  x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA,
  null = "null", na = "null")
write_tsv <- function(x, path) data.table::fwrite(
  x, path, sep = "\t", na = "NA")
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 1e-300, "<1e-300",
    formatC(x, format = "e", digits = 2)))
}

assignment_path <- file.path(peak2_report,
  "all_four_samples_il7_peak2_pt_cell_assignments.tsv.gz")
stopifnot(file.exists(assignment_path))
cells <- fread(cmd = paste("zcat", shQuote(assignment_path)))
needed <- c("source_cell_id", "sample_id", "genotype", "knockout_target",
  "support_status", "mapped_pt")
stopifnot(all(needed %in% names(cells)), !anyDuplicated(cells$source_cell_id),
  nrow(cells) == 11263L, uniqueN(cells$sample_id) == 4L,
  setequal(cells$genotype, c("WT", "KO")),
  setequal(cells$knockout_target, c("Flt3l", "Il7")),
  all(cells$support_status %in% c("supported", "borderline")),
  all(is.finite(cells$mapped_pt)))

bundle <- wb_open_bundle(ko_bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_ko_fp))

session <- wb_session_resume(bundle, session_workspace)
stopifnot(is.null(session$durable$promotion_receipt))
find_run <- function(label, operation) {
  marker <- paste0("[", label, "] ")
  ids <- names(session$runs)[vapply(session$runs, function(x) {
    note <- if (is.null(x$notes)) "" else as.character(x$notes)
    identical(x$op_name, operation) && startsWith(note, marker)
  }, logical(1))]
  stopifnot(length(ids) == 1L)
  session$runs[[ids[[1L]]]]
}

normal_run <- find_run("normalize_rna", "matrix_log_normalize")
matrix_id <- normal_run$outputs$matrix$matrix_id
stopifnot(identical(matrix_id, "RNA_log1p_cp10k_il7_peak2_il7r_v1"))

contrast_specs <- data.table(
  contrast_id = c("flt3l_wt_vs_ko", "il7_wt_vs_ko"),
  run_label = c("de_flt3l", "de_il7"),
  target = c("Flt3l", "Il7"),
  comparison = c("FLT3L WT versus FLT3LKO", "IL-7 WT versus IL7KO"))

all_de <- rbindlist(lapply(seq_len(nrow(contrast_specs)), function(i) {
  spec <- contrast_specs[i]
  de <- as.data.table(find_run(spec$run_label, "de_wilcox")$outputs$de_table$table)
  stopifnot(nrow(de) == 33696L, !anyDuplicated(de$feature_id))
  de[, `:=`(contrast_id = spec$contrast_id,
    knockout_target = spec$target, comparison = spec$comparison)]
  de
}), use.names = TRUE)

# Group A was WT and group B was KO in both saved scWorkbench tests.
# Rank independently by BH-adjusted P, then raw P and absolute effect.
rankable <- all_de[is.finite(avg_log1p_fc) & avg_log1p_fc != 0 &
  is.finite(p_val) & is.finite(p_val_adj)]
rankable[, direction := fifelse(avg_log1p_fc > 0,
  "WT higher", "KO higher")]
rankable[, abs_effect := abs(avg_log1p_fc)]
setorder(rankable, contrast_id, direction, p_val_adj, p_val,
  -abs_effect, feature_id)
top <- rankable[, head(.SD, 10L), by = .(contrast_id, knockout_target,
  comparison, direction)]
top[, abs_effect := NULL]
top[, rank := seq_len(.N), by = .(contrast_id, direction)]
stopifnot(nrow(top) == 40L,
  all(top[, .N, by = .(contrast_id, direction)]$N == 10L))
setcolorder(top, c("contrast_id", "knockout_target", "comparison",
  "direction", "rank", "feature_id", "pct_a", "pct_b", "avg_expr_a",
  "avg_expr_b", "avg_log1p_fc", "p_val", "p_val_adj"))
write_tsv(top, file.path(stage_dir, "top10_up_down_genes.tsv"))

staged <- wb_session_view(session)
on.exit(try(wb_close_bundle(staged), silent = TRUE), add = TRUE)
cell_order <- as.character(cells$source_cell_id)
features <- unique(top$feature_id)
expr <- wb_matrix_read_features(staged, matrix_id, features, cell_order,
  missing = "NA")
stopifnot(nrow(expr) == length(features), ncol(expr) == length(cell_order),
  identical(rownames(expr), features), identical(colnames(expr), cell_order))
expr_long <- as.data.table(as.table(expr))
setnames(expr_long, c("feature_id", "source_cell_id", "log1p_cp10k"))
expr_long[, `:=`(feature_id = as.character(feature_id),
  source_cell_id = as.character(source_cell_id),
  log1p_cp10k = as.numeric(log1p_cp10k))]
expr_long <- merge(expr_long, cells[, .(source_cell_id, sample_id,
  genotype, knockout_target)], by = "source_cell_id", all.x = TRUE,
  sort = FALSE)
stopifnot(nrow(expr_long) == length(features) * nrow(cells),
  !anyNA(expr_long$sample_id), !anyNA(expr_long$log1p_cp10k))

plot_manifest <- list()
sample_colours <- c("WT" = "#0072B2", "KO" = "#CC79A7")
for (i in seq_len(nrow(contrast_specs))) {
  spec <- contrast_specs[i]
  for (dir in c("WT higher", "KO higher")) {
    selected <- top[contrast_id == spec$contrast_id & direction == dir]
    pd <- expr_long[knockout_target == spec$target &
      feature_id %in% selected$feature_id]
    gene_levels <- rev(selected[order(rank), feature_id])
    pd[, feature_id := factor(feature_id, levels = gene_levels)]
    pd[, genotype := factor(genotype, levels = c("WT", "KO"))]
    subtitle <- if (dir == "WT higher")
      "Positive effect: expression is higher in WT (group A)" else
      "Negative effect: expression is higher in KO (group B)"
    p <- ggplot(pd, aes(genotype, log1p_cp10k, fill = genotype)) +
      geom_violin(scale = "width", trim = TRUE, linewidth = 0.25,
        colour = "grey30") +
      stat_summary(fun = median, geom = "point", shape = 21, size = 1.7,
        fill = "white", colour = "black", stroke = 0.35) +
      facet_wrap(~feature_id, scales = "free_y", ncol = 5) +
      scale_fill_manual(values = sample_colours) +
      labs(x = NULL, y = "Expression (log1p CP10K)",
        title = paste0(spec$comparison, ": top 10 ", dir, " genes"),
        subtitle = paste0(subtitle,
          "; ranked by genome-wide BH-adjusted Wilcoxon P"),
        caption = paste0("Exact IL-7 peak-2 pseudotime interval; ",
          "cell-level descriptive comparison (one sample per group).")) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))
    stem <- paste0(spec$contrast_id, "_",
      ifelse(dir == "WT higher", "wt_higher", "ko_higher"), "_top10")
    ggsave(file.path(stage_dir, "figures", paste0(stem, ".pdf")), p,
      width = 13, height = 6.8, device = cairo_pdf)
    ggsave(file.path(stage_dir, "figures", paste0(stem, ".png")), p,
      width = 13, height = 6.8, dpi = 300)
    plot_manifest[[length(plot_manifest) + 1L]] <- data.table(
      contrast_id = spec$contrast_id, comparison = spec$comparison,
      direction = dir, png = file.path("figures", paste0(stem, ".png")),
      pdf = file.path("figures", paste0(stem, ".pdf")))
  }
}
plot_manifest <- rbindlist(plot_manifest)
write_tsv(plot_manifest, file.path(stage_dir, "plot_manifest.tsv"))

table_rows <- vapply(seq_len(nrow(top)), function(i) {
  x <- top[i]
  paste0("<tr><td>", x$comparison, "</td><td>", x$direction,
    "</td><td>", x$rank, "</td><td><code>", x$feature_id,
    "</code></td><td>", sprintf("%+.4f", x$avg_log1p_fc),
    "</td><td>", fmt_p(x$p_val_adj), "</td></tr>")
}, character(1))
figure_blocks <- vapply(seq_len(nrow(plot_manifest)), function(i) {
  x <- plot_manifest[i]
  paste0("<h2>", x$comparison, ": ", x$direction, "</h2>",
    "<p><a href='", x$pdf, "'>Vector PDF</a></p>",
    "<img src='", x$png, "'>")
}, character(1))
html <- c("<!doctype html><html><head><meta charset='utf-8'>",
  "<title>IL-7 peak-2 top differential genes</title>",
  "<style>body{font-family:sans-serif;max-width:1400px;margin:30px auto;line-height:1.45}img{max-width:100%}table{border-collapse:collapse}th,td{border:1px solid #bbb;padding:5px 8px;text-align:right}th:first-child,td:first-child,th:nth-child(2),td:nth-child(2),th:nth-child(4),td:nth-child(4){text-align:left}.note{background:#fff3cd;padding:12px}</style></head><body>",
  "<h1>Top differential genes in the IL-7 peak-2 interval</h1>",
  "<p>Four separate plots show the ten strongest WT-higher and ten strongest KO-higher genes for the FLT3L WT versus FLT3LKO and IL-7 WT versus IL7KO comparisons. Genes were ranked within direction by genome-wide BH-adjusted P, then raw P, absolute effect, and feature ID. White points mark medians.</p>",
  "<p class='note'><strong>Inference limitation:</strong> each group contains one biological sample. The native scWorkbench Wilcoxon results describe cell-level separation and are not replicate-level significance tests.</p>",
  figure_blocks,
  "<h2>Selected genes and statistics</h2>",
  "<p><a href='top10_up_down_genes.tsv'>Download complete selection table</a></p>",
  "<table><tr><th>Comparison</th><th>Direction</th><th>Rank</th><th>Gene</th><th>Average log1p effect</th><th>BH P</th></tr>",
  table_rows, "</table></body></html>")
writeLines(html, file.path(stage_dir, "REPORT.html"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
checks <- list(
  exact_peak2_assignments = nrow(cells) == 11263L,
  saved_tests_cover_all_genes = all(all_de[, .N, by = contrast_id]$N == 33696L),
  ten_genes_per_direction_and_comparison = nrow(top) == 40L,
  selected_expression_complete = nrow(expr_long) ==
    length(features) * nrow(cells) && !anyNA(expr_long$log1p_cp10k),
  four_plots_complete = nrow(plot_manifest) == 4L &&
    all(file.exists(file.path(stage_dir, plot_manifest$png))) &&
    all(file.exists(file.path(stage_dir, plot_manifest$pdf))),
  parent_bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after),
  result_promotions_zero = is.null(session$durable$promotion_receipt))
stopifnot(all(unlist(checks)))

completion <- list(status = "complete_unpromoted_report_only",
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script_version = script_version, report = file.path(output_dir, "REPORT.html"),
  source_bundle = ko_bundle_dir,
  source_bundle_fingerprint = fingerprint_before,
  source_peak2_report = peak2_report,
  source_session_workspace = session_workspace,
  normalized_matrix = matrix_id,
  reused_operations = c("matrix_log_normalize", "de_wilcox:de_flt3l",
    "de_wilcox:de_il7"),
  new_analysis_operations = 0L,
  ranking = "within contrast and effect direction: p_val_adj, p_val, descending absolute avg_log1p_fc, feature_id",
  direction_definition = "group A is WT; positive is WT higher; negative is KO higher",
  inference_scope = "descriptive_cell_level_only_one_biological_sample_per_group",
  selected_genes = lapply(seq_len(nrow(top)), function(i) as.list(top[i])),
  checks = checks, promotions = 0L)
write_json(completion, file.path(stage_dir, "completion.json"))

stopifnot(file.rename(stage_dir, output_dir))
cat("Completed:", file.path(output_dir, "REPORT.html"), "\n")
