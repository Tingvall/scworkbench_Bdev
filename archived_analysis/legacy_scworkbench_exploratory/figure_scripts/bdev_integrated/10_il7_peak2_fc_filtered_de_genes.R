#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 10_il7_peak2_fc_filtered_de_genes.R KO_BUNDLE",
  "PEAK2_TRACK_REPORT SESSION_WORKSPACE OUTPUT_DIR EXPECTED_KO_FP",
  "SCRIPT_VERSION"
), call. = FALSE)

ko_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
peak2_report <- normalizePath(args[[2L]], mustWork = TRUE)
session_workspace <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
expected_ko_fp <- args[[5L]]
script_version <- args[[6L]]
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
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
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")

assignment_path <- file.path(peak2_report,
  "all_four_samples_il7_peak2_pt_cell_assignments.tsv.gz")
stopifnot(file.exists(assignment_path))
cells <- fread(cmd = paste("zcat", shQuote(assignment_path)))
stopifnot(nrow(cells) == 11263L, !anyDuplicated(cells$source_cell_id),
  uniqueN(cells$sample_id) == 4L,
  setequal(cells$knockout_target, c("Flt3l", "Il7")),
  setequal(cells$genotype, c("WT", "KO")),
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
  comparison = c("FLT3L WT versus FLT3LKO", "IL-7 WT versus IL7KO"),
  ko_label = c("FLT3LKO", "IL7KO"))
de_tables <- lapply(seq_len(nrow(contrast_specs)), function(i) {
  d <- as.data.table(find_run(contrast_specs$run_label[[i]],
    "de_wilcox")$outputs$de_table$table)
  stopifnot(nrow(d) == 33696L, !anyDuplicated(d$feature_id))
  d
})
names(de_tables) <- contrast_specs$contrast_id
features <- as.character(de_tables[[1L]]$feature_id)
stopifnot(setequal(features, de_tables[[2L]]$feature_id))

staged <- wb_session_view(session)
on.exit(try(wb_close_bundle(staged), silent = TRUE), add = TRUE)
cell_order <- as.character(cells$source_cell_id)
group_indices <- list(
  flt3l_wt = which(cells$knockout_target == "Flt3l" & cells$genotype == "WT"),
  flt3l_ko = which(cells$knockout_target == "Flt3l" & cells$genotype == "KO"),
  il7_wt = which(cells$knockout_target == "Il7" & cells$genotype == "WT"),
  il7_ko = which(cells$knockout_target == "Il7" & cells$genotype == "KO"))
stopifnot(all(vapply(group_indices, length, integer(1)) >= 100L))

blocks <- split(features, ceiling(seq_along(features) / 128L))
means <- rbindlist(lapply(blocks, function(feats) {
  m <- wb_matrix_read_features(staged, matrix_id, feats, cell_order,
    missing = "NA")
  stopifnot(!anyNA(m), identical(rownames(m), feats),
    identical(colnames(m), cell_order))
  linear_cp10k <- expm1(as.matrix(m))
  data.table(feature_id = feats,
    flt3l_wt_mean_cp10k = rowMeans(linear_cp10k[, group_indices$flt3l_wt,
      drop = FALSE]),
    flt3l_ko_mean_cp10k = rowMeans(linear_cp10k[, group_indices$flt3l_ko,
      drop = FALSE]),
    il7_wt_mean_cp10k = rowMeans(linear_cp10k[, group_indices$il7_wt,
      drop = FALSE]),
    il7_ko_mean_cp10k = rowMeans(linear_cp10k[, group_indices$il7_ko,
      drop = FALSE]))
}))
stopifnot(nrow(means) == 33696L, !anyDuplicated(means$feature_id))

log2fc_cutoff <- log2(1.5)
selected <- rbindlist(lapply(seq_len(nrow(contrast_specs)), function(i) {
  spec <- contrast_specs[i]
  d <- merge(copy(de_tables[[spec$contrast_id]]), means,
    by = "feature_id", all.x = TRUE)
  if (spec$target == "Flt3l") {
    d[, `:=`(mean_cp10k_wt = flt3l_wt_mean_cp10k,
      mean_cp10k_ko = flt3l_ko_mean_cp10k)]
  } else {
    d[, `:=`(mean_cp10k_wt = il7_wt_mean_cp10k,
      mean_cp10k_ko = il7_ko_mean_cp10k)]
  }
  d[, standard_log2fc := log2((mean_cp10k_wt + 1) /
    (mean_cp10k_ko + 1))]
  d <- d[is.finite(p_val_adj) & p_val_adj < 0.05 &
    is.finite(standard_log2fc) & abs(standard_log2fc) >= log2fc_cutoff]
  d[, `:=`(contrast_id = spec$contrast_id,
    knockout_target = spec$target, comparison = spec$comparison,
    direction = fifelse(standard_log2fc > 0, "WT higher",
      paste0(spec$ko_label, " higher")))]
  d[, abs_standard_log2fc := abs(standard_log2fc)]
  setorder(d, direction, p_val_adj, p_val, -abs_standard_log2fc, feature_id)
  d[, rank_within_direction := seq_len(.N), by = direction]
  d[, c("flt3l_wt_mean_cp10k", "flt3l_ko_mean_cp10k",
    "il7_wt_mean_cp10k", "il7_ko_mean_cp10k",
    "abs_standard_log2fc") := NULL]
  d
}), use.names = TRUE)
stopifnot(nrow(selected) > 0L,
  all(selected$p_val_adj < 0.05),
  all(abs(selected$standard_log2fc) >= log2fc_cutoff))
setcolorder(selected, c("contrast_id", "knockout_target", "comparison",
  "direction", "rank_within_direction", "feature_id", "mean_cp10k_wt",
  "mean_cp10k_ko", "standard_log2fc", "pct_a", "pct_b", "avg_expr_a",
  "avg_expr_b", "avg_log1p_fc", "p_val", "p_val_adj"))
write_tsv(selected, file.path(stage_dir,
  "fc1p5_bh05_selected_genes.tsv"))
summary <- selected[, .(n_genes = .N,
  minimum_abs_log2fc = min(abs(standard_log2fc)),
  maximum_abs_log2fc = max(abs(standard_log2fc))),
  by = .(contrast_id, comparison, direction)]
setorder(summary, contrast_id, direction)
write_tsv(summary, file.path(stage_dir, "selection_summary.tsv"))

selected_features <- unique(selected$feature_id)
expr <- wb_matrix_read_features(staged, matrix_id, selected_features,
  cell_order, missing = "NA")
stopifnot(!anyNA(expr), identical(rownames(expr), selected_features),
  identical(colnames(expr), cell_order))
expr_long <- as.data.table(as.table(expr))
setnames(expr_long, c("feature_id", "source_cell_id", "log1p_cp10k"))
expr_long[, `:=`(feature_id = as.character(feature_id),
  source_cell_id = as.character(source_cell_id),
  log1p_cp10k = as.numeric(log1p_cp10k))]
expr_long <- merge(expr_long, cells[, .(source_cell_id, genotype,
  knockout_target)], by = "source_cell_id", all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(expr_long$genotype), !anyNA(expr_long$log1p_cp10k))

plot_manifest <- list()
page_size <- 20L
for (i in seq_len(nrow(contrast_specs))) {
  spec <- contrast_specs[i]
  directions <- unique(selected[contrast_id == spec$contrast_id,
    direction])
  for (dir in directions) {
    genes <- selected[contrast_id == spec$contrast_id & direction == dir][
      order(rank_within_direction), feature_id]
    pages <- split(genes, ceiling(seq_along(genes) / page_size))
    for (page in seq_along(pages)) {
      page_genes <- pages[[page]]
      pd <- expr_long[knockout_target == spec$target &
        feature_id %in% page_genes]
      pd[, feature_id := factor(feature_id, levels = rev(page_genes))]
      pd[, genotype := factor(genotype, levels = c("WT", "KO"),
        labels = c("WT", spec$ko_label))]
      p <- ggplot(pd, aes(genotype, log1p_cp10k, fill = genotype)) +
        geom_violin(scale = "width", trim = TRUE, linewidth = 0.25,
          colour = "grey30") +
        stat_summary(fun = median, geom = "point", shape = 21, size = 1.6,
          fill = "white", colour = "black", stroke = 0.35) +
        facet_wrap(~feature_id, scales = "free_y", ncol = 5) +
        scale_fill_manual(values = c("WT" = "#0072B2",
          setNames("#CC79A7", spec$ko_label))) +
        labs(x = NULL, y = "Expression (log1p CP10K)",
          title = paste0(spec$comparison, ": ", dir,
            " genes passing 1.5-fold and BH cutoffs"),
          subtitle = paste0("Standard |log2FC| ≥ ",
            sprintf("%.3f", log2fc_cutoff),
            " and genome-wide BH-adjusted P < 0.05; page ", page,
            " of ", length(pages)),
          caption = paste0("Exact IL-7 peak-2 pseudotime interval; ",
            "cell-level descriptive comparison (one sample per group).")) +
        theme_bw(base_size = 10) +
        theme(legend.position = "none", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))
      direction_id <- if (startsWith(dir, "WT")) "wt_higher" else "ko_higher"
      stem <- paste0(spec$contrast_id, "_", direction_id,
        "_fc1p5_bh05_page", page)
      plot_height <- 3.6 + 2.8 * ceiling(length(page_genes) / 5)
      ggsave(file.path(stage_dir, "figures", paste0(stem, ".pdf")), p,
        width = 13, height = plot_height, device = cairo_pdf)
      ggsave(file.path(stage_dir, "figures", paste0(stem, ".png")), p,
        width = 13, height = plot_height, dpi = 300)
      plot_manifest[[length(plot_manifest) + 1L]] <- data.table(
        contrast_id = spec$contrast_id, comparison = spec$comparison,
        direction = dir, page = page, n_genes = length(page_genes),
        png = file.path("figures", paste0(stem, ".png")),
        pdf = file.path("figures", paste0(stem, ".pdf")))
    }
  }
}
plot_manifest <- rbindlist(plot_manifest)
write_tsv(plot_manifest, file.path(stage_dir, "plot_manifest.tsv"))

summary_rows <- vapply(seq_len(nrow(summary)), function(i) {
  x <- summary[i]
  paste0("<tr><td>", x$comparison, "</td><td>", x$direction,
    "</td><td>", x$n_genes, "</td><td>",
    sprintf("%.3f", x$minimum_abs_log2fc), "</td><td>",
    sprintf("%.3f", x$maximum_abs_log2fc), "</td></tr>")
}, character(1))
figure_blocks <- vapply(seq_len(nrow(plot_manifest)), function(i) {
  x <- plot_manifest[i]
  paste0("<h2>", x$comparison, ": ", x$direction,
    if (max(plot_manifest[contrast_id == x$contrast_id &
      direction == x$direction, page]) > 1L) paste0(" (page ", x$page, ")")
      else "", "</h2><p><a href='", x$pdf, "'>Vector PDF</a></p>",
    "<img src='", x$png, "'>")
}, character(1))
html <- c("<!doctype html><html><head><meta charset='utf-8'>",
  "<title>IL-7 peak-2 genes passing 1.5-fold and BH cutoffs</title>",
  "<style>body{font-family:sans-serif;max-width:1400px;margin:30px auto;line-height:1.45}img{max-width:100%}table{border-collapse:collapse}th,td{border:1px solid #bbb;padding:6px 9px;text-align:right}th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left}.note{background:#fff3cd;padding:12px}</style></head><body>",
  "<h1>IL-7 peak-2 differential genes with at least 1.5-fold change</h1>",
  "<p>Genes pass both genome-wide BH-adjusted Wilcoxon P &lt; 0.05 and standard |log2FC| &ge; log2(1.5). Standard log2FC is log2((mean linear CP10K in WT + 1)/(mean linear CP10K in KO + 1)). No normalization or differential test was repeated.</p>",
  "<p class='note'><strong>Inference limitation:</strong> each group contains one biological sample. These are descriptive cell-level comparisons, not replicate-level genotype tests.</p>",
  "<h2>Selection summary</h2><table><tr><th>Comparison</th><th>Direction</th><th>Genes</th><th>Minimum |log2FC|</th><th>Maximum |log2FC|</th></tr>",
  summary_rows, "</table><p><a href='fc1p5_bh05_selected_genes.tsv'>Download selected genes and statistics</a></p>",
  figure_blocks, "</body></html>")
writeLines(html, file.path(stage_dir, "REPORT.html"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
checks <- list(exact_peak2_assignments = nrow(cells) == 11263L,
  saved_tests_cover_all_genes = all(vapply(de_tables, nrow, integer(1)) ==
    33696L), selected_genes_meet_both_cutoffs = all(selected$p_val_adj < 0.05) &&
    all(abs(selected$standard_log2fc) >= log2fc_cutoff),
  all_selected_expression_complete = !anyNA(expr_long$log1p_cp10k),
  all_plot_files_complete = nrow(plot_manifest) >= 4L &&
    all(file.exists(file.path(stage_dir, plot_manifest$png))) &&
    all(file.exists(file.path(stage_dir, plot_manifest$pdf))),
  parent_bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after), result_promotions_zero =
    is.null(session$durable$promotion_receipt))
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
    "de_wilcox:de_il7"), new_analysis_operations = 0L,
  standard_log2fc = "log2((mean linear CP10K WT + 1)/(mean linear CP10K KO + 1))",
  log2fc_cutoff = log2fc_cutoff, adjusted_p_cutoff = 0.05,
  selection_summary = lapply(seq_len(nrow(summary)),
    function(i) as.list(summary[i])),
  inference_scope = "descriptive_cell_level_only_one_biological_sample_per_group",
  checks = checks, promotions = 0L)
write_json(completion, file.path(stage_dir, "completion.json"))
stopifnot(file.rename(stage_dir, output_dir))
cat("Completed:", file.path(output_dir, "REPORT.html"), "\n")
