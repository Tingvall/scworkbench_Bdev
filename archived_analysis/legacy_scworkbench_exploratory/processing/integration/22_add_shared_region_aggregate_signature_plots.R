#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 22_add_shared_region_aggregate_signature_plots.R",
  "REPORT_DIR SCRIPT_VERSION execute"
), call. = FALSE)
report_dir <- normalizePath(args[[1L]], mustWork = TRUE)
script_version <- args[[2L]]
stopifnot(identical(args[[3L]], "execute"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

atomic_replace <- function(path, writer) {
  tmp <- paste0(path, ".tmp")
  if (file.exists(tmp)) unlink(tmp)
  writer(tmp)
  if (!file.rename(tmp, path)) stop("Could not replace: ", path,
    call. = FALSE)
}
write_tsv_atomic <- function(x, path, compress = FALSE) atomic_replace(
  path, function(tmp) fwrite(x, tmp, sep = "\t", na = "NA",
    compress = if (compress) "gzip" else "none"))
write_json_atomic <- function(x, path) atomic_replace(path, function(tmp)
  jsonlite::write_json(x, tmp, pretty = TRUE, auto_unbox = TRUE,
    digits = NA, null = "null", na = "null"))
save_plot_atomic <- function(plot, path, type) atomic_replace(path,
  function(tmp) {
    if (type == "png") ggsave(tmp, plot, width = 8.5, height = 5.1,
      dpi = 300, device = "png") else ggsave(tmp, plot, width = 8.5,
      height = 5.1, device = cairo_pdf)
  })
format_p <- function(x) ifelse(x < 1e-300, "<1e-300",
  formatC(x, format = "e", digits = 2))

completion_path <- file.path(report_dir, "completion.json")
completion <- read_json(completion_path, simplifyVector = FALSE)
source_report <- normalizePath(completion$source_original_region_report,
  mustWork = TRUE)
parse_path <- file.path(source_report,
  "flt3l_il7_parse_signal_aggregate_cell_scores.tsv.gz")
lineage_path <- file.path(source_report,
  "lineage_aggregate_cell_scores.tsv.gz")
stopifnot(file.exists(parse_path), file.exists(lineage_path))
assignments <- fread(cmd = paste("zcat", shQuote(file.path(report_dir,
  "shared_region_cell_assignments.tsv.gz"))))
summary_table <- fread(file.path(report_dir,
  "differential_expression_summary.tsv"))
stopifnot(!anyDuplicated(assignments$source_cell_id),
  all(c("source_cell_id", "region_id", "genotype") %in%
    names(assignments)))

parse_scores <- fread(cmd = paste("zcat", shQuote(parse_path)))[
  direction == "up" & cytokine %in% c("Flt3L", "IL-7"),
  .(cell_id, family = "Parse cytokine UP",
    signature = paste(cytokine, "UP"), aggregate_score, mapped_genes)]
lineage_scores <- fread(cmd = paste("zcat", shQuote(lineage_path)))[,
  .(cell_id, family = "Lineage", signature = lineage,
    aggregate_score, mapped_genes)]
scores <- rbindlist(list(parse_scores, lineage_scores), use.names = TRUE)
scores <- merge(scores, assignments[, .(source_cell_id, region_id,
  genotype, cytokine_comparison = cytokine)],
  by.x = "cell_id", by.y = "source_cell_id", all = FALSE, sort = FALSE)
scores[, `:=`(family = factor(family,
  levels = c("Parse cytokine UP", "Lineage")),
  signature = factor(signature,
    levels = c("Flt3L UP", "IL-7 UP", "B lineage", "T lineage")),
  genotype = factor(genotype, levels = c("WT", "KO")))]
stopifnot(nrow(scores) == 4L * nrow(assignments),
  !anyNA(scores$aggregate_score), all(scores[, uniqueN(genotype),
    by = .(region_id, family, signature)]$V1 == 2L),
  all(scores[, uniqueN(mapped_genes),
    by = .(family, signature)]$V1 == 1L))

stats <- scores[, {
  wt <- aggregate_score[genotype == "WT"]
  ko <- aggregate_score[genotype == "KO"]
  test <- wilcox.test(wt, ko, alternative = "two.sided", exact = FALSE,
    correct = FALSE)
  list(n_wt = length(wt), n_ko = length(ko),
    mean_wt = mean(wt), mean_ko = mean(ko),
    mean_difference_wt_minus_ko = mean(wt) - mean(ko),
    median_wt = median(wt), median_ko = median(ko),
    median_difference_wt_minus_ko = median(wt) - median(ko),
    wilcoxon_w = unname(test$statistic), p_val = test$p.value,
    mapped_genes = unique(mapped_genes))
}, by = .(region_id, cytokine_comparison, family, signature)]
stats[, p_val_adj := p.adjust(p_val, method = "BH")]
setorder(stats, region_id, family, signature)
stopifnot(nrow(stats) == 8L, all(is.finite(stats$p_val)),
  all(is.finite(stats$p_val_adj)))
write_tsv_atomic(scores, file.path(report_dir,
  "aggregate_signature_scores_shared_regions.tsv.gz"), compress = TRUE)
write_tsv_atomic(stats, file.path(report_dir,
  "aggregate_signature_significance.tsv"))

palette <- c("WT" = "#D55E00", "FLT3LKO" = "#0072B2",
  "IL7KO" = "#009E73")
manifest <- list()
manifest_index <- 0L
for (i in seq_len(nrow(summary_table))) {
  reg <- summary_table[i]
  for (family_value in c("Parse cytokine UP", "Lineage")) {
    pd <- scores[region_id == reg$region_id & family == family_value]
    ps <- stats[region_id == reg$region_id & family == family_value]
    signature_levels <- if (family_value == "Parse cytokine UP")
      c("Flt3L UP", "IL-7 UP") else c("B lineage", "T lineage")
    pd[, signature := factor(as.character(signature),
      levels = signature_levels)]
    ps[, p_text := format_p(p_val)]
    ps[, padj_text := format_p(p_val_adj)]
    facet_labels <- setNames(paste0(as.character(ps$signature),
      "\nmean WT-KO: ", sprintf("%.3f", ps$mean_difference_wt_minus_ko),
      "; P: ", ps$p_text, "\nBH padj: ", ps$padj_text),
      as.character(ps$signature))
    ko_name <- if (reg$cytokine == "FLT3L") "FLT3LKO" else "IL7KO"
    genotype_colours <- c("WT" = palette[["WT"]],
      "KO" = palette[[ko_name]])
    family_title <- if (family_value == "Parse cytokine UP")
      "Parse cytokine UP signals" else "B- and T-lineage gene sets"
    p <- ggplot(pd, aes(genotype, aggregate_score, fill = genotype)) +
      geom_violin(scale = "width", trim = TRUE, linewidth = 0.25,
        colour = "grey30") +
      geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.28,
        fill = "white", alpha = 0.82) +
      facet_wrap(~signature, scales = "free_y", nrow = 1,
        labeller = as_labeller(facet_labels)) +
      scale_fill_manual(values = genotype_colours) +
      labs(x = NULL, y = "Mean log1p CP10K expression",
        title = paste0(reg$cytokine, " shared-enriched region: ",
          family_title),
        subtitle = paste("Two-sided cell-level Wilcoxon; BH adjustment across",
          "all 8 requested tests"),
        caption = paste("One biological sample per group; significance is",
          "descriptive cell-level evidence.")) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        plot.margin = margin(5.5, 5.5, 5.5, 18))
    family_stem <- if (family_value == "Parse cytokine UP")
      "parse_cytokine_up_signals" else "b_t_lineage_signals"
    stem <- paste0(tolower(reg$region_id), "_", family_stem, "_violin")
    pdf_rel <- file.path("figures", paste0(stem, ".pdf"))
    png_rel <- file.path("figures", paste0(stem, ".png"))
    save_plot_atomic(p, file.path(report_dir, pdf_rel), "pdf")
    save_plot_atomic(p, file.path(report_dir, png_rel), "png")
    manifest_index <- manifest_index + 1L
    manifest[[manifest_index]] <- data.table(region_id = reg$region_id,
      cytokine = reg$cytokine, family = family_value,
      signatures = paste(signature_levels, collapse = ";"),
      png = png_rel, pdf = pdf_rel)
  }
}
manifest <- rbindlist(manifest)
write_tsv_atomic(manifest, file.path(report_dir,
  "aggregate_signature_violin_manifest.tsv"))

blocks <- paste(vapply(seq_len(nrow(manifest)), function(i) {
  x <- manifest[i]
  family_title <- if (x$family == "Parse cytokine UP")
    "Parse cytokine UP signals" else "B- and T-lineage gene sets"
  paste0("<h3>", x$cytokine, ": ", family_title,
    "</h3><p>Facet strips show mean WT-minus-KO score difference, raw ",
    "two-sided Wilcoxon P, and BH-adjusted P across all eight tests. ",
    "<a href='", x$pdf, "'>Vector PDF</a></p><img src='", x$png, "'>")
}, character(1)), collapse = "\n")
report_path <- file.path(report_dir, "REPORT.html")
html <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
if (!grepl("Aggregate cytokine and lineage signatures", html, fixed = TRUE)) {
  marker <- "<h2>FLT3L: WT-higher violin plots</h2>"
  insertion <- paste0(
    "<h2>Aggregate cytokine and lineage signatures</h2>",
    "<p class='note'><strong>Inference limitation:</strong> these Wilcoxon ",
    "tests compare cells, with one biological sample per group; they are not ",
    "replicate-level genotype tests. <a href='aggregate_signature_significance.tsv'>",
    "Download complete significance metrics</a>.</p>\n", blocks, "\n", marker)
  html <- sub(marker, insertion, html, fixed = TRUE)
}
atomic_replace(report_path, function(tmp) writeLines(html, tmp))

completion$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
completion$update_script_version <- script_version
completion$aggregate_signature_plots <- list(
  parse_signatures = list("Flt3L UP", "IL-7 UP"),
  lineage_signatures = list("B lineage", "T lineage"),
  source_report = source_report,
  score = "arithmetic mean log1p-CP10K across established mapped genes",
  test = "two-sided cell-level Wilcoxon rank-sum",
  p_adjustment = "BH jointly across 8 requested tests",
  inference = "descriptive_cell_level_only_one_biological_sample_per_group")
completion$checks$aggregate_score_rows_complete <-
  nrow(scores) == 4L * nrow(assignments)
completion$checks$aggregate_significance_complete <- nrow(stats) == 8L
completion$checks$aggregate_signature_figures_complete <-
  all(file.exists(c(file.path(report_dir, manifest$png),
    file.path(report_dir, manifest$pdf))))
stopifnot(all(unlist(completion$checks)))
write_json_atomic(completion, completion_path)
cat("Added aggregate signature plots in place:", report_path, "\n")
print(stats[, .(region_id, family, signature,
  mean_difference_wt_minus_ko, p_val, p_val_adj)])
