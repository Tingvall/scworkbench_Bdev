#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 18_update_shared_enriched_violin_top2_per_direction.R",
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

required <- c("genes_passing_fc1p5_bh05.tsv",
  "violin_plot_normalized_expression.tsv.gz",
  "differential_expression_summary.tsv", "completion.json",
  file.path("figures", "shared_enriched_regions.png"),
  file.path("figures", "shared_region_de_volcano.png"))
stopifnot(all(file.exists(file.path(report_dir, required))))

selected <- fread(file.path(report_dir, "genes_passing_fc1p5_bh05.tsv"))
stopifnot(all(c("region_id", "cytokine", "feature_id", "direction",
  "p_val_adj", "p_val", "standard_log2fc") %in% names(selected)))
selected[, abs_standard_log2fc := abs(standard_log2fc)]
setorder(selected, region_id, direction, p_val_adj, p_val,
  -abs_standard_log2fc, feature_id)
plot_genes <- selected[, head(.SD, 2L), by = .(region_id, direction)]
plot_genes[, abs_standard_log2fc := NULL]
plot_genes[, direction_order := fifelse(direction == "WT_higher", 1L, 2L)]
setorder(plot_genes, region_id, direction_order, p_val_adj, p_val, feature_id)
plot_genes[, direction_order := NULL]
plot_genes[, plot_rank := seq_len(.N), by = region_id]
counts <- plot_genes[, .N, by = .(region_id, direction)]
stopifnot(nrow(counts) == 4L, all(counts$N == 2L),
  uniqueN(plot_genes[, .(region_id, feature_id)]) == 8L)

expr_path <- file.path(report_dir,
  "violin_plot_normalized_expression.tsv.gz")
expr <- fread(cmd = paste("zcat", shQuote(expr_path)))
stopifnot(all(c("source_cell_id", "feature_id", "log1p_cp10k",
  "region_id", "genotype") %in% names(expr)))
expr <- expr[plot_genes[, .(region_id, feature_id)],
  on = .(region_id, feature_id), nomatch = 0L]
stopifnot(uniqueN(expr[, .(region_id, feature_id)]) == 8L,
  !anyNA(expr$log1p_cp10k), all(expr[, uniqueN(genotype),
    by = .(region_id, feature_id)]$V1 == 2L))

write_tsv_atomic(plot_genes,
  file.path(report_dir, "violin_plot_gene_selection.tsv"))
write_tsv_atomic(expr, expr_path, compress = TRUE)

summary_table <- fread(file.path(report_dir,
  "differential_expression_summary.tsv"))
violin_manifest <- vector("list", nrow(summary_table))
for (i in seq_len(nrow(summary_table))) {
  reg <- summary_table[i]
  genes <- plot_genes[region_id == reg$region_id][order(plot_rank)]
  pd <- expr[region_id == reg$region_id & feature_id %in% genes$feature_id]
  facet_labels <- setNames(paste0(genes$feature_id, " (",
    fifelse(genes$direction == "WT_higher", "WT up", "KO up"), ")"),
    genes$feature_id)
  pd[, feature_id := factor(feature_id, levels = genes$feature_id)]
  pd[, genotype := factor(genotype, levels = c("WT", "KO"))]
  p <- ggplot(pd, aes(genotype, log1p_cp10k, fill = genotype)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.25,
      colour = "grey30") +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.6,
      fill = "white", colour = "black", stroke = 0.35) +
    facet_wrap(~feature_id, scales = "free_y", ncol = 2,
      labeller = as_labeller(facet_labels)) +
    scale_fill_manual(values = c("WT" = "#0072B2", "KO" = "#CC79A7")) +
    labs(x = NULL, y = "Expression (log1p CP10K)",
      title = paste0(reg$cytokine,
        " shared-enriched region: cutoff-passing genes"),
      subtitle = paste0("Top 2 per direction (4 total) among ",
        reg$passing_fc1p5_bh05,
        " passing genes; BH P < 0.05 and |log2FC| >= log2(1.5)"),
      caption = "White points mark medians; one biological sample per group.") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 5.5, 5.5, 18))
  stem <- paste0(tolower(reg$region_id),
    "_top2_per_direction_cutoff_genes_violin")
  pdf_path <- file.path(report_dir, "figures", paste0(stem, ".pdf"))
  png_path <- file.path(report_dir, "figures", paste0(stem, ".png"))
  atomic_replace(pdf_path, function(tmp) ggsave(tmp, p, width = 8.5,
    height = 7, device = cairo_pdf))
  atomic_replace(png_path, function(tmp) ggsave(tmp, p, width = 8.5,
    height = 7, dpi = 300, device = "png"))
  violin_manifest[[i]] <- data.table(region_id = reg$region_id,
    cytokine = reg$cytokine, passing_genes = reg$passing_fc1p5_bh05,
    plotted_genes = nrow(genes),
    png = file.path("figures", paste0(stem, ".png")),
    pdf = file.path("figures", paste0(stem, ".pdf")))
}
violin_manifest <- rbindlist(violin_manifest)
write_tsv_atomic(violin_manifest,
  file.path(report_dir, "violin_plot_manifest.tsv"))
stale_violin_files <- list.files(file.path(report_dir, "figures"),
  pattern = "_top20_cutoff_genes_violin[.](png|pdf)$", full.names = TRUE)
if (length(stale_violin_files)) unlink(stale_violin_files)

summary_rows <- vapply(seq_len(nrow(summary_table)), function(i) {
  x <- summary_table[i]
  paste0("<tr><td>", x$cytokine, "</td><td>", x$first_bin, "-",
    x$last_bin, "</td><td>", x$n_wt_cells, "</td><td>",
    x$n_ko_cells, "</td><td>", x$genes_tested, "</td><td>",
    x$significant_bh05, "</td><td>", x$passing_fc1p5_bh05,
    "</td><td>", x$passing_wt_higher, "</td><td>",
    x$passing_ko_higher, "</td></tr>")
}, character(1))
violin_blocks <- vapply(seq_len(nrow(violin_manifest)), function(i) {
  x <- violin_manifest[i]
  paste0("<h2>", x$cytokine, " violin plots</h2><p>",
    x$passing_genes, " genes pass both cutoffs; the two most significant in ",
    "each direction are shown (", x$plotted_genes,
    " total). <a href='", x$pdf,
    "'>Vector PDF</a></p><img src='", x$png, "'>")
}, character(1))
html <- c("<!doctype html><html><head><meta charset='utf-8'>",
  "<title>Shared-enriched trajectory-region RNA DE</title>",
  "<style>body{font-family:sans-serif;max-width:1400px;margin:30px auto;line-height:1.45}img{max-width:100%}table{border-collapse:collapse}th,td{border:1px solid #bbb;padding:6px 9px;text-align:right}th:first-child,td:first-child{text-align:left}.note{background:#fff3cd;padding:12px}</style></head><body>",
  "<h1>Shared WT-and-KO enriched trajectory-region RNA differential expression</h1>",
  "<p>Regions are maximal contiguous bins whose Bdev-normalized rolling cell fractions are strictly above 1 in both WT and KO. IL-7 is restricted to the original peak 2. The saved native scWorkbench log1p-CP10K matrix was reused; new native Wilcoxon tests were required for the changed cell sets.</p>",
  "<p class='note'><strong>Inference limitation:</strong> each genotype contains one biological sample per cytokine. Results describe cell-level separation and are not replicate-level genotype tests.</p>",
  "<img src='figures/shared_enriched_regions.png'>",
  "<h2>DE summary</h2><table><tr><th>Cytokine</th><th>Bins</th><th>WT cells</th><th>KO cells</th><th>Genes tested</th><th>BH &lt; 0.05</th><th>Both cutoffs</th><th>WT higher</th><th>KO higher</th></tr>",
  summary_rows, "</table>",
  "<p><a href='all_genes_shared_regions_normalized_counts_and_de.tsv'>Download complete normalized-count and DEG table</a> | <a href='genes_passing_fc1p5_bh05.tsv'>Download cutoff-passing genes</a></p>",
  "<h2>Full DEG volcano plots</h2><p><a href='figures/shared_region_de_volcano.pdf'>Vector PDF</a></p><img src='figures/shared_region_de_volcano.png'>",
  violin_blocks, "</body></html>")
atomic_replace(file.path(report_dir, "REPORT.html"),
  function(tmp) writeLines(html, tmp))

completion_path <- file.path(report_dir, "completion.json")
completion <- jsonlite::read_json(completion_path, simplifyVector = FALSE)
completion$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
completion$update_script_version <- script_version
completion$violin_rule <- paste("two most significant cutoff-passing genes",
  "per direction and comparison")
completion$checks$violin_cap_respected <- NULL
completion$checks$violin_exact_two_per_direction <- TRUE
completion$checks$violin_expression_complete <- TRUE
completion$checks$figures_complete <- all(file.exists(file.path(report_dir,
  violin_manifest$png))) && all(file.exists(file.path(report_dir,
  violin_manifest$pdf)))
stopifnot(all(unlist(completion$checks)))
write_json_atomic(completion, completion_path)

cat("Updated in place:", file.path(report_dir, "REPORT.html"), "\n")
print(plot_genes[, .(region_id, direction, feature_id, standard_log2fc,
  p_val_adj)])
