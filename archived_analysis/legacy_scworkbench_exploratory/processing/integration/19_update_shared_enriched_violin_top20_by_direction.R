#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 19_update_shared_enriched_violin_top20_by_direction.R",
  "REPORT_DIR SCRIPT_VERSION execute"
), call. = FALSE)

report_dir <- normalizePath(args[[1L]], mustWork = TRUE)
script_version <- args[[2L]]
stopifnot(identical(args[[3L]], "execute"))

suppressPackageStartupMessages({
  library(scworkbench)
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
  "differential_expression_summary.tsv", "completion.json",
  "shared_region_cell_assignments.tsv.gz",
  file.path("figures", "shared_enriched_regions.png"),
  file.path("figures", "shared_region_de_volcano.png"))
stopifnot(all(file.exists(file.path(report_dir, required))))

completion_path <- file.path(report_dir, "completion.json")
completion <- jsonlite::read_json(completion_path, simplifyVector = FALSE)
selected <- fread(file.path(report_dir, "genes_passing_fc1p5_bh05.tsv"))
assignments <- fread(cmd = paste("zcat", shQuote(file.path(report_dir,
  "shared_region_cell_assignments.tsv.gz"))))
summary_table <- fread(file.path(report_dir,
  "differential_expression_summary.tsv"))

stopifnot(all(c("region_id", "cytokine", "feature_id", "direction",
  "p_val_adj", "p_val", "standard_log2fc") %in% names(selected)),
  all(c("source_cell_id", "region_id", "genotype") %in% names(assignments)))
selected[, abs_standard_log2fc := abs(standard_log2fc)]
setorder(selected, region_id, direction, p_val_adj, p_val,
  -abs_standard_log2fc, feature_id)
plot_genes <- selected[, head(.SD, 20L), by = .(region_id, direction)]
plot_genes[, abs_standard_log2fc := NULL]
setorder(plot_genes, region_id, direction, p_val_adj, p_val, feature_id)
plot_genes[, plot_rank := seq_len(.N), by = .(region_id, direction)]

expected_counts <- selected[, .(expected = min(.N, 20L)),
  by = .(region_id, direction)]
actual_counts <- plot_genes[, .(actual = .N), by = .(region_id, direction)]
count_check <- merge(expected_counts, actual_counts,
  by = c("region_id", "direction"), all = TRUE)
stopifnot(nrow(count_check) == 4L,
  all(count_check$expected == count_check$actual))

bundle <- wb_open_bundle(completion$source_bundle, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
stopifnot(identical(wb_manifest(bundle)$current_state$fingerprint,
  completion$source_bundle_fingerprint),
  all(assignments$source_cell_id %in%
    wb_cells(bundle, "set:qc_use_downstream")))
session <- wb_session_resume(bundle, completion$source_session_workspace)
stopifnot(is.null(session$durable$promotion_receipt))
staged <- wb_session_view(session)
on.exit(try(wb_close_bundle(staged), silent = TRUE), add = TRUE)

plot_features <- unique(plot_genes$feature_id)
plot_cells <- as.character(assignments$source_cell_id)
expr_matrix <- wb_matrix_read_features(staged, completion$normalized_matrix,
  plot_features, plot_cells, missing = "NA")
stopifnot(!anyNA(expr_matrix),
  identical(rownames(expr_matrix), plot_features),
  identical(colnames(expr_matrix), plot_cells))
expr <- as.data.table(as.table(expr_matrix))
setnames(expr, c("feature_id", "source_cell_id", "log1p_cp10k"))
expr[, `:=`(feature_id = as.character(feature_id),
  source_cell_id = as.character(source_cell_id),
  log1p_cp10k = as.numeric(log1p_cp10k))]
expr <- merge(expr, assignments[, .(source_cell_id, region_id, genotype)],
  by = "source_cell_id", all.x = TRUE, sort = FALSE)
expr <- expr[plot_genes[, .(region_id, feature_id)],
  on = .(region_id, feature_id), nomatch = 0L]
stopifnot(uniqueN(expr[, .(region_id, feature_id)]) ==
    nrow(unique(plot_genes[, .(region_id, feature_id)])),
  !anyNA(expr$log1p_cp10k), all(expr[, uniqueN(genotype),
    by = .(region_id, feature_id)]$V1 == 2L))

write_tsv_atomic(plot_genes,
  file.path(report_dir, "violin_plot_gene_selection.tsv"))
write_tsv_atomic(expr, file.path(report_dir,
  "violin_plot_normalized_expression.tsv.gz"), compress = TRUE)

violin_manifest <- list()
manifest_index <- 0L
for (i in seq_len(nrow(summary_table))) {
  reg <- summary_table[i]
  for (direction_value in c("WT_higher", "KO_higher")) {
    genes <- plot_genes[region_id == reg$region_id &
      direction == direction_value][order(plot_rank)]
    pd <- expr[region_id == reg$region_id &
      feature_id %in% genes$feature_id]
    pd[, feature_id := factor(feature_id, levels = genes$feature_id)]
    pd[, genotype := factor(genotype, levels = c("WT", "KO"))]
    genotype_colours <- c("WT" = "#D55E00", "KO" =
      if (reg$cytokine == "FLT3L") "#0072B2" else "#009E73")
    direction_label <- if (direction_value == "WT_higher")
      "WT-higher" else "KO-higher"
    p <- ggplot(pd, aes(genotype, log1p_cp10k, fill = genotype)) +
      geom_violin(scale = "width", trim = TRUE, linewidth = 0.25,
        colour = "grey30") +
      stat_summary(fun = median, geom = "point", shape = 21, size = 1.4,
        fill = "white", colour = "black", stroke = 0.35) +
      facet_wrap(~feature_id, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = genotype_colours) +
      labs(x = NULL, y = "Expression (log1p CP10K)",
        title = paste0(reg$cytokine, " shared-enriched region: ",
          direction_label, " genes"),
        subtitle = paste0("Top ", nrow(genes), " of ",
          selected[region_id == reg$region_id &
            direction == direction_value, .N],
          " passing genes in this direction; BH P < 0.05 and ",
          "|log2FC| >= log2(1.5)"),
        caption = paste("White points mark medians; one biological sample",
          "per group.")) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        plot.margin = margin(5.5, 5.5, 5.5, 18))
    direction_stem <- if (direction_value == "WT_higher")
      "wt_higher" else "ko_higher"
    stem <- paste0(tolower(reg$region_id), "_", direction_stem,
      "_top20_cutoff_genes_violin")
    plot_height <- max(6.5, 2.35 * ceiling(nrow(genes) / 4))
    pdf_path <- file.path(report_dir, "figures", paste0(stem, ".pdf"))
    png_path <- file.path(report_dir, "figures", paste0(stem, ".png"))
    atomic_replace(pdf_path, function(tmp) ggsave(tmp, p, width = 12,
      height = plot_height, device = cairo_pdf))
    atomic_replace(png_path, function(tmp) ggsave(tmp, p, width = 12,
      height = plot_height, dpi = 300, device = "png"))
    manifest_index <- manifest_index + 1L
    violin_manifest[[manifest_index]] <- data.table(
      region_id = reg$region_id, cytokine = reg$cytokine,
      direction = direction_value,
      passing_direction_genes = selected[region_id == reg$region_id &
        direction == direction_value, .N], plotted_genes = nrow(genes),
      png = file.path("figures", paste0(stem, ".png")),
      pdf = file.path("figures", paste0(stem, ".pdf")))
  }
}
violin_manifest <- rbindlist(violin_manifest)
write_tsv_atomic(violin_manifest,
  file.path(report_dir, "violin_plot_manifest.tsv"))

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
  direction_label <- if (x$direction == "WT_higher")
    "WT-higher" else "KO-higher"
  paste0("<h2>", x$cytokine, ": ", direction_label,
    " violin plots</h2><p>", x$passing_direction_genes,
    " genes pass both cutoffs in this direction; the top ",
    x$plotted_genes, " are shown. <a href='", x$pdf,
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

stale_violin_files <- list.files(file.path(report_dir, "figures"),
  pattern = "_top2_per_direction_cutoff_genes_violin[.](png|pdf)$",
  full.names = TRUE)
if (length(stale_violin_files)) unlink(stale_violin_files)

completion$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
completion$update_script_version <- script_version
completion$violin_rule <- paste("up to 20 most significant cutoff-passing",
  "genes per direction and comparison, with a separate figure for each",
  "direction")
completion$violin_colour_scale <- list(WT = "#D55E00",
  FLT3L_KO = "#0072B2", IL7_KO = "#009E73")
completion$checks$violin_exact_two_per_direction <- NULL
completion$checks$violin_top20_per_direction <-
  all(count_check$expected == count_check$actual)
completion$checks$violin_expression_complete <- TRUE
completion$checks$figures_complete <- all(file.exists(file.path(report_dir,
  violin_manifest$png))) && all(file.exists(file.path(report_dir,
  violin_manifest$pdf)))
stopifnot(all(unlist(completion$checks)))
write_json_atomic(completion, completion_path)

cat("Updated in place:", file.path(report_dir, "REPORT.html"), "\n")
print(violin_manifest)
