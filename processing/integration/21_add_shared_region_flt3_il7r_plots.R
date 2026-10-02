#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 21_add_shared_region_flt3_il7r_plots.R",
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
save_plot_atomic <- function(plot, path, type) atomic_replace(path,
  function(tmp) {
    if (type == "png") ggsave(tmp, plot, width = 8.5, height = 4.8,
      dpi = 300, device = "png") else ggsave(tmp, plot, width = 8.5,
      height = 4.8, device = cairo_pdf)
  })

completion_path <- file.path(report_dir, "completion.json")
completion <- read_json(completion_path, simplifyVector = FALSE)
assignments <- fread(cmd = paste("zcat", shQuote(file.path(report_dir,
  "shared_region_cell_assignments.tsv.gz"))))
de <- fread(file.path(report_dir,
  "all_genes_shared_regions_normalized_counts_and_de.tsv"))
summary_table <- fread(file.path(report_dir,
  "differential_expression_summary.tsv"))
target_features <- c("Flt3", "Il7r")
target_stats <- de[feature_id %in% target_features]
stopifnot(nrow(target_stats) == 4L,
  all(target_stats[, .N, by = region_id]$N == 2L))

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
plot_cells <- as.character(assignments$source_cell_id)
expr_matrix <- wb_matrix_read_features(staged, completion$normalized_matrix,
  target_features, plot_cells, missing = "NA")
stopifnot(!anyNA(expr_matrix), identical(rownames(expr_matrix),
  target_features), identical(colnames(expr_matrix), plot_cells))
expr <- as.data.table(as.table(expr_matrix))
setnames(expr, c("feature_id", "source_cell_id", "log1p_cp10k"))
expr[, `:=`(feature_id = as.character(feature_id),
  source_cell_id = as.character(source_cell_id),
  log1p_cp10k = as.numeric(log1p_cp10k))]
expr <- merge(expr, assignments[, .(source_cell_id, region_id, genotype)],
  by = "source_cell_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(expr) == 2L * nrow(assignments), !anyNA(expr$region_id),
  all(expr[, uniqueN(genotype), by = .(region_id, feature_id)]$V1 == 2L))
write_tsv_atomic(expr, file.path(report_dir,
  "target_genes_flt3_il7r_normalized_expression.tsv.gz"), compress = TRUE)

palette <- c("WT" = "#D55E00", "FLT3LKO" = "#0072B2",
  "IL7KO" = "#009E73")
manifest <- list()
for (i in seq_len(nrow(summary_table))) {
  reg <- summary_table[i]
  stats <- target_stats[region_id == reg$region_id]
  pd <- expr[region_id == reg$region_id]
  pd[, feature_id := factor(feature_id, levels = target_features)]
  pd[, genotype := factor(genotype, levels = c("WT", "KO"))]
  ko_name <- if (reg$cytokine == "FLT3L") "FLT3LKO" else "IL7KO"
  genotype_colours <- c("WT" = palette[["WT"]],
    "KO" = palette[[ko_name]])
  stats[, padj_text := fifelse(p_val_adj < 1e-300, "<1e-300",
    formatC(p_val_adj, format = "e", digits = 2))]
  facet_labels <- setNames(paste0(stats$feature_id, "\nlog2FC: ",
    sprintf("%.3f", stats$standard_log2fc), "; BH padj: ",
    stats$padj_text), stats$feature_id)
  p <- ggplot(pd, aes(genotype, log1p_cp10k, fill = genotype)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.25,
      colour = "grey30") +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.6,
      fill = "white", colour = "black", stroke = 0.35) +
    facet_wrap(~feature_id, scales = "free_y", nrow = 1,
      labeller = as_labeller(facet_labels)) +
    scale_fill_manual(values = genotype_colours) +
    labs(x = NULL, y = "Expression (log1p CP10K)",
      title = paste0(reg$cytokine,
        " shared-enriched region: Flt3 and Il7r"),
      subtitle = "Standard log2FC is WT / KO; P values are genome-wide BH adjusted",
      caption = "White points mark medians; one biological sample per group.") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 5.5, 5.5, 18))
  stem <- paste0(tolower(reg$region_id),
    "_flt3_il7r_target_gene_violin")
  pdf_rel <- file.path("figures", paste0(stem, ".pdf"))
  png_rel <- file.path("figures", paste0(stem, ".png"))
  save_plot_atomic(p, file.path(report_dir, pdf_rel), "pdf")
  save_plot_atomic(p, file.path(report_dir, png_rel), "png")
  manifest[[i]] <- data.table(region_id = reg$region_id,
    cytokine = reg$cytokine, genes = "Flt3;Il7r", png = png_rel,
    pdf = pdf_rel)
}
manifest <- rbindlist(manifest)
write_tsv_atomic(manifest, file.path(report_dir,
  "target_gene_violin_manifest.tsv"))

blocks <- paste(vapply(seq_len(nrow(manifest)), function(i) {
  x <- manifest[i]
  paste0("<h3>", x$cytokine,
    " comparison</h3><p>Standard log2FC and genome-wide BH-adjusted P are ",
    "shown in each facet. <a href='", x$pdf,
    "'>Vector PDF</a></p><img src='", x$png, "'>")
}, character(1)), collapse = "\n")
report_path <- file.path(report_dir, "REPORT.html")
html <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
if (!grepl("Targeted Flt3 and Il7r expression", html, fixed = TRUE)) {
  marker <- "<h2>FLT3L: WT-higher violin plots</h2>"
  insertion <- paste0("<h2>Targeted Flt3 and Il7r expression</h2>\n",
    blocks, "\n", marker)
  html <- sub(marker, insertion, html, fixed = TRUE)
}
atomic_replace(report_path, function(tmp) writeLines(html, tmp))

completion$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
completion$update_script_version <- script_version
completion$targeted_gene_plots <- list(genes = as.list(target_features),
  comparisons = as.list(manifest$region_id),
  normalized_expression = file.path(report_dir,
    "target_genes_flt3_il7r_normalized_expression.tsv.gz"))
completion$checks$target_gene_expression_complete <-
  nrow(expr) == 2L * nrow(assignments) && !anyNA(expr$log1p_cp10k)
completion$checks$target_gene_figures_complete <- all(file.exists(c(
  file.path(report_dir, manifest$png), file.path(report_dir, manifest$pdf))))
stopifnot(all(unlist(completion$checks)))
write_json_atomic(completion, completion_path)
cat("Added Flt3 and Il7r plots in place:", report_path, "\n")
print(target_stats[, .(region_id, feature_id, standard_log2fc, p_val_adj)])
