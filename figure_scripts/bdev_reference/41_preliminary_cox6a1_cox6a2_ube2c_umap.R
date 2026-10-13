#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 41_preliminary_cox6a1_cox6a2_ube2c_umap.R",
  "<bundle> <output_dir> <expected_fingerprint>"))
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
expected_fingerprint <- args[[3L]]
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con))
  write.table(x, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}

genes <- c("Cox6a1", "Cox6a2", "Ube2c")
bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "bdev_reference"),
  identical(manifest$current_state$fingerprint, expected_fingerprint),
  all(genes %in% wb_feature_ids(bundle, "SCT_data")))
cells <- as.character(wb_cells(bundle, "all"))
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", genes,
  cells, missing = "NA"))
stopifnot(identical(rownames(expression), genes),
  identical(colnames(expression), cells), all(is.finite(expression)))
space <- wb_space(bundle, "umap.atac", cells = cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells),
    all(c("x", "y") %in% names(space)))
  umap <- data.frame(cell_id = cells, x = as.numeric(space$x),
    y = as.numeric(space$y), stringsAsFactors = FALSE)
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  umap <- data.frame(cell_id = cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]), stringsAsFactors = FALSE)
}
stopifnot(all(is.finite(umap$x)), all(is.finite(umap$y)))

quantile_rows <- list()
value_rows <- list()
plots <- list()
for (i in seq_along(genes)) {
  gene <- genes[[i]]
  value <- as.numeric(expression[gene, ])
  positive <- value[value > 0]
  stopifnot(length(positive) > 0L)
  q <- as.numeric(quantile(positive, c(0, 0.75, 0.90, 0.99),
    names = FALSE, type = 7))
  stopifnot(all(is.finite(q)), q[[4L]] > 0)
  d <- umap
  d$expression <- value
  d$relative_to_p99 <- pmin(pmax(value, 0) / q[[4L]], 1)
  d <- d[order(d$relative_to_p99), ]
  value_rows[[i]] <- data.frame(cell_id = cells, gene = gene,
    sct_data = value, display_p99_clipped = pmin(pmax(value, 0), q[[4L]]))
  quantile_rows[[i]] <- data.frame(gene = gene,
    detected_fraction = mean(value > 0), positive_min = q[[1L]],
    positive_p75 = q[[2L]], positive_p90 = q[[3L]],
    positive_p99 = q[[4L]], observed_max = max(value))
  plots[[i]] <- ggplot(d, aes(x, y, color = relative_to_p99)) +
    geom_point(size = 0.31, alpha = 1, stroke = 0) +
    scale_color_gradient(low = "#FFFFFF", high = "#C51B5A",
      limits = c(0, 1), oob = scales::squish,
      breaks = c(0, 0.5, 0.75, 1),
      labels = c("0", "0.5", "0.75", "p99+"),
      name = "Expression /\npositive-cell p99") +
    coord_equal() +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9), legend.position = "bottom",
      legend.key.width = grid::unit(17, "mm")) +
    labs(title = gene,
      subtitle = sprintf("Positive-cell p99 = %.3g; detected in %.1f%%",
        q[[4L]], 100 * mean(value > 0)))
}
quantiles <- do.call(rbind, quantile_rows)
write.table(quantiles, file.path(stage_dir, "expression_color_scales.tsv"),
  sep = "\t", quote = TRUE, row.names = FALSE)
write_tsv_gz(do.call(rbind, value_rows), file.path(stage_dir,
  "cox6a1_cox6a2_ube2c_umap_expression.tsv.gz"))

figure <- wrap_plots(plots, nrow = 1, guides = "collect") +
  plot_annotation(title = "Cox6a1, Cox6a2, and Ube2c expression on Bdev umap.atac",
    subtitle = paste("White-to-crimson scales are gene-specific; values above",
      "each positive-cell p99 are clipped for display; high cells plotted last"),
    caption = "SCT_data values; bundle opened read-only.",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 8.5), legend.position = "bottom"))
ggsave(file.path(stage_dir, "cox6a1_cox6a2_ube2c_expression_umap.png"),
  figure, width = 13.5, height = 5.2, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "cox6a1_cox6a2_ube2c_expression_umap.pdf"),
  figure, width = 13.5, height = 5.2, device = cairo_pdf, bg = "white")

diagnostics <- list(status = "complete_preliminary_unregistered",
  bundle_id = manifest$bundle_id, bundle_fingerprint = expected_fingerprint,
  genes = as.list(genes), matrix = "SCT_data", space = "umap.atac",
  scale = "gene-specific white-to-crimson, zero to positive-cell p99",
  high_expression_drawn_last = TRUE, bundle_mutated = FALSE,
  capsule_registered = FALSE)
jsonlite::write_json(diagnostics, file.path(stage_dir, "diagnostics.json"),
  pretty = TRUE, auto_unbox = TRUE)
writeLines(c("# Cox6a1, Cox6a2, and Ube2c expression UMAP", "",
  "The three panels show SCT_data on the existing Bdev umap.atac.",
  "Each gene uses its own white-to-crimson zero-to-positive-p99 scale, and high-expression cells are drawn last."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
stopifnot(identical(wb_manifest(bundle)$current_state$fingerprint,
  expected_fingerprint))
wb_close_bundle(bundle)
bundle <- NULL
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed Cox6a1/Cox6a2/Ube2c expression UMAP at", output_dir, "\n")
