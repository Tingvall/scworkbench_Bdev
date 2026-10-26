#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "Usage: 20_update_shared_enriched_plot_annotations.R",
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
write_tsv_atomic <- function(x, path) atomic_replace(path,
  function(tmp) fwrite(x, tmp, sep = "\t", na = "NA"))
write_json_atomic <- function(x, path) atomic_replace(path, function(tmp)
  jsonlite::write_json(x, tmp, pretty = TRUE, auto_unbox = TRUE,
    digits = NA, null = "null", na = "null"))
save_plot_atomic <- function(plot, path, width, height, type) {
  atomic_replace(path, function(tmp) {
    if (type == "png") ggsave(tmp, plot, width = width, height = height,
      dpi = 300, device = "png") else ggsave(tmp, plot, width = width,
      height = height, device = cairo_pdf)
  })
}

completion_path <- file.path(report_dir, "completion.json")
completion <- read_json(completion_path, simplifyVector = FALSE)
regions <- fread(file.path(report_dir, "shared_enriched_pt_regions.tsv"))
combined <- fread(file.path(report_dir,
  "all_genes_shared_regions_normalized_counts_and_de.tsv"))
plot_genes <- fread(file.path(report_dir,
  "violin_plot_gene_selection.tsv"))
expr <- fread(cmd = paste("zcat", shQuote(file.path(report_dir,
  "violin_plot_normalized_expression.tsv.gz"))))
violin_manifest <- fread(file.path(report_dir, "violin_plot_manifest.tsv"))
distribution <- fread(file.path(completion$source_distribution_report,
  "query_mapped_pt_cell_distribution.tsv"))

palette <- c("WT" = "#D55E00", "FLT3LKO" = "#0072B2",
  "IL7KO" = "#009E73")

# Apply the same biological-group palette to the selected-region curves.
region_plot_data <- copy(distribution[target %in% c("Flt3l", "Il7")])
region_plot_data[, target_label := factor(target,
  levels = c("Flt3l", "Il7"), labels = c("FLT3L", "IL-7"))]
region_plot_data[, colour_group := fifelse(genotype == "WT", "WT",
  fifelse(target == "Flt3l", "FLT3LKO", "IL7KO"))]
region_plot_regions <- copy(regions)
region_plot_regions[, target_label := factor(target,
  levels = c("Flt3l", "Il7"), labels = c("FLT3L", "IL-7"))]
p_regions <- ggplot(region_plot_data,
    aes(bin_midpoint, bdev_normalized_rolling_fraction,
      colour = colour_group, group = colour_group)) +
  geom_rect(data = region_plot_regions,
    aes(xmin = pt_left, xmax = pt_right, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#7B61A8", alpha = 0.18) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey25") +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_text(data = region_plot_regions,
    aes(x = (pt_left + pt_right) / 2, y = Inf,
      label = paste0("bins ", first_bin, "-", last_bin)),
    inherit.aes = FALSE, vjust = 1.3, fontface = "bold", size = 3.2) +
  facet_wrap(~target_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = palette) +
  labs(x = "Mapped Bdev pseudotime",
    y = "Bdev-normalized rolling cell fraction", colour = NULL,
    title = "Shared WT-and-KO enriched pseudotime regions",
    subtitle = paste0("Shading requires both genotype curves > 1; ",
      "IL-7 is restricted to the original peak 2")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank())
save_plot_atomic(p_regions, file.path(report_dir, "figures",
  "shared_enriched_regions.pdf"), 11, 5.2, "pdf")
save_plot_atomic(p_regions, file.path(report_dir, "figures",
  "shared_enriched_regions.png"), 11, 5.2, "png")

# Label every genome-wide significant gene with an absolute standard log2FC
# of at least one. Fixed label columns keep all labels visible without relying
# on a non-container dependency such as ggrepel.
volcano_data <- copy(combined)
volcano_data[, neglog10_padj := -log10(pmax(p_val_adj, 1e-300))]
volcano_data[, status := fifelse(passes_fc1p5_bh05 &
  standard_log2fc > 0, "WT higher: passes both", fifelse(
  passes_fc1p5_bh05 & standard_log2fc < 0 & target == "Flt3l",
  "FLT3LKO higher: passes both", fifelse(passes_fc1p5_bh05 &
  standard_log2fc < 0 & target == "Il7", "IL7KO higher: passes both",
  "Does not pass both")))]
volcano_data[, region_label := factor(region_id,
  levels = regions$region_id, labels = paste0(regions$cytokine,
    " shared region"))]
volcano_labels <- volcano_data[is.finite(p_val_adj) & p_val_adj < 0.05 &
  is.finite(standard_log2fc) & abs(standard_log2fc) >= 1]
volcano_labels[, side := fifelse(standard_log2fc > 0, "positive",
  "negative")]
volcano_labels[, label_y := seq(max(neglog10_padj),
  max(neglog10_padj) * 0.32, length.out = .N),
  by = .(region_id, side)]
volcano_ranges <- volcano_data[is.finite(standard_log2fc),
  .(xmax = max(abs(standard_log2fc))), by = region_id]
volcano_labels <- merge(volcano_labels, volcano_ranges, by = "region_id")
volcano_labels[, `:=`(label_x = fifelse(standard_log2fc > 0,
  0.98 * xmax, -0.98 * xmax), label_hjust = fifelse(
  standard_log2fc > 0, 1, 0))]
write_tsv_atomic(volcano_labels[, .(region_id, cytokine, feature_id,
  direction, standard_log2fc, p_val, p_val_adj)],
  file.path(report_dir, "volcano_labeled_genes.tsv"))
p_volcano <- ggplot(volcano_data[is.finite(neglog10_padj)],
    aes(standard_log2fc, neglog10_padj, colour = status)) +
  geom_point(size = 0.7, alpha = 0.55) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
    linetype = "dashed", colour = "grey35") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted",
    colour = "grey35") +
  geom_segment(data = volcano_labels,
    aes(x = standard_log2fc, y = neglog10_padj, xend = label_x,
      yend = label_y), inherit.aes = FALSE, linewidth = 0.25,
    colour = "grey35") +
  geom_text(data = volcano_labels,
    aes(x = label_x, y = label_y, label = feature_id,
      hjust = label_hjust), inherit.aes = FALSE, size = 2.7) +
  facet_wrap(~region_label, ncol = 2) +
  scale_colour_manual(values = c("WT higher: passes both" = palette[["WT"]],
    "FLT3LKO higher: passes both" = palette[["FLT3LKO"]],
    "IL7KO higher: passes both" = palette[["IL7KO"]],
    "Does not pass both" = "grey75")) +
  labs(x = "Standard log2FC (WT / KO)", y = "-log10(BH-adjusted P)",
    colour = NULL, title = "Full shared-region differential expression",
    subtitle = paste("All 33,696 genes tested in each comparison; labels:",
      "BH padj < 0.05 and |log2FC| >= 1")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"))
save_plot_atomic(p_volcano, file.path(report_dir, "figures",
  "shared_region_de_volcano.pdf"), 11, 6.5, "pdf")
save_plot_atomic(p_volcano, file.path(report_dir, "figures",
  "shared_region_de_volcano.png"), 11, 6.5, "png")

# Re-render the existing four direction-specific violin figures with the new
# palette and a BH-adjusted P value in every gene strip.
for (i in seq_len(nrow(violin_manifest))) {
  item <- violin_manifest[i]
  genes <- plot_genes[region_id == item$region_id &
    direction == item$direction][order(plot_rank)]
  pd <- expr[region_id == item$region_id &
    feature_id %in% genes$feature_id]
  pd[, feature_id := factor(feature_id, levels = genes$feature_id)]
  pd[, genotype := factor(genotype, levels = c("WT", "KO"))]
  ko_name <- if (item$cytokine == "FLT3L") "FLT3LKO" else "IL7KO"
  genotype_colours <- c("WT" = palette[["WT"]],
    "KO" = palette[[ko_name]])
  padj_labels <- setNames(ifelse(genes$p_val_adj < 1e-300,
    "<1e-300", formatC(genes$p_val_adj, format = "e", digits = 2)),
    genes$feature_id)
  facet_labels <- setNames(paste0(genes$feature_id, "\nBH padj: ",
    padj_labels[genes$feature_id]), genes$feature_id)
  direction_label <- if (item$direction == "WT_higher")
    "WT-higher" else "KO-higher"
  p <- ggplot(pd, aes(genotype, log1p_cp10k, fill = genotype)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.25,
      colour = "grey30") +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.4,
      fill = "white", colour = "black", stroke = 0.35) +
    facet_wrap(~feature_id, scales = "free_y", ncol = 4,
      labeller = as_labeller(facet_labels)) +
    scale_fill_manual(values = genotype_colours) +
    labs(x = NULL, y = "Expression (log1p CP10K)",
      title = paste0(item$cytokine, " shared-enriched region: ",
        direction_label, " genes"),
      subtitle = paste0("Top ", nrow(genes), " of ",
        item$passing_direction_genes,
        " passing genes in this direction; BH P < 0.05 and ",
        "|log2FC| >= log2(1.5)"),
      caption = paste("White points mark medians; one biological sample",
        "per group.")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 8.5),
      axis.text.x = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 5.5, 5.5, 18))
  plot_height <- max(6.5, 2.55 * ceiling(nrow(genes) / 4))
  save_plot_atomic(p, file.path(report_dir, item$pdf), 12, plot_height, "pdf")
  save_plot_atomic(p, file.path(report_dir, item$png), 12, plot_height, "png")
}

report_path <- file.path(report_dir, "REPORT.html")
html <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
if (!grepl("volcano_labeled_genes.tsv", html, fixed = TRUE)) html <- sub(
  "<h2>Full DEG volcano plots</h2><p>",
  paste0("<h2>Full DEG volcano plots</h2><p>Genes with BH-adjusted P &lt; ",
    "0.05 and |log2FC| &gt;= 1 are labeled. ",
    "<a href='volcano_labeled_genes.tsv'>Download labeled genes</a> | "),
  html, fixed = TRUE)
html <- gsub(" are shown. <a href=", paste0(
  " are shown; facet strips report BH-adjusted P. <a href="),
  html, fixed = TRUE)
atomic_replace(report_path, function(tmp) writeLines(html, tmp))

completion$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
completion$update_script_version <- script_version
completion$figure_colour_scale <- as.list(palette)
completion$violin_colour_scale <- as.list(palette)
completion$violin_facet_annotation <- "BH-adjusted P"
completion$volcano_label_rule <- paste("BH-adjusted P < 0.05 and absolute",
  "standard log2FC >= 1")
completion$volcano_labeled_genes <- nrow(volcano_labels)
completion$checks$requested_colour_scale_applied <- TRUE
completion$checks$violin_padj_annotations_complete <-
  all(is.finite(plot_genes$p_val_adj))
completion$checks$volcano_label_selection_valid <-
  all(volcano_labels$p_val_adj < 0.05 &
    abs(volcano_labels$standard_log2fc) >= 1)
completion$checks$updated_figures_complete <- all(file.exists(c(
  file.path(report_dir, "figures", c("shared_enriched_regions.png",
    "shared_enriched_regions.pdf", "shared_region_de_volcano.png",
    "shared_region_de_volcano.pdf")),
  file.path(report_dir, violin_manifest$png),
  file.path(report_dir, violin_manifest$pdf))))
stopifnot(all(unlist(completion$checks)))
write_json_atomic(completion, completion_path)
cat("Updated plot annotations in place:", report_path, "\n")
cat("Volcano labels:", nrow(volcano_labels), "\n")
