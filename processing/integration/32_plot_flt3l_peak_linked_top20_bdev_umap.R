#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "Usage: 32_plot_flt3l_peak_linked_top20_bdev_umap.R LOCAL_DEG_REPORT",
  "BDEV_BUNDLE OUTPUT_DIR EXPECTED_LOCAL_COMPLETION_SHA256",
  "EXPECTED_BDEV_FINGERPRINT SCRIPT_VERSION TOP_N"
), call. = FALSE)

local_report <- normalizePath(args[[1L]], mustWork = TRUE)
bdev_bundle_dir <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
expected_local_sha <- args[[4L]]
expected_bdev_fingerprint <- args[[5L]]
script_version <- args[[6L]]
top_n <- as.integer(args[[7L]])
stopifnot(is.finite(top_n), top_n > 0L)

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}

local_completion <- file.path(local_report, "completion.json")
stopifnot(identical(sha256_file(local_completion), expected_local_sha))
gene_path <- file.path(local_report,
  "all_unique_peak_linked_significant_genes.tsv")
genes <- fread(gene_path)
required <- c("gene", "p_adj_bh", "p_value", "log2fc_flt3lko_over_wt",
  "direction", "linked_peak_count", "linked_peak_ids")
stopifnot(all(required %in% names(genes)), nrow(genes) > 0L,
  !anyDuplicated(genes$gene), all(genes$p_adj_bh < 0.05),
  all(genes$direction %chin% c("WT_higher", "FLT3LKO_higher")))
genes[, abs_log2fc := abs(log2fc_flt3lko_over_wt)]

rank_direction <- function(x, direction_name) {
  y <- copy(x[direction == direction_name])
  setorder(y, -abs_log2fc, p_adj_bh, p_value, gene)
  y[, direction_rank := seq_len(.N)]
  y
}
ranked <- rbindlist(list(
  rank_direction(genes, "WT_higher"),
  rank_direction(genes, "FLT3LKO_higher")
), use.names = TRUE, fill = TRUE)
stopifnot(nrow(ranked) == nrow(genes), !anyDuplicated(ranked$gene))

bundle <- wb_open_bundle(bdev_bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_bdev_fingerprint))

cells <- as.character(wb_cells(bundle, "all"))
stopifnot(length(cells) == 26723L, !anyDuplicated(cells))
space <- wb_space(bundle, "umap.atac", cells = cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells))
  coordinates <- data.table(x = as.numeric(space$x), y = as.numeric(space$y))
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  coordinates <- data.table(x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]))
}
stopifnot(all(is.finite(coordinates$x)), all(is.finite(coordinates$y)))

sct_features <- as.character(wb_feature_ids(bundle, "SCT_data"))
ranked[, present_in_bdev_sct_data := gene %chin% sct_features]
selected <- ranked[present_in_bdev_sct_data == TRUE,
  head(.SD, top_n), by = direction]
selected[, bdev_plot_rank := seq_len(.N), by = direction]
selected <- selected[order(match(direction,
  c("WT_higher", "FLT3LKO_higher")), bdev_plot_rank)]
stopifnot(nrow(selected) == 2L * top_n,
  selected[, all(.N == top_n), by = direction]$V1,
  !anyDuplicated(selected$gene))
selection_cutoffs <- selected[, .(maximum_source_direction_rank =
  max(direction_rank)), by = direction]
skipped <- ranked[present_in_bdev_sct_data == FALSE][selection_cutoffs,
  on = "direction"][direction_rank <= maximum_source_direction_rank]
skipped <- skipped[order(match(direction,
  c("WT_higher", "FLT3LKO_higher")), direction_rank)]
write_tsv(selected, file.path(stage_dir, "top20_per_direction_selection.tsv"))
write_tsv(skipped, file.path(stage_dir,
  "higher_ranked_genes_absent_from_bdev_sct_data.tsv"))

selected_genes <- selected$gene
log_event("Reading", length(selected_genes), "selected genes from Bdev SCT_data")
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data",
  selected_genes, cells, missing = "NA"))
stopifnot(identical(rownames(expression), selected_genes),
  identical(colnames(expression), cells), all(is.finite(expression)))

expression_summaries <- vector("list", nrow(selected))
plot_values <- vector("list", nrow(selected))
for (i in seq_len(nrow(selected))) {
  values <- as.numeric(expression[selected$gene[[i]], ])
  detected <- values > 0
  positive <- values[detected]
  cap <- if (length(positive)) as.numeric(quantile(positive, 0.99,
    names = FALSE, type = 7)) else NA_real_
  scaled <- rep(NA_real_, length(values))
  if (is.finite(cap) && cap > 0)
    scaled[detected] <- pmin(values[detected], cap) / cap
  plot_values[[selected$gene[[i]]]] <- scaled
  expression_summaries[[i]] <- data.table(
    gene = selected$gene[[i]], direction = selected$direction[[i]],
    source_direction_rank = selected$direction_rank[[i]],
    bdev_plot_rank = selected$bdev_plot_rank[[i]], cells = length(values),
    detected_cells = sum(detected), detected_percent = 100 * mean(detected),
    positive_p99_display_cap = cap, minimum_sct_data = min(values),
    maximum_sct_data = max(values))
}
expression_summary <- rbindlist(expression_summaries)
write_tsv(expression_summary,
  file.path(stage_dir, "bdev_sct_expression_summary.tsv"))

direction_specs <- list(
  WT_higher = list(label = "WT-higher", low = "#F4D3C3", high = "#D55E00"),
  FLT3LKO_higher = list(label = "FLT3LKO-higher", low = "#C7DFF0", high = "#0072B2")
)

make_plot <- function(direction_name) {
  local <- selected[direction == direction_name][order(bdev_plot_rank)]
  labels <- paste0(local$bdev_plot_rank, ". ", local$gene, "\nlog2FC ",
    sprintf("%+.2f", local$log2fc_flt3lko_over_wt))
  rows <- lapply(seq_len(nrow(local)), function(i) {
    data.table(x = coordinates$x, y = coordinates$y,
      relative_expression = plot_values[[local$gene[[i]]]],
      gene_panel = labels[[i]], panel_rank = i)
  })
  plot_data <- rbindlist(rows)
  setorder(plot_data, panel_rank, relative_expression, na.last = FALSE)
  plot_data[, gene_panel := factor(gene_panel, levels = labels)]
  spec <- direction_specs[[direction_name]]
  ggplot(plot_data, aes(x, y)) +
    geom_point(data = plot_data[is.na(relative_expression)],
      colour = "#CDD2DA", size = 0.10, alpha = 0.65, stroke = 0) +
    geom_point(data = plot_data[!is.na(relative_expression)],
      aes(colour = relative_expression), size = 0.26, alpha = 1, stroke = 0) +
    facet_wrap(~gene_panel, ncol = 4) +
    scale_colour_gradient(low = spec$low, high = spec$high,
      limits = c(0, 1), oob = scales::squish,
      breaks = c(0, 0.5, 1), labels = c("low", "mid", "p99+"),
      name = "Positive-cell\nSCT_data") +
    coord_equal() +
    labs(title = paste("Top", top_n, spec$label,
      "FLT3L-responsive peak-linked genes"),
      subtitle = paste("Bdev SCT_data on the existing ATAC-derived UMAP;",
        "ranked by absolute FLT3LKO/WT log2FC"),
      caption = paste("Positive values are clipped at each gene's positive-cell",
        "99th percentile and scaled within gene. Gray cells have SCT_data <= 0.",
        "Facet log2FC is FLT3LKO/WT.")) +
    theme_void(base_size = 9) +
    theme(strip.text = element_text(size = 8.5, face = "bold",
      lineheight = 0.95, margin = margin(2, 2, 2, 2)),
      plot.title = element_text(size = 17, face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      plot.caption = element_text(size = 8, colour = "grey25", hjust = 0),
      plot.margin = margin(8, 8, 8, 8), legend.position = "bottom",
      legend.key.width = grid::unit(24, "mm"))
}

figure_rows <- list()
for (direction_name in names(direction_specs)) {
  stem <- paste0(tolower(direction_name), "_top", top_n,
    "_bdev_sct_data_atac_umap")
  plot <- make_plot(direction_name)
  png_path <- file.path(stage_dir, "figures", paste0(stem, ".png"))
  pdf_path <- file.path(stage_dir, "figures", paste0(stem, ".pdf"))
  log_event("Rendering", basename(png_path))
  ggsave(png_path, plot, width = 14, height = 16, dpi = 240,
    units = "in", bg = "white", limitsize = FALSE)
  ggsave(pdf_path, plot, width = 14, height = 16, units = "in",
    bg = "white", useDingbats = FALSE, limitsize = FALSE)
  figure_rows[[length(figure_rows) + 1L]] <- data.table(
    direction = direction_name, genes = top_n,
    png = file.path("figures", basename(png_path)),
    pdf = file.path("figures", basename(pdf_path)))
}
figure_inventory <- rbindlist(figure_rows)
write_tsv(figure_inventory, file.path(stage_dir, "figure_inventory.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after),
  all(file.exists(file.path(stage_dir, figure_inventory$png))),
  all(file.exists(file.path(stage_dir, figure_inventory$pdf))))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
input_paths <- c(local_completion, gene_path,
  file.path(bdev_bundle_dir, "manifest.json"), script_path)
input_manifest <- data.table(
  input = c("local DEG completion", "unique linked-gene table",
    "Bdev bundle manifest", "figure script"),
  path = input_paths,
  sha256 = vapply(input_paths, sha256_file, character(1)))
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))

completion <- list(
  schema_version = "flt3l_peak_linked_top20_bdev_umap.v1",
  script_version = script_version,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  inputs = list(local_deg_report = local_report,
    local_deg_completion_sha256 = expected_local_sha,
    bdev_bundle = bdev_bundle_dir,
    bdev_bundle_fingerprint = fingerprint_before),
  selection = list(top_per_direction = top_n,
    rank_rule = "descending absolute FLT3LKO/WT log2FC among genes present in Bdev SCT_data; ties by BH-adjusted P, raw P, gene",
    directions = c("WT_higher", "FLT3LKO_higher"),
    genes = nrow(selected), higher_ranked_genes_absent_from_bdev_sct_data =
      nrow(skipped)),
  display = list(expression_matrix = "SCT_data", embedding = "umap.atac",
    cells = length(cells), positive_value_cap = "per-gene positive-cell p99",
    nonpositive_cells = "gray", colours = list(
      WT_higher = "#D55E00", FLT3LKO_higher = "#0072B2")),
  figures = split(figure_inventory[, .(png, pdf)], figure_inventory$direction),
  bdev_bundle_unchanged = TRUE, result_promotions = 0L,
  scworkbench_version = as.character(packageVersion("scworkbench")))
write_json(completion, file.path(stage_dir, "completion.json"))
writeLines(c("# Top peak-linked genes on the Bdev UMAP", "",
  paste0("Two figures show the top ", top_n,
    " largest-magnitude genes in each FLT3LKO/WT direction."), "",
  paste0("- Genes come from ", nrow(genes),
    " unique significant genes within ±3 Mb of the validated FLT3L-responsive peaks."),
  "- Ranking is performed separately by direction using absolute FLT3LKO/WT log2FC among genes available in Bdev `SCT_data`.",
  paste0("- ", nrow(skipped),
    " higher-ranked gene(s) lacked Bdev `SCT_data` coverage and are recorded separately."),
  "- Expression is read directly from the accepted Bdev `SCT_data` matrix.",
  "- Coordinates are the existing Bdev ATAC-derived `umap.atac` embedding.",
  "- Positive expression is clipped at each gene's positive-cell p99 and scaled within gene; non-positive cells are gray.",
  "- WT-higher panels use orange and FLT3LKO-higher panels use blue.",
  "- The Bdev bundle was opened read-only and its fingerprint was unchanged.",
  "", "Outputs:", "",
  "- `figures/wt_higher_top20_bdev_sct_data_atac_umap.png` and `.pdf`",
  "- `figures/flt3lko_higher_top20_bdev_sct_data_atac_umap.png` and `.pdf`",
  "- `top20_per_direction_selection.tsv`",
  "- `higher_ranked_genes_absent_from_bdev_sct_data.tsv`",
  "- `bdev_sct_expression_summary.tsv`"),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

output_files <- list.files(stage_dir, recursive = TRUE, full.names = TRUE)
output_files <- output_files[basename(output_files) != "output_manifest.tsv"]
output_manifest <- data.table(
  file = substring(output_files, nchar(stage_dir) + 2L),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, sha256_file, character(1)))
setorder(output_manifest, file)
write_tsv(output_manifest, file.path(stage_dir, "output_manifest.tsv"))

stopifnot(file.rename(stage_dir, output_dir))
log_event("Completed top-per-direction Bdev UMAP report at", output_dir)
