#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste("Usage: 04_flt3_il7_refined_aggregate_sct_umap.R",
    "<bdev_bundle> <final_membership.tsv> <output_dir>",
    "<expected_bundle_fingerprint>"))
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
membership_path <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = TRUE)
expected_fingerprint <- args[[4L]]

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(jsonlite)
})

read_tsv <- function(path) read.delim(path, check.names = FALSE,
  na.strings = c("NA", ""), stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
log_event <- function(...) {
  cat(format(Sys.time(), tz = "UTC"), ..., "\n")
  flush.console()
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)

members <- read_tsv(membership_path)
required <- c("signature_id", "cytokine", "direction", "final_rank",
  "feature_id", "present_in_bdev_sct_data")
stopifnot(all(required %in% names(members)))
members <- members[members$cytokine %in% c("Flt3L", "IL-7") &
  members$direction %in% c("up", "down"), required, drop = FALSE]
members$final_rank <- as.integer(members$final_rank)
members$present_in_bdev_sct_data <- as.logical(members$present_in_bdev_sct_data)
members <- members[order(match(members$cytokine, c("Flt3L", "IL-7")),
  match(members$direction, c("up", "down")), members$final_rank), ]
rownames(members) <- NULL
expected_counts <- c("Flt3L.up" = 45L, "Flt3L.down" = 20L,
  "IL-7.up" = 20L, "IL-7.down" = 20L)
observed <- table(paste(members$cytokine, members$direction, sep = "."))
stopifnot(identical(as.integer(observed[names(expected_counts)]),
  as.integer(expected_counts)), all(members$present_in_bdev_sct_data))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint))
cells <- as.character(wb_cells(bundle, "all"))
stopifnot(length(cells) == 26723L, !anyDuplicated(cells))

space <- wb_space(bundle, "umap.atac", cells = cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells))
  coordinates <- data.frame(cell_id = cells, x = as.numeric(space$x),
    y = as.numeric(space$y), stringsAsFactors = FALSE)
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  coordinates <- data.frame(cell_id = cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]), stringsAsFactors = FALSE)
}

genes <- unique(members$feature_id)
sct_features <- as.character(wb_feature_ids(bundle, "SCT_data"))
stopifnot(all(genes %in% sct_features))
log_event("Reading", length(genes), "unique final genes from SCT_data")
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", genes,
  cells, missing = "NA"))
stopifnot(identical(rownames(expression), genes),
  identical(colnames(expression), cells), all(is.finite(expression)))

panel_order <- c("Flt3L up", "Flt3L down", "IL-7 up", "IL-7 down")
score_rows <- list()
summary_rows <- list()
index <- 0L
for (cytokine in c("Flt3L", "IL-7")) for (direction in c("up", "down")) {
  index <- index + 1L
  local <- members[members$cytokine == cytokine &
    members$direction == direction, , drop = FALSE]
  genes_local <- local$feature_id
  raw_score <- colMeans(expression[genes_local, , drop = FALSE])
  positive <- is.finite(raw_score) & raw_score > 0
  cap <- if (any(positive)) as.numeric(quantile(raw_score[positive], 0.99,
    names = FALSE, type = 7)) else NA_real_
  display_score <- rep(NA_real_, length(raw_score))
  if (is.finite(cap) && cap > 0) {
    display_score[positive] <- pmin(raw_score[positive], cap) / cap
  }
  panel <- paste(cytokine, direction)
  score_rows[[index]] <- data.frame(cell_id = cells, cytokine = cytokine,
    direction = direction, panel = panel, genes = length(genes_local),
    mean_sct_data = raw_score, relative_display_score = display_score,
    stringsAsFactors = FALSE)
  summary_rows[[index]] <- data.frame(cytokine = cytokine,
    direction = direction, genes = length(genes_local), cells = length(cells),
    positive_cells = sum(positive), positive_percent = 100 * mean(positive),
    score_minimum = min(raw_score), score_median = median(raw_score),
    score_mean = mean(raw_score), score_maximum = max(raw_score),
    positive_p99_display_cap = cap, stringsAsFactors = FALSE)
}
scores <- do.call(rbind, score_rows)
summaries <- do.call(rbind, summary_rows)
scores$panel <- factor(scores$panel, levels = panel_order)
plot_data <- merge(scores, coordinates, by = "cell_id", sort = FALSE)
plot_data <- plot_data[order(plot_data$panel,
  ifelse(is.na(plot_data$relative_display_score), -Inf,
    plot_data$relative_display_score)), ]
stopifnot(nrow(plot_data) == length(cells) * 4L)

panel_labels <- setNames(vapply(seq_len(nrow(summaries)), function(i) {
  x <- summaries[i, ]
  paste0(x$cytokine, " ", x$direction, " — ", x$genes, " genes\n",
    sprintf("%.1f%% cells with mean SCT_data > 0", x$positive_percent))
}, character(1)), panel_order)
plot_data$panel_label <- factor(panel_labels[as.character(plot_data$panel)],
  levels = panel_labels[panel_order])

plot <- ggplot(plot_data, aes(x, y)) +
  geom_point(data = plot_data[is.na(plot_data$relative_display_score), ],
    color = "#CDD2DA", size = 0.16, alpha = 0.7, stroke = 0) +
  geom_point(data = plot_data[!is.na(plot_data$relative_display_score), ],
    aes(color = relative_display_score), size = 0.34, alpha = 1, stroke = 0) +
  facet_wrap(~panel_label, ncol = 2) +
  scale_color_gradientn(
    colours = c("#0066FF", "#00D9FF", "#00F060", "#FFF000", "#FF5A00", "#FF00A8"),
    values = c(0, 0.12, 0.30, 0.52, 0.76, 1), limits = c(0, 1),
    oob = scales::squish, na.value = "#CDD2DA",
    breaks = c(0, 0.5, 1), labels = c("low", "mid", "p99+"),
    name = "Mean gene-set\nSCT_data") +
  coord_equal() +
  labs(title = "Refined Flt3L and IL-7 aggregate expression scores",
    subtitle = paste("Mean SCT_data expression across each final gene set on",
      "the existing Bdev ATAC-derived UMAP"),
    caption = paste("Scores are clipped separately at each panel's positive-cell",
      "99th percentile for display. Gray cells have a mean score <= 0.")) +
  theme_void(base_size = 10) +
  theme(strip.text = element_text(size = 11, face = "bold", lineheight = 1.05,
    margin = margin(5, 3, 5, 3)), plot.title = element_text(size = 20,
    face = "bold"), plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 9, color = "grey25", hjust = 0),
    plot.margin = margin(10, 10, 10, 10), legend.position = "bottom",
    legend.key.width = grid::unit(30, "mm"))

figure_file <- "flt3_il7_refined_aggregate_sct_data_atac_umap_bright.png"
ggsave(file.path(output_dir, figure_file), plot, width = 13.5, height = 12,
  dpi = 260, units = "in", bg = "white")
write_tsv(summaries, file.path(output_dir, "aggregate_score_summary.tsv"))
write_tsv_gz(scores[, c("cell_id", "cytokine", "direction", "genes",
  "mean_sct_data", "relative_display_score")], file.path(output_dir,
    "aggregate_cell_scores.tsv.gz"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after),
  file.exists(file.path(output_dir, figure_file)))

input_manifest <- data.frame(
  input = c("Bdev bundle manifest", "final membership", "figure script"),
  path = c(file.path(bundle_dir, "manifest.json"), membership_path, script_path),
  sha256 = vapply(c(file.path(bundle_dir, "manifest.json"), membership_path,
    script_path), wb_source_sha256, character(1)), stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(output_dir, "input_manifest.tsv"))

report <- c("# Aggregate refined-set SCT-data UMAP", "",
  "The figure shows one panel for each refined Flt3L/IL-7 direction.",
  "Each cell's score is the arithmetic mean of its Bdev `SCT_data` values",
  "across every gene in that final set. Scores are not UCell scores.", "",
  "| Cytokine | Direction | Genes | Positive cells | Mean score | Display p99 |",
  "|---|---|---:|---:|---:|---:|")
for (i in seq_len(nrow(summaries))) {
  x <- summaries[i, ]
  report <- c(report, paste0("| ", x$cytokine, " | ", x$direction,
    " | ", x$genes, " | ", sprintf("%.1f%%", x$positive_percent),
    " | ", sprintf("%.6f", x$score_mean), " | ",
    sprintf("%.6f", x$positive_p99_display_cap), " |"))
}
report <- c(report, "", paste0("Figure: `", figure_file, "`."),
  "Cell-level raw and display scores are in `aggregate_cell_scores.tsv.gz`.",
  "No bundle state or feature set was changed or promoted.")
writeLines(report, file.path(output_dir, "REPORT.md"))

completion <- list(status = "complete",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  bdev_bundle = bundle_dir, bdev_bundle_fingerprint = fingerprint_before,
  bdev_bundle_unchanged = TRUE, source_membership = membership_path,
  expression_matrix = "SCT_data",
  aggregate_score = "arithmetic mean SCT_data across all genes in each final set",
  embedding = "umap.atac", embedding_source_modality = "ATAC/peaks",
  cells = length(cells), memberships = nrow(members), unique_genes = length(genes),
  figure = figure_file, result_promotions = 0L,
  scworkbench_version = as.character(packageVersion("scworkbench")))
write_json(completion, file.path(output_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))

output_files <- list.files(output_dir, full.names = TRUE)
output_files <- output_files[basename(output_files) != "output_manifest.tsv"]
output_manifest <- data.frame(file = basename(output_files),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE)
write_tsv(output_manifest[order(output_manifest$file), ],
  file.path(output_dir, "output_manifest.tsv"))
log_event("Completed aggregate SCT-data UMAP at", output_dir)
