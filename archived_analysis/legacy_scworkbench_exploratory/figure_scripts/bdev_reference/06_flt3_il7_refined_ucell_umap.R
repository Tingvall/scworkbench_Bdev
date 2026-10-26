#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste("Usage: 06_flt3_il7_refined_ucell_umap.R",
    "<bdev_bundle> <trajectory_report_dir> <output_dir>",
    "<expected_bundle_fingerprint>"))
}
bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
score_dir <- normalizePath(args[[2L]], mustWork = TRUE)
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
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)

score_path <- file.path(score_dir, "ucell_cell_scores.tsv.gz")
native_summary_path <- file.path(score_dir, "native_ucell_summary.tsv")
source_completion_path <- file.path(score_dir, "completion.json")
stopifnot(file.exists(score_path), file.exists(native_summary_path),
  file.exists(source_completion_path))
source_completion <- read_json(source_completion_path, simplifyVector = TRUE)
stopifnot(identical(source_completion$status, "complete"),
  identical(source_completion$score, "native scWorkbench UCell on RNA_counts"),
  identical(as.integer(source_completion$max_rank), 1500L),
  identical(source_completion$bdev_bundle_fingerprint, expected_fingerprint))

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

scores <- read_tsv(score_path)
required <- c("key", "cell_id", "value")
stopifnot(all(required %in% names(scores)))
panel_order <- c("Flt3L up", "Flt3L down", "IL-7 up", "IL-7 down")
stopifnot(setequal(unique(scores$key), panel_order),
  nrow(scores) == length(cells) * length(panel_order),
  all(is.finite(scores$value)), all(scores$value >= 0 & scores$value <= 1))
for (key in panel_order) {
  local <- scores[scores$key == key, , drop = FALSE]
  stopifnot(identical(as.character(local$cell_id), cells),
    !anyDuplicated(local$cell_id))
}

native_summary <- read_tsv(native_summary_path)
stopifnot(all(c("signature", "n_features_input", "n_features_matched",
  "n_features_missing") %in% names(native_summary)),
  setequal(native_summary$signature, panel_order),
  all(native_summary$n_features_missing == 0L))
gene_counts <- setNames(as.integer(native_summary$n_features_matched),
  native_summary$signature)

summaries <- do.call(rbind, lapply(panel_order, function(key) {
  values <- scores$value[scores$key == key]
  positive <- values > 0
  cap <- as.numeric(quantile(values[positive], 0.99, names = FALSE, type = 7))
  data.frame(key = key, genes = gene_counts[[key]], cells = length(values),
    positive_cells = sum(positive), positive_percent = 100 * mean(positive),
    score_minimum = min(values), score_q25 = quantile(values, 0.25,
      names = FALSE, type = 7), score_median = median(values),
    score_mean = mean(values), score_q75 = quantile(values, 0.75,
      names = FALSE, type = 7), score_maximum = max(values),
    positive_p99_display_cap = cap, stringsAsFactors = FALSE)
}))
scores$display_score <- NA_real_
for (key in panel_order) {
  ix <- scores$key == key & scores$value > 0
  cap <- summaries$positive_p99_display_cap[summaries$key == key]
  scores$display_score[ix] <- pmin(scores$value[ix], cap) / cap
}
scores$panel <- factor(scores$key, levels = panel_order)
plot_data <- merge(scores, coordinates, by = "cell_id", sort = FALSE)
plot_data <- plot_data[order(plot_data$panel,
  ifelse(is.na(plot_data$display_score), -Inf, plot_data$display_score)), ]
stopifnot(nrow(plot_data) == length(cells) * 4L)

panel_labels <- setNames(vapply(panel_order, function(key) {
  x <- summaries[summaries$key == key, ]
  paste0(key, " — ", x$genes, " genes\n",
    sprintf("%.1f%% cells with UCell > 0", x$positive_percent))
}, character(1)), panel_order)
plot_data$panel_label <- factor(panel_labels[as.character(plot_data$panel)],
  levels = panel_labels[panel_order])

add_points_and_scale <- function(p, data) {
  p +
    geom_point(data = data[is.na(data$display_score), ],
      color = "#CDD2DA", size = 0.16, alpha = 0.7, stroke = 0) +
    geom_point(data = data[!is.na(data$display_score), ],
      aes(color = display_score), size = 0.34, alpha = 1, stroke = 0) +
    scale_color_gradientn(
      colours = c("#0066FF", "#00D9FF", "#00F060", "#FFF000", "#FF5A00", "#FF00A8"),
      values = c(0, 0.12, 0.30, 0.52, 0.76, 1), limits = c(0, 1),
      oob = scales::squish, na.value = "#CDD2DA",
      breaks = c(0, 0.5, 1), labels = c("low", "mid", "p99+"),
      name = "UCell score\n(panel scale)") +
    coord_equal() + theme_void(base_size = 10) +
    theme(strip.text = element_text(size = 11, face = "bold", lineheight = 1.05,
      margin = margin(5, 3, 5, 3)), plot.title = element_text(size = 20,
      face = "bold"), plot.subtitle = element_text(size = 12),
      plot.caption = element_text(size = 9, color = "grey25", hjust = 0),
      plot.margin = margin(10, 10, 10, 10), legend.position = "bottom",
      legend.key.width = grid::unit(30, "mm"))
}
caption <- paste("Native scWorkbench UCell on Bdev RNA_counts (max_rank = 1500).",
  "Colors are scaled separately to each panel's positive-score 99th percentile;",
  "gray cells have score 0.")
overview <- add_points_and_scale(ggplot(plot_data, aes(x, y)), plot_data) +
  facet_wrap(~panel_label, ncol = 2) +
  labs(title = "Refined Flt3L and IL-7 UCell scores",
    subtitle = "Four final gene sets on the existing Bdev ATAC-derived UMAP",
    caption = caption)
overview_file <- "flt3_il7_refined_ucell_atac_umap_bright.png"
ggsave(file.path(output_dir, overview_file), overview, width = 13.5,
  height = 12, dpi = 260, units = "in", bg = "white")

figure_rows <- list(data.frame(key = "overview", file = overview_file,
  stringsAsFactors = FALSE))
for (key in panel_order) {
  local <- plot_data[plot_data$key == key, , drop = FALSE]
  p <- add_points_and_scale(ggplot(local, aes(x, y)), local) +
    labs(title = paste(key, "UCell score"),
      subtitle = paste(gene_counts[[key]],
        "genes on the existing Bdev ATAC-derived UMAP"), caption = caption)
  slug <- tolower(gsub("[^A-Za-z0-9]+", "_", key))
  file <- paste0(slug, "_refined_ucell_atac_umap_bright.png")
  ggsave(file.path(output_dir, file), p, width = 8.5, height = 7.5,
    dpi = 260, units = "in", bg = "white")
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(key = key,
    file = file, stringsAsFactors = FALSE)
}
figures <- do.call(rbind, figure_rows)
write_tsv(figures, file.path(output_dir, "figure_inventory.tsv"))
write_tsv(summaries, file.path(output_dir, "ucell_score_summary.tsv"))
write_tsv_gz(scores[, c("cell_id", "key", "value", "display_score")],
  file.path(output_dir, "ucell_cell_scores_with_display.tsv.gz"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after),
  all(file.exists(file.path(output_dir, figures$file))))
input_paths <- c(file.path(bundle_dir, "manifest.json"), score_path,
  native_summary_path, source_completion_path, script_path)
input_manifest <- data.frame(input = c("Bdev bundle manifest",
  "native UCell cell scores", "native UCell summary",
  "trajectory completion record", "figure script"), path = input_paths,
  sha256 = vapply(input_paths, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(output_dir, "input_manifest.tsv"))

report <- c("# Refined-set UCell UMAP", "",
  "The overview and four separate figures show the saved native scWorkbench",
  "UCell scores for the refined Flt3L and IL-7 sets on Bdev `umap.atac`.",
  "Raw scores are unchanged. Display colors are normalized separately within",
  "each panel to its positive-score 99th percentile.", "",
  "The score source is Bdev `RNA_counts`, with `max_rank = 1500`.",
  "No UCell score was recalculated, and no bundle object was changed or promoted.")
writeLines(report, file.path(output_dir, "REPORT.md"))
completion <- list(status = "complete",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  bdev_bundle = bundle_dir, bdev_bundle_fingerprint = fingerprint_before,
  bdev_bundle_unchanged = TRUE, source_trajectory_report = score_dir,
  score = "reused native scWorkbench UCell on RNA_counts", max_rank = 1500L,
  embedding = "umap.atac", embedding_source_modality = "ATAC/peaks",
  cells = length(cells), gene_sets = length(panel_order),
  display_scale = "within-panel positive-score p99", figure = overview_file,
  separate_figures = figures$file[figures$key != "overview"],
  result_promotions = 0L,
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
cat("Completed refined UCell UMAP report at", output_dir, "\n")
