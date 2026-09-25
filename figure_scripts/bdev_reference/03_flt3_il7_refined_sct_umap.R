#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste("Usage: 03_flt3_il7_refined_sct_umap.R",
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
observed_counts <- table(paste(members$cytokine, members$direction, sep = "."))
stopifnot(identical(as.integer(observed_counts[names(expected_counts)]),
  as.integer(expected_counts)), !anyDuplicated(members[, c("cytokine",
    "direction", "feature_id")]), all(members$present_in_bdev_sct_data))
write_tsv(members, file.path(output_dir, "selected_final_membership.tsv"))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint))

cells <- as.character(wb_cells(bundle, "all"))
stopifnot(length(cells) == 26723L, !anyDuplicated(cells))
space <- wb_space(bundle, "umap.atac", cells = cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells))
  coordinates <- data.frame(x = as.numeric(space$x), y = as.numeric(space$y))
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  coordinates <- data.frame(x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]))
}
stopifnot(all(is.finite(coordinates$x)), all(is.finite(coordinates$y)))

genes <- unique(members$feature_id)
sct_features <- as.character(wb_feature_ids(bundle, "SCT_data"))
stopifnot(all(genes %in% sct_features))
log_event("Reading", length(genes), "final genes from SCT_data")
expression <- as.matrix(wb_matrix_read_features(bundle, "SCT_data", genes,
  cells, missing = "NA"))
stopifnot(identical(rownames(expression), genes),
  identical(colnames(expression), cells), all(is.finite(expression)))

summary_rows <- vector("list", nrow(members))
plot_values <- vector("list", nrow(members))
for (i in seq_len(nrow(members))) {
  gene <- members$feature_id[[i]]
  values <- as.numeric(expression[gene, ])
  detected <- values > 0
  positive <- values[detected]
  cap <- if (length(positive)) as.numeric(quantile(positive, 0.99,
    names = FALSE, type = 7)) else NA_real_
  scaled <- rep(NA_real_, length(values))
  if (is.finite(cap) && cap > 0) {
    scaled[detected] <- pmin(values[detected], cap) / cap
  }
  key <- paste(members$cytokine[[i]], members$direction[[i]], gene, sep = "::")
  plot_values[[key]] <- scaled
  summary_rows[[i]] <- data.frame(signature_id = members$signature_id[[i]],
    cytokine = members$cytokine[[i]], direction = members$direction[[i]],
    final_rank = members$final_rank[[i]], feature_id = gene,
    cells = length(values), detected_cells = sum(detected),
    detected_percent = 100 * mean(detected),
    positive_p99_display_cap = cap, minimum_sct_data = min(values),
    maximum_sct_data = max(values), stringsAsFactors = FALSE)
}
expression_summary <- do.call(rbind, summary_rows)
write_tsv(expression_summary, file.path(output_dir, "expression_summary.tsv"))

format_detection <- function(x) {
  ifelse(x < 0.1, sprintf("%.3f%%", x),
    ifelse(x < 1, sprintf("%.2f%%", x), sprintf("%.1f%%", x)))
}

make_plot <- function(cytokine, direction, columns = 5L) {
  local <- members[members$cytokine == cytokine &
    members$direction == direction, , drop = FALSE]
  summaries <- expression_summary[match(paste(local$signature_id,
    local$feature_id), paste(expression_summary$signature_id,
      expression_summary$feature_id)), ]
  stopifnot(nrow(local) == expected_counts[[paste(cytokine, direction,
    sep = ".")]], identical(local$feature_id, summaries$feature_id))
  labels <- paste0(summaries$feature_id, "\n",
    format_detection(summaries$detected_percent), " detected")
  rows <- lapply(seq_len(nrow(local)), function(i) {
    key <- paste(cytokine, direction, local$feature_id[[i]], sep = "::")
    data.frame(x = coordinates$x, y = coordinates$y,
      relative_expression = plot_values[[key]], gene_panel = labels[[i]],
      panel_rank = i, stringsAsFactors = FALSE)
  })
  plot_data <- do.call(rbind, rows)
  plot_data <- plot_data[order(plot_data$panel_rank,
    ifelse(is.na(plot_data$relative_expression), -Inf,
      plot_data$relative_expression)), ]
  plot_data$gene_panel <- factor(plot_data$gene_panel, levels = labels)
  ggplot(plot_data, aes(x, y)) +
    geom_point(data = plot_data[is.na(plot_data$relative_expression), ],
      color = "#CDD2DA", size = 0.10, alpha = 0.65, stroke = 0) +
    geom_point(data = plot_data[!is.na(plot_data$relative_expression), ],
      aes(color = relative_expression), size = 0.26, alpha = 1, stroke = 0) +
    facet_wrap(~gene_panel, ncol = columns) +
    scale_color_gradientn(
      colours = c("#0066FF", "#00D9FF", "#00F060", "#FFF000", "#FF5A00", "#FF00A8"),
      values = c(0, 0.12, 0.30, 0.52, 0.76, 1), limits = c(0, 1),
      oob = scales::squish, na.value = "#CDD2DA",
      breaks = c(0, 0.5, 1), labels = c("low", "mid", "p99+"),
      name = "Positive-cell\nSCT_data") +
    coord_equal() +
    labs(title = paste(cytokine, direction, "refined final genes"),
      subtitle = paste("Bdev SCT_data on the existing ATAC-derived UMAP;",
        "panels follow final set order"),
      caption = paste("Positive values are clipped at each gene's positive-cell",
        "99th percentile, scaled within gene, and overlaid as larger saturated",
        "points. Gray cells have SCT_data <= 0.")) +
    theme_void(base_size = 9) +
    theme(strip.text = element_text(size = 8.2, face = "bold",
      lineheight = 0.95, margin = margin(2, 2, 2, 2)),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 11),
      plot.caption = element_text(size = 8, color = "grey25", hjust = 0),
      plot.margin = margin(8, 8, 8, 8), legend.position = "bottom",
      legend.key.width = grid::unit(24, "mm"))
}

figure_rows <- list()
for (cytokine in c("Flt3L", "IL-7")) for (direction in c("up", "down")) {
  n <- expected_counts[[paste(cytokine, direction, sep = ".")]]
  panel_rows <- ceiling(n / 5L)
  filename <- paste0(tolower(gsub("[^A-Za-z0-9]+", "", cytokine)), "_",
    direction, "_refined_sct_data_atac_umap_bright.png")
  path <- file.path(output_dir, filename)
  log_event("Rendering", filename)
  ggsave(path, make_plot(cytokine, direction), width = 15,
    height = 4.0 + 2.25 * panel_rows, dpi = 240, units = "in",
    bg = "white", limitsize = FALSE)
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(
    cytokine = cytokine, direction = direction, genes = n,
    file = filename, stringsAsFactors = FALSE)
}
figure_inventory <- do.call(rbind, figure_rows)
write_tsv(figure_inventory, file.path(output_dir, "figure_inventory.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after),
  all(file.exists(file.path(output_dir, figure_inventory$file))))

input_manifest <- data.frame(
  input = c("Bdev bundle manifest", "final membership", "figure script"),
  path = c(file.path(bundle_dir, "manifest.json"), membership_path, script_path),
  sha256 = vapply(c(file.path(bundle_dir, "manifest.json"), membership_path,
    script_path), wb_source_sha256, character(1)), stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(output_dir, "input_manifest.tsv"))

report <- c("# Separate bright SCT-data UMAP figures", "",
  "Four figures display the refined Flt3L and IL-7 final genes separately by cytokine and direction.", "",
  "- Expression values are read directly from Bdev `SCT_data`.",
  "- Coordinates are the existing ATAC-derived `umap.atac` embedding.",
  "- Cells with `SCT_data <= 0` are light gray.",
  "- Every positive value uses a bright cyan-to-magenta scale after within-gene p99 clipping.",
  "- No bundle state or feature set was changed or promoted.", "",
  "| Cytokine | Direction | Genes | Figure |", "|---|---|---:|---|")
for (i in seq_len(nrow(figure_inventory))) {
  x <- figure_inventory[i, ]
  report <- c(report, paste0("| ", x$cytokine, " | ", x$direction,
    " | ", x$genes, " | `", x$file, "` |"))
}
writeLines(report, file.path(output_dir, "REPORT.md"))

completion <- list(status = "complete",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  bdev_bundle = bundle_dir, bdev_bundle_fingerprint = fingerprint_before,
  bdev_bundle_unchanged = TRUE, source_membership = membership_path,
  expression_matrix = "SCT_data", embedding = "umap.atac",
  embedding_source_modality = "ATAC/peaks", cells = length(cells),
  genes = nrow(members), figures = as.list(setNames(figure_inventory$file,
    paste(figure_inventory$cytokine, figure_inventory$direction, sep = "."))),
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
log_event("Completed four separate bright UMAP figures at", output_dir)
