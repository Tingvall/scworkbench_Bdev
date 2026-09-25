#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: 02_il7_flt3_top25_core_atac_umap.R",
    "<bdev_bundle> <top25_core_membership.tsv> <output_dir>",
    "<expected_bundle_fingerprint>"
  ))
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
core_path <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = TRUE)
expected_fingerprint <- args[[4L]]

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(patchwork)
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

core <- read_tsv(core_path)
required_columns <- c("signature_id", "cytokine", "direction", "variant",
  "core_rank", "feature_id")
stopifnot(all(required_columns %in% names(core)))
selected <- core[core$variant == "top25" &
  core$cytokine %in% c("IL-7", "Flt3L") &
  core$direction %in% c("up", "down"), required_columns, drop = FALSE]
selected$core_rank <- as.integer(selected$core_rank)
selected <- selected[order(match(selected$direction, c("up", "down")),
  match(selected$cytokine, c("IL-7", "Flt3L")), selected$core_rank), ]
rownames(selected) <- NULL

expected_counts <- c("IL-7.up" = 17L, "Flt3L.up" = 17L,
  "IL-7.down" = 7L, "Flt3L.down" = 7L)
observed_counts <- table(paste(selected$cytokine, selected$direction, sep = "."))
stopifnot(identical(as.integer(observed_counts[names(expected_counts)]),
  as.integer(expected_counts)), !anyDuplicated(selected[, c("cytokine",
    "direction", "feature_id")]))
write_tsv(selected, file.path(output_dir, "selected_top25_core_membership.tsv"))

bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint))

cells <- as.character(wb_cells(bundle, "all"))
stopifnot(length(cells) == 26723L, !anyDuplicated(cells))
space_info <- wb_space_info(bundle, "umap.atac")
stopifnot(nrow(space_info) == 1L, as.integer(space_info$dims[[1L]]) == 2L)
space <- wb_space(bundle, "umap.atac", cells = cells)
if (is.data.frame(space)) {
  stopifnot(identical(as.character(space$cell_id), cells),
    all(c("x", "y") %in% names(space)))
  coordinates <- data.frame(cell_id = cells, x = as.numeric(space$x),
    y = as.numeric(space$y), stringsAsFactors = FALSE)
} else {
  stopifnot(identical(rownames(space), cells), ncol(space) == 2L)
  coordinates <- data.frame(cell_id = cells, x = as.numeric(space[, 1L]),
    y = as.numeric(space[, 2L]), stringsAsFactors = FALSE)
}
stopifnot(all(is.finite(coordinates$x)), all(is.finite(coordinates$y)))

genes <- unique(as.character(selected$feature_id))
matrix_features <- as.character(wb_feature_ids(bundle, "RNA_counts"))
stopifnot(all(genes %in% matrix_features))
log_event("Reading", length(genes), "unique core genes from RNA_counts")
counts <- wb_matrix_read_features(bundle, "RNA_counts", genes, cells,
  missing = "NA")
stopifnot(identical(rownames(counts), genes),
  identical(colnames(counts), cells), all(is.finite(counts)),
  all(counts >= 0))
depth_metadata <- wb_meta(bundle, cells = cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "nCount_RNA"))
stopifnot(identical(as.character(depth_metadata$cell_id), cells))
rna_depth <- as.numeric(depth_metadata$nCount_RNA)
stopifnot(length(rna_depth) == length(cells), all(is.finite(rna_depth)),
  all(rna_depth > 0), all(rna_depth >= colSums(counts)))
expression <- log1p(sweep(counts, 2L, rna_depth, "/") * 10000)
stopifnot(all(is.finite(expression)), all(expression >= 0))

summary_rows <- vector("list", nrow(selected))
plot_values <- vector("list", nrow(selected))
for (i in seq_len(nrow(selected))) {
  gene <- selected$feature_id[[i]]
  values <- as.numeric(expression[gene, ])
  detected <- values > 0
  positive <- values[detected]
  cap <- if (length(positive)) {
    as.numeric(quantile(positive, 0.99, names = FALSE, type = 7))
  } else 1
  if (!is.finite(cap) || cap <= 0) cap <- max(c(positive, 1))
  relative <- pmin(values, cap) / cap
  relative[!detected] <- NA_real_
  summary_rows[[i]] <- data.frame(
    signature_id = selected$signature_id[[i]],
    cytokine = selected$cytokine[[i]],
    direction = selected$direction[[i]],
    core_rank = selected$core_rank[[i]],
    feature_id = gene,
    cells = length(values),
    detected_cells = sum(detected),
    detected_percent = 100 * mean(detected),
    nonzero_p99_display_cap = cap,
    maximum_expression = max(values),
    stringsAsFactors = FALSE
  )
  plot_values[[i]] <- relative
}
expression_summary <- do.call(rbind, summary_rows)
names(plot_values) <- paste(selected$cytokine, selected$direction,
  selected$feature_id, sep = "::")
write_tsv(expression_summary, file.path(output_dir, "expression_summary.tsv"))

make_group_plot <- function(cytokine, direction, columns) {
  membership <- selected[selected$cytokine == cytokine &
    selected$direction == direction, , drop = FALSE]
  summaries <- expression_summary[expression_summary$cytokine == cytokine &
    expression_summary$direction == direction, , drop = FALSE]
  stopifnot(identical(membership$feature_id, summaries$feature_id))
  format_detection <- function(x) {
    ifelse(x < 0.1, sprintf("%.3f%%", x),
      ifelse(x < 1, sprintf("%.2f%%", x), sprintf("%.1f%%", x)))
  }
  labels <- paste0(summaries$feature_id, "\n",
    format_detection(summaries$detected_percent), " detected")
  rows <- lapply(seq_len(nrow(membership)), function(i) {
    key <- paste(cytokine, direction, membership$feature_id[[i]], sep = "::")
    data.frame(
      x = coordinates$x,
      y = coordinates$y,
      relative_expression = plot_values[[key]],
      gene_panel = labels[[i]],
      panel_rank = i,
      stringsAsFactors = FALSE
    )
  })
  plot_data <- do.call(rbind, rows)
  plot_data <- plot_data[order(plot_data$panel_rank,
    ifelse(is.na(plot_data$relative_expression), -Inf,
      plot_data$relative_expression)), ]
  plot_data$gene_panel <- factor(plot_data$gene_panel, levels = labels)
  ggplot(plot_data, aes(x, y, color = relative_expression)) +
    geom_point(size = 0.10, alpha = 0.90, stroke = 0) +
    facet_wrap(~gene_panel, ncol = columns) +
    scale_color_viridis_c(
      option = "magma", limits = c(0, 1), oob = scales::squish,
      na.value = "#E5E7EB", breaks = c(0, 0.5, 1),
      labels = c("low", "mid", "p99+"),
      name = "Within-gene\nexpression"
    ) +
    coord_equal() +
    labs(title = paste(cytokine, paste0("top25 ", direction, " core"))) +
    theme_void(base_size = 8) +
    theme(
      strip.text = element_text(size = 7.5, face = "bold", lineheight = 0.95,
        margin = margin(1.5, 1.5, 1.5, 1.5)),
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      plot.margin = margin(4, 4, 4, 4),
      legend.position = "bottom",
      legend.key.width = grid::unit(18, "mm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
}

render_direction <- function(direction, width, height, columns) {
  il7 <- make_group_plot("IL-7", direction, columns)
  flt3 <- make_group_plot("Flt3L", direction, columns)
  combined <- (il7 / flt3) +
    plot_layout(guides = "collect", heights = c(1, 1)) +
    plot_annotation(
      title = paste0("IL-7 and Flt3L top25 core genes — ", direction),
      subtitle = paste(
        "Bdev log1p(CP10k) RNA expression on the existing ATAC-derived UMAP",
        "(umap.atac); panels follow core rank"
      ),
      caption = paste(
        "Log1p(CP10k) expression is scaled separately per gene to its 99th",
        "percentile among detected cells. Gray indicates zero expression;",
        "the strip gives the fraction of Bdev cells with expression > 0."
      ),
      theme = theme(
        plot.title = element_text(size = 17, face = "bold"),
        plot.subtitle = element_text(size = 11),
        plot.caption = element_text(size = 8, color = "grey30", hjust = 0)
      )
    ) & theme(legend.position = "bottom")
  filename <- file.path(output_dir,
    paste0("il7_flt3_top25_core_", direction, "_atac_umap.png"))
  log_event("Rendering", filename)
  ggsave(filename, combined, width = width, height = height, dpi = 220,
    units = "in", bg = "white", limitsize = FALSE)
  filename
}

up_figure <- render_direction("up", width = 14.5, height = 19.5, columns = 5L)
down_figure <- render_direction("down", width = 12.5, height = 10.5, columns = 4L)

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after),
  file.exists(up_figure), file.exists(down_figure))

input_manifest <- data.frame(
  input = c("Bdev bundle manifest", "top25 core membership", "figure script"),
  path = c(file.path(bundle_dir, "manifest.json"), core_path, script_path),
  sha256 = vapply(c(file.path(bundle_dir, "manifest.json"), core_path,
    script_path), wb_source_sha256, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(input_manifest, file.path(output_dir, "input_manifest.tsv"))

report <- c(
  "# IL-7 and Flt3L top25 core-gene expression on Bdev ATAC UMAP",
  "",
  "This read-only sanity check maps normalized Bdev RNA expression onto the",
  "existing `umap.atac` embedding. It uses the top25 three-fold core memberships",
  "from the unpromoted Bdev-compatible Parse review.",
  "",
  sprintf("- Up figure: `%s` (%d IL-7 and %d Flt3L genes).",
    basename(up_figure), expected_counts[["IL-7.up"]],
    expected_counts[["Flt3L.up"]]),
  sprintf("- Down figure: `%s` (%d IL-7 and %d Flt3L genes).",
    basename(down_figure), expected_counts[["IL-7.down"]],
    expected_counts[["Flt3L.down"]]),
  "- Expression: `log1p(RNA_counts / nCount_RNA * 10,000)`.",
  "- Embedding: `umap.atac`, imported from the peaks assay.",
  "- Display: zero-expression cells are gray; positive values are clipped at",
  "  each gene's non-zero 99th percentile and scaled from 0 to 1.",
  "- `expression_summary.tsv` records detection and unscaled expression ranges.",
  "- No bundle state or feature set was changed or promoted.",
  "",
  paste0("Bdev bundle fingerprint: `", fingerprint_before, "`.")
)
writeLines(report, file.path(output_dir, "REPORT.md"))

completion <- list(
  status = "complete",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  bdev_bundle = bundle_dir,
  bdev_bundle_fingerprint = fingerprint_before,
  bdev_bundle_unchanged = TRUE,
  expression_matrix = "RNA_counts",
  expression_transform = "log1p(count / nCount_RNA * 10000)",
  embedding = "umap.atac",
  embedding_source_modality = "ATAC/peaks",
  cells = length(cells),
  selected_core_memberships = nrow(selected),
  unique_genes = length(genes),
  gene_counts = as.list(expected_counts),
  figures = c(up = basename(up_figure), down = basename(down_figure)),
  result_promotions = 0L,
  scworkbench_version = as.character(packageVersion("scworkbench"))
)
write_json(completion, file.path(output_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))

output_files <- list.files(output_dir, full.names = TRUE)
output_files <- output_files[basename(output_files) != "output_manifest.tsv"]
output_manifest <- data.frame(
  file = basename(output_files),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(output_manifest[order(output_manifest$file), ],
  file.path(output_dir, "output_manifest.tsv"))
log_event("Completed Bdev ATAC-UMAP core-gene figures at", output_dir)
