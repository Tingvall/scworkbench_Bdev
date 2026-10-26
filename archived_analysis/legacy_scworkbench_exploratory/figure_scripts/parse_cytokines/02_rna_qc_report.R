#!/usr/bin/env Rscript

# Render factual RNA-QC distributions and neutral empirical-tail sensitivity
# previews from promoted scWorkbench artifacts. This script is read-only with
# respect to the bundle and never creates a cell set.
#
# Usage:
#   Rscript --vanilla 02_rna_qc_report.R BUNDLE_DIR OUTPUT_DIR

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 02_rna_qc_report.R BUNDLE_DIR OUTPUT_DIR", call. = FALSE)
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
if (file.exists(output_dir) || dir.exists(output_dir)) {
  stop("Output destination already exists; refusing to overwrite: ", output_dir, call. = FALSE)
}

suppressPackageStartupMessages({
  library(scworkbench)
  library(ggplot2)
  library(jsonlite)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Could not create output directory: ", output_dir, call. = FALSE)

write_tsv <- function(value, filename) {
  utils::write.table(
    value, file.path(output_dir, filename), sep = "\t", quote = TRUE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
}

markdown_table <- function(value, digits = 8L) {
  if (!nrow(value)) return("_No rows._")
  format_value <- function(x) {
    if (is.numeric(x)) {
      out <- ifelse(is.na(x), "NA", format(signif(x, digits), trim = TRUE, scientific = FALSE, big.mark = ","))
    } else {
      out <- ifelse(is.na(x), "NA", as.character(x))
    }
    gsub("|", "\\\\|", out, fixed = TRUE)
  }
  formatted <- lapply(value, format_value)
  header <- paste0("| ", paste(names(value), collapse = " | "), " |")
  divider <- paste0("| ", paste(rep("---", ncol(value)), collapse = " | "), " |")
  rows <- vapply(seq_len(nrow(value)), function(i) {
    paste0("| ", paste(vapply(formatted, `[[`, character(1L), i), collapse = " | "), " |")
  }, character(1L))
  paste(c(header, divider, rows), collapse = "\n")
}

artifact_row <- function(artifacts, type) {
  rows <- artifacts[
    artifacts$source_op == "rna_qc_metrics" &
      artifacts$type == type &
      grepl("parse_rna_qc", artifacts$name, fixed = TRUE),
    , drop = FALSE
  ]
  if (nrow(rows) != 1L) {
    stop("Expected exactly one current parse_rna_qc artifact of type ", type, "; found ", nrow(rows), call. = FALSE)
  }
  rows
}

metric_labels <- c(
  total_counts = "Total RNA counts",
  detected_features = "Detected RNA features",
  mitochondrial_counts = "Mitochondrial counts",
  mitochondrial_percent = "Mitochondrial percent",
  ribosomal_counts = "Rp[sl]-prefix counts",
  ribosomal_percent = "Rp[sl]-prefix percent"
)

# The processing stage records a successful recipient validation and finalized
# fingerprint. This read-only renderer verifies those identities and artifacts
# directly without repeating the full, expensive recipient audit.
bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest <- wb_manifest(bundle)
stopifnot(identical(manifest$bundle_id, "parse_cytokines"))

all_cells <- as.character(wb_cells(bundle, cell_set = "all"))
imported_cells <- as.character(wb_cells(bundle, cell_set = "set:parse_cytokines"))
stopifnot(identical(all_cells, imported_cells), length(all_cells) == 251599L)
cell_sets <- sort(as.character(wb_cell_sets(bundle)))

artifacts <- wb_art_list(bundle)
metrics_row <- artifact_row(artifacts, "cell_vector_family")
features_row <- artifact_row(artifacts, "feature_table")
summary_row <- artifact_row(artifacts, "obs_table")
artifact_ids <- c(
  metrics = as.character(metrics_row$artifact_id[[1L]]),
  qc_features = as.character(features_row$artifact_id[[1L]]),
  summary = as.character(summary_row$artifact_id[[1L]])
)

qc_features_primitive <- wb_art_load(bundle, artifact_ids[["qc_features"]], expected_kind = "feature_table")
qc_summary_primitive <- wb_art_load(bundle, artifact_ids[["summary"]], expected_kind = "obs_table")
qc_features <- qc_features_primitive$table
qc_summary <- qc_summary_primitive$table
stopifnot(identical(qc_summary_primitive$meta$threshold_selection, "not_performed"))
write_tsv(qc_features, "qc_feature_membership.tsv")

metrics <- wb_art_load(bundle, artifact_ids[["metrics"]], expected_kind = "cell_vector_family")
metric_values <- lapply(names(metric_labels), function(key) {
  rows <- metrics$key == key
  values <- as.numeric(metrics$value[rows])
  names(values) <- as.character(metrics$cell_id[rows])
  stopifnot(
    identical(names(values), all_cells),
    length(values) == length(all_cells),
    all(is.finite(values)),
    all(values >= 0)
  )
  as.numeric(values)
})
names(metric_values) <- names(metric_labels)
stopifnot(all(metric_values$total_counts >= 500))

metadata <- wb_meta(
  bundle,
  cells = all_cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "sample_id", "condition", "replicate")
)
stopifnot(
  identical(as.character(metadata$cell_id), all_cells),
  !anyNA(metadata$sample_id),
  length(unique(metadata$sample_id)) == 96L
)

quantile_probabilities <- c(0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975, 0.99, 0.995, 0.999, 1)
quantile_labels <- c("min", "p00.1", "p00.5", "p01", "p02.5", "p05", "p10", "p25", "p50", "p75", "p90", "p95", "p97.5", "p99", "p99.5", "p99.9", "max")

overall_quantiles <- do.call(rbind, lapply(names(metric_values), function(key) {
  data.frame(
    metric = key,
    metric_label = unname(metric_labels[[key]]),
    quantile = quantile_labels,
    probability = quantile_probabilities,
    value = as.numeric(stats::quantile(metric_values[[key]], quantile_probabilities, names = FALSE, type = 7)),
    stringsAsFactors = FALSE
  )
}))
write_tsv(overall_quantiles, "rna_qc_overall_quantiles.tsv")

sample_ids <- unique(as.character(metadata$sample_id))
sample_indices <- split(seq_along(all_cells), factor(metadata$sample_id, levels = sample_ids))
sample_probabilities <- c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99)
sample_labels <- c("p01", "p05", "p25", "p50", "p75", "p95", "p99")
sample_summary <- do.call(rbind, lapply(names(metric_values), function(key) {
  do.call(rbind, lapply(seq_along(sample_indices), function(i) {
    indices <- sample_indices[[i]]
    values <- metric_values[[key]][indices]
    sample_rows <- unique(metadata[indices, c("sample_id", "condition", "replicate"), drop = FALSE])
    stopifnot(nrow(sample_rows) == 1L)
    quantiles <- as.numeric(stats::quantile(values, sample_probabilities, names = FALSE, type = 7))
    data.frame(
      sample_id = as.character(sample_rows$sample_id),
      condition = as.character(sample_rows$condition),
      replicate = as.character(sample_rows$replicate),
      cells = length(indices),
      metric = key,
      metric_label = unname(metric_labels[[key]]),
      mean = mean(values),
      stats::setNames(as.list(quantiles), sample_labels),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
}))
rownames(sample_summary) <- NULL
write_tsv(sample_summary, "rna_qc_by_sample_summary.tsv")

# These are empirical tail-sensitivity calculations, not recommended or chosen
# biological thresholds. Strict inequalities mirror a future min/max rule and
# make the effect of ties explicit.
tail_probabilities <- c(0.001, 0.005, 0.01, 0.025, 0.05)
tail_sensitivity <- do.call(rbind, lapply(names(metric_values), function(key) {
  values <- metric_values[[key]]
  lower <- lapply(tail_probabilities, function(probability) {
    threshold <- as.numeric(stats::quantile(values, probability, names = FALSE, type = 1))
    excluded <- sum(values < threshold)
    data.frame(
      metric = key, direction = "lower", empirical_tail_probability = probability,
      candidate_threshold = threshold, candidate_rule = paste0(key, " < ", format(threshold, scientific = FALSE, trim = TRUE)),
      cells_excluded = excluded, cells_retained = length(values) - excluded,
      percent_excluded = 100 * excluded / length(values), stringsAsFactors = FALSE
    )
  })
  upper <- lapply(tail_probabilities, function(probability) {
    threshold <- as.numeric(stats::quantile(values, 1 - probability, names = FALSE, type = 1))
    excluded <- sum(values > threshold)
    data.frame(
      metric = key, direction = "upper", empirical_tail_probability = probability,
      candidate_threshold = threshold, candidate_rule = paste0(key, " > ", format(threshold, scientific = FALSE, trim = TRUE)),
      cells_excluded = excluded, cells_retained = length(values) - excluded,
      percent_excluded = 100 * excluded / length(values), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, c(lower, upper))
}))
tail_sensitivity$status <- "preview_only_not_selected_or_applied"
tail_sensitivity$cell_ids_materialized <- FALSE
write_tsv(tail_sensitivity, "rna_qc_tail_sensitivity.tsv")

histogram_rows <- function(values, metric, transform = c("log1p10", "identity"), bins = 100L) {
  transform <- match.arg(transform)
  plotted <- if (transform == "log1p10") log10(1 + values) else values
  breaks <- if (transform == "identity" && all(plotted >= 0) && all(plotted <= 100)) {
    seq(0, 100, length.out = bins + 1L)
  } else {
    seq(min(plotted), max(plotted), length.out = bins + 1L)
  }
  h <- graphics::hist(plotted, breaks = breaks, plot = FALSE, include.lowest = TRUE, right = TRUE)
  data.frame(
    metric = metric,
    metric_label = unname(metric_labels[[metric]]),
    midpoint = h$mids,
    width = diff(h$breaks),
    count = h$counts,
    stringsAsFactors = FALSE
  )
}

count_metrics <- c("total_counts", "detected_features", "mitochondrial_counts", "ribosomal_counts")
count_histograms <- do.call(rbind, lapply(count_metrics, function(key) {
  histogram_rows(metric_values[[key]], key, transform = "log1p10")
}))
percent_metrics <- c("mitochondrial_percent", "ribosomal_percent")
percent_histograms <- do.call(rbind, lapply(percent_metrics, function(key) {
  histogram_rows(metric_values[[key]], key, transform = "identity")
}))

plot_theme <- theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8, hjust = 0),
    axis.text.x = element_text(size = 8)
  )

count_plot <- ggplot(count_histograms, aes(x = midpoint, y = count, width = width)) +
  geom_col(fill = "#2B6CB0", colour = NA) +
  facet_wrap(~metric_label, scales = "free", ncol = 2) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    title = "Parse RNA QC count distributions",
    subtitle = "All 251,599 imported bundle cells; x-axis is log10(1 + value)",
    x = "log10(1 + metric value)", y = "Cells",
    caption = "Factual promoted scWorkbench metrics. No threshold or exclusion is shown or applied."
  ) +
  plot_theme

percent_plot <- ggplot(percent_histograms, aes(x = midpoint, y = count, width = width)) +
  geom_col(fill = "#2F855A", colour = NA) +
  facet_wrap(~metric_label, scales = "free_y", ncol = 1) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    title = "Parse RNA QC percentage distributions",
    subtitle = "All 251,599 imported bundle cells; exact feature membership is reported separately",
    x = "Percent of total RNA counts", y = "Cells",
    caption = "The Rp[sl] metric is prefix-defined, not a curated structural-ribosome annotation. No threshold is applied."
  ) +
  plot_theme

sample_plot_data <- sample_summary[sample_summary$metric %in% c(
  "total_counts", "detected_features", "mitochondrial_percent", "ribosomal_percent"
), , drop = FALSE]
sample_plot_data$sample_id <- factor(sample_plot_data$sample_id, levels = sample_ids)
sample_plot_data$plot_p05 <- sample_plot_data$p05
sample_plot_data$plot_p50 <- sample_plot_data$p50
sample_plot_data$plot_p95 <- sample_plot_data$p95
count_rows <- sample_plot_data$metric %in% c("total_counts", "detected_features")
sample_plot_data$plot_p05[count_rows] <- log10(1 + sample_plot_data$plot_p05[count_rows])
sample_plot_data$plot_p50[count_rows] <- log10(1 + sample_plot_data$plot_p50[count_rows])
sample_plot_data$plot_p95[count_rows] <- log10(1 + sample_plot_data$plot_p95[count_rows])

sample_plot <- ggplot(sample_plot_data, aes(x = sample_id, y = plot_p50)) +
  geom_linerange(aes(ymin = plot_p05, ymax = plot_p95), colour = "grey55", linewidth = 0.3) +
  geom_point(aes(colour = condition), size = 0.9, show.legend = FALSE) +
  facet_wrap(~metric_label, scales = "free_y", ncol = 1) +
  labs(
    title = "RNA QC distributions by source sample",
    subtitle = "Point: median; interval: p05-p95. Count metrics are log10(1 + value); percentages are untransformed.",
    x = "Source sample (source order)", y = "Metric summary",
    caption = "Per-sample distributions are descriptive and do not define sample- or cell-level exclusions."
  ) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))

ggsave(file.path(output_dir, "rna_qc_count_distributions.png"), count_plot, width = 10, height = 7.5, dpi = 240, bg = "white")
ggsave(file.path(output_dir, "rna_qc_percentage_distributions.png"), percent_plot, width = 10, height = 7.5, dpi = 240, bg = "white")
ggsave(file.path(output_dir, "rna_qc_by_sample.png"), sample_plot, width = 22, height = 12, dpi = 240, bg = "white")

qc_feature_counts <- as.data.frame(table(qc_features$metric_group), stringsAsFactors = FALSE)
names(qc_feature_counts) <- c("metric_group", "matched_features")
universe_table <- data.frame(
  population = c("Imported cells", "QC exclusion population", "Accepted QC analysis population"),
  durable_cell_set = c("set:parse_cytokines", "not defined", "not defined"),
  status = c("retained", "not proposed", "awaiting user-selected thresholds"),
  cells = c(length(imported_cells), NA_integer_, NA_integer_),
  stringsAsFactors = FALSE
)

display_quantiles <- overall_quantiles[overall_quantiles$quantile %in% c(
  "min", "p01", "p05", "p25", "p50", "p75", "p95", "p99", "max"
), c("metric_label", "quantile", "value"), drop = FALSE]
display_quantiles <- reshape(
  display_quantiles, idvar = "metric_label", timevar = "quantile", direction = "wide"
)
names(display_quantiles) <- sub("^value\\.", "", names(display_quantiles))
display_quantiles <- display_quantiles[, c("metric_label", "min", "p01", "p05", "p25", "p50", "p75", "p95", "p99", "max")]

report_path <- file.path(output_dir, "rna_qc_report.md")
report <- paste0(
  "# Parse cytokine RNA QC distributions\n\n",
  "Generated: `", format(Sys.time(), tz = "UTC", usetz = TRUE), "`\n\n",
  "## Decision boundary\n\n",
  "This report describes factual RNA QC metrics for every cell already imported into the `parse_cytokines` bundle. ",
  "It does **not** choose a filtering threshold, propose a final exclusion rule, remove a cell, or derive an accepted analysis population. ",
  "The empirical tail table is a neutral sensitivity aid: each row reports what a strict one-metric tail rule would exclude, and none is selected or applied.\n\n",
  markdown_table(universe_table), "\n\n",
  "The imported bundle universe is the previously accepted >=500-total-count source object (251,599 cells). ",
  "That historical source boundary remains part of provenance; this QC stage retains all 251,599 imported cells.\n\n",
  "## Promoted scWorkbench evidence\n\n",
  "- Metrics artifact: `", artifact_ids[["metrics"]], "`\n",
  "- Feature-membership artifact: `", artifact_ids[["qc_features"]], "`\n",
  "- Summary artifact: `", artifact_ids[["summary"]], "`\n",
  "- Bundle state fingerprint: `", manifest$current_state$fingerprint, "`\n",
  "- Threshold selection: `not_performed`\n\n",
  "Feature matching is case-sensitive: `^mt-` for mitochondrial features and `^Rp[sl]` for the ribosomal-prefix metric. ",
  "The latter is an auditable symbol-prefix definition and is not represented as a curated structural-ribosome gene set.\n\n",
  markdown_table(qc_feature_counts), "\n\n",
  "Exact membership is in `qc_feature_membership.tsv`.\n\n",
  "## Overall distributions\n\n",
  markdown_table(display_quantiles), "\n\n",
  "![Count distributions](rna_qc_count_distributions.png)\n\n",
  "![Percentage distributions](rna_qc_percentage_distributions.png)\n\n",
  "## Per-sample distributions\n\n",
  "The per-sample figure shows medians and p05-p95 intervals for all 96 source samples. ",
  "The complete numeric table retains sample, condition, replicate, cell count, mean, and seven quantiles.\n\n",
  "![Per-sample distributions](rna_qc_by_sample.png)\n\n",
  "## Exclusion preview status\n\n",
  "No exclusion rule is proposed as selected, and no cell IDs are materialized as exclusions. ",
  "`rna_qc_tail_sensitivity.tsv` lists strict lower- and upper-tail effects at empirical 0.1%, 0.5%, 1%, 2.5%, and 5% points for review. ",
  "Because ties are retained at the boundary, the observed excluded percentage can be smaller than the named tail probability. ",
  "After explicit thresholds are supplied, the next step is to show the resulting union of proposed exclusions (including overlap and per-sample effects) before creating any candidate or accepted cell set.\n"
)
writeLines(report, report_path, useBytes = TRUE)

completion <- list(
  status = "complete",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  bundle_id = "parse_cytokines",
  bundle_dir = bundle_dir,
  bundle_state_fingerprint = manifest$current_state$fingerprint,
  imported_cell_set = "set:parse_cytokines",
  imported_cells = length(imported_cells),
  cell_sets = cell_sets,
  artifact_ids = as.list(artifact_ids),
  metric_keys = names(metric_labels),
  mitochondrial_features = sum(qc_features$metric_group == "mitochondrial"),
  ribosomal_prefix_features = sum(qc_features$metric_group == "ribosomal"),
  samples = length(sample_ids),
  threshold_selection = "not_performed",
  exclusion_set_created = FALSE,
  accepted_analysis_set_created = FALSE,
  tail_sensitivity_status = "preview_only_not_selected_or_applied",
  outputs = list(
    report = report_path,
    feature_membership = file.path(output_dir, "qc_feature_membership.tsv"),
    overall_quantiles = file.path(output_dir, "rna_qc_overall_quantiles.tsv"),
    sample_summary = file.path(output_dir, "rna_qc_by_sample_summary.tsv"),
    tail_sensitivity = file.path(output_dir, "rna_qc_tail_sensitivity.tsv")
  )
)
jsonlite::write_json(
  completion, file.path(output_dir, "completion.json"), pretty = TRUE,
  auto_unbox = TRUE, null = "null", na = "null", digits = NA
)
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
wb_close_bundle(bundle)
bundle <- NULL

cat(report_path, "\n")
