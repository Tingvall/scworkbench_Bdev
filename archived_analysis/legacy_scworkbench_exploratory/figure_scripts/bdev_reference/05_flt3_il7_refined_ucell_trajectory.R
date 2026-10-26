#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste("Usage: 05_flt3_il7_refined_ucell_trajectory.R",
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
source("src/sliding_window_signature_summary.R")

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
  "feature_id")
stopifnot(all(required %in% names(members)))
members <- members[members$cytokine %in% c("Flt3L", "IL-7") &
  members$direction %in% c("up", "down"), required, drop = FALSE]
members$final_rank <- as.integer(members$final_rank)
members <- members[order(match(members$cytokine, c("Flt3L", "IL-7")),
  match(members$direction, c("up", "down")), members$final_rank), ]
expected_keys <- c("Flt3L up", "Flt3L down", "IL-7 up", "IL-7 down")
expected_counts <- c("Flt3L up" = 45L, "Flt3L down" = 20L,
  "IL-7 up" = 20L, "IL-7 down" = 20L)
members$key <- paste(members$cytokine, members$direction)
stopifnot(identical(as.integer(table(factor(members$key,
  levels = expected_keys))), as.integer(expected_counts)))
signatures <- setNames(lapply(expected_keys, function(key) {
  members$feature_id[members$key == key]
}), expected_keys)

bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_fingerprint))
cells <- as.character(wb_cells(bundle, "all"))
metadata <- wb_meta(bundle, cells = cells, cols = list(
  wb_ref_metadata_field("core:cell_id"), "pt", "bb_up"))
stopifnot(length(cells) == 26723L,
  identical(as.character(metadata$cell_id), cells),
  sum(is.finite(metadata$pt)) == 12642L, sum(is.na(metadata$pt)) == 14081L)
rna_features <- as.character(wb_feature_ids(bundle, "RNA_counts"))
stopifnot(all(unlist(signatures, use.names = FALSE) %in% rna_features))
pt <- data.frame(cell_id = cells[is.finite(metadata$pt)],
  value = as.numeric(metadata$pt[is.finite(metadata$pt)]),
  stringsAsFactors = FALSE)
change_labels <- sort(unique(na.omit(as.character(metadata$bb_up))))
stopifnot(identical(change_labels, c("up1", "up2")))
change_points <- do.call(rbind, lapply(change_labels, function(label) {
  selected <- !is.na(metadata$bb_up) & metadata$bb_up == label &
    is.finite(metadata$pt)
  values <- as.numeric(metadata$pt[selected])
  stopifnot(length(values) > 0L)
  data.frame(metadata_field = "bb_up", label = label,
    n_cells = length(values), pt_min = min(values),
    change_point_up_location = median(values), pt_max = max(values),
    location_summary = "median original pt among cells with this bb_up label",
    stringsAsFactors = FALSE)
}))

# Match the earlier trajectory workflow: native raw-count UCell in an
# unpromoted durable session, followed by native distribution summaries.
session_id <- "flt3_il7_refined_top100_ucell_v1"
session_root <- file.path(output_dir, "sessions")
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) wb_session_resume(bundle, workspace) else
  wb_session_create(bundle, workspace_root = session_root, session_id = session_id)
if (!is.null(session$durable$promotion_receipt)) stop("Review session was promoted")
score_runs <- names(session$runs)[vapply(session$runs, function(x)
  identical(x$op_name, "ucell"), logical(1))]
matching <- score_runs[vapply(score_runs, function(id) {
  tab <- session$runs[[id]]$outputs$signature_summary$table
  setequal(as.character(tab$signature), names(signatures))
}, logical(1))]
if (!length(matching)) {
  if (length(score_runs)) stop("Unexpected UCell run in review session")
  log_event("Scoring four final sets with native scWorkbench UCell")
  scored <- wb_session_run_op(session, bundle, "ucell", inputs = list(
    matrix = "RNA_counts", cell_set = "all", cell_ids = cells,
    min_cells = 10L, features = signatures, name = session_id,
    max_rank = 1500L, ncores = 4L, parallel_backend = "auto",
    chunk_size = 400L, missing_features_policy = "error"),
    notes = "Unpromoted native Bdev UCell scores for four refined Flt3L/IL-7 trajectory plots")
  session <- scored$session
  run <- scored$run
} else {
  stopifnot(length(matching) == 1L)
  run <- session$runs[[matching[[1L]]]]
  log_event("Reusing completed native UCell score run")
}
native_summary <- run$outputs$signature_summary$table
native_scores <- run$outputs$signature_scores
scores <- data.frame(key = as.character(native_scores$key),
  cell_id = as.character(native_scores$cell_id),
  value = as.numeric(native_scores$value), stringsAsFactors = FALSE)
stopifnot(length(scores$value) == length(cells) * 4L,
  setequal(unique(scores$key), expected_keys), all(is.finite(scores$value)),
  all(scores$value >= 0 & scores$value <= 1),
  all(native_summary$n_features_missing == 0L))
for (key in expected_keys) {
  local <- scores[scores$key == key, ]
  stopifnot(identical(local$cell_id, cells), !anyDuplicated(local$cell_id))
}
saveRDS(native_scores, file.path(output_dir, "ucell_scores.rds"))
write_tsv_gz(scores, file.path(output_dir, "ucell_cell_scores.tsv.gz"))
write_tsv(native_summary, file.path(output_dir, "native_ucell_summary.tsv"))

width_fraction <- 0.05
step_fraction <- 0.01
min_cells <- 30L
log_event("Calculating native 5% sliding-window median and IQR summaries")
result <- sliding_window_signature_summary(scores, pt,
  width_fractions = width_fraction, step_fraction = step_fraction,
  min_cells = min_cells)
summaries <- result$summaries
windows <- result$windows
label_info <- unique(members[, c("key", "cytokine", "direction")])
label_info$genes <- expected_counts[label_info$key]
label_info <- label_info[match(expected_keys, label_info$key), ]
label_map <- setNames(paste0(label_info$cytokine, " ",
  toupper(label_info$direction), " — ", label_info$genes, " genes"),
  label_info$key)
summaries$direction <- factor(label_info$direction[
  match(summaries$feature_set, label_info$key)], levels = c("up", "down"))
summaries$cytokine <- label_info$cytokine[
  match(summaries$feature_set, label_info$key)]
summaries$panel_label <- factor(label_map[summaries$feature_set],
  levels = label_map[expected_keys])
summaries$window_midpoint <- (summaries$window_left +
  summaries$window_right) / 2
windows$window_midpoint <- (windows$window_left + windows$window_right) / 2

low <- windows[!windows$eligible_window, , drop = FALSE]
limits <- range(pt$value)
half_step <- result$parameters$step_pt_units / 2
low$left <- pmax(limits[[1L]], low$window_midpoint - half_step)
low$right <- pmin(limits[[2L]], low$window_midpoint + half_step)
palette <- c(up = "#D73027", down = "#2166AC")
caption <- sprintf(paste0(
  "Window: 5%% of pt range (%.3f units); step: 1%% (%.3f units). ",
  "Line: median; band: 25th–75th percentiles (IQR).\n",
  "At least %d cells per window; endpoint windows are clipped. ",
  "Each x value is (window left + window right)/2. ",
  "Dotted lines: median pt for Bdev bb_up groups up1 and up2. ",
  "Bands show cell-to-cell spread, not confidence intervals."),
  width_fraction * diff(limits), result$parameters$step_pt_units, min_cells)

plot_curves <- function(data, title, facet) {
  p <- ggplot(data[data$eligible_window, ],
      aes(window_midpoint, score_median, group = feature_set,
        color = direction, fill = direction)) +
    geom_rect(data = low,
      aes(xmin = left, xmax = right, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "#64748B", alpha = 0.13) +
    geom_vline(data = change_points,
      aes(xintercept = change_point_up_location), inherit.aes = FALSE,
      linewidth = 0.6, linetype = "dotted", color = "#111827") +
    geom_ribbon(aes(ymin = score_q25, ymax = score_q75),
      alpha = 0.22, colour = NA) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = palette) + scale_fill_manual(values = palette) +
    labs(title = title,
      subtitle = paste("Pooled", length(pt$value),
        "cells with defined original pt; 5% sliding-window median and IQR"),
      x = "Existing Bdev pseudotime (pt)", y = "UCell score",
      color = "Direction", fill = "Direction", caption = caption) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
      legend.position = "bottom", strip.text = element_text(face = "bold"),
      plot.title.position = "plot")
  if (facet == "panel") p + facet_wrap(~panel_label, scales = "free_y", ncol = 2)
  else p + facet_wrap(~direction, scales = "free_y", nrow = 1,
    labeller = labeller(direction = c(up = "UP", down = "DOWN")))
}

overview <- plot_curves(summaries,
  "Refined Flt3L and IL-7 UCell scores along Bdev pseudotime", "panel")
overview_png <- file.path(output_dir,
  "flt3_il7_refined_ucell_trajectory_5pct_iqr.png")
overview_pdf <- file.path(output_dir,
  "flt3_il7_refined_ucell_trajectory_5pct_iqr.pdf")
ggsave(overview_png, overview, width = 13, height = 9.5, dpi = 240,
  units = "in", bg = "white")
ggsave(overview_pdf, overview, width = 13, height = 9.5,
  units = "in", bg = "white")

figure_rows <- list()
for (cytokine in c("Flt3L", "IL-7")) {
  p <- plot_curves(summaries[summaries$cytokine == cytokine, ],
    paste(cytokine, "refined UCell scores along Bdev pseudotime"), "direction")
  slug <- tolower(gsub("[^A-Za-z0-9]+", "", cytokine))
  filename <- paste0(slug, "_refined_ucell_trajectory_5pct_iqr.png")
  ggsave(file.path(output_dir, filename), p, width = 11.5, height = 6.5,
    dpi = 240, units = "in", bg = "white")
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(cytokine = cytokine,
    file = filename, stringsAsFactors = FALSE)
}
figure_inventory <- do.call(rbind, figure_rows)

write_tsv(summaries, file.path(output_dir, "sliding_window_summary.tsv"))
write_tsv(windows, file.path(output_dir, "window_support.tsv"))
write_tsv(change_points, file.path(output_dir,
  "change_point_up_locations.tsv"))
write_tsv(data.frame(cell_id = pt$cell_id, pt = pt$value),
  file.path(output_dir, "source_pseudotime.tsv"))
write_tsv(figure_inventory, file.path(output_dir, "figure_inventory.tsv"))
saveRDS(result$membership, file.path(output_dir, "window_membership.rds"))
saveRDS(result$native_summaries,
  file.path(output_dir, "native_window_summaries.rds"))
write_json(result$parameters, file.path(output_dir, "parameters.json"))
support_summary <- data.frame(width_percent = 5,
  width_pt_units = width_fraction * diff(limits), centers = nrow(windows),
  retained_centers = sum(windows$eligible_window),
  unsupported_centers = sum(!windows$eligible_window),
  minimum_cells = min(windows$n_cells), maximum_cells = max(windows$n_cells))
write_tsv(support_summary, file.path(output_dir, "window_support_summary.tsv"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, fingerprint_after),
  file.exists(overview_png), file.exists(overview_pdf),
  all(file.exists(file.path(output_dir, figure_inventory$file))))
source_files <- c(script_path, "src/sliding_window_signature_summary.R")
input_manifest <- data.frame(input = c("Bdev bundle manifest",
  "final membership", "trajectory script", "sliding-window helper"),
  path = c(file.path(bundle_dir, "manifest.json"), membership_path, source_files),
  sha256 = vapply(c(file.path(bundle_dir, "manifest.json"), membership_path,
    source_files), wb_source_sha256, character(1)), stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(output_dir, "input_manifest.tsv"))

report <- c("# Refined Flt3L/IL-7 UCell scores along Bdev pseudotime", "",
  paste0("The analysis uses all **", length(pt$value),
    "** Bdev cells with defined existing `pt`; the **",
    sum(is.na(metadata$pt)), "** missing source values remain missing."), "",
  "- Score: native scWorkbench UCell on raw Bdev `RNA_counts`, `max_rank = 1500`.",
  "- Window: 5% of the observed pseudotime range; step: 1%.",
  "- Summary: median line and type-7 25th–75th percentile band.",
  "- Minimum support: 30 cells; endpoint windows are clipped.",
  "- X coordinate: the arithmetic midpoint of each actual window (`(window_left + window_right) / 2`), including clipped endpoint windows.",
  "- Vertical dotted lines: median original `pt` for cells labeled `up1` or `up2` in the source Bdev `bb_up` metadata column.",
  "- No graph smoothing, extra curve smoother, pseudotime refitting, imputation, bootstrap, or promotion.", "",
  paste0("All **", support_summary$retained_centers, "/",
    support_summary$centers, "** window centers met the support requirement."), "",
  "## Figures", "",
  "- `flt3_il7_refined_ucell_trajectory_5pct_iqr.png`: four-panel overview.",
  "- `flt3l_refined_ucell_trajectory_5pct_iqr.png`: Flt3L up/down.",
  "- `il7_refined_ucell_trajectory_5pct_iqr.png`: IL-7 up/down.",
  "- A PDF copy of the overview is included.", "",
  "The IQR describes cell-to-cell spread within overlapping windows and is not a confidence interval.")
writeLines(report, file.path(output_dir, "REPORT.md"))

completion <- list(status = "complete", promoted = FALSE,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  bdev_bundle = bundle_dir, bdev_bundle_fingerprint = fingerprint_before,
  bdev_bundle_unchanged = TRUE, source_membership = membership_path,
  score = "native scWorkbench UCell on RNA_counts", max_rank = 1500L,
  scoring_session_workspace = workspace,
  pseudotime_source = "Bdev@meta.data$pt", defined_pt_cells = length(pt$value),
  missing_pt_cells = sum(is.na(metadata$pt)), width_fraction = width_fraction,
  step_fraction = step_fraction, min_cells = min_cells,
  center = "median", band = "IQR (type-7 q25 to q75)",
  window_x_coordinate = "(window_left + window_right) / 2",
  change_point_up_source = "Bdev@meta.data$bb_up",
  change_point_up_location = "median original pt within each non-missing bb_up group",
  overview_figure = basename(overview_png), result_promotions = 0L,
  scworkbench_version = as.character(packageVersion("scworkbench")))
write_json(completion, file.path(output_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))

output_files <- list.files(output_dir, full.names = TRUE)
output_files <- output_files[basename(output_files) != "output_manifest.tsv" &
  basename(output_files) != "sessions"]
output_manifest <- data.frame(file = basename(output_files),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE)
write_tsv(output_manifest[order(output_manifest$file), ],
  file.path(output_dir, "output_manifest.tsv"))
log_event("Completed 5% native-UCell trajectory report at", output_dir)
