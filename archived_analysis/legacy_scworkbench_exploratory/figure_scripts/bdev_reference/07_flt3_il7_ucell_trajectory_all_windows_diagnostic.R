#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(paste("Usage: 07_flt3_il7_ucell_trajectory_all_windows_diagnostic.R",
    "<source_trajectory_report> <output_dir>"))
}
source_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = TRUE)

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
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)

summary_path <- file.path(source_dir, "sliding_window_summary.tsv")
change_path <- file.path(source_dir, "change_point_up_locations.tsv")
completion_path <- file.path(source_dir, "completion.json")
stopifnot(file.exists(summary_path), file.exists(change_path),
  file.exists(completion_path))
source_completion <- read_json(completion_path, simplifyVector = TRUE)
stopifnot(identical(source_completion$status, "complete"),
  identical(source_completion$window_x_coordinate,
    "(window_left + window_right) / 2"))

d <- read_tsv(summary_path)
change_points <- read_tsv(change_path)
required <- c("feature_set", "cytokine", "direction", "panel_label",
  "window_index", "window_left", "window_right", "window_midpoint",
  "n_cells", "eligible_window", "score_q25", "score_median", "score_q75")
stopifnot(all(required %in% names(d)), nrow(d) == 404L,
  length(unique(d$feature_set)) == 4L, all(is.finite(d$score_median)),
  all(abs(d$window_midpoint - (d$window_left + d$window_right) / 2) < 1e-12),
  all(c("up1", "up2") %in% change_points$label))
d$direction <- factor(d$direction, levels = c("up", "down"))
panel_order <- c("Flt3L up", "Flt3L down", "IL-7 up", "IL-7 down")
label_map <- setNames(vapply(panel_order, function(key) {
  as.character(unique(d$panel_label[d$feature_set == key]))
}, character(1)), panel_order)
d$panel_label <- factor(d$panel_label, levels = label_map)
low <- d[!d$eligible_window, , drop = FALSE]
stopifnot(nrow(low) == 24L,
  identical(sort(unique(low$window_index)), c(41L, 42L, 98L, 99L, 100L, 101L)),
  identical(sort(unique(low$n_cells)), c(2L, 3L, 23L, 25L)))

palette <- c(up = "#D73027", down = "#2166AC")
caption <- paste0(
  "All windows are retained. Line: median; band: cell-level IQR. ",
  "Hollow points mark windows below the usual 30-cell threshold (n = 2–25).\n",
  "Window width: 5% of pt range; step: 1%; x is the actual window midpoint. ",
  "Dotted lines: median pt of Bdev bb_up groups up1 and up2.")
plot_curves <- function(data, title, facet) {
  p <- ggplot(data, aes(window_midpoint, score_median,
      group = feature_set, color = direction, fill = direction)) +
    geom_vline(data = change_points,
      aes(xintercept = change_point_up_location), inherit.aes = FALSE,
      linewidth = 0.6, linetype = "dotted", color = "#111827") +
    geom_ribbon(aes(ymin = score_q25, ymax = score_q75),
      alpha = 0.22, colour = NA) +
    geom_line(linewidth = 1.0) +
    geom_point(data = low[low$feature_set %in% unique(data$feature_set), ],
      aes(window_midpoint, score_median, color = direction),
      inherit.aes = FALSE, shape = 21, fill = "white", stroke = 0.8,
      size = 2.1) +
    scale_color_manual(values = palette) + scale_fill_manual(values = palette) +
    labs(title = title,
      subtitle = "Diagnostic view retaining windows with fewer than 30 cells",
      x = "Existing Bdev pseudotime (pt)", y = "UCell score",
      color = "Direction", fill = "Direction", caption = caption) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
      legend.position = "bottom", strip.text = element_text(face = "bold"),
      plot.title.position = "plot")
  if (facet == "panel") p + facet_wrap(~panel_label, scales = "free_y", ncol = 2)
  else p + facet_wrap(~direction, scales = "free_y", nrow = 1,
    labeller = labeller(direction = c(up = "UP", down = "DOWN")))
}

overview <- plot_curves(d,
  "Refined Flt3L and IL-7 UCell trajectories — all windows retained", "panel")
overview_png <- file.path(output_dir,
  "flt3_il7_refined_ucell_trajectory_all_windows_diagnostic.png")
overview_pdf <- file.path(output_dir,
  "flt3_il7_refined_ucell_trajectory_all_windows_diagnostic.pdf")
ggsave(overview_png, overview, width = 13, height = 9.5, dpi = 240,
  units = "in", bg = "white")
ggsave(overview_pdf, overview, width = 13, height = 9.5,
  units = "in", bg = "white")

figure_rows <- list(data.frame(cytokine = "overview",
  file = basename(overview_png), stringsAsFactors = FALSE))
for (cytokine in c("Flt3L", "IL-7")) {
  local <- d[d$cytokine == cytokine, , drop = FALSE]
  p <- plot_curves(local,
    paste(cytokine, "refined UCell trajectories — all windows retained"),
    "direction")
  slug <- tolower(gsub("[^A-Za-z0-9]+", "", cytokine))
  file <- paste0(slug,
    "_refined_ucell_trajectory_all_windows_diagnostic.png")
  ggsave(file.path(output_dir, file), p, width = 11.5, height = 6.5,
    dpi = 240, units = "in", bg = "white")
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(cytokine = cytokine,
    file = file, stringsAsFactors = FALSE)
}
figures <- do.call(rbind, figure_rows)
write_tsv(figures, file.path(output_dir, "figure_inventory.tsv"))
write_tsv(low[, required], file.path(output_dir,
  "retained_low_support_windows.tsv"))

input_paths <- c(summary_path, change_path, completion_path, script_path)
input_manifest <- data.frame(input = c("sliding-window summary",
  "bb_up change-point locations", "source completion record", "figure script"),
  path = input_paths,
  sha256 = vapply(input_paths, wb_source_sha256, character(1)),
  stringsAsFactors = FALSE)
write_tsv(input_manifest, file.path(output_dir, "input_manifest.tsv"))
report <- c("# All-window trajectory diagnostic", "",
  "This diagnostic redraws the Flt3L and IL-7 native-UCell trajectories while",
  "retaining the six windows that did not meet the usual 30-cell threshold.",
  "Hollow points identify the low-support estimates; the original native",
  "window medians and IQRs are unchanged.", "",
  "- Windows 41 and 42 contain 23 and 25 cells.",
  "- Windows 98–101 contain 3, 2, 2 and 2 cells.",
  "- This is a sensitivity display, not a replacement for the minimum-support primary figure.")
writeLines(report, file.path(output_dir, "REPORT.md"))
completion <- list(status = "complete", promoted = FALSE,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  source_trajectory_report = source_dir,
  source_bdev_bundle_fingerprint = source_completion$bdev_bundle_fingerprint,
  score_recalculated = FALSE, window_summaries_recalculated = FALSE,
  low_support_windows_retained = 6L,
  low_support_cell_counts = sort(unique(low$n_cells)),
  low_support_marker = "hollow point", result_promotions = 0L,
  overview_figure = basename(overview_png))
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
cat("Completed all-window diagnostic at", output_dir, "\n")
