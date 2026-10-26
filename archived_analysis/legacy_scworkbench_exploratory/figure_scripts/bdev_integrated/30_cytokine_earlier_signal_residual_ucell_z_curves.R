#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 30_cytokine_earlier_signal_residual_ucell_z_curves.R",
  "<three_method_curves.tsv.gz> <earlier_window_statistics.tsv.gz>",
  "<row_membership.tsv> <bb.tsv> <earlier_completion.json> <output_dir>"))
comparison_path <- normalizePath(args[[1L]], mustWork = TRUE)
earlier_path <- normalizePath(args[[2L]], mustWork = TRUE)
membership_path <- normalizePath(args[[3L]], mustWork = TRUE)
bb_path <- normalizePath(args[[4L]], mustWork = TRUE)
earlier_completion_path <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  as.data.table(read.delim(con, check.names = FALSE,
    stringsAsFactors = FALSE))
}
short_label <- function(x) {
  x <- sub("^parse_cytokine_rest_only_refinement_all_top100_v1\\.", "", x)
  toupper(gsub("\\.", " ", x))
}

current <- read_tsv_gz(comparison_path)
earlier <- read_tsv_gz(earlier_path)
membership <- fread(membership_path)
bb <- fread(bb_path)
earlier_completion <- read_json(earlier_completion_path,
  simplifyVector = TRUE)
stopifnot(earlier_completion$tested_sets == 38L,
  identical(earlier_completion$score_method, "depth-adjusted mean"),
  uniqueN(earlier$feature_set) == 38L,
  setequal(unique(earlier$feature_set), unique(current$feature_set)))

curves <- rbindlist(list(
  earlier[, .(feature_set, window_index, window_center,
    method = "Earlier depth-adjusted cytokine signal",
    original_value = score_mean)],
  current[method == "RNA-feature residualized UCell",
    .(feature_set, window_index, window_center, method,
      original_value)],
  current[method == "Original expression z-score",
    .(feature_set, window_index, window_center, method,
      original_value)]))
curves[, scaled_0_1 := {
  limits <- range(original_value)
  if (diff(limits) > 0) (original_value - limits[[1L]]) / diff(limits) else 0
}, by = .(feature_set, method)]
stopifnot(min(curves$scaled_0_1) >= 0, max(curves$scaled_0_1) <= 1)

setorder(membership, row_position_top_to_bottom)
display_order <- membership$feature_set
stopifnot(setequal(display_order, unique(curves$feature_set)))
curves[, set_label := factor(short_label(feature_set),
  levels = short_label(display_order))]
curves[, method := factor(method, levels = c(
  "Earlier depth-adjusted cytokine signal",
  "RNA-feature residualized UCell", "Original expression z-score"))]

p <- ggplot(curves, aes(window_center, scaled_0_1, color = method,
    linetype = method)) +
  geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
    linetype = "dotted", linewidth = .28, color = "grey35") +
  geom_line(linewidth = .45) +
  facet_wrap(~set_label, ncol = 5) +
  scale_color_manual(values = c(
    "Earlier depth-adjusted cytokine signal" = "#252525",
    "RNA-feature residualized UCell" = "#D7301F",
    "Original expression z-score" = "#2166AC")) +
  scale_linetype_manual(values = c(
    "Earlier depth-adjusted cytokine signal" = "solid",
    "RNA-feature residualized UCell" = "solid",
    "Original expression z-score" = "longdash")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1),
    expand = expansion(mult = c(.01, .03))) +
  labs(title = "Cytokine signals: earlier depth-adjusted score, residualized UCell, and expression z-score",
    subtitle = paste("Each cytokine/method curve independently min-max scaled",
      "to 0-1; dotted lines mark BB1 and BB2"),
    x = "Existing Bdev pseudotime (pt)",
    y = "Within-curve scaled value (0-1)", color = NULL,
    linetype = NULL,
    caption = paste("Black: earlier mean log1p(CP10K) cytokine score with",
      "conditional RNA-depth component removed; scaling is visualization-only.")) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
    strip.text = element_text(size = 6.5),
    plot.title = element_text(face = "bold"))

stem <- "cytokine_earlier_depth_adjusted_residual_ucell_original_z_scaled_0_1_curves"
ggsave(file.path(stage_dir, paste0(stem, ".png")), p,
  width = 17, height = 21, dpi = 240, bg = "white")
ggsave(file.path(stage_dir, paste0(stem, ".pdf")), p,
  width = 17, height = 21, bg = "white")
fwrite(curves, file.path(stage_dir, paste0(stem, ".tsv.gz")),
  sep = "\t", quote = FALSE, compress = "gzip")

writeLines(c("# Earlier cytokine signal / residual-UCell / z-score curves", "",
  "- Black: the earlier depth-adjusted mean cytokine score from `reviewed_cytokine_38_positive_peak_windows_depth_adjusted_v2_20261011`.",
  "- The earlier score is mean log1p(CP10K) expression across signature genes after subtracting only the fitted conditional RNA-depth component from `score ~ ns(pt, df=8) + scaled log1p(nCount_RNA)`.",
  "- Red: cell-level UCell residualized using `UCell ~ ns(log1p(nFeature_RNA), df=4)`, summarized by window median.",
  "- Blue dashed: the later per-gene depth-adjusted, per-gene standardized expression trajectory z-score.",
  "- Each cytokine/method curve is independently min-max scaled to 0-1 for visualization only.",
  paste0("- Earlier analysis has ", uniqueN(earlier$window_index),
    " eligible windows (minimum 30 cells); the other two have ",
    uniqueN(current$window_index), " windows (minimum 20 cells).")),
  file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", cytokine_sets = 38L,
  methods = levels(curves$method),
  earlier_signal_source = dirname(earlier_path),
  earlier_score_method = earlier_completion$score_method,
  scaling = "independent min-max 0-1 per cytokine and method; visualization only",
  earlier_windows = uniqueN(earlier$window_index),
  current_windows = uniqueN(current$window_index)),
  file.path(stage_dir, "completion.json"), pretty = TRUE,
  auto_unbox = TRUE, digits = NA)
stopifnot(file.rename(stage_dir, output_dir))
message("Completed earlier-signal cytokine comparison at: ", output_dir)
