#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: 28_epigenetic_raw_residual_ucell_z_curves.R",
  "<residual_diagnostic_dir> <z_report> <selection.tsv> <bb.tsv>",
  "<output_dir>"))
diagnostic_dir <- normalizePath(args[[1L]], mustWork = TRUE)
z_root <- normalizePath(args[[2L]], mustWork = TRUE)
selection_path <- normalizePath(args[[3L]], mustWork = TRUE)
bb_path <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = TRUE)

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

selection <- fread(selection_path)
bb <- fread(bb_path)
ucell <- read_tsv_gz(file.path(diagnostic_dir, "window_curves.tsv.gz"))
z <- rbindlist(lapply(c("bp", "mf", "cc"), function(g)
  read_tsv_gz(file.path(z_root, g, "window_curves.tsv.gz"))[
    feature_set %chin% selection$feature_set]))
stopifnot(nrow(selection) == 9L,
  setequal(unique(ucell$feature_set), selection$feature_set),
  setequal(unique(z$feature_set), selection$feature_set))

curves <- rbindlist(list(
  ucell[method == "raw_ucell", .(analysis_group, feature_set,
    window_index, window_center, method = "Raw UCell",
    original_value = window_median)],
  ucell[method == "complexity_residual", .(analysis_group, feature_set,
    window_index, window_center, method = "RNA-feature residualized UCell",
    original_value = window_median)],
  z[, .(analysis_group, feature_set, window_index, window_center,
    method = "Original expression z-score", original_value = z_score)]))
stopifnot(curves[, .N, by = .(feature_set, method)][, uniqueN(N)] == 1L)

# Independent 0--1 scaling is intentionally restricted to visualization: the
# three methods have incompatible raw units, and the requested comparison is
# of trajectory shape and peak timing.
curves[, scaled_0_1 := {
  limits <- range(original_value)
  if (diff(limits) > 0) (original_value - limits[[1L]]) / diff(limits) else 0
}, by = .(feature_set, method)]
curves[, set_label := selection$set_name[match(feature_set,
  selection$feature_set)]]
curves[, set_label := factor(set_label,
  levels = selection$set_name)]
curves[, method := factor(method, levels = c("Raw UCell",
  "RNA-feature residualized UCell", "Original expression z-score"))]

p <- ggplot(curves, aes(window_center, scaled_0_1, color = method,
    linetype = method)) +
  geom_vline(data = bb, aes(xintercept = pt_median), inherit.aes = FALSE,
    linetype = "dotted", linewidth = .35, color = "grey35") +
  geom_line(linewidth = .62) +
  facet_wrap(~set_label, ncol = 3) +
  scale_color_manual(values = c("Raw UCell" = "#252525",
    "RNA-feature residualized UCell" = "#D7301F",
    "Original expression z-score" = "#2166AC")) +
  scale_linetype_manual(values = c("Raw UCell" = "solid",
    "RNA-feature residualized UCell" = "solid",
    "Original expression z-score" = "longdash")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1),
    expand = expansion(mult = c(.01, .03))) +
  labs(title = "Raw UCell, RNA-complexity residuals, and expression z-score",
    subtitle = paste("Each gene-set/method curve independently min–max scaled",
      "to 0–1; dotted lines mark BB1 and BB2"),
    x = "Existing Bdev pseudotime (pt)",
    y = "Within-curve scaled value (0–1)", color = NULL,
    linetype = NULL,
    caption = "Scaling is visualization-only; permutation tests used unscaled values.") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
    strip.text = element_text(size = 7),
    plot.title = element_text(face = "bold"))

png_path <- file.path(output_dir,
  "raw_residual_ucell_original_z_scaled_0_1_curves.png")
pdf_path <- file.path(output_dir,
  "raw_residual_ucell_original_z_scaled_0_1_curves.pdf")
table_path <- file.path(output_dir,
  "raw_residual_ucell_original_z_scaled_0_1_curves.tsv.gz")
for (path in c(png_path, pdf_path, table_path))
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
ggsave(png_path, p, width = 12, height = 9, dpi = 240)
ggsave(pdf_path, p, width = 12, height = 9)
fwrite(curves, table_path, sep = "\t", quote = FALSE, compress = "gzip")

write_json(list(status = "complete", selected_sets = nrow(selection),
  methods = levels(curves$method), scaling = paste("independent min-max 0-1",
    "within each gene-set/method curve; visualization only"),
  png = basename(png_path), pdf = basename(pdf_path),
  plotted_values = basename(table_path)),
  file.path(output_dir,
    "raw_residual_ucell_original_z_scaled_0_1_completion.json"),
  pretty = TRUE, auto_unbox = TRUE)
message("Completed three-method trajectory plot at: ", png_path)
