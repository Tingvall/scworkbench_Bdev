#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 39_mapped_density_enriched_vs_bdev_residual_regions.R",
  "<mapped_density.tsv> <shared_regions.tsv>",
  "<residual_cytokine_regions.tsv> <output_dir>"))
density_path <- normalizePath(args[[1L]], mustWork = TRUE)
shared_path <- normalizePath(args[[2L]], mustWork = TRUE)
residual_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(stage_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

density <- fread(density_path)[target %chin% c("Flt3l", "Il7")]
stopifnot(nrow(density) > 0L,
  all(c("target", "genotype", "bin_midpoint",
    "bdev_normalized_rolling_fraction") %in% names(density)),
  setequal(unique(density$target), c("Flt3l", "Il7")),
  setequal(unique(density$genotype), c("WT", "KO")))
density[, panel := factor(target, levels = c("Flt3l", "Il7"),
  labels = c("FLT3L", "IL-7"))]
density[, curve := fifelse(genotype == "WT", "WT",
  fifelse(target == "Flt3l", "FLT3LKO", "IL7KO"))]

shared <- fread(shared_path)
stopifnot(nrow(shared) == 2L,
  setequal(shared$target, c("Flt3l", "Il7")))
shared_regions <- shared[, .(panel = factor(target,
    levels = c("Flt3l", "Il7"), labels = c("FLT3L", "IL-7")),
  region_id, region_type = "Mapped WT+KO enriched",
  xmin = pt_left, xmax = pt_right)]

residual <- fread(residual_path)
flt_key <- "parse_cytokine_rest_only_refinement_all_top100_v1.flt3l.up"
il7_key <- "parse_cytokine_rest_only_refinement_all_top100_v1.il_7.up"
residual <- residual[feature_set %chin% c(flt_key, il7_key)]
stopifnot(nrow(residual) == 3L,
  sum(residual$feature_set == flt_key) == 2L,
  sum(residual$feature_set == il7_key) == 1L,
  all(residual$region_supported == TRUE))
residual_regions <- residual[, .(
  panel = factor(fifelse(feature_set == flt_key, "Flt3l", "Il7"),
    levels = c("Flt3l", "Il7"), labels = c("FLT3L", "IL-7")),
  region_id, region_type = "Bdev residual-UCell activity",
  xmin = region_left, xmax = region_right)]
regions <- rbind(shared_regions, residual_regions)
regions[, region_type := factor(region_type,
  levels = c("Mapped WT+KO enriched", "Bdev residual-UCell activity"))]

overlap <- merge(shared_regions[, .(panel, enriched_id = region_id,
    enriched_left = xmin, enriched_right = xmax)],
  residual_regions[, .(panel, bdev_id = region_id,
    bdev_left = xmin, bdev_right = xmax)], by = "panel", allow.cartesian = TRUE)
overlap[, `:=`(overlap_left = pmax(enriched_left, bdev_left),
  overlap_right = pmin(enriched_right, bdev_right))]
overlap[, overlap_width := pmax(0, overlap_right - overlap_left)]

labels <- copy(regions)
labels[, label := ifelse(region_type == "Mapped WT+KO enriched",
  "mapped enrichment", "Bdev residual-UCell")]
labels[, x := (xmin + xmax) / 2]

p <- ggplot(density, aes(bin_midpoint,
    bdev_normalized_rolling_fraction, colour = curve, group = curve)) +
  geom_rect(data = regions,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
      fill = region_type), inherit.aes = FALSE, alpha = 0.16) +
  geom_vline(data = regions, aes(xintercept = xmin, colour = region_type),
    inherit.aes = FALSE, linewidth = 0.35, linetype = "dashed",
    show.legend = FALSE) +
  geom_vline(data = regions, aes(xintercept = xmax, colour = region_type),
    inherit.aes = FALSE, linewidth = 0.35, linetype = "dashed",
    show.legend = FALSE) +
  geom_hline(yintercept = 1, colour = "grey25", linewidth = 0.45,
    linetype = "dotted") +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_text(data = labels,
    aes(x = x, y = Inf, label = label, colour = region_type),
    inherit.aes = FALSE, angle = 90, hjust = 1.05, vjust = 0.4,
    size = 2.8, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~panel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("WT" = "#D55E00",
    "FLT3LKO" = "#0072B2", "IL7KO" = "#009E73",
    "Mapped WT+KO enriched" = "#6A3D9A",
    "Bdev residual-UCell activity" = "#C77C00"),
    breaks = c("WT", "FLT3LKO", "IL7KO"), name = NULL) +
  scale_fill_manual(values = c("Mapped WT+KO enriched" = "#6A3D9A",
    "Bdev residual-UCell activity" = "#F2B134"), name = "Highlighted interval") +
  labs(
    title = "Mapped WT/KO density relative to Bdev cytokine-activity regions",
    subtitle = paste("Purple: both mapped WT and KO enriched over Bdev;",
      "gold: independently retrieved Bdev residual-UCell positive region"),
    x = "Mapped Bdev pseudotime",
    y = "Bdev-normalized rolling cell density",
    caption = "Dotted horizontal line at 1 denotes the Bdev reference density.") +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom", legend.box = "vertical",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(8, 8, 8, 8))

png_path <- file.path(stage_dir,
  "mapped_density_enriched_vs_bdev_residual_regions.png")
pdf_path <- file.path(stage_dir,
  "mapped_density_enriched_vs_bdev_residual_regions.pdf")
ggsave(png_path, p, width = 12, height = 5.8, dpi = 300, bg = "white")
ggsave(pdf_path, p, width = 12, height = 5.8, device = cairo_pdf)
fwrite(regions, file.path(stage_dir, "highlighted_regions.tsv"), sep = "\t")
fwrite(overlap, file.path(stage_dir, "region_overlap.tsv"), sep = "\t")

writeLines(c(
  "# Mapped density versus Bdev residual-UCell cytokine regions", "",
  "The mapped-density curves are centered three-bin rolling cell fractions divided by the corresponding Bdev reference fraction; 1 denotes equal density.",
  "Purple intervals require both mapped WT and KO density ratios to exceed 1. Gold intervals are independently defined positive regions from the completed Bdev residual-UCell cytokine permutation analysis.",
  "FLT3L has two supported Bdev residual-UCell regions; IL-7 has one. Exact interval bounds and pairwise overlaps are exported."),
  file.path(stage_dir, "REPORT.md"))
write_json(list(status = "complete", promoted = FALSE,
  density_source = density_path, enriched_region_source = shared_path,
  bdev_residual_region_source = residual_path,
  mapped_enriched_regions = nrow(shared_regions),
  bdev_residual_regions = nrow(residual_regions),
  output_png = basename(png_path), output_pdf = basename(pdf_path)),
  file.path(stage_dir, "completion.json"), pretty = TRUE, auto_unbox = TRUE)

stopifnot(file.rename(stage_dir, output_dir))
message("Completed mapped-density region comparison at: ", output_dir)
