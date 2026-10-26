#!/usr/bin/env Rscript

# Diagnose why IL7 query cells form an apparent enrichment peak at the start
# of transferred Bdev pseudotime. This is report-only and leaves all bundles
# and prior reports unchanged.

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) stop(paste(
  "Usage: 80_diagnose_il7_early_mapping_peak.R BDEV_BUNDLE PROJECTION_TSV_GZ",
  "DENSITY_TSV REGION_TSV ALL_NEIGHBORS_TSV_GZ PT_NEIGHBORS_TSV_GZ",
  "LINEAGE_SCORES_TSV_GZ MAPPING_METRICS_JSON OUTPUT_DIR SCRIPT_VERSION"
), call. = FALSE)

bdev_bundle <- normalizePath(args[[1L]], mustWork = TRUE)
projection_path <- normalizePath(args[[2L]], mustWork = TRUE)
density_path <- normalizePath(args[[3L]], mustWork = TRUE)
region_path <- normalizePath(args[[4L]], mustWork = TRUE)
all_neighbor_path <- normalizePath(args[[5L]], mustWork = TRUE)
pt_neighbor_path <- normalizePath(args[[6L]], mustWork = TRUE)
lineage_path <- normalizePath(args[[7L]], mustWork = TRUE)
metrics_path <- normalizePath(args[[8L]], mustWork = TRUE)
output_dir <- normalizePath(args[[9L]], mustWork = FALSE)
script_version <- args[[10L]]

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output already exists: ", stage_dir, call. = FALSE)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage_dir, "figures"), showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
})

write_json <- function(x, path) jsonlite::write_json(
  x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA,
  null = "null", na = "null")
write_tsv <- function(x, path) data.table::fwrite(x, path, sep = "\t", na = "NA")
qtext <- function(x) paste(format(round(quantile(x, c(.25, .5, .75),
  na.rm = TRUE), 3), trim = TRUE), collapse = ";")

regions <- fread(region_path)[target == "Il7"]
stopifnot(nrow(regions) == 2L,
  identical(as.character(regions$region_id), c("Il7_peak1", "Il7_peak2")))
setorder(regions, peak_number)

projection <- fread(cmd = paste("zcat", shQuote(projection_path)))
required_projection <- c("cell_id", "source_cell_id", "sample_id", "genotype",
  "knockout_target", "support_status", "mapped_pt", "predicted_cluster_k30",
  "median30_distance", "projected_x", "projected_y")
stopifnot(all(required_projection %in% names(projection)))
query <- projection[knockout_target == "Il7" &
  support_status %in% c("supported", "borderline") & is.finite(mapped_pt)]
query[, region_id := "Other mapped IL-7 cells"]
for (i in seq_len(nrow(regions))) {
  query[mapped_pt >= regions$pt_left[[i]] & mapped_pt < regions$pt_right[[i]],
    region_id := regions$region_id[[i]]]
}
query[, region_label := factor(fcase(
  region_id == "Il7_peak1", "Peak 1 (early)",
  region_id == "Il7_peak2", "Peak 2",
  default = "Other"), levels = c("Other", "Peak 1 (early)", "Peak 2"))]
stopifnot(query[region_id == "Il7_peak1", .N] == 6025L,
  query[region_id == "Il7_peak2", .N] == 3381L)

# The original support call used all Bdev cells. Recalculate the corresponding
# distance summaries for the separate neighbor pool restricted to cells with pt.
all_neighbors <- fread(cmd = paste("zcat", shQuote(all_neighbor_path)))[rank <= 30L]
pt_neighbors <- fread(cmd = paste("zcat", shQuote(pt_neighbor_path)))[rank <= 30L]
all_distance <- all_neighbors[, .(
  all_nearest_distance = min(distance),
  all_median30_distance = median(distance),
  all_k30_distance = max(distance)), by = query_cell_id]
pt_distance <- pt_neighbors[, .(
  pt_nearest_distance = min(distance),
  pt_median30_distance = median(distance),
  pt_k30_distance = max(distance)), by = query_cell_id]
distance <- merge(all_distance, pt_distance, by = "query_cell_id")
distance[, pt_to_all_median30_ratio := pt_median30_distance / all_median30_distance]
query <- merge(query, distance, by.x = "cell_id", by.y = "query_cell_id",
  all.x = TRUE, sort = FALSE)
stopifnot(nrow(query) == 11256L, !anyNA(query$pt_median30_distance))

# Recover immutable reference metadata and native UMAP. The reference pt field
# is only read and is checked before/after opening the bundle.
bdev <- wb_open_bundle(bdev_bundle, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(bdev), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bdev)$current_state$fingerprint
bdev_cells <- as.character(wb_cells(bdev, "all"))
bdev_md <- as.data.table(wb_meta(bdev, cells = bdev_cells))
bdev_umap <- as.data.table(wb_space(bdev, "umap.atac", cells = bdev_cells))
stopifnot(identical(as.character(bdev_md$cell_id), bdev_cells),
  identical(as.character(bdev_umap$cell_id), bdev_cells))
reference_pt_copy <- bdev_md$pt
reference <- data.table(
  reference_cell_id = paste0("bdev::", bdev_cells),
  source_cell_id = bdev_cells,
  x = bdev_umap$x, y = bdev_umap$y,
  reference_cluster = as.character(bdev_md[["peaks_snn_res.1.2"]]),
  reference_pt = as.numeric(bdev_md$pt),
  reference_phase = as.character(bdev_md$Phase),
  BT_peak_group1 = as.character(bdev_md$BT_peak_group1),
  B_gene1 = as.numeric(bdev_md$B_gene1),
  T_gene1 = as.numeric(bdev_md$T_gene1),
  B_specific_peakscore1 = as.numeric(bdev_md$B_specific_peakscore1),
  T_specific_peakscore1 = as.numeric(bdev_md$T_specific_peakscore1))
stopifnot(identical(reference_pt_copy, bdev_md$pt))

# Attach reference identity to both neighbor searches to reveal whether the
# cell-support pool and pt-transfer pool represent the same biological state.
query_key <- query[, .(cell_id, genotype, region_id, region_label)]
all_neighbor_annotated <- merge(all_neighbors, query_key,
  by.x = "query_cell_id", by.y = "cell_id")
all_neighbor_annotated <- merge(all_neighbor_annotated,
  reference[, .(reference_cell_id, reference_cluster, reference_pt)],
  by = "reference_cell_id")
all_neighbor_annotated[, neighbor_pool := "All reference cells"]
pt_neighbor_annotated <- merge(pt_neighbors, query_key,
  by.x = "query_cell_id", by.y = "cell_id")
pt_neighbor_annotated <- merge(pt_neighbor_annotated,
  reference[, .(reference_cell_id, reference_cluster, reference_pt)],
  by = "reference_cell_id")
pt_neighbor_annotated[, neighbor_pool := "Reference cells with defined pt"]
stopifnot(all(is.finite(pt_neighbor_annotated$reference_pt)))

cluster18_reference <- reference[reference_cluster == "18"]
stopifnot(nrow(cluster18_reference) == 447L,
  sum(is.finite(cluster18_reference$reference_pt)) == 0L)

region_counts <- query[, .(n_cells = .N,
  supported_fraction = mean(support_status == "supported")),
  by = .(region_id, genotype)]
region_counts[, fraction_of_mapped_il7 := n_cells / sum(n_cells), by = genotype]
setorder(region_counts, region_id, genotype)

cluster_composition <- query[region_id %in% regions$region_id,
  .(n_cells = .N), by = .(region_id, genotype,
    predicted_cluster = as.character(predicted_cluster_k30))]
cluster_composition[, fraction := n_cells / sum(n_cells),
  by = .(region_id, genotype)]
setorder(cluster_composition, region_id, genotype, -n_cells)

distance_summary <- query[region_id %in% regions$region_id, .(
  n_cells = .N,
  all_median30_q25_median_q75 = qtext(all_median30_distance),
  pt_median30_q25_median_q75 = qtext(pt_median30_distance),
  ratio_q25_median_q75 = qtext(pt_to_all_median30_ratio),
  pt_distance_above_all_reference_q99_fraction =
    mean(pt_median30_distance > read_json(metrics_path,
      simplifyVector = TRUE)$distance_calibration$q99),
  cluster18_fraction = mean(as.character(predicted_cluster_k30) == "18")),
  by = .(region_id, genotype)]
setorder(distance_summary, region_id, genotype)

lineage <- fread(cmd = paste("zcat", shQuote(lineage_path)))[
  knockout_target == "Il7" & lineage %in% c("B lineage", "T lineage"),
  .(source_cell_id = cell_id, lineage, aggregate_score)]
lineage <- merge(lineage,
  query[region_id %in% regions$region_id,
    .(source_cell_id, region_id, genotype)], by = "source_cell_id")
lineage_wide <- dcast(lineage,
  source_cell_id + region_id + genotype ~ lineage,
  value.var = "aggregate_score")
lineage_wide[, T_minus_B := `T lineage` - `B lineage`]
lineage_summary <- lineage_wide[, .(
  n_cells = .N, mean_T_minus_B = mean(T_minus_B),
  median_T_minus_B = median(T_minus_B),
  q25_T_minus_B = quantile(T_minus_B, .25),
  q75_T_minus_B = quantile(T_minus_B, .75),
  fraction_T_score_above_B_score = mean(T_minus_B > 0)),
  by = .(region_id, genotype)]
setorder(lineage_summary, region_id, genotype)

reference_cluster_summary <- reference[, .(
  n_reference_cells = .N,
  n_defined_pt = sum(is.finite(reference_pt)),
  fraction_defined_pt = mean(is.finite(reference_pt)),
  fraction_T_peak_group = mean(BT_peak_group1 == "T", na.rm = TRUE),
  mean_T_minus_B_gene_score = mean(T_gene1 - B_gene1, na.rm = TRUE),
  mean_T_minus_B_peak_score = mean(
    T_specific_peakscore1 - B_specific_peakscore1, na.rm = TRUE)),
  by = reference_cluster]
setorder(reference_cluster_summary, reference_cluster)

neighbor_composition <- rbindlist(list(
  all_neighbor_annotated[region_id %in% regions$region_id,
    .(neighbor_records = .N, median_distance = median(distance),
      median_reference_pt = median(reference_pt, na.rm = TRUE),
      fraction_reference_pt_defined = mean(is.finite(reference_pt))),
    by = .(neighbor_pool, region_id, genotype, reference_cluster)],
  pt_neighbor_annotated[region_id %in% regions$region_id,
    .(neighbor_records = .N, median_distance = median(distance),
      median_reference_pt = median(reference_pt),
      fraction_reference_pt_defined = 1),
    by = .(neighbor_pool, region_id, genotype, reference_cluster)]),
  use.names = TRUE)
setorder(neighbor_composition, neighbor_pool, region_id, genotype,
  -neighbor_records)

write_tsv(region_counts, file.path(stage_dir, "region_cell_counts.tsv"))
write_tsv(cluster_composition,
  file.path(stage_dir, "query_cluster_composition.tsv"))
write_tsv(distance_summary,
  file.path(stage_dir, "mapping_distance_summary.tsv"))
write_tsv(lineage_summary,
  file.path(stage_dir, "lineage_score_summary.tsv"))
write_tsv(reference_cluster_summary,
  file.path(stage_dir, "reference_cluster_pt_and_lineage_summary.tsv"))
write_tsv(neighbor_composition,
  file.path(stage_dir, "neighbor_pool_cluster_composition.tsv"))

base_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    strip.text = element_text(face = "bold"))
peak_colours <- c("Other" = "grey75", "Peak 1 (early)" = "#009E73",
  "Peak 2" = "#CC79A7")
genotype_colours <- c("WT" = "#C84545", "KO" = "#3B77A8")

density <- fread(density_path)[target == "Il7"]
density_long <- rbindlist(list(
  density[, .(genotype, bin_midpoint,
    metric = "Raw query rolling fraction (%)", value = 100 * rolling_fraction)],
  density[, .(genotype, bin_midpoint,
    metric = "Query / Bdev rolling fraction", value = bdev_normalized_rolling_fraction)]))
density_long[, metric := factor(metric,
  levels = c("Raw query rolling fraction (%)", "Query / Bdev rolling fraction"))]
region_rect <- rbindlist(lapply(seq_len(nrow(regions)), function(i) data.table(
  metric = levels(density_long$metric),
  xmin = regions$pt_left[[i]], xmax = regions$pt_right[[i]],
  region_label = factor(if (i == 1L) "Peak 1 (early)" else "Peak 2",
    levels = levels(query$region_label)))))
p_density <- ggplot(density_long, aes(bin_midpoint, value, colour = genotype)) +
  geom_rect(data = region_rect,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
      fill = region_label), inherit.aes = FALSE, alpha = .12) +
  geom_hline(data = data.table(metric = factor(
    "Query / Bdev rolling fraction", levels = levels(density_long$metric)), y = 1),
    aes(yintercept = y), linetype = "dotted", colour = "grey30") +
  geom_line(linewidth = .8, na.rm = TRUE) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = genotype_colours) +
  scale_fill_manual(values = peak_colours, drop = FALSE) +
  labs(x = "Mapped Bdev pseudotime", y = NULL, colour = "Genotype",
    fill = "Mapped interval", title = "A. The early peak exists before Bdev normalization",
    subtitle = "Peak intervals were originally defined from WT query/Bdev ratio > 1") +
  base_theme

reference_plot <- reference[, .(x, y)]
query_plot <- query[is.finite(projected_x) & is.finite(projected_y)]
query_plot[, draw_order := region_label != "Other"]
setorder(query_plot, draw_order)
p_umap <- ggplot() +
  geom_point(data = reference_plot, aes(x, y), colour = "grey90", size = .09) +
  geom_point(data = query_plot, aes(projected_x, projected_y,
    colour = region_label), size = .22, alpha = .62) +
  facet_wrap(~genotype) + coord_equal() +
  scale_colour_manual(values = peak_colours, drop = FALSE) +
  labs(x = "Bdev umap.atac 1", y = "Bdev umap.atac 2",
    colour = "Mapped interval",
    title = "B. Peak 1 is the lower off-trajectory island",
    subtitle = "Query coordinates are Bdev-neighbor barycentres; grey points are the native Bdev reference") +
  base_theme

cluster_plot <- copy(cluster_composition)
cluster_plot[, region_label := factor(fcase(
  region_id == "Il7_peak1", "Peak 1 (early)",
  default = "Peak 2"), levels = c("Peak 1 (early)", "Peak 2"))]
cluster_plot[, cluster_display := predicted_cluster]
cluster_plot[!predicted_cluster %in% c("18", "6", "17", "23", "21", "8", "7", "3"),
  cluster_display := "Other"]
cluster_plot <- cluster_plot[, .(n_cells = sum(n_cells)),
  by = .(region_label, genotype, cluster_display)]
cluster_plot[, fraction := n_cells / sum(n_cells), by = .(region_label, genotype)]
p_cluster <- ggplot(cluster_plot,
    aes(interaction(region_label, genotype, sep = "\n"), fraction,
      fill = cluster_display)) +
  geom_col(width = .75) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_viridis_d(option = "turbo", name = "Transferred\npeak cluster") +
  labs(x = NULL, y = "Cell fraction",
    title = "C. Peak 1 is almost entirely cluster 18") + base_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

lineage_wide[, region_label := factor(fcase(
  region_id == "Il7_peak1", "Peak 1 (early)", default = "Peak 2"),
  levels = c("Peak 1 (early)", "Peak 2"))]
p_lineage <- ggplot(lineage_wide,
    aes(interaction(region_label, genotype, sep = "\n"), T_minus_B,
      fill = genotype)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey35") +
  geom_violin(scale = "width", trim = TRUE, linewidth = .25) +
  geom_boxplot(width = .13, outlier.shape = NA, fill = "white", linewidth = .3) +
  scale_fill_manual(values = genotype_colours) +
  labs(x = NULL, y = "T-lineage score minus B-lineage score",
    title = "D. Peak 1 is T-lineage-like") + base_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")

distance_long <- melt(query[region_id %in% regions$region_id,
  .(cell_id, region_id, genotype, all_median30_distance, pt_median30_distance)],
  id.vars = c("cell_id", "region_id", "genotype"),
  variable.name = "neighbor_pool", value.name = "median30_distance")
distance_long[, region_label := factor(fcase(
  region_id == "Il7_peak1", "Peak 1 (early)", default = "Peak 2"),
  levels = c("Peak 1 (early)", "Peak 2"))]
distance_long[, neighbor_pool := factor(neighbor_pool,
  levels = c("all_median30_distance", "pt_median30_distance"),
  labels = c("All Bdev cells", "Only Bdev cells with pt"))]
p_distance <- ggplot(distance_long,
    aes(neighbor_pool, median30_distance, fill = neighbor_pool)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = .25) +
  geom_boxplot(width = .12, outlier.shape = NA, fill = "white", linewidth = .3) +
  geom_hline(yintercept = read_json(metrics_path,
    simplifyVector = TRUE)$distance_calibration$q99,
    linetype = "dotted", colour = "#D55E00", linewidth = .6) +
  facet_grid(genotype ~ region_label) +
  scale_fill_manual(values = c("All Bdev cells" = "#56B4E9",
    "Only Bdev cells with pt" = "#E69F00")) +
  labs(x = NULL, y = "Median distance to 30 neighbors",
    fill = "Neighbor pool",
    title = "E. Peak-1 cells are close to off-trajectory reference cells but far from cells with pseudotime",
    subtitle = "Dotted line: q99 support threshold calibrated for the all-reference search") +
  base_theme

combined <- p_density / p_umap / (p_cluster | p_lineage) / p_distance +
  plot_annotation(
    title = "Diagnostic of the apparent early IL-7 mapped-pseudotime peak",
    subtitle = paste(
      "Cluster 18 has no reference pt; its T-like query cells are assigned to the nearest",
      "available defined-pt cells at the trajectory start")) +
  plot_layout(heights = c(1.25, 1.35, 1, 1.15))
ggsave(file.path(stage_dir, "figures", "il7_early_peak_diagnostic.pdf"),
  combined, width = 12, height = 17, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "il7_early_peak_diagnostic.png"),
  combined, width = 12, height = 17, dpi = 220, bg = "white")

p1_wt <- region_counts[region_id == "Il7_peak1" & genotype == "WT"]
p1_ko <- region_counts[region_id == "Il7_peak1" & genotype == "KO"]
p2_wt <- region_counts[region_id == "Il7_peak2" & genotype == "WT"]
p2_ko <- region_counts[region_id == "Il7_peak2" & genotype == "KO"]
p1_cluster_wt <- distance_summary[region_id == "Il7_peak1" & genotype == "WT"]
p1_cluster_ko <- distance_summary[region_id == "Il7_peak1" & genotype == "KO"]
p1_lin_wt <- lineage_summary[region_id == "Il7_peak1" & genotype == "WT"]
p1_lin_ko <- lineage_summary[region_id == "Il7_peak1" & genotype == "KO"]

report_lines <- c(
  "# Diagnostic of the apparent early IL-7 mapping peak", "",
  "## Conclusion", "",
  paste(
    "The strong first peak is primarily an **off-trajectory T-lineage population",
    "that is forced onto the beginning of Bdev pseudotime**, not a second",
    "well-supported B-cell developmental accumulation."), "",
  "The mechanism is structural:", "",
  "1. Mapping support and transferred cluster identity were calculated against all Bdev reference cells.",
  "2. Reference peak cluster 18 contains 447 cells, but none has a defined `pt` value.",
  "3. Pseudotime was transferred in a separate neighbor search that excluded all reference cells without `pt`.",
  "4. Cluster-18-like query cells therefore retain strong global mapping support but receive pseudotime from much more distant, defined-pt cells at the earliest available trajectory states (mainly clusters 8, 3, and 7).", "",
  "## Evidence", "",
  sprintf("- Peak 1 contains %s WT cells (%.1f%% of mapped IL-7 WT) and %s KO cells (%.1f%% of mapped IL-7 KO), so it is present in raw counts and is not created by Bdev normalization.",
    format(p1_wt$n_cells, big.mark = ","), 100 * p1_wt$fraction_of_mapped_il7,
    format(p1_ko$n_cells, big.mark = ","), 100 * p1_ko$fraction_of_mapped_il7),
  sprintf("- %.1f%% of peak-1 WT and %.1f%% of peak-1 KO cells transfer to cluster 18.",
    100 * p1_cluster_wt$cluster18_fraction,
    100 * p1_cluster_ko$cluster18_fraction),
  sprintf("- T-lineage score exceeds B-lineage score in %.1f%% of peak-1 WT and %.1f%% of peak-1 KO cells. For peak 2 the corresponding fractions are %.1f%% and %.1f%%.",
    100 * p1_lin_wt$fraction_T_score_above_B_score,
    100 * p1_lin_ko$fraction_T_score_above_B_score,
    100 * lineage_summary[region_id == "Il7_peak2" & genotype == "WT",
      fraction_T_score_above_B_score],
    100 * lineage_summary[region_id == "Il7_peak2" & genotype == "KO",
      fraction_T_score_above_B_score]),
  sprintf("- Peak-1 median distance to 30 neighbors rises from %s to %s in WT and from %s to %s in KO when the pool is restricted to reference cells with defined pt.",
    strsplit(p1_cluster_wt$all_median30_q25_median_q75, ";", fixed = TRUE)[[1]][2],
    strsplit(p1_cluster_wt$pt_median30_q25_median_q75, ";", fixed = TRUE)[[1]][2],
    strsplit(p1_cluster_ko$all_median30_q25_median_q75, ";", fixed = TRUE)[[1]][2],
    strsplit(p1_cluster_ko$pt_median30_q25_median_q75, ";", fixed = TRUE)[[1]][2]),
  sprintf("- Peak 2 contains %s WT and %s KO cells, is B-lineage-skewed, and its defined-pt neighbor distance is substantially smaller; it is therefore much more compatible with a true mapped Bdev region, although part of peak 2 also maps to clusters lacking pt and should retain a caveat.",
    format(p2_wt$n_cells, big.mark = ","), format(p2_ko$n_cells, big.mark = ",")), "",
  "## Interpretation", "",
  paste(
    "For Bdev trajectory analyses, peak 1 should not be treated as an IL-7 B-cell",
    "developmental peak. The appropriate correction is to require trajectory support",
    "in the pseudotime-defined reference pool (or exclude transferred off-trajectory/T-like",
    "clusters) before defining mapped enrichment regions. Peak 2 remains the biologically",
    "interpretable IL-7 B-lineage region under this diagnostic."), "",
  "The analysis is descriptive because there is one sample per genotype. It diagnoses mapping geometry and composition; it does not test biological replication.", "",
  "## Outputs", "",
  "- `figures/il7_early_peak_diagnostic.pdf` and `.png`: combined diagnostic figure.",
  "- `region_cell_counts.tsv`: peak abundance in each genotype.",
  "- `query_cluster_composition.tsv`: transferred cluster composition.",
  "- `mapping_distance_summary.tsv`: all-reference versus defined-pt neighbor distances.",
  "- `lineage_score_summary.tsv`: B/T lineage score contrast.",
  "- `reference_cluster_pt_and_lineage_summary.tsv`: reference cluster support for pt and lineage.",
  "- `neighbor_pool_cluster_composition.tsv`: identities of neighbors used by each search.", "",
  sprintf("Script version: `%s`.", script_version))
writeLines(report_lines, file.path(stage_dir, "REPORT.md"))

input_manifest <- data.table(
  role = c("bdev_bundle", "query_projection", "mapped_density", "peak_regions",
    "all_reference_neighbors", "defined_pt_neighbors", "lineage_scores",
    "mapping_metrics"),
  path = c(bdev_bundle, projection_path, density_path, region_path,
    all_neighbor_path, pt_neighbor_path, lineage_path, metrics_path))
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))
write_json(list(
  status = "complete_unpromoted_report_only",
  script_version = script_version,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  target = "Il7",
  query_cells = nrow(query),
  peak1_cells = query[region_id == "Il7_peak1", .N],
  peak2_cells = query[region_id == "Il7_peak2", .N],
  reference_cluster18_cells = nrow(cluster18_reference),
  reference_cluster18_defined_pt = sum(is.finite(cluster18_reference$reference_pt)),
  source_bundle_fingerprint_before = fingerprint_before,
  source_bundle_fingerprint_after = wb_manifest(bdev)$current_state$fingerprint,
  reference_pt_preserved = identical(reference_pt_copy, bdev_md$pt),
  promoted = FALSE), file.path(stage_dir, "completion.json"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

stopifnot(identical(fingerprint_before,
  wb_manifest(bdev)$current_state$fingerprint))
if (!file.rename(stage_dir, output_dir))
  stop("Could not atomically finalize output directory", call. = FALSE)
cat("Diagnostic report:", file.path(output_dir, "REPORT.md"), "\n")
