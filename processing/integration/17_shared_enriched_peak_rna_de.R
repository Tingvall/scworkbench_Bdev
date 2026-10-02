#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "Usage: 17_shared_enriched_peak_rna_de.R KO_BUNDLE DISTRIBUTION_REPORT",
  "PROJECTION_REPORT ORIGINAL_REGION_REPORT SESSION_WORKSPACE OUTPUT_DIR",
  "EXPECTED_KO_FP SCRIPT_VERSION"
), call. = FALSE)

ko_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
distribution_report <- normalizePath(args[[2L]], mustWork = TRUE)
projection_report <- normalizePath(args[[3L]], mustWork = TRUE)
original_region_report <- normalizePath(args[[4L]], mustWork = TRUE)
session_workspace <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
expected_ko_fp <- args[[7L]]
script_version <- args[[8L]]
stopifnot(identical(args[[9L]], "execute"))

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output already exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

write_json <- function(x, path) jsonlite::write_json(
  x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA,
  null = "null", na = "null")
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
log_event <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(..., collapse = " "), "\n")
html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}

distribution_path <- file.path(distribution_report,
  "query_mapped_pt_cell_distribution.tsv")
projection_path <- file.path(projection_report,
  "query_reference_umap_projection.tsv.gz")
original_regions_path <- file.path(original_region_report,
  "wt_enriched_pt_regions.tsv")
stopifnot(file.exists(distribution_path), file.exists(projection_path),
  file.exists(original_regions_path))
distribution <- fread(distribution_path)
projection <- fread(cmd = paste("zcat", shQuote(projection_path)))
original_regions <- fread(original_regions_path)
stopifnot(all(c("target", "genotype", "bin", "bin_left", "bin_right",
  "bin_midpoint", "bdev_normalized_rolling_fraction") %in%
  names(distribution)),
  all(c("source_cell_id", "sample_id", "genotype", "knockout_target",
    "support_status", "mapped_pt") %in% names(projection)))

# A shared-enriched bin must be above the Bdev-normalized ratio-1 line in both
# genotypes. IL-7 candidates are additionally restricted to the original
# WT-defined Il7_peak2 interval; the shared portion of Il7_peak1 is excluded.
wide <- dcast(distribution[target %in% c("Flt3l", "Il7")],
  target + bin + bin_left + bin_right + bin_midpoint ~ genotype,
  value.var = "bdev_normalized_rolling_fraction")
stopifnot(all(c("WT", "KO") %in% names(wide)))
il7_peak2 <- original_regions[region_id == "Il7_peak2"]
stopifnot(nrow(il7_peak2) == 1L)
shared_bins <- wide[is.finite(WT) & is.finite(KO) & WT > 1 & KO > 1 &
  (target == "Flt3l" | (target == "Il7" &
    bin >= il7_peak2$first_bin & bin <= il7_peak2$last_bin))]
setorder(shared_bins, target, bin)
shared_bins[, run := cumsum(seq_len(.N) == 1L |
  target != shift(target, fill = target[[1L]]) |
  bin != shift(bin, fill = bin[[1L]] - 1L) + 1L)]
regions <- shared_bins[, .(
  first_bin = min(bin), last_bin = max(bin),
  pt_left = min(bin_left), pt_right = max(bin_right), n_bins = .N,
  maximum_wt_ratio = max(WT), maximum_ko_ratio = max(KO),
  maximum_minimum_genotype_ratio = max(pmin(WT, KO))),
  by = .(target, run)]
setorder(regions, target, pt_left)
stopifnot(nrow(regions) == 2L,
  identical(regions$target, c("Flt3l", "Il7")),
  identical(regions$first_bin, c(19L, 21L)),
  identical(regions$last_bin, c(26L, 31L)))
regions[, `:=`(
  cytokine = fifelse(target == "Flt3l", "FLT3L", "IL-7"),
  source_peak = fifelse(target == "Flt3l", "Flt3l_peak1", "Il7_peak2"),
  region_id = fifelse(target == "Flt3l", "Flt3l_shared_peak1",
    "Il7_shared_peak2"))]
setcolorder(regions, c("region_id", "target", "cytokine", "source_peak",
  "first_bin", "last_bin", "pt_left", "pt_right", "n_bins",
  "maximum_wt_ratio", "maximum_ko_ratio",
  "maximum_minimum_genotype_ratio", "run"))
write_tsv(regions[, !"run"], file.path(stage_dir,
  "shared_enriched_pt_regions.tsv"))
write_tsv(shared_bins[, !"run"], file.path(stage_dir,
  "shared_enriched_bins.tsv"))

eligible <- projection[support_status %in% c("supported", "borderline") &
  is.finite(mapped_pt)]
assignments <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
  reg <- regions[i]
  eligible[knockout_target == reg$target & mapped_pt >= reg$pt_left &
    mapped_pt < reg$pt_right,
    .(source_cell_id, sample_id, genotype, knockout_target, support_status,
      mapped_pt, region_id = reg$region_id, cytokine = reg$cytokine,
      source_peak = reg$source_peak, pt_left = reg$pt_left,
      pt_right = reg$pt_right)]
}))
stopifnot(!anyDuplicated(assignments$source_cell_id),
  setequal(assignments$genotype, c("WT", "KO")))
setorder(assignments, knockout_target, genotype, mapped_pt)
fwrite(assignments, file.path(stage_dir,
  "shared_region_cell_assignments.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")
region_counts <- assignments[, .(n_cells = .N, n_samples = uniqueN(sample_id),
  samples = paste(sort(unique(sample_id)), collapse = ";")),
  by = .(region_id, knockout_target, cytokine, genotype, pt_left, pt_right)]
setorder(region_counts, knockout_target, genotype)
stopifnot(nrow(region_counts) == 4L, all(region_counts$n_samples == 1L),
  all(region_counts$n_cells >= 100L))
write_tsv(region_counts, file.path(stage_dir, "region_cell_counts.tsv"))

region_plot_data <- copy(distribution[target %in% c("Flt3l", "Il7")])
region_plot_data[, target_label := factor(target,
  levels = c("Flt3l", "Il7"), labels = c("FLT3L", "IL-7"))]
region_plot_data[, colour_group := fifelse(genotype == "WT", "WT",
  fifelse(target == "Flt3l", "FLT3LKO", "IL7KO"))]
region_plot_regions <- copy(regions)
region_plot_regions[, target_label := factor(target,
  levels = c("Flt3l", "Il7"), labels = c("FLT3L", "IL-7"))]
p_regions <- ggplot(region_plot_data,
    aes(bin_midpoint, bdev_normalized_rolling_fraction,
      colour = colour_group, group = colour_group)) +
  geom_rect(data = region_plot_regions,
    aes(xmin = pt_left, xmax = pt_right, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#7B61A8", alpha = 0.18) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey25") +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_text(data = region_plot_regions,
    aes(x = (pt_left + pt_right) / 2, y = Inf,
      label = paste0("bins ", first_bin, "-", last_bin)),
    inherit.aes = FALSE, vjust = 1.3, fontface = "bold", size = 3.2) +
  facet_wrap(~target_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("WT" = "#D55E00",
    "FLT3LKO" = "#0072B2", "IL7KO" = "#009E73")) +
  labs(x = "Mapped Bdev pseudotime",
    y = "Bdev-normalized rolling cell fraction", colour = NULL,
    title = "Shared WT-and-KO enriched pseudotime regions",
    subtitle = paste0("Shading requires both genotype curves > 1; ",
      "IL-7 is restricted to the original peak 2")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank())
ggsave(file.path(stage_dir, "figures", "shared_enriched_regions.pdf"),
  p_regions, width = 11, height = 5.2, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "shared_enriched_regions.png"),
  p_regions, width = 11, height = 5.2, dpi = 300)

bundle <- wb_open_bundle(ko_bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_ko_fp),
  all(assignments$source_cell_id %in%
    wb_cells(bundle, "set:qc_use_downstream")))
session <- wb_session_resume(bundle, session_workspace)
stopifnot(is.null(session$durable$promotion_receipt))
find_run <- function(label, operation) {
  marker <- paste0("[", label, "] ")
  ids <- names(session$runs)[vapply(session$runs, function(x) {
    note <- if (is.null(x$notes)) "" else as.character(x$notes)
    identical(x$op_name, operation) && startsWith(note, marker)
  }, logical(1))]
  stopifnot(length(ids) <= 1L)
  if (length(ids)) session$runs[[ids[[1L]]]] else NULL
}
normal_run <- find_run("normalize_downstream_rna", "matrix_log_normalize")
stopifnot(!is.null(normal_run), "matrix" %in% names(normal_run$outputs))
matrix_id <- normal_run$outputs$matrix$matrix_id
stopifnot(identical(matrix_id, "RNA_log1p_cp10k_wt_enriched_peak_de_v1"))
write_json(wb_op_interface("de_wilcox", bundle = bundle),
  file.path(stage_dir, "de_wilcox_interface.json"))

de_tables <- list()
for (i in seq_len(nrow(regions))) {
  reg <- regions[i]
  label <- if (reg$target == "Flt3l") "de_shared_flt3l" else
    "de_shared_il7_peak2"
  region_cells <- assignments[region_id == reg$region_id]
  wt_cells <- as.character(region_cells[genotype == "WT", source_cell_id])
  ko_cells <- as.character(region_cells[genotype == "KO", source_cell_id])
  stopifnot(length(wt_cells) >= 100L, length(ko_cells) >= 100L)
  run <- find_run(label, "de_wilcox")
  if (is.null(run)) {
    log_event("Running full-gene shared-region DE for", reg$region_id,
      "WT", length(wt_cells), "KO", length(ko_cells))
    result <- wb_session_run_op(session, bundle, "de_wilcox",
      inputs = list(matrix = matrix_id,
        cell_set = "set:qc_use_downstream", min_cells = 10L,
        group_a = wt_cells, group_b = ko_cells, min_pct = 0,
        p_adjust_method = "BH"),
      notes = paste0("[", label, "] Group A=WT; group B=KO; exact ",
        reg$region_id, " bins where both genotype Bdev-normalized rolling ",
        "fractions are >1; all genes; descriptive cell-level test; no ",
        "biological replicates."))
    session <- result$session
    run <- result$run
  } else log_event("Reusing completed shared-region DE for", reg$region_id)
  d <- as.data.table(run$outputs$de_table$table)
  stopifnot(nrow(d) == 33696L, !anyDuplicated(d$feature_id))
  d[, `:=`(region_id = reg$region_id, target = reg$target,
    cytokine = reg$cytokine, source_peak = reg$source_peak,
    first_bin = reg$first_bin, last_bin = reg$last_bin,
    pt_left = reg$pt_left, pt_right = reg$pt_right,
    n_wt_cells = length(wt_cells), n_ko_cells = length(ko_cells))]
  de_tables[[reg$region_id]] <- d
}

staged <- wb_session_view(session)
on.exit(try(wb_close_bundle(staged), silent = TRUE), add = TRUE)
all_features <- as.character(de_tables[[1L]]$feature_id)
cell_order <- as.character(assignments$source_cell_id)
blocks <- split(all_features, ceiling(seq_along(all_features) / 128L))
mean_tables <- rbindlist(lapply(blocks, function(feats) {
  m <- wb_matrix_read_features(staged, matrix_id, feats, cell_order,
    missing = "NA")
  stopifnot(!anyNA(m), identical(rownames(m), feats),
    identical(colnames(m), cell_order))
  linear <- expm1(as.matrix(m))
  out <- list(feature_id = feats)
  for (i in seq_len(nrow(regions))) {
    reg <- regions[i]
    prefix <- if (reg$target == "Flt3l") "flt3l" else "il7"
    idx_wt <- which(assignments$region_id == reg$region_id &
      assignments$genotype == "WT")
    idx_ko <- which(assignments$region_id == reg$region_id &
      assignments$genotype == "KO")
    out[[paste0(prefix, "_mean_cp10k_wt")]] <- rowMeans(
      linear[, idx_wt, drop = FALSE])
    out[[paste0(prefix, "_mean_cp10k_ko")]] <- rowMeans(
      linear[, idx_ko, drop = FALSE])
  }
  as.data.table(out)
}))
stopifnot(nrow(mean_tables) == 33696L,
  !anyDuplicated(mean_tables$feature_id))

combined <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
  reg <- regions[i]
  d <- merge(copy(de_tables[[reg$region_id]]), mean_tables,
    by = "feature_id", all.x = TRUE)
  prefix <- if (reg$target == "Flt3l") "flt3l" else "il7"
  d[, `:=`(mean_cp10k_wt = get(paste0(prefix, "_mean_cp10k_wt")),
    mean_cp10k_ko = get(paste0(prefix, "_mean_cp10k_ko")))]
  d[, standard_log2fc := log2((mean_cp10k_wt + 1) /
    (mean_cp10k_ko + 1))]
  d[, `:=`(direction = fifelse(standard_log2fc > 0, "WT_higher",
    fifelse(standard_log2fc < 0, "KO_higher", "equal")),
    significant_bh05 = is.finite(p_val_adj) & p_val_adj < 0.05,
    passes_fc1p5_bh05 = is.finite(p_val_adj) & p_val_adj < 0.05 &
      is.finite(standard_log2fc) & abs(standard_log2fc) >= log2(1.5))]
  d[, c("flt3l_mean_cp10k_wt", "flt3l_mean_cp10k_ko",
    "il7_mean_cp10k_wt", "il7_mean_cp10k_ko") := NULL]
  setcolorder(d, c("region_id", "target", "cytokine", "source_peak",
    "first_bin", "last_bin", "pt_left", "pt_right", "n_wt_cells",
    "n_ko_cells", "feature_id", "direction", "mean_cp10k_wt",
    "mean_cp10k_ko", "standard_log2fc", "significant_bh05",
    "passes_fc1p5_bh05", "pct_a", "pct_b", "avg_expr_a", "avg_expr_b",
    "avg_log1p_fc", "p_val", "p_val_adj"))
  d
}), use.names = TRUE)
combined[, abs_standard_log2fc := abs(standard_log2fc)]
setorder(combined, region_id, p_val_adj, p_val, -abs_standard_log2fc,
  feature_id)
combined[, abs_standard_log2fc := NULL]
write_tsv(combined, file.path(stage_dir,
  "all_genes_shared_regions_normalized_counts_and_de.tsv"))
for (reg in regions$region_id) write_tsv(combined[region_id == reg],
  file.path(stage_dir, paste0(tolower(reg),
    "_all_genes_normalized_counts_and_de.tsv")))

selected <- combined[passes_fc1p5_bh05 == TRUE]
selected[, abs_standard_log2fc := abs(standard_log2fc)]
setorder(selected, region_id, p_val_adj, p_val, -abs_standard_log2fc,
  feature_id)
selected[, abs_standard_log2fc := NULL]
selected[, significance_rank := seq_len(.N), by = region_id]
write_tsv(selected, file.path(stage_dir,
  "genes_passing_fc1p5_bh05.tsv"))
selection_summary <- combined[, .(
  genes_tested = .N,
  significant_bh05 = sum(significant_bh05),
  passing_fc1p5_bh05 = sum(passes_fc1p5_bh05),
  passing_wt_higher = sum(passes_fc1p5_bh05 & direction == "WT_higher"),
  passing_ko_higher = sum(passes_fc1p5_bh05 & direction == "KO_higher")),
  by = .(region_id, cytokine, first_bin, last_bin, pt_left, pt_right,
    n_wt_cells, n_ko_cells)]
write_tsv(selection_summary, file.path(stage_dir,
  "differential_expression_summary.tsv"))

plot_genes <- selected[, head(.SD, 20L), by = .(region_id, direction)]
setorder(plot_genes, region_id, direction, p_val_adj, p_val)
plot_genes[, plot_rank := seq_len(.N), by = .(region_id, direction)]
write_tsv(plot_genes, file.path(stage_dir,
  "violin_plot_gene_selection.tsv"))
target_features <- c("Flt3", "Il7r")
plot_features <- unique(c(plot_genes$feature_id, target_features))
plot_cells <- as.character(assignments$source_cell_id)
expr <- wb_matrix_read_features(staged, matrix_id, plot_features,
  plot_cells, missing = "NA")
stopifnot(!anyNA(expr))
expr_long <- as.data.table(as.table(expr))
setnames(expr_long, c("feature_id", "source_cell_id", "log1p_cp10k"))
expr_long[, `:=`(feature_id = as.character(feature_id),
  source_cell_id = as.character(source_cell_id),
  log1p_cp10k = as.numeric(log1p_cp10k))]
expr_long <- merge(expr_long, assignments[, .(source_cell_id, region_id,
  genotype)], by = "source_cell_id", all.x = TRUE, sort = FALSE)
fwrite(expr_long[feature_id %in% plot_features], file.path(stage_dir,
  "violin_plot_normalized_expression.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")
fwrite(expr_long[feature_id %in% target_features], file.path(stage_dir,
  "target_genes_flt3_il7r_normalized_expression.tsv.gz"), sep = "\t",
  na = "NA", compress = "gzip")

violin_manifest <- list()
manifest_index <- 0L
for (i in seq_len(nrow(regions))) {
  reg <- regions[i]
  for (direction_value in c("WT_higher", "KO_higher")) {
  genes <- plot_genes[region_id == reg$region_id &
    direction == direction_value][order(plot_rank)]
  pd <- expr_long[region_id == reg$region_id &
    feature_id %in% genes$feature_id]
  pd[, feature_id := factor(feature_id, levels = genes$feature_id)]
  pd[, genotype := factor(genotype, levels = c("WT", "KO"))]
  genotype_colours <- c("WT" = "#D55E00", "KO" =
    if (reg$target == "Flt3l") "#0072B2" else "#009E73")
  padj_labels <- setNames(ifelse(genes$p_val_adj < 1e-300,
    "<1e-300", formatC(genes$p_val_adj, format = "e", digits = 2)),
    genes$feature_id)
  facet_labels <- setNames(paste0(genes$feature_id, "\nBH padj: ",
    padj_labels[genes$feature_id]), genes$feature_id)
  p <- ggplot(pd, aes(genotype, log1p_cp10k, fill = genotype)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.25,
      colour = "grey30") +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.6,
      fill = "white", colour = "black", stroke = 0.35) +
    facet_wrap(~feature_id, scales = "free_y", ncol = 4,
      labeller = as_labeller(facet_labels)) +
    scale_fill_manual(values = genotype_colours) +
    labs(x = NULL, y = "Expression (log1p CP10K)",
      title = paste0(reg$cytokine, " shared-enriched region: ",
        fifelse(direction_value == "WT_higher", "WT-higher", "KO-higher"),
        " genes"),
      subtitle = paste0("Top ", nrow(genes), " of ",
        selected[region_id == reg$region_id & direction == direction_value,
          .N], " passing genes in this direction; BH P < 0.05 and ",
        "|log2FC| >= log2(1.5)"),
      caption = "White points mark medians; one biological sample per group.") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 5.5, 5.5, 18))
  direction_stem <- if (direction_value == "WT_higher") "wt_higher" else
    "ko_higher"
  stem <- paste0(tolower(reg$region_id), "_", direction_stem,
    "_top20_cutoff_genes_violin")
  plot_height <- max(6.5, 2.35 * ceiling(nrow(genes) / 4))
  ggsave(file.path(stage_dir, "figures", paste0(stem, ".pdf")), p,
    width = 12, height = plot_height, device = cairo_pdf)
  ggsave(file.path(stage_dir, "figures", paste0(stem, ".png")), p,
    width = 12, height = plot_height, dpi = 300)
  manifest_index <- manifest_index + 1L
  violin_manifest[[manifest_index]] <- data.table(region_id = reg$region_id,
    cytokine = reg$cytokine, direction = direction_value,
    passing_direction_genes = selected[region_id == reg$region_id &
      direction == direction_value, .N], plotted_genes = nrow(genes),
    png = file.path("figures", paste0(stem, ".png")),
    pdf = file.path("figures", paste0(stem, ".pdf")))
  }
}
violin_manifest <- rbindlist(violin_manifest)
write_tsv(violin_manifest, file.path(stage_dir, "violin_plot_manifest.tsv"))

target_manifest <- list()
for (i in seq_len(nrow(regions))) {
  reg <- regions[i]
  stats <- combined[region_id == reg$region_id &
    feature_id %in% target_features]
  stopifnot(nrow(stats) == 2L)
  pd <- expr_long[region_id == reg$region_id &
    feature_id %in% target_features]
  pd[, feature_id := factor(feature_id, levels = target_features)]
  pd[, genotype := factor(genotype, levels = c("WT", "KO"))]
  genotype_colours <- c("WT" = "#D55E00", "KO" =
    if (reg$target == "Flt3l") "#0072B2" else "#009E73")
  stats[, padj_text := fifelse(p_val_adj < 1e-300, "<1e-300",
    formatC(p_val_adj, format = "e", digits = 2))]
  facet_labels <- setNames(paste0(stats$feature_id, "\nlog2FC: ",
    sprintf("%.3f", stats$standard_log2fc), "; BH padj: ",
    stats$padj_text), stats$feature_id)
  p <- ggplot(pd, aes(genotype, log1p_cp10k, fill = genotype)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.25,
      colour = "grey30") +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.6,
      fill = "white", colour = "black", stroke = 0.35) +
    facet_wrap(~feature_id, scales = "free_y", nrow = 1,
      labeller = as_labeller(facet_labels)) +
    scale_fill_manual(values = genotype_colours) +
    labs(x = NULL, y = "Expression (log1p CP10K)",
      title = paste0(reg$cytokine,
        " shared-enriched region: Flt3 and Il7r"),
      subtitle = "Standard log2FC is WT / KO; P values are genome-wide BH adjusted",
      caption = "White points mark medians; one biological sample per group.") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 5.5, 5.5, 18))
  stem <- paste0(tolower(reg$region_id),
    "_flt3_il7r_target_gene_violin")
  ggsave(file.path(stage_dir, "figures", paste0(stem, ".pdf")), p,
    width = 8.5, height = 4.8, device = cairo_pdf)
  ggsave(file.path(stage_dir, "figures", paste0(stem, ".png")), p,
    width = 8.5, height = 4.8, dpi = 300)
  target_manifest[[i]] <- data.table(region_id = reg$region_id,
    cytokine = reg$cytokine,
    png = file.path("figures", paste0(stem, ".png")),
    pdf = file.path("figures", paste0(stem, ".pdf")))
}
target_manifest <- rbindlist(target_manifest)
write_tsv(target_manifest, file.path(stage_dir,
  "target_gene_violin_manifest.tsv"))

volcano_data <- copy(combined)
volcano_data[, neglog10_padj := -log10(pmax(p_val_adj, 1e-300))]
volcano_data[, status := fifelse(passes_fc1p5_bh05 &
  standard_log2fc > 0, "WT higher: passes both", fifelse(
  passes_fc1p5_bh05 & standard_log2fc < 0 & target == "Flt3l",
  "FLT3LKO higher: passes both", fifelse(passes_fc1p5_bh05 &
  standard_log2fc < 0 & target == "Il7", "IL7KO higher: passes both",
  "Does not pass both")))]
volcano_data[, region_label := factor(region_id,
  levels = regions$region_id, labels = paste0(regions$cytokine,
    " shared region"))]
volcano_labels <- volcano_data[is.finite(p_val_adj) & p_val_adj < 0.05 &
  is.finite(standard_log2fc) & abs(standard_log2fc) >= 1]
volcano_labels[, side := fifelse(standard_log2fc > 0, "positive",
  "negative")]
volcano_labels[, label_y := seq(max(neglog10_padj),
  max(neglog10_padj) * 0.32, length.out = .N),
  by = .(region_id, side)]
volcano_ranges <- volcano_data[is.finite(standard_log2fc),
  .(xmax = max(abs(standard_log2fc))), by = region_id]
volcano_labels <- merge(volcano_labels, volcano_ranges, by = "region_id")
volcano_labels[, `:=`(label_x = fifelse(standard_log2fc > 0,
  0.98 * xmax, -0.98 * xmax), label_hjust = fifelse(
  standard_log2fc > 0, 1, 0))]
p_volcano <- ggplot(volcano_data[is.finite(neglog10_padj)],
    aes(standard_log2fc, neglog10_padj, colour = status)) +
  geom_point(size = 0.7, alpha = 0.55) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
    linetype = "dashed", colour = "grey35") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted",
    colour = "grey35") +
  geom_segment(data = volcano_labels,
    aes(x = standard_log2fc, y = neglog10_padj, xend = label_x,
      yend = label_y), inherit.aes = FALSE, linewidth = 0.25,
    colour = "grey35") +
  geom_text(data = volcano_labels,
    aes(x = label_x, y = label_y, label = feature_id,
      hjust = label_hjust), inherit.aes = FALSE, size = 2.7) +
  facet_wrap(~region_label, ncol = 2) +
  scale_colour_manual(values = c("WT higher: passes both" = "#D55E00",
    "FLT3LKO higher: passes both" = "#0072B2",
    "IL7KO higher: passes both" = "#009E73",
    "Does not pass both" = "grey75")) +
  labs(x = "Standard log2FC (WT / KO)", y = "-log10(BH-adjusted P)",
    colour = NULL, title = "Full shared-region differential expression",
    subtitle = paste("All 33,696 genes tested in each comparison; labels:",
      "BH padj < 0.05 and |log2FC| >= 1")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"))
ggsave(file.path(stage_dir, "figures", "shared_region_de_volcano.pdf"),
  p_volcano, width = 11, height = 5.8, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "shared_region_de_volcano.png"),
  p_volcano, width = 11, height = 5.8, dpi = 300)

summary_rows <- vapply(seq_len(nrow(selection_summary)), function(i) {
  x <- selection_summary[i]
  paste0("<tr><td>", x$cytokine, "</td><td>", x$first_bin, "-",
    x$last_bin, "</td><td>", x$n_wt_cells, "</td><td>",
    x$n_ko_cells, "</td><td>", x$genes_tested, "</td><td>",
    x$significant_bh05, "</td><td>", x$passing_fc1p5_bh05,
    "</td><td>", x$passing_wt_higher, "</td><td>",
    x$passing_ko_higher, "</td></tr>")
}, character(1))
violin_blocks <- vapply(seq_len(nrow(violin_manifest)), function(i) {
  x <- violin_manifest[i]
  direction_label <- if (x$direction == "WT_higher") "WT-higher" else
    "KO-higher"
  paste0("<h2>", x$cytokine, ": ", direction_label,
    " violin plots</h2><p>", x$passing_direction_genes,
    " genes pass both cutoffs in this direction; the top ", x$plotted_genes,
    " are shown. <a href='", x$pdf,
    "'>Vector PDF</a></p><img src='", x$png, "'>")
}, character(1))
target_blocks <- vapply(seq_len(nrow(target_manifest)), function(i) {
  x <- target_manifest[i]
  paste0("<h2>", x$cytokine,
    ": Flt3 and Il7r</h2><p>Targeted normalized-expression violins with ",
    "standard log2FC and BH-adjusted P. <a href='", x$pdf,
    "'>Vector PDF</a></p><img src='", x$png, "'>")
}, character(1))
html <- c("<!doctype html><html><head><meta charset='utf-8'>",
  "<title>Shared-enriched trajectory-region RNA DE</title>",
  "<style>body{font-family:sans-serif;max-width:1400px;margin:30px auto;line-height:1.45}img{max-width:100%}table{border-collapse:collapse}th,td{border:1px solid #bbb;padding:6px 9px;text-align:right}th:first-child,td:first-child{text-align:left}.note{background:#fff3cd;padding:12px}</style></head><body>",
  "<h1>Shared WT-and-KO enriched trajectory-region RNA differential expression</h1>",
  "<p>Regions are maximal contiguous bins whose Bdev-normalized rolling cell fractions are strictly above 1 in both WT and KO. IL-7 is restricted to the original peak 2. The saved native scWorkbench log1p-CP10K matrix was reused; new native Wilcoxon tests were required for the changed cell sets.</p>",
  "<p class='note'><strong>Inference limitation:</strong> each genotype contains one biological sample per cytokine. Results describe cell-level separation and are not replicate-level genotype tests.</p>",
  "<img src='figures/shared_enriched_regions.png'>",
  "<h2>DE summary</h2><table><tr><th>Cytokine</th><th>Bins</th><th>WT cells</th><th>KO cells</th><th>Genes tested</th><th>BH &lt; 0.05</th><th>Both cutoffs</th><th>WT higher</th><th>KO higher</th></tr>",
  summary_rows, "</table>",
  "<p><a href='all_genes_shared_regions_normalized_counts_and_de.tsv'>Download complete normalized-count and DEG table</a> | <a href='genes_passing_fc1p5_bh05.tsv'>Download cutoff-passing genes</a></p>",
  "<h2>Full DEG volcano plots</h2><p><a href='figures/shared_region_de_volcano.pdf'>Vector PDF</a></p><img src='figures/shared_region_de_volcano.png'>",
  "<h2>Targeted Flt3 and Il7r expression</h2>", target_blocks,
  violin_blocks, "</body></html>")
writeLines(html, file.path(stage_dir, "REPORT.html"))

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
violin_expected <- selected[, .(expected = min(.N, 20L)),
  by = .(region_id, direction)]
violin_actual <- plot_genes[, .(actual = .N),
  by = .(region_id, direction)]
violin_count_check <- merge(violin_expected, violin_actual,
  by = c("region_id", "direction"), all = TRUE)
checks <- list(
  exact_shared_regions = identical(regions$first_bin, c(19L, 21L)) &&
    identical(regions$last_bin, c(26L, 31L)),
  both_genotypes_above_one_in_every_region_bin = all(shared_bins$WT > 1 &
    shared_bins$KO > 1),
  il7_restricted_to_original_peak2 = regions[target == "Il7",
    first_bin] >= il7_peak2$first_bin && regions[target == "Il7",
    last_bin] <= il7_peak2$last_bin,
  both_genotypes_in_each_region = all(assignments[, uniqueN(genotype),
    by = region_id]$V1 == 2L),
  all_genes_tested_in_both_comparisons = nrow(combined) == 2L * 33696L,
  normalized_means_complete = all(is.finite(combined$mean_cp10k_wt)) &&
    all(is.finite(combined$mean_cp10k_ko)),
  cutoff_selection_valid = all(selected$p_val_adj < 0.05) &&
    all(abs(selected$standard_log2fc) >= log2(1.5)),
  violin_top20_per_direction = nrow(violin_count_check) == 4L &&
    all(violin_count_check$expected == violin_count_check$actual),
  figures_complete = file.exists(file.path(stage_dir, "figures",
    "shared_enriched_regions.png")) && file.exists(file.path(stage_dir,
    "figures", "shared_region_de_volcano.png")) &&
    all(file.exists(file.path(stage_dir, violin_manifest$png))) &&
    all(file.exists(file.path(stage_dir, target_manifest$png))),
  parent_bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after), result_promotions_zero =
    is.null(session$durable$promotion_receipt))
stopifnot(all(unlist(checks)))

completion <- list(status = "complete_unpromoted_report_only",
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script_version = script_version, report = file.path(output_dir, "REPORT.html"),
  source_bundle = ko_bundle_dir,
  source_bundle_fingerprint = fingerprint_before,
  source_distribution_report = distribution_report,
  source_projection_report = projection_report,
  source_original_region_report = original_region_report,
  source_session_workspace = session_workspace,
  normalized_matrix = matrix_id,
  normalization_recomputed = FALSE,
  new_operations = c("de_wilcox:de_shared_flt3l",
    "de_wilcox:de_shared_il7_peak2"),
  region_rule = "maximal contiguous bins with both WT and KO Bdev-normalized rolling fraction strictly above 1; IL-7 restricted to original Il7_peak2",
  regions = lapply(seq_len(nrow(regions)), function(i) as.list(regions[i, !"run"])),
  standard_log2fc = "log2((mean linear CP10K WT + 1)/(mean linear CP10K KO + 1))",
  selection_cutoffs = list(p_val_adj = 0.05,
    absolute_log2fc = log2(1.5), absolute_fold_change = 1.5),
  violin_rule = "up to 20 most significant cutoff-passing genes per direction and comparison, with a separate figure for each direction",
  violin_colour_scale = list(WT = "#D55E00", FLT3L_KO = "#0072B2",
    IL7_KO = "#009E73"),
  volcano_label_rule = "BH-adjusted P < 0.05 and absolute standard log2FC >= 1",
  targeted_gene_plots = list(genes = target_features,
    comparisons = as.list(target_manifest$region_id)),
  differential_summary = lapply(seq_len(nrow(selection_summary)),
    function(i) as.list(selection_summary[i])),
  inference_scope = "descriptive_cell_level_only_one_biological_sample_per_group",
  checks = checks, promotions = 0L)
write_json(completion, file.path(stage_dir, "completion.json"))
stopifnot(file.rename(stage_dir, output_dir))
cat("Completed:", file.path(output_dir, "REPORT.html"), "\n")
