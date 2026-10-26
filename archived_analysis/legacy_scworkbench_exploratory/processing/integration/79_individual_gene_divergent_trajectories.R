#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "Usage: 79_individual_gene_divergent_trajectories.R",
  "<dynamic_gsea_dir> <session_stage_bundle> <projection.tsv.gz>",
  "<output_dir> <top_report_n> <top_curve_n>"))
input_dir <- normalizePath(args[[1L]], mustWork = TRUE)
stage_bundle_dir <- normalizePath(args[[2L]], mustWork = TRUE)
projection_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
top_report_n <- as.integer(args[[5L]])
top_curve_n <- as.integer(args[[6L]])
stopifnot(is.finite(top_report_n), top_report_n >= 1L,
  is.finite(top_curve_n), top_curve_n >= 1L,
  top_curve_n <= top_report_n)
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

ranking_path <- file.path(input_dir, "gene_dynamic_change_rankings.tsv.gz")
endpoints_path <- file.path(input_dir, "endpoint_bin_definitions.tsv")
diagnostics_path <- file.path(input_dir, "diagnostics.json")
stopifnot(file.exists(ranking_path), file.exists(endpoints_path),
  file.exists(diagnostics_path))
genes <- fread(cmd = paste("gzip -cd", shQuote(ranking_path)))

# A gene must satisfy the biological sign definition under both endpoint
# policies. Ranking by the smaller absolute divergence prevents a large result
# under only one policy from dominating the report.
genes[, pattern := fcase(
  wt_change > 0 & ko_change < wt_change &
    sensitivity_wt_change > 0 &
    sensitivity_ko_change < sensitivity_wt_change,
  "WT induction attenuated in KO",
  wt_change < 0 & ko_change > wt_change &
    sensitivity_wt_change < 0 &
    sensitivity_ko_change > sensitivity_wt_change,
  "Failure to silence in KO",
  default = NA_character_)]
divergent <- genes[!is.na(pattern)]
divergent[, robust_divergence := pmin(abs(differential_change),
  abs(sensitivity_differential_change))]
divergent[, primary_direction_consistent :=
  sign(differential_change) == sign(sensitivity_differential_change)]
stopifnot(all(divergent$primary_direction_consistent),
  all(divergent$robust_divergence >= 0))
setorder(divergent, comparison_id, pattern, -robust_divergence,
  feature_id)
divergent[, pattern_rank := seq_len(.N),
  by = .(comparison_id, pattern)]
divergent[, robust_percentile := 1 - (pattern_rank - 1) / .N,
  by = .(comparison_id, pattern)]
fwrite(divergent, file.path(stage_dir,
  "all_direction_consistent_divergent_genes.tsv.gz"), sep = "\t",
  compress = "gzip")

inventory <- divergent[, .(direction_consistent_genes = .N,
  median_robust_divergence = median(robust_divergence),
  robust_divergence_at_rank_10 = robust_divergence[min(10L, .N)],
  robust_divergence_at_rank_100 = robust_divergence[min(100L, .N)]),
  by = .(comparison_id, region_id, pattern)]
fwrite(inventory, file.path(stage_dir, "divergent_gene_inventory.tsv"),
  sep = "\t")

top <- divergent[, head(.SD, top_report_n),
  by = .(comparison_id, region_id, pattern)]
fwrite(top, file.path(stage_dir, "top_divergent_genes.tsv"),
  sep = "\t", na = "NA")
curve_selection <- top[pattern_rank <= top_curve_n]

parent_diagnostics <- jsonlite::read_json(diagnostics_path,
  simplifyVector = TRUE)
source_dir <- normalizePath(parent_diagnostics$source_run, mustWork = TRUE)
source_diagnostics_path <- file.path(source_dir, "diagnostics.json")
stopifnot(file.exists(source_diagnostics_path))
source_diagnostics <- jsonlite::read_json(source_diagnostics_path,
  simplifyVector = TRUE)
means_path <- file.path(source_dir,
  "pseudotime_bin_gene_expression_means.tsv.gz")
stopifnot(file.exists(means_path))
means <- fread(cmd = paste("gzip -cd", shQuote(means_path)))
endpoints <- fread(endpoints_path)
means <- merge(means,
  endpoints[, .(region_id, comparison_id, bin, primary_window)],
  by = c("region_id", "comparison_id", "bin"), all.x = TRUE)
means[, `:=`(log2_wt = log2(mean_cp10k_wt + 1),
  log2_ko = log2(mean_cp10k_ko + 1))]

curve_data <- merge(curve_selection[, .(comparison_id, region_id, pattern,
    feature_id, pattern_rank, wt_change, ko_change, differential_change,
    sensitivity_wt_change, sensitivity_ko_change,
    sensitivity_differential_change, robust_divergence,
    maximum_mean_cp10k)],
  means[, .(comparison_id, region_id, bin, bin_midpoint, bdev_weight,
    reliable_20, primary_window, feature_id, log2_wt, log2_ko)],
  by = c("comparison_id", "region_id", "feature_id"))
stopifnot(uniqueN(curve_data[, .(comparison_id, pattern, feature_id)]) ==
  nrow(curve_selection))
baselines <- curve_data[primary_window == "early", .(
  wt_baseline = weighted.mean(log2_wt, bdev_weight),
  ko_baseline = weighted.mean(log2_ko, bdev_weight)),
  by = .(comparison_id, region_id, pattern, feature_id)]
curve_data <- merge(curve_data, baselines,
  by = c("comparison_id", "region_id", "pattern", "feature_id"))
curve_data[, `:=`(WT = log2_wt - wt_baseline,
  KO = log2_ko - ko_baseline)]
curves <- melt(curve_data,
  id.vars = c("comparison_id", "region_id", "pattern", "feature_id",
    "pattern_rank", "bin", "bin_midpoint", "reliable_20",
    "primary_window", "wt_change", "ko_change", "differential_change",
    "sensitivity_wt_change", "sensitivity_ko_change",
    "sensitivity_differential_change", "robust_divergence",
    "maximum_mean_cp10k"), measure.vars = c("WT", "KO"),
  variable.name = "genotype",
  value.name = "baseline_centered_log2_expression")
fwrite(curves, file.path(stage_dir, "top_divergent_gene_curves.tsv.gz"),
  sep = "\t", compress = "gzip")
mean_curves <- copy(curves)

# Cell-level IQR display. Selection and endpoint ranking remain based on the
# saved bin means above. For visualization, read only the selected curve genes
# from the retained scWorkbench log1p-CP10K matrix and summarize the exact
# mapped cells in each bin by median and 25th/75th percentiles.
counts_path <- file.path(source_dir, "pseudotime_bin_cell_counts.tsv")
stopifnot(file.exists(counts_path))
bin_counts <- fread(counts_path)
projection <- fread(cmd = paste("gzip -cd", shQuote(projection_path)))[
  support_status %chin% c("supported", "borderline") & is.finite(mapped_pt),
  .(source_cell_id, genotype, knockout_target, mapped_pt)]
stopifnot(!anyDuplicated(projection$source_cell_id))
assignments <- rbindlist(lapply(seq_len(nrow(bin_counts)), function(i) {
  b <- bin_counts[i]
  target <- if (b$comparison_id == "Flt3L_KO_vs_WT") "Flt3l" else "Il7"
  projection[knockout_target == target & mapped_pt >= b$bin_left &
      mapped_pt < b$bin_right,
    .(source_cell_id, genotype, mapped_pt, region_id = b$region_id,
      comparison_id = b$comparison_id, bin = b$bin,
      bin_midpoint = b$bin_midpoint, bdev_weight = b$bdev_weight,
      reliable_20 = b$reliable_20)]
}))
stopifnot(!anyDuplicated(assignments$source_cell_id))
assignment_counts <- assignments[, .N,
  by = .(comparison_id, region_id, bin, genotype)]
expected_counts <- melt(bin_counts,
  id.vars = c("comparison_id", "region_id", "bin"),
  measure.vars = c("n_wt", "n_ko"), variable.name = "genotype",
  value.name = "expected_n")
expected_counts[, genotype := fifelse(genotype == "n_wt", "WT", "KO")]
count_audit <- merge(expected_counts, assignment_counts,
  by = c("comparison_id", "region_id", "bin", "genotype"), all = TRUE)
stopifnot(all(count_audit$expected_n == count_audit$N))

staged <- wb_open_bundle(stage_bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(try(wb_close_bundle(staged), silent = TRUE), add = TRUE)
matrix_id <- source_diagnostics$normalized_matrix
stopifnot(matrix_id %in% wb_matrices(staged)$matrix_id)
curve_features <- unique(curve_selection$feature_id)
cell_order <- assignments$source_cell_id
matrix <- wb_matrix_read_features(staged, matrix_id, curve_features,
  cell_order, missing = "NA")
stopifnot(!anyNA(matrix), identical(rownames(matrix), curve_features),
  identical(colnames(matrix), cell_order))
# scWorkbench log1p-CP10K uses natural log; divide by log(2) to obtain
# log2(CP10K + 1) on the same display scale as the endpoint analysis.
cell_log2 <- as.matrix(matrix) / log(2)
cell_iqr <- rbindlist(lapply(seq_along(curve_features), function(gene_i) {
  gene <- curve_features[[gene_i]]
  rbindlist(lapply(seq_len(nrow(bin_counts)), function(bin_i) {
    b <- bin_counts[bin_i]
    rbindlist(lapply(c("WT", "KO"), function(group) {
      idx <- which(assignments$comparison_id == b$comparison_id &
        assignments$bin == b$bin & assignments$genotype == group)
      values <- cell_log2[gene_i, idx]
      q <- unname(quantile(values, c(0.25, 0.5, 0.75), names = FALSE,
        type = 7))
      data.table(comparison_id = b$comparison_id, region_id = b$region_id,
        feature_id = gene, bin = b$bin, bin_midpoint = b$bin_midpoint,
        bdev_weight = b$bdev_weight, reliable_20 = b$reliable_20,
        genotype = group, n_cells = length(idx), q25 = q[[1L]],
        median = q[[2L]], q75 = q[[3L]])
    }))
  }))
}))
cell_iqr <- merge(cell_iqr,
  endpoints[, .(comparison_id, region_id, bin, primary_window)],
  by = c("comparison_id", "region_id", "bin"))
iqr_baselines <- cell_iqr[primary_window == "early", .(
  baseline_median = weighted.mean(median, bdev_weight)),
  by = .(comparison_id, region_id, feature_id, genotype)]
cell_iqr <- merge(cell_iqr, iqr_baselines,
  by = c("comparison_id", "region_id", "feature_id", "genotype"))
cell_iqr[, `:=`(centered_q25 = q25 - baseline_median,
  centered_median = median - baseline_median,
  centered_q75 = q75 - baseline_median)]
cell_iqr <- merge(cell_iqr,
  curve_selection[, .(comparison_id, region_id, pattern, feature_id,
    pattern_rank, wt_change, ko_change, differential_change,
    sensitivity_differential_change, robust_divergence,
    maximum_mean_cp10k)],
  by = c("comparison_id", "region_id", "feature_id"))
cell_iqr <- merge(cell_iqr,
  mean_curves[, .(comparison_id, region_id, feature_id, bin, genotype,
    mean_curve = baseline_centered_log2_expression)],
  by = c("comparison_id", "region_id", "feature_id", "bin", "genotype"))
fwrite(cell_iqr, file.path(stage_dir,
  "top_divergent_gene_cell_iqr_curves.tsv.gz"), sep = "\t",
  compress = "gzip")
curves <- cell_iqr
wb_close_bundle(staged)

comparison_label <- c(Flt3L_KO_vs_WT = "FLT3L", IL7_KO_vs_WT = "IL-7")
pattern_stub <- c(`WT induction attenuated in KO` = "induction_failure",
  `Failure to silence in KO` = "silencing_failure")
for (current_comparison in unique(curves$comparison_id)) {
  for (current_pattern in unique(curves[comparison_id == current_comparison]$pattern)) {
    d <- copy(curves[comparison_id == current_comparison &
      pattern == current_pattern])
    label_data <- unique(d[, .(feature_id, pattern_rank, wt_change,
      ko_change, robust_divergence)])
    setorder(label_data, pattern_rank)
    label_data[, facet := sprintf(
      "%s\nWT %+.2f, KO %+.2f; robust |Δ| %.2f",
      feature_id, wt_change, ko_change, robust_divergence)]
    d <- merge(d, label_data[, .(feature_id, facet)], by = "feature_id")
    d[, facet := factor(facet, levels = label_data$facet)]
    p <- ggplot(d, aes(bin_midpoint, mean_curve,
        colour = genotype, fill = genotype, group = genotype)) +
      geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
      geom_ribbon(aes(ymin = centered_q25, ymax = centered_q75),
        alpha = 0.16, colour = NA) +
      geom_line(linewidth = 0.72) +
      geom_point(aes(shape = reliable_20), size = 1.6) +
      facet_wrap(~facet, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c(WT = "#B2182B", KO = "#2166AC")) +
      scale_fill_manual(values = c(WT = "#B2182B", KO = "#2166AC"),
        guide = "none") +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
        labels = c(`TRUE` = "≥20 cells/genotype",
          `FALSE` = "Sparse terminal"), name = "Bin coverage") +
      theme_bw(base_size = 8.5) +
      theme(strip.text = element_text(size = 7), legend.position = "bottom",
        panel.grid.minor = element_blank()) +
      labs(title = paste(comparison_label[[current_comparison]],
          current_pattern, "genes"),
        subtitle = paste("Mean-expression line with cell-level IQR; top", top_curve_n,
          "genes ranked by the smaller primary/sensitivity mean divergence"),
        x = "Mapped Bdev pseudotime",
        y = "Mean-expression change with cell-level IQR",
        colour = "Genotype")
    stub <- paste0(tolower(sub("_KO_vs_WT", "", current_comparison)), "_",
      pattern_stub[[current_pattern]], "_gene_curves")
    height <- max(8, 2.5 + 2.45 * ceiling(uniqueN(d$feature_id) / 3))
    ggsave(file.path(stage_dir, "figures", paste0(stub, ".png")), p,
      width = 12, height = height, dpi = 300, bg = "white")
    ggsave(file.path(stage_dir, "figures", paste0(stub, ".pdf")), p,
      width = 12, height = height, device = cairo_pdf, bg = "white")
  }
}

lines <- c("# Individual-gene divergent trajectory screen", "",
  "This descriptive screen applies the same early-to-late trajectory logic used for gene sets to individual genes.", "",
  "- WT induction attenuated in KO: WT change > 0 and KO change < WT change.",
  "- Failure to silence in KO: WT change < 0 and KO change > WT change.",
  "- A gene must meet the same sign definition in both the primary reliable-bin endpoints and the sensitivity endpoints that include the sparse terminal bin.",
  "- Ranking score: the smaller absolute WT-minus-KO differential change across the two endpoint policies.",
  "- Curves retain the selected log2(mean CP10K + 1) trajectory as the line and add the cell-level interquartile range of log2(CP10K + 1) as a ribbon. Lines are centered to the mean-expression early baseline; ribbons are centered to the cell-median early baseline. Their separation reveals zero inflation or skew rather than sampling uncertainty.",
  "- No gene-level P value or FDR is assigned: there is one biological sample per genotype, and cells/bins are not biological replicates.", "")
for (current_comparison in unique(inventory$comparison_id)) {
  lines <- c(lines, paste0("## ", current_comparison), "")
  for (current_pattern in c("WT induction attenuated in KO",
      "Failure to silence in KO")) {
    inv <- inventory[comparison_id == current_comparison &
      pattern == current_pattern]
    d <- top[comparison_id == current_comparison & pattern == current_pattern]
    lines <- c(lines, paste0("### ", current_pattern), "",
      sprintf("%s genes retain this direction under both endpoint policies. The top %d are shown.",
        format(inv$direction_consistent_genes, big.mark = ","), nrow(d)), "",
      "| Rank | Gene | WT change | KO change | Differential | Sensitivity differential | Robust score | Max mean CP10K |",
      "|---:|---|---:|---:|---:|---:|---:|---:|",
      vapply(seq_len(nrow(d)), function(i) sprintf(
        "| %d | %s | %+.3f | %+.3f | %+.3f | %+.3f | %.3f | %.3f |",
        d$pattern_rank[[i]], d$feature_id[[i]], d$wt_change[[i]],
        d$ko_change[[i]], d$differential_change[[i]],
        d$sensitivity_differential_change[[i]],
        d$robust_divergence[[i]], d$maximum_mean_cp10k[[i]]),
        character(1)), "")
  }
}
lines <- c(lines, "## Interpretation", "",
  "The ranking identifies reproducible endpoint divergence within these measured samples, not gene-level population inference. Curves should be used to distinguish gradual divergence from a single-bin transition, and gene identity should be interpreted alongside expression magnitude and known biology.")
writeLines(lines, file.path(stage_dir, "REPORT.md"))

diagnostics <- list(status = "complete_report_only",
  source_dynamic_gsea = input_dir,
  primary_endpoints = "first/last two bins with >=20 cells per genotype",
  sensitivity_endpoints = "literal first/last two mapped-region bins",
  direction_consistency_required = TRUE,
  ranking = paste("minimum absolute differential change across primary and",
    "sensitivity endpoint policies"),
  top_report_n = top_report_n, top_curve_n = top_curve_n,
  curve_summary = paste("log2(mean CP10K + 1) line plus cell-level IQR of",
    "log2(CP10K + 1); not a confidence interval"),
  normalized_matrix = matrix_id,
  stage_bundle = stage_bundle_dir,
  projection = projection_path,
  inferential_gene_level_p_values = FALSE,
  reason = "one biological sample per genotype",
  inventory = inventory)
jsonlite::write_json(diagnostics, file.path(stage_dir, "diagnostics.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
message("Completed individual-gene divergence report at: ", output_dir)
