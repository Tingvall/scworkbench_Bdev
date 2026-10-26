#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) stop(paste(
  "Usage: 12_wt_enriched_peak_de.R KO_BUNDLE DISTRIBUTION_REPORT",
  "PROJECTION_REPORT SESSION_ROOT OUTPUT_DIR EXPECTED_KO_FP",
  "PARSE_SIGNAL_REPORT LINEAGE_WORKBOOK LINEAGE_TSV",
  "EXPECTED_LINEAGE_SHA SCRIPT_VERSION"
), call. = FALSE)

ko_bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
distribution_report <- normalizePath(args[[2L]], mustWork = TRUE)
projection_report <- normalizePath(args[[3L]], mustWork = TRUE)
session_root <- normalizePath(args[[4L]], mustWork = FALSE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
expected_ko_fp <- args[[6L]]
parse_signal_report <- normalizePath(args[[7L]], mustWork = TRUE)
lineage_workbook <- normalizePath(args[[8L]], mustWork = TRUE)
lineage_tsv <- normalizePath(args[[9L]], mustWork = TRUE)
expected_lineage_sha <- args[[10L]]
script_version <- args[[11L]]

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output already exists: ", stage_dir, call. = FALSE)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage_dir, "figures"), showWarnings = FALSE)
dir.create(session_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(scworkbench)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

write_json <- function(x, path) jsonlite::write_json(
  x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA,
  null = "null", na = "null")
write_tsv <- function(x, path) data.table::fwrite(x, path, sep = "\t", na = "NA")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}
log_event <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paste(..., collapse = " "), "\n")

distribution_path <- file.path(distribution_report,
  "query_mapped_pt_cell_distribution.tsv")
projection_path <- file.path(projection_report,
  "query_reference_umap_projection.tsv.gz")
parse_signal_membership_path <- file.path(parse_signal_report,
  "passing_gene_membership.tsv")
parse_signal_completion_path <- file.path(parse_signal_report,
  "completion.json")
bb_boundary_path <- file.path(parse_signal_report, "bb_up_change_points.tsv")
required <- c(distribution_path, projection_path,
  file.path(distribution_report, "completion.json"),
  file.path(projection_report, "completion.json"), lineage_workbook,
  lineage_tsv, parse_signal_membership_path, parse_signal_completion_path)
required <- c(required, bb_boundary_path)
stopifnot(all(file.exists(required)))
stopifnot(identical(sha256_file(lineage_workbook), expected_lineage_sha))

lineage_source_rows <- fread(lineage_tsv)
stopifnot(identical(names(lineage_source_rows), c("source_sheet",
  "source_row", "lineage", "feature_id")),
  all(lineage_source_rows$source_sheet == "Table S1B"),
  identical(lineage_source_rows[, .N, by = lineage][order(lineage)]$N,
    c(254L, 234L)))
lineage_source_rows[, duplicate_within_lineage := duplicated(feature_id) |
  duplicated(feature_id, fromLast = TRUE), by = lineage]

parse_signal_completion <- read_json(parse_signal_completion_path,
  simplifyVector = TRUE)
stopifnot(identical(parse_signal_completion$status, "complete"),
  identical(parse_signal_completion$promoted, FALSE),
  identical(as.integer(parse_signal_completion$passing_signatures), 38L))
parse_signal_membership_all <- fread(parse_signal_membership_path)
bb_boundaries <- fread(bb_boundary_path)
stopifnot(nrow(bb_boundaries) == 2L,
  all(is.finite(bb_boundaries$pt_median)))
bb_boundaries[, boundary_label := paste0("BB", seq_len(.N))]
parse_signal_membership <- parse_signal_membership_all[
  cytokine %in% c("Flt3L", "IL-7") & direction %in% c("up", "down")]
setorder(parse_signal_membership, cytokine, direction, final_rank)
stopifnot(uniqueN(parse_signal_membership$signature_id) == 4L,
  nrow(parse_signal_membership) == 105L,
  identical(parse_signal_membership[, .N,
    by = .(cytokine, direction)][order(cytokine, direction)]$N,
    c(20L, 45L, 20L, 20L)))

distribution <- fread(distribution_path)
projection <- fread(cmd = paste("zcat", shQuote(projection_path)))
stopifnot(all(c("target", "genotype", "bin", "bin_left", "bin_right",
  "bdev_normalized_rolling_fraction") %in% names(distribution)),
  all(c("source_cell_id", "sample_id", "genotype", "knockout_target",
    "support_status", "mapped_pt") %in% names(projection)))

# Panel F defines enriched regions from WT alone. A region is a maximal run of
# adjacent bins whose Bdev-normalized rolling cell fraction is strictly > 1.
wt_above <- distribution[genotype == "WT" &
  is.finite(bdev_normalized_rolling_fraction) &
  bdev_normalized_rolling_fraction > 1]
setorder(wt_above, target, bin)
wt_above[, run := cumsum(seq_len(.N) == 1L |
  target != shift(target, fill = target[[1L]]) |
  bin != shift(bin, fill = bin[[1L]] - 1L) + 1L)]
regions <- wt_above[, .(
  first_bin = min(bin), last_bin = max(bin),
  pt_left = min(bin_left), pt_right = max(bin_right),
  n_bins = .N,
  maximum_wt_bdev_normalized_fraction = max(bdev_normalized_rolling_fraction)
), by = .(target, run)]
setorder(regions, target, pt_left)
regions[, peak_number := seq_len(.N), by = target]
regions[, region_id := paste0(target, "_peak", peak_number)]
regions[, cytokine := fifelse(target == "Flt3l", "Flt3L", "IL-7")]
setcolorder(regions, c("region_id", "target", "cytokine", "peak_number",
  "first_bin", "last_bin", "pt_left", "pt_right", "n_bins",
  "maximum_wt_bdev_normalized_fraction", "run"))
stopifnot(identical(regions$target, c("Flt3l", "Il7", "Il7")),
  identical(regions$peak_number, c(1L, 1L, 2L)))
write_tsv(regions[, !"run"], file.path(stage_dir,
  "wt_enriched_pt_regions.tsv"))

eligible <- projection[support_status %in% c("supported", "borderline") &
  is.finite(mapped_pt)]
assignments <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
  reg <- regions[i]
  eligible[knockout_target == reg$target & mapped_pt >= reg$pt_left &
    mapped_pt < reg$pt_right,
    .(cell_id, source_cell_id, sample_id, genotype, knockout_target,
      support_status, mapped_pt,
      region_id = reg$region_id, cytokine = reg$cytokine,
      peak_number = reg$peak_number, pt_left = reg$pt_left,
      pt_right = reg$pt_right)]
}))
stopifnot(!anyDuplicated(assignments[, .(source_cell_id, knockout_target)]),
  setequal(unique(assignments$genotype), c("WT", "KO")))
setorder(assignments, knockout_target, peak_number, genotype, mapped_pt)
fwrite(assignments, file.path(stage_dir, "cell_region_assignments.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip")

region_counts <- assignments[, .(n_cells = .N,
  n_samples = uniqueN(sample_id), samples = paste(sort(unique(sample_id)),
    collapse = ";")),
  by = .(region_id, knockout_target, cytokine, peak_number, genotype,
    pt_left, pt_right)]
setorder(region_counts, knockout_target, peak_number, genotype)
stopifnot(all(region_counts$n_samples == 1L), nrow(region_counts) == 6L,
  all(region_counts$n_cells >= 100L))
write_tsv(region_counts, file.path(stage_dir, "region_cell_counts.tsv"))

# Visual audit of the exact pseudotime regions used by the three DE contrasts.
region_plot_data <- copy(distribution[target %in% c("Flt3l", "Il7")])
region_plot_data[, target_label := factor(target,
  levels = c("Flt3l", "Il7"), labels = c("Flt3L", "IL-7"))]
region_plot_regions <- copy(regions)
region_plot_regions[, target_label := factor(target,
  levels = c("Flt3l", "Il7"), labels = c("Flt3L", "IL-7"))]
region_plot_regions[, region_label := fifelse(target == "Flt3l",
  "Flt3L peak 1", paste0("IL-7 peak ", peak_number))]
region_fill_values <- c("Flt3l_peak1" = "#E69F00",
  "Il7_peak1" = "#009E73", "Il7_peak2" = "#CC79A7")
enrichment_region_plot <- ggplot(region_plot_data,
    aes(bin_midpoint, bdev_normalized_rolling_fraction,
      colour = genotype, group = genotype)) +
  geom_rect(data = region_plot_regions,
    aes(xmin = pt_left, xmax = pt_right, ymin = -Inf, ymax = Inf,
      fill = region_id), inherit.aes = FALSE, alpha = 0.18,
    colour = NA) +
  geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.55,
    colour = "grey25") +
  geom_vline(data = bb_boundaries, aes(xintercept = pt_median),
    inherit.aes = FALSE, linetype = "dashed", linewidth = 0.45,
    colour = "grey35") +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_text(data = region_plot_regions,
    aes(x = (pt_left + pt_right) / 2, y = Inf, label = region_label),
    inherit.aes = FALSE, vjust = 1.3, size = 3.2, fontface = "bold") +
  geom_text(data = bb_boundaries,
    aes(x = pt_median, y = -Inf, label = boundary_label),
    inherit.aes = FALSE, vjust = -0.5, hjust = -0.08, size = 3,
    colour = "grey20") +
  facet_wrap(~target_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("WT" = "#C84545", "KO" = "#3B77A8")) +
  scale_fill_manual(values = region_fill_values) +
  scale_x_continuous(limits = range(c(region_plot_data$bin_left,
    region_plot_data$bin_right)), expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Mapped Bdev pseudotime",
    y = "Bdev-normalized rolling cell fraction", colour = "Genotype",
    fill = "DE region",
    title = "Pseudotime regions used for WT-versus-KO differential expression",
    subtitle = "Shading marks maximal contiguous WT bins above the ratio-1 line; identical bounds were applied to WT and KO cells") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold",
    size = 11), panel.grid.minor = element_blank())
ggsave(file.path(stage_dir, "figures",
  "wt_enriched_de_pseudotime_regions.pdf"), enrichment_region_plot,
  width = 11, height = 5.2, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures",
  "wt_enriched_de_pseudotime_regions.png"), enrichment_region_plot,
  width = 11, height = 5.2, dpi = 300)

bundle <- wb_open_bundle(ko_bundle_dir, validate = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
fingerprint_before <- wb_manifest(bundle)$current_state$fingerprint
stopifnot(identical(fingerprint_before, expected_ko_fp))
downstream_cells <- as.character(wb_cells(bundle, "set:qc_use_downstream"))
stopifnot(all(assignments$source_cell_id %in% downstream_cells))
metadata <- as.data.table(wb_meta(bundle, cells = downstream_cells))
stopifnot(identical(as.character(metadata$cell_id), downstream_cells))
metadata_key <- metadata[, .(source_cell_id = cell_id, library_id,
  genotype, knockout_target)]
alignment <- merge(assignments[, .(source_cell_id, sample_id, genotype,
  knockout_target)], metadata_key, by = "source_cell_id", all.x = TRUE,
  suffixes = c("_projection", "_bundle"))
stopifnot(nrow(alignment) == nrow(assignments),
  !anyNA(alignment$library_id),
  identical(as.character(alignment$sample_id),
    as.character(alignment$library_id)),
  identical(as.character(alignment$genotype_projection),
    as.character(alignment$genotype_bundle)),
  identical(as.character(alignment$knockout_target_projection),
    as.character(alignment$knockout_target_bundle)))

session_id <- "ko_wt_enriched_pt_peak_de_v1"
workspace <- file.path(session_root, session_id)
session <- if (dir.exists(workspace)) wb_session_resume(bundle, workspace) else
  wb_session_create(bundle, workspace_root = session_root,
    session_id = session_id)
if (!is.null(session$durable$promotion_receipt))
  stop("The report-only DE session was unexpectedly promoted")
find_run <- function(session, label, op) {
  marker <- paste0("[", label, "] ")
  ids <- names(session$runs)[vapply(session$runs, function(x) {
    note <- if (is.null(x$notes)) "" else as.character(x$notes)
    identical(x$op_name, op) && startsWith(note, marker)
  }, logical(1))]
  stopifnot(length(ids) <= 1L)
  if (length(ids)) session$runs[[ids[[1L]]]] else NULL
}

write_json(wb_op_interface("matrix_log_normalize", bundle = bundle),
  file.path(stage_dir, "matrix_log_normalize_interface.json"))
write_json(wb_op_interface("de_wilcox", bundle = bundle),
  file.path(stage_dir, "de_wilcox_interface.json"))

normalized_matrix_id <- "RNA_log1p_cp10k_wt_enriched_peak_de_v1"
normal_run <- find_run(session, "normalize_downstream_rna", "matrix_log_normalize")
if (is.null(normal_run)) {
  log_event("Staging full-bundle log1p-CP10K RNA matrix")
  result <- wb_session_run_op(session, bundle, "matrix_log_normalize",
    inputs = list(matrix = "RNA_counts", cell_set = "all",
      min_cells = 10L, output_matrix_id = normalized_matrix_id,
      scale_factor = 10000, overwrite = FALSE),
    notes = paste0("[normalize_downstream_rna] Full-bundle ",
      "log1p(count/library size*10000) matrix (operation requires cell_set=all); ",
      "DE groups remain restricted to marked downstream cells; staged only, ",
      "no promotion."))
  session <- result$session
  normal_run <- result$run
} else log_event("Reusing staged normalized RNA matrix")
stopifnot("matrix" %in% names(normal_run$outputs))

de_tables <- list()
for (i in seq_len(nrow(regions))) {
  reg <- regions[i]
  label <- paste0("de_", tolower(reg$region_id))
  de_run <- find_run(session, label, "de_wilcox")
  region_cells <- assignments[region_id == reg$region_id]
  wt_cells <- as.character(region_cells[genotype == "WT", source_cell_id])
  ko_cells <- as.character(region_cells[genotype == "KO", source_cell_id])
  stopifnot(length(wt_cells) >= 100L, length(ko_cells) >= 100L,
    !length(intersect(wt_cells, ko_cells)))
  if (is.null(de_run)) {
    log_event("Running cell-level Wilcoxon DE for", reg$region_id,
      "WT", length(wt_cells), "KO", length(ko_cells))
    result <- wb_session_run_op(session, bundle, "de_wilcox",
      inputs = list(matrix = normalized_matrix_id,
        cell_set = "set:qc_use_downstream", min_cells = 10L,
        group_a = wt_cells, group_b = ko_cells, min_pct = 0.01,
        p_adjust_method = "BH"),
      notes = paste0("[", label, "] Group A=WT; group B=KO; region ",
        reg$region_id, " defined only from WT Bdev-normalized rolling ",
        "fraction >1; descriptive cell-level test; no biological replicates."))
    session <- result$session
    de_run <- result$run
  } else log_event("Reusing completed DE run for", reg$region_id)
  de <- as.data.table(de_run$outputs$de_table$table)
  stopifnot(nrow(de) > 1000L,
    all(c("feature_id", "pct_a", "pct_b", "avg_expr_a", "avg_expr_b",
      "avg_log1p_fc", "p_val", "p_val_adj") %in% names(de)))
  de[, `:=`(region_id = reg$region_id, target = reg$target,
    cytokine = reg$cytokine, peak_number = reg$peak_number,
    pt_left = reg$pt_left, pt_right = reg$pt_right,
    n_wt_cells = length(wt_cells), n_ko_cells = length(ko_cells),
    direction = fifelse(avg_log1p_fc > 0, "WT_higher",
      fifelse(avg_log1p_fc < 0, "KO_higher", "equal")),
    significant_cell_level_bh05 = is.finite(p_val_adj) & p_val_adj < 0.05)]
  setcolorder(de, c("region_id", "target", "cytokine", "peak_number",
    "pt_left", "pt_right", "n_wt_cells", "n_ko_cells", "feature_id",
    "direction", "significant_cell_level_bh05", "pct_a", "pct_b",
    "avg_expr_a", "avg_expr_b", "avg_log1p_fc", "p_val", "p_val_adj"))
  write_tsv(de, file.path(stage_dir, paste0(reg$region_id,
    "_wt_vs_ko_de.tsv")))
  de_tables[[reg$region_id]] <- de
}
combined_de <- rbindlist(de_tables, use.names = TRUE)
write_tsv(combined_de, file.path(stage_dir, "all_regions_wt_vs_ko_de.tsv"))

ranked <- combined_de[order(region_id, p_val_adj, -abs(avg_log1p_fc))]
top_genes <- ranked[direction != "equal", head(.SD, 20L),
  by = .(region_id, direction)]
write_tsv(top_genes, file.path(stage_dir, "top20_genes_per_region_direction.tsv"))

# Violin selection is performed independently for each direction within each
# region: significant genes with absolute effect >0.2, capped at the 10
# smallest BH adjusted P values per direction when necessary.
violin_candidates <- combined_de[is.finite(p_val_adj) & p_val_adj < 0.05 &
  is.finite(avg_log1p_fc) & abs(avg_log1p_fc) > 0.2 &
  direction %in% c("WT_higher", "KO_higher")]
violin_candidates[, absolute_effect := abs(avg_log1p_fc)]
setorder(violin_candidates, region_id, direction, p_val_adj,
  -absolute_effect, feature_id)
violin_genes <- violin_candidates[, head(.SD, 10L),
  by = .(region_id, direction)]
violin_genes[, selection_rank := seq_len(.N),
  by = .(region_id, direction)]
violin_totals <- violin_candidates[, .(
  qualifying_genes_in_direction = .N), by = .(region_id, direction)]
violin_genes <- merge(violin_genes, violin_totals,
  by = c("region_id", "direction"), all.x = TRUE)
setorder(violin_genes, region_id, direction, selection_rank)
stopifnot(all(violin_genes$qualifying_genes_in_direction >=
  violin_genes$selection_rank), all(violin_genes$selection_rank <= 10L))
write_tsv(violin_genes, file.path(stage_dir,
  "violin_gene_selection_abs_fc_gt_0.2.tsv"))

staged_bundle <- wb_session_view(session)
on.exit(try(wb_close_bundle(staged_bundle), silent = TRUE), add = TRUE)

# scWorkbench does not expose per-cell arithmetic-mean gene-set scores. With
# explicit user approval, calculate them reproducibly from the staged native
# log1p-CP10K matrix returned by scWorkbench.
lineage_membership <- unique(lineage_source_rows[, .(lineage, feature_id)])
normalized_features <- rownames(wb_matrix(staged_bundle,
  normalized_matrix_id))
lineage_membership[, mapped_to_matrix := feature_id %in% normalized_features]
lineage_coverage <- lineage_membership[, .(
  unique_source_genes = .N,
  mapped_genes = sum(mapped_to_matrix),
  missing_genes = sum(!mapped_to_matrix),
  mapped_fraction = mean(mapped_to_matrix),
  missing_feature_ids = paste(feature_id[!mapped_to_matrix], collapse = ";")
), by = lineage]
setorder(lineage_coverage, lineage)
stopifnot(identical(lineage_coverage$unique_source_genes, c(251L, 234L)),
  identical(lineage_coverage$mapped_genes, c(246L, 233L)))
lineage_source_rows <- merge(lineage_source_rows,
  lineage_membership, by = c("lineage", "feature_id"), all.x = TRUE,
  sort = FALSE)
setorder(lineage_source_rows, source_row, lineage)
write_tsv(lineage_source_rows, file.path(stage_dir,
  "lineage_gene_source_rows_and_mapping.tsv"))
write_tsv(lineage_coverage, file.path(stage_dir,
  "lineage_gene_coverage.tsv"))

lineage_scores <- rbindlist(lapply(c("B lineage", "T lineage"),
  function(lineage_value) {
    genes <- lineage_membership[lineage == lineage_value &
      mapped_to_matrix, feature_id]
    expression <- as.matrix(wb_matrix_read_features(staged_bundle,
      normalized_matrix_id, genes, downstream_cells, missing = "NA"))
    stopifnot(identical(rownames(expression), genes),
      identical(colnames(expression), downstream_cells), !anyNA(expression))
    data.table(cell_id = downstream_cells, lineage = lineage_value,
      aggregate_score = colMeans(expression), mapped_genes = length(genes))
  }))
lineage_scores <- merge(lineage_scores,
  metadata[, .(cell_id, sample_id = library_id, genotype, knockout_target)],
  by = "cell_id", all.x = TRUE, sort = FALSE)
sample_order <- c("FLT3LKO_BM_WT", "FLT3LKO_BM_KO", "IL7KO_BM_WT",
  "IL7KO_BM_KO")
stopifnot(setequal(unique(lineage_scores$sample_id), sample_order),
  nrow(lineage_scores) == length(downstream_cells) * 2L)
sample_labels <- c("FLT3LKO_BM_WT" = "Flt3L WT",
  "FLT3LKO_BM_KO" = "Flt3L KO", "IL7KO_BM_WT" = "IL-7 WT",
  "IL7KO_BM_KO" = "IL-7 KO")
lineage_scores[, sample_id := factor(sample_id, levels = sample_order)]
lineage_scores[, lineage := factor(lineage,
  levels = c("B lineage", "T lineage"))]
lineage_score_summary <- lineage_scores[, .(
  n_cells = .N, mean_score = mean(aggregate_score),
  median_score = median(aggregate_score),
  q25 = quantile(aggregate_score, 0.25),
  q75 = quantile(aggregate_score, 0.75)),
  by = .(lineage, sample_id, genotype, knockout_target, mapped_genes)]
setorder(lineage_score_summary, lineage, sample_id)
fwrite(lineage_scores, file.path(stage_dir,
  "lineage_aggregate_cell_scores.tsv.gz"), sep = "\t", na = "NA",
  compress = "gzip")
write_tsv(lineage_score_summary, file.path(stage_dir,
  "lineage_aggregate_score_summary.tsv"))

lineage_plot <- ggplot(lineage_scores, aes(sample_id, aggregate_score,
    fill = sample_id)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.3,
    colour = "grey25") +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.3,
    fill = "white", alpha = 0.82) +
  facet_wrap(~lineage, nrow = 1, scales = "free_y") +
  scale_x_discrete(labels = sample_labels) +
  scale_fill_manual(values = c("FLT3LKO_BM_WT" = "#D95F5F",
    "FLT3LKO_BM_KO" = "#4C78A8", "IL7KO_BM_WT" = "#F2A65A",
    "IL7KO_BM_KO" = "#6F4C9B"), labels = sample_labels) +
  labs(x = NULL, y = "Mean log1p CP10K expression",
    fill = "Sample",
    title = "B- and T-lineage aggregate expression across all four samples",
    subtitle = "Arithmetic mean across exactly mapped Table S1B genes: B=246, T=233") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold",
    size = 11), panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"))
ggsave(file.path(stage_dir, "figures",
  "b_t_lineage_aggregate_score_all_four_samples.pdf"), lineage_plot,
  width = 11, height = 5.2, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures",
  "b_t_lineage_aggregate_score_all_four_samples.png"), lineage_plot,
  width = 11, height = 5.2, dpi = 300)

# Apply the same explicitly approved arithmetic-mean score to the four accepted
# Parse-derived Flt3L/IL-7 UP and DOWN signatures.
parse_signal_membership[, mapped_to_matrix :=
  feature_id %in% normalized_features]
parse_signal_coverage <- parse_signal_membership[, .(
  source_genes = .N, mapped_genes = sum(mapped_to_matrix),
  missing_genes = sum(!mapped_to_matrix),
  mapped_fraction = mean(mapped_to_matrix),
  missing_feature_ids = paste(feature_id[!mapped_to_matrix], collapse = ";")),
  by = .(signature_id, cytokine, direction)]
setorder(parse_signal_coverage, cytokine, direction)
stopifnot(all(parse_signal_coverage$source_genes ==
    parse_signal_coverage$mapped_genes),
  identical(parse_signal_coverage$mapped_genes, c(20L, 45L, 20L, 20L)))
write_tsv(parse_signal_membership, file.path(stage_dir,
  "flt3l_il7_parse_signal_gene_membership_and_mapping.tsv"))
write_tsv(parse_signal_coverage, file.path(stage_dir,
  "flt3l_il7_parse_signal_gene_coverage.tsv"))

parse_signal_scores <- rbindlist(lapply(seq_len(nrow(parse_signal_coverage)),
  function(i) {
    spec <- parse_signal_coverage[i]
    genes <- parse_signal_membership[signature_id == spec$signature_id &
      mapped_to_matrix, feature_id]
    expression <- as.matrix(wb_matrix_read_features(staged_bundle,
      normalized_matrix_id, genes, downstream_cells, missing = "NA"))
    stopifnot(identical(rownames(expression), genes),
      identical(colnames(expression), downstream_cells), !anyNA(expression))
    data.table(cell_id = downstream_cells,
      signature_id = spec$signature_id, cytokine = spec$cytokine,
      direction = spec$direction, aggregate_score = colMeans(expression),
      mapped_genes = length(genes))
  }))
parse_signal_scores <- merge(parse_signal_scores,
  metadata[, .(cell_id, sample_id = library_id, genotype, knockout_target)],
  by = "cell_id", all.x = TRUE, sort = FALSE)
parse_signal_scores[, panel_label := paste(cytokine, toupper(direction))]
signal_panel_order <- c("Flt3L UP", "Flt3L DOWN", "IL-7 UP", "IL-7 DOWN")
parse_signal_scores[, panel_label := factor(panel_label,
  levels = signal_panel_order)]
parse_signal_scores[, sample_id := factor(sample_id, levels = sample_order)]
stopifnot(nrow(parse_signal_scores) == length(downstream_cells) * 4L,
  setequal(as.character(unique(parse_signal_scores$sample_id)), sample_order),
  !anyNA(parse_signal_scores$panel_label))
parse_signal_summary <- parse_signal_scores[, .(
  n_cells = .N, mean_score = mean(aggregate_score),
  median_score = median(aggregate_score),
  q25 = quantile(aggregate_score, 0.25),
  q75 = quantile(aggregate_score, 0.75)),
  by = .(panel_label, cytokine, direction, sample_id, genotype,
    knockout_target, mapped_genes)]
setorder(parse_signal_summary, panel_label, sample_id)
fwrite(parse_signal_scores, file.path(stage_dir,
  "flt3l_il7_parse_signal_aggregate_cell_scores.tsv.gz"), sep = "\t",
  na = "NA", compress = "gzip")
write_tsv(parse_signal_summary, file.path(stage_dir,
  "flt3l_il7_parse_signal_aggregate_score_summary.tsv"))

parse_signal_plot <- ggplot(parse_signal_scores,
    aes(sample_id, aggregate_score, fill = sample_id)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.3,
    colour = "grey25") +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.3,
    fill = "white", alpha = 0.82) +
  facet_wrap(~panel_label, nrow = 1, scales = "free_y") +
  scale_x_discrete(labels = sample_labels) +
  scale_fill_manual(values = c("FLT3LKO_BM_WT" = "#D95F5F",
    "FLT3LKO_BM_KO" = "#4C78A8", "IL7KO_BM_WT" = "#F2A65A",
    "IL7KO_BM_KO" = "#6F4C9B"), labels = sample_labels) +
  labs(x = NULL, y = "Mean log1p CP10K expression", fill = "Sample",
    title = "Parse-derived Flt3L and IL-7 signal scores across all four samples",
    subtitle = "Arithmetic mean across accepted UP/DOWN signature genes; exact RNA mapping") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold",
    size = 10), panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 32, hjust = 1, face = "bold"))
ggsave(file.path(stage_dir, "figures",
  "flt3l_il7_parse_signal_aggregate_score_all_four_samples.pdf"),
  parse_signal_plot, width = 16, height = 5.3, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures",
  "flt3l_il7_parse_signal_aggregate_score_all_four_samples.png"),
  parse_signal_plot, width = 16, height = 5.3, dpi = 300)

# Sample-specific profiles on the mapped Bdev pseudotime axis. The underlying
# grid is the same 100-bin grid used by panel F; a centered five-bin weighted
# mean is shown only when at least 20 mapped cells contribute to the window.
trajectory_bins <- unique(distribution[, .(bin, bin_left, bin_right,
  bin_midpoint)])
setorder(trajectory_bins, bin)
stopifnot(nrow(trajectory_bins) == 100L,
  identical(trajectory_bins$bin, seq_len(100L)))
trajectory_breaks <- c(trajectory_bins$bin_left[[1L]],
  trajectory_bins$bin_right)
make_score_trajectory <- function(scores, panel_order, output_stem,
    title, subtitle_prefix) {
  mapped <- merge(scores, eligible[, .(source_cell_id, mapped_pt)],
    by.x = "cell_id", by.y = "source_cell_id", all = FALSE, sort = FALSE)
  stopifnot(all(is.finite(mapped$mapped_pt)),
    setequal(as.character(unique(mapped$sample_id)), sample_order))
  mapped[, bin := findInterval(mapped_pt, trajectory_breaks,
    rightmost.closed = TRUE, all.inside = TRUE)]
  binned <- mapped[, .(n_cells = .N, score_sum = sum(aggregate_score),
    bin_mean = mean(aggregate_score)),
    by = .(panel_label = as.character(panel_label),
      sample_id = as.character(sample_id), bin)]
  grid <- CJ(panel_label = panel_order, sample_id = sample_order,
    bin = seq_len(100L), unique = TRUE)
  binned <- merge(grid, binned,
    by = c("panel_label", "sample_id", "bin"), all.x = TRUE, sort = FALSE)
  binned <- merge(binned, trajectory_bins, by = "bin", all.x = TRUE,
    sort = FALSE)
  setorder(binned, panel_label, sample_id, bin)
  binned[is.na(n_cells), `:=`(n_cells = 0L, score_sum = 0)]
  binned[, `:=`(
    rolling_n_cells = frollsum(n_cells, n = 5L, align = "center",
      fill = NA_real_),
    rolling_score_sum = frollsum(score_sum, n = 5L, align = "center",
      fill = NA_real_)), by = .(panel_label, sample_id)]
  binned[, rolling_mean := fifelse(rolling_n_cells >= 20L,
    rolling_score_sum / rolling_n_cells, NA_real_)]
  binned[, panel_label := factor(panel_label, levels = panel_order)]
  binned[, sample_id := factor(sample_id, levels = sample_order)]
  plot <- ggplot(binned,
      aes(bin_midpoint, rolling_mean, colour = sample_id,
        group = sample_id)) +
    geom_vline(data = bb_boundaries, aes(xintercept = pt_median),
      inherit.aes = FALSE, linetype = "dashed", linewidth = 0.45,
      colour = "grey30") +
    geom_line(linewidth = 0.9) +
    geom_text(data = bb_boundaries,
      aes(x = pt_median, y = Inf, label = boundary_label),
      inherit.aes = FALSE, vjust = 1.25, hjust = -0.08, size = 3,
      colour = "grey20") +
    facet_wrap(~panel_label, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = c("FLT3LKO_BM_WT" = "#D95F5F",
      "FLT3LKO_BM_KO" = "#4C78A8", "IL7KO_BM_WT" = "#F2A65A",
      "IL7KO_BM_KO" = "#6F4C9B"), labels = sample_labels) +
    scale_x_continuous(limits = range(trajectory_breaks),
      expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = "Mapped Bdev pseudotime",
      y = "Rolling mean aggregate score", colour = "Sample",
      title = title,
      subtitle = paste0(subtitle_prefix,
        "; centered 5-bin cell-weighted mean on 100 bins; minimum 20 cells")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold",
      size = 10), panel.grid.minor = element_blank())
  ggsave(file.path(stage_dir, "figures", paste0(output_stem, ".pdf")),
    plot, width = 12, height = 5.1, device = cairo_pdf)
  ggsave(file.path(stage_dir, "figures", paste0(output_stem, ".png")),
    plot, width = 12, height = 5.1, dpi = 300)
  binned[]
}

lineage_trajectory_scores <- copy(lineage_scores)
lineage_trajectory_scores[, panel_label := as.character(lineage)]
lineage_trajectory <- make_score_trajectory(lineage_trajectory_scores,
  c("B lineage", "T lineage"),
  "b_t_lineage_aggregate_score_trajectory_all_four_samples",
  "B- and T-lineage aggregate expression along mapped Bdev pseudotime",
  "All four samples; dashed lines mark BB1 and BB2")
write_tsv(lineage_trajectory, file.path(stage_dir,
  "lineage_aggregate_score_trajectory.tsv"))

parse_up_scores <- copy(parse_signal_scores[direction == "up"])
parse_up_scores[, panel_label := factor(as.character(panel_label),
  levels = c("Flt3L UP", "IL-7 UP"))]
parse_up_trajectory <- make_score_trajectory(parse_up_scores,
  c("Flt3L UP", "IL-7 UP"),
  "flt3l_il7_parse_up_signal_aggregate_score_trajectory_all_four_samples",
  "Parse-derived cytokine UP signals along mapped Bdev pseudotime",
  "All four samples; dashed lines mark BB1 and BB2")
write_tsv(parse_up_trajectory, file.path(stage_dir,
  "flt3l_il7_parse_up_signal_aggregate_score_trajectory.tsv"))

violin_files <- data.table(region_id = character(), direction = character(),
  file = character())
for (i in seq_len(nrow(regions))) {
  reg <- regions[i]
  region_cells <- assignments[region_id == reg$region_id]
  cell_order <- as.character(region_cells$source_cell_id)
  for (direction_value in c("WT_higher", "KO_higher")) {
    selected <- violin_genes[region_id == reg$region_id &
      direction == direction_value]
    stopifnot(nrow(selected) >= 1L, nrow(selected) <= 10L)
    genes <- as.character(selected$feature_id)
    expression <- as.matrix(wb_matrix_read_features(staged_bundle,
      normalized_matrix_id, genes, cell_order, missing = "NA"))
    stopifnot(identical(rownames(expression), genes),
      identical(colnames(expression), cell_order), !anyNA(expression))
    expression_long <- as.data.table(t(expression),
      keep.rownames = "source_cell_id")
    expression_long <- melt(expression_long, id.vars = "source_cell_id",
      variable.name = "feature_id", value.name = "log1p_cp10k")
    expression_long <- merge(expression_long,
      region_cells[, .(source_cell_id, genotype)], by = "source_cell_id",
      all.x = TRUE, sort = FALSE)
    expression_long <- merge(expression_long,
      selected[, .(feature_id, avg_log1p_fc, p_val_adj, selection_rank)],
      by = "feature_id", all.x = TRUE, sort = FALSE)
    expression_long[, genotype := factor(genotype, levels = c("WT", "KO"))]
    expression_long[, facet_label := sprintf("%02d. %s\neffect %+.3f",
      selection_rank, feature_id, avg_log1p_fc)]
    facet_levels <- unique(expression_long[order(selection_rank)]$facet_label)
    expression_long[, facet_label := factor(facet_label,
      levels = facet_levels)]
    direction_label <- if (direction_value == "WT_higher")
      "WT-higher genes" else "KO-higher genes"
    subtitle <- if (selected$qualifying_genes_in_direction[[1L]] > 10L)
      sprintf("Top 10 by BH-adjusted P among %d qualifying %s",
        selected$qualifying_genes_in_direction[[1L]], direction_label) else
      sprintf("All %d qualifying %s",
        selected$qualifying_genes_in_direction[[1L]], direction_label)
    violin_plot <- ggplot(expression_long, aes(genotype, log1p_cp10k,
        fill = genotype)) +
      geom_violin(scale = "width", trim = TRUE, linewidth = 0.25,
        colour = "grey25") +
      geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.25,
        fill = "white", alpha = 0.8) +
      facet_wrap(~facet_label, ncol = 5, scales = "free_y") +
      scale_fill_manual(values = c("WT" = "#C84545", "KO" = "#3B77A8")) +
      labs(x = NULL, y = "log1p CP10K expression",
        title = sprintf("%s peak %d — %s", reg$cytokine,
          reg$peak_number, direction_label), subtitle = subtitle) +
      theme_bw(base_size = 9) +
      theme(legend.position = "none", strip.text = element_text(size = 8,
        face = "bold"), panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold"))
    direction_slug <- if (direction_value == "WT_higher")
      "wt_higher" else "ko_higher"
    stem <- paste0(reg$region_id, "_", direction_slug,
      "_significant_abs_fc_gt_0.2_violins")
    height <- if (nrow(selected) > 5L) 6.2 else 3.8
    ggsave(file.path(stage_dir, "figures", paste0(stem, ".pdf")),
      violin_plot, width = 12, height = height, device = cairo_pdf)
    ggsave(file.path(stage_dir, "figures", paste0(stem, ".png")),
      violin_plot, width = 12, height = height, dpi = 300)
    violin_files <- rbind(violin_files, data.table(region_id = reg$region_id,
      direction = direction_value, file = paste0(stem, ".png")))
  }
}
wb_close_bundle(staged_bundle)

plot_data <- copy(combined_de)
plot_data[, neglog10_padj := -log10(pmax(p_val_adj, 1e-300))]
plot_data[, status := fifelse(p_val_adj < 0.05 & avg_log1p_fc > 0,
  "WT higher", fifelse(p_val_adj < 0.05 & avg_log1p_fc < 0,
    "KO higher", "BH >= 0.05"))]
labels <- plot_data[status != "BH >= 0.05"][order(p_val_adj,
  -abs(avg_log1p_fc)), head(.SD, 10L), by = region_id]
region_labels <- regions[, .(region_id,
  facet_label = sprintf("%s peak %d\npt [%.3f, %.3f)", cytokine,
    peak_number, pt_left, pt_right))]
plot_data <- merge(plot_data, region_labels, by = "region_id")
labels <- merge(labels, region_labels, by = "region_id")
volcano <- ggplot(plot_data, aes(avg_log1p_fc, neglog10_padj,
    colour = status)) +
  geom_point(size = 0.55, alpha = 0.55) +
  geom_text(data = labels, aes(label = feature_id), colour = "black",
    size = 2.5, check_overlap = TRUE, vjust = -0.5) +
  facet_wrap(~facet_label, ncol = 3, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  scale_colour_manual(values = c("KO higher" = "#3B77A8",
    "WT higher" = "#C84545", "BH >= 0.05" = "grey75")) +
  labs(x = "scWorkbench avg_log1p_fc (positive = WT higher)",
    y = expression(-log[10](BH~adjusted~italic(P))), colour = NULL,
    title = "WT versus KO expression within WT-enriched trajectory regions",
    subtitle = "Cell-level Wilcoxon tests on staged log1p-CP10K RNA; descriptive only (one sample per genotype)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(),
    plot.margin = margin(8, 14, 8, 14))
ggsave(file.path(stage_dir, "figures", "region_wt_vs_ko_volcano.pdf"),
  volcano, width = 12, height = 4.7, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "region_wt_vs_ko_volcano.png"),
  volcano, width = 12, height = 4.7, dpi = 300)

dag <- wb_session_dag(session, bundle = bundle, include_artifacts = TRUE,
  include_active = TRUE)
write_json(dag, file.path(stage_dir, "native_session_dag.json"))
if (is.list(dag) && is.data.frame(dag$nodes))
  fwrite(dag$nodes, file.path(stage_dir, "dag_nodes.csv"), na = "NA")
if (is.list(dag) && is.data.frame(dag$edges))
  fwrite(dag$edges, file.path(stage_dir, "dag_edges.csv"), na = "NA")

fingerprint_after <- wb_manifest(bundle)$current_state$fingerprint
checks <- list(
  exact_expected_region_pattern = identical(regions$target,
    c("Flt3l", "Il7", "Il7")),
  regions_derived_from_wt_only = TRUE,
  threshold_strictly_above_one = TRUE,
  enrichment_region_figure_exists = file.exists(file.path(stage_dir,
    "figures", "wt_enriched_de_pseudotime_regions.png")),
  supported_or_borderline_only = all(assignments$support_status %in%
    c("supported", "borderline")),
  finite_mapped_pt_only = all(is.finite(assignments$mapped_pt)),
  both_genotypes_in_every_region = all(region_counts[, uniqueN(genotype),
    by = region_id]$V1 == 2L),
  one_sample_per_genotype_per_region = all(region_counts$n_samples == 1L),
  source_ids_align_to_bundle = nrow(alignment) == nrow(assignments),
  source_metadata_matches_projection = TRUE,
  normalized_matrix_staged = TRUE,
  de_regions_complete = length(de_tables) == 3L,
  violin_thresholds_applied = all(violin_genes$p_val_adj < 0.05 &
    abs(violin_genes$avg_log1p_fc) > 0.2),
  violin_direction_cap_respected = all(violin_genes[, .N,
    by = .(region_id, direction)]$N <= 10L),
  both_violin_directions_present = nrow(violin_files) == 6L,
  violin_figures_complete = all(file.exists(file.path(stage_dir, "figures",
    violin_files$file))),
  lineage_source_exact_sheet_s1b = all(lineage_source_rows$source_sheet ==
    "Table S1B"),
  lineage_expected_unique_gene_counts = identical(
    lineage_coverage$unique_source_genes, c(251L, 234L)),
  lineage_expected_mapped_gene_counts = identical(
    lineage_coverage$mapped_genes, c(246L, 233L)),
  lineage_scores_cover_all_downstream_cells = nrow(lineage_scores) ==
    length(downstream_cells) * 2L,
  lineage_scores_include_all_four_samples =
    setequal(as.character(unique(lineage_scores$sample_id)), sample_order),
  lineage_aggregate_figure_exists = file.exists(file.path(stage_dir,
    "figures", "b_t_lineage_aggregate_score_all_four_samples.png")),
  lineage_trajectory_all_panel_sample_combinations =
    lineage_trajectory[is.finite(rolling_mean), uniqueN(paste(panel_label,
      sample_id))] == 8L,
  lineage_trajectory_figure_exists = file.exists(file.path(stage_dir,
    "figures", "b_t_lineage_aggregate_score_trajectory_all_four_samples.png")),
  parse_signal_four_profiles_exact = uniqueN(
    parse_signal_membership$signature_id) == 4L,
  parse_signal_all_genes_mapped = all(parse_signal_coverage$source_genes ==
    parse_signal_coverage$mapped_genes),
  parse_signal_scores_cover_all_downstream_cells =
    nrow(parse_signal_scores) == length(downstream_cells) * 4L,
  parse_signal_scores_include_all_four_samples =
    setequal(as.character(unique(parse_signal_scores$sample_id)), sample_order),
  parse_signal_aggregate_figure_exists = file.exists(file.path(stage_dir,
    "figures", "flt3l_il7_parse_signal_aggregate_score_all_four_samples.png")),
  parse_up_trajectory_all_panel_sample_combinations =
    parse_up_trajectory[is.finite(rolling_mean), uniqueN(paste(panel_label,
      sample_id))] == 8L,
  parse_up_trajectory_figure_exists = file.exists(file.path(stage_dir,
    "figures",
    "flt3l_il7_parse_up_signal_aggregate_score_trajectory_all_four_samples.png")),
  parent_bundle_fingerprint_unchanged = identical(fingerprint_before,
    fingerprint_after),
  result_promotions = 0L)
stopifnot(all(unlist(checks[names(checks) != "result_promotions"])),
  checks$result_promotions == 0L,
  is.null(session$durable$promotion_receipt))
write_json(list(status = "passed", checks = checks),
  file.path(stage_dir, "checks.json"))

region_lines <- apply(regions, 1L, function(x) sprintf(
  "- `%s`: bins %s-%s, pt `[%.6f, %.6f)`, %s bins.",
  x[["region_id"]], x[["first_bin"]], x[["last_bin"]],
  as.numeric(x[["pt_left"]]), as.numeric(x[["pt_right"]]), x[["n_bins"]]))
count_lines <- apply(region_counts, 1L, function(x) sprintf(
  "- `%s`, %s: %s cells from `%s`.", x[["region_id"]],
  x[["genotype"]], x[["n_cells"]], x[["samples"]]))
report_md <- c(
  "# WT-enriched trajectory-region RNA differential expression", "",
  "## Interpretation warning", "",
  paste("Each cytokine contrast contains **one WT sample and one KO sample**.",
    "The scWorkbench Wilcoxon tests treat cells as observations and therefore",
    "do **not** provide biological-replicate-level inference. Adjusted P values",
    "may be extremely small because of cell number and must be interpreted as",
    "descriptive cell-level evidence, not genotype-level significance."), "",
  "## Region definition", "",
  paste("Regions were derived exclusively from the WT curve in panel F of the",
    "accepted v4 report. A region is a maximal contiguous run of trajectory",
    "bins where the WT Bdev-normalized rolling cell fraction is strictly > 1.",
    "The same pseudotime bounds were then applied to WT and KO query cells.",
    "Only supported or borderline mappings with finite mapped pseudotime were used."), "",
  region_lines, "",
  "![Highlighted DE pseudotime regions](figures/wt_enriched_de_pseudotime_regions.png)", "",
  "The dotted horizontal line is the ratio-1 reference; dashed vertical lines mark BB1 and BB2.", "",
  "## Included cells", "", count_lines, "",
  "## Differential-expression method", "",
  paste("RNA counts were library-size normalized within a non-promoted",
    "scWorkbench session as `log1p(counts / library_size * 10,000)`.",
    "For each region, native `de_wilcox` compared group A (WT) with group B (KO),",
    "using `min_pct=0.01` and Benjamini-Hochberg adjustment.",
    "Positive `avg_log1p_fc` denotes higher expression in WT."), "",
  "![Region-wise volcano plots](figures/region_wt_vs_ko_volcano.png)", "",
  "## B- and T-lineage aggregate expression", "",
  paste("Gene lists were extracted from `Table S1B` of the cited supplementary",
    "workbook. Duplicate source rows are preserved in the source table but",
    "deduplicated for scoring. Exact bundle-feature mapping retained 246 of",
    "251 unique B-lineage genes and 233 of 234 T-lineage genes."), "",
  paste("Because scWorkbench does not export per-cell arithmetic-mean gene-set",
    "scores, this explicitly approved report-only calculation takes the mean",
    "log1p-CP10K expression across mapped genes for every downstream cell."), "",
  "![B- and T-lineage aggregate scores](figures/b_t_lineage_aggregate_score_all_four_samples.png)", "",
  paste("Trajectory profiles use supported and borderline query mappings with",
    "finite mapped pseudotime. Lines are centered five-bin, cell-weighted means",
    "on the same 100-bin axis used by panel F, requiring at least 20 cells per",
    "window. Dashed lines mark BB1 and BB2."), "",
  "![B- and T-lineage trajectory profiles](figures/b_t_lineage_aggregate_score_trajectory_all_four_samples.png)", "",
  "## Parse-derived Flt3L and IL-7 UP/DOWN aggregate expression", "",
  paste("The four accepted Parse-derived profiles contain 45 Flt3L-UP genes",
    "and 20 genes each for Flt3L-DOWN, IL-7-UP, and IL-7-DOWN. All 105",
    "members map exactly to the KO RNA matrix."), "",
  paste("The same approved calculation takes each cell's arithmetic mean",
    "log1p-CP10K expression across the mapped genes. All downstream cells",
    "from all four samples are shown in every panel."), "",
  "![Parse-derived cytokine aggregate scores](figures/flt3l_il7_parse_signal_aggregate_score_all_four_samples.png)", "",
  paste("The Flt3L-UP and IL-7-UP scores are also shown along mapped Bdev",
    "pseudotime for all four samples using the identical rolling profile",
    "definition and BB boundaries."), "",
  "![Parse-derived UP-signal trajectory profiles](figures/flt3l_il7_parse_up_signal_aggregate_score_trajectory_all_four_samples.png)", "",
  "## Expression violins for genes with |effect| > 0.2", "",
  paste("Within each region, genes require BH-adjusted P < 0.05 and",
    "absolute `avg_log1p_fc` > 0.2. WT-higher and KO-higher genes are",
    "selected and displayed separately, with up to 10 genes per direction",
    "ranked by adjusted P value."), "",
  unlist(lapply(seq_len(nrow(regions)), function(i) {
    reg_id <- regions$region_id[[i]]
    wt_file <- violin_files[region_id == reg_id &
      direction == "WT_higher", file]
    ko_file <- violin_files[region_id == reg_id &
      direction == "KO_higher", file]
    c(paste0("### ", regions$cytokine[[i]], " peak ",
        regions$peak_number[[i]]), "", "#### WT-higher genes", "",
      paste0("![", reg_id, " WT-higher violins](figures/", wt_file, ")"),
      "", "#### KO-higher genes", "",
      paste0("![", reg_id, " KO-higher violins](figures/", ko_file, ")"),
      "")
  })),
  "## Outputs", "",
  "- `wt_enriched_pt_regions.tsv`: exact bin and pseudotime definitions.",
  "- `figures/wt_enriched_de_pseudotime_regions.*`: highlighted DE intervals.",
  "- `cell_region_assignments.tsv.gz`: marked WT and KO cells.",
  "- `region_cell_counts.tsv`: cells and samples per contrast.",
  "- `*_wt_vs_ko_de.tsv`: complete scWorkbench result for each region.",
  "- `all_regions_wt_vs_ko_de.tsv`: combined DE results.",
  "- `top20_genes_per_region_direction.tsv`: convenient ranked summary.",
  "- `violin_gene_selection_abs_fc_gt_0.2.tsv`: exact violin selection.",
  "- `lineage_gene_source_rows_and_mapping.tsv`: exact Table S1B rows and mapping.",
  "- `lineage_gene_coverage.tsv`: unique-gene coverage by lineage.",
  "- `lineage_aggregate_cell_scores.tsv.gz`: per-cell aggregate scores.",
  "- `lineage_aggregate_score_summary.tsv`: four-sample summaries.",
  "- `lineage_aggregate_score_trajectory.tsv`: B/T trajectory profiles.",
  "- `flt3l_il7_parse_signal_gene_membership_and_mapping.tsv`: exact four-profile memberships.",
  "- `flt3l_il7_parse_signal_gene_coverage.tsv`: RNA mapping coverage.",
  "- `flt3l_il7_parse_signal_aggregate_cell_scores.tsv.gz`: per-cell scores.",
  "- `flt3l_il7_parse_signal_aggregate_score_summary.tsv`: four-sample summaries.",
  "- `flt3l_il7_parse_up_signal_aggregate_score_trajectory.tsv`: UP-signal trajectory profiles.",
  "- `native_session_dag.json`: unpromoted scWorkbench analysis DAG.", "",
  "## Provenance", "",
  paste0("- Source bundle: `", ko_bundle_dir, "`"),
  paste0("- Source bundle fingerprint: `", fingerprint_before, "`"),
  paste0("- Distribution report: `", distribution_report, "`"),
  paste0("- Projection report: `", projection_report, "`"),
  paste0("- Accepted Parse-derived signal report: `", parse_signal_report, "`"),
  paste0("- Durable session: `", workspace, "`"),
  paste0("- Lineage workbook: `", lineage_workbook, "`"),
  paste0("- Lineage workbook SHA-256: `", expected_lineage_sha, "`"),
  "- Source URL: <https://ars.els-cdn.com/content/image/1-s2.0-S2589004225011162-mmc2.xlsx>",
  paste0("- Script version: `", script_version, "`"),
  "- Promotions: none.")
writeLines(report_md, file.path(stage_dir, "REPORT.md"))

report_html <- c("<!doctype html><html><head><meta charset='utf-8'>",
  "<title>WT-enriched trajectory-region DE</title>",
  "<style>body{font-family:Arial,sans-serif;max-width:1200px;margin:30px auto;line-height:1.45;color:#222}code{background:#f1f1f1;padding:2px 4px}.warning{background:#fff3cd;border-left:5px solid #d79b00;padding:12px}table{border-collapse:collapse}th,td{border:1px solid #ccc;padding:5px 8px}img{max-width:100%}</style></head><body>",
  "<h1>WT-enriched trajectory-region RNA differential expression</h1>",
  paste0("<div class='warning'><b>Interpretation warning:</b> Each cytokine ",
    "contrast has one WT sample and one KO sample. These are descriptive ",
    "cell-level Wilcoxon results, not biological-replicate-level inference.</div>"),
  "<h2>Region definition</h2><p>Maximal contiguous WT bins with Bdev-normalized rolling cell fraction &gt; 1; identical bounds applied to WT and KO. Supported and borderline mappings only.</p>",
  "<ul>", paste0("<li>", html_escape(region_lines), "</li>"), "</ul>",
  "<img src='figures/wt_enriched_de_pseudotime_regions.png'>",
  "<p>The dotted horizontal line marks ratio 1. Dashed vertical lines mark BB1 and BB2. Shaded intervals are the exact regions used for DE.</p>",
  "<h2>Included cells</h2><ul>", paste0("<li>", html_escape(count_lines), "</li>"), "</ul>",
  "<h2>Results</h2><img src='figures/region_wt_vs_ko_volcano.png'>",
  "<p>Positive <code>avg_log1p_fc</code> denotes higher expression in WT. Complete and combined result tables are included beside this report.</p>",
  "<h2>B- and T-lineage aggregate expression</h2>",
  "<p>Gene lists were extracted from <code>Table S1B</code> of the <a href='https://ars.els-cdn.com/content/image/1-s2.0-S2589004225011162-mmc2.xlsx'>cited supplementary workbook</a>. Duplicate source rows were deduplicated for scoring. Exact mapping retained 246/251 unique B-lineage genes and 233/234 T-lineage genes.</p>",
  "<p>With explicit user approval, per-cell scores were calculated outside a native terminal operation as the arithmetic mean of staged log1p-CP10K expression across mapped genes. All downstream cells from all four samples are shown.</p>",
  "<img src='figures/b_t_lineage_aggregate_score_all_four_samples.png'>",
  "<p>Trajectory profiles use supported and borderline mappings with finite mapped pseudotime, a centered five-bin cell-weighted mean on the same 100-bin axis as panel F, and at least 20 cells per window. Dashed lines mark BB1 and BB2.</p>",
  "<img src='figures/b_t_lineage_aggregate_score_trajectory_all_four_samples.png'>",
  "<h2>Parse-derived Flt3L and IL-7 UP/DOWN aggregate expression</h2>",
  "<p>The accepted profiles contain 45 Flt3L-UP genes and 20 genes each for Flt3L-DOWN, IL-7-UP, and IL-7-DOWN. All 105 genes map exactly to the KO RNA matrix. Scores use the same approved arithmetic mean of staged log1p-CP10K expression, and all four samples are shown in every panel.</p>",
  "<img src='figures/flt3l_il7_parse_signal_aggregate_score_all_four_samples.png'>",
  "<p>Flt3L-UP and IL-7-UP scores are also shown along mapped Bdev pseudotime for all four samples with the same rolling definition and BB boundaries.</p>",
  "<img src='figures/flt3l_il7_parse_up_signal_aggregate_score_trajectory_all_four_samples.png'>",
  "<h2>Expression violins for significant genes with |effect| &gt; 0.2</h2>",
  "<p>Genes require BH-adjusted P &lt; 0.05. WT-higher and KO-higher genes are selected and displayed separately, with up to 10 genes per direction ranked by adjusted P value.</p>",
  unlist(lapply(seq_len(nrow(regions)), function(i) {
    reg_id <- regions$region_id[[i]]
    wt_file <- violin_files[region_id == reg_id &
      direction == "WT_higher", file]
    ko_file <- violin_files[region_id == reg_id &
      direction == "KO_higher", file]
    c(paste0("<h3>", html_escape(regions$cytokine[[i]]), " peak ",
        regions$peak_number[[i]], "</h3>"),
      "<h4 style='color:#a52f2f'>WT-higher genes</h4>",
      paste0("<img src='figures/", wt_file, "'>"),
      "<h4 style='color:#245f8d'>KO-higher genes</h4>",
      paste0("<img src='figures/", ko_file, "'>"))
  })),
  "<h2>Method</h2><p>Native scWorkbench <code>matrix_log_normalize</code> followed by <code>de_wilcox</code> on each region; min_pct=0.01, BH adjustment. The durable session was not promoted and the accepted parent fingerprint was unchanged.</p>",
  "</body></html>")
writeLines(report_html, file.path(stage_dir, "REPORT.html"))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
extractor_path <- normalizePath("processing/integration/extract_lineage_gene_sets.py",
  mustWork = TRUE)
input_manifest <- data.table(role = c("distribution", "projection",
  "distribution_completion", "projection_completion", "source_manifest",
  "parse_signal_membership", "parse_signal_completion",
  "bb_boundaries",
  "lineage_workbook", "lineage_extracted_tsv", "lineage_extractor",
  "analysis_script"), path = normalizePath(c(distribution_path,
  projection_path, file.path(distribution_report, "completion.json"),
  file.path(projection_report, "completion.json"),
  file.path(ko_bundle_dir, "manifest.json"), parse_signal_membership_path,
  parse_signal_completion_path, bb_boundary_path, lineage_workbook, lineage_tsv,
  extractor_path, script_path), mustWork = TRUE))
input_manifest[, sha256 := vapply(path, sha256_file, character(1))]
write_tsv(input_manifest, file.path(stage_dir, "input_manifest.tsv"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

completion <- list(status = "complete_unpromoted_report_only",
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  report = file.path(output_dir, "REPORT.html"),
  source_bundle = ko_bundle_dir,
  source_bundle_fingerprint = fingerprint_before,
  source_bundle_fingerprint_unchanged = TRUE,
  distribution_report = distribution_report,
  projection_report = projection_report,
  region_ids = regions$region_id,
  region_count = nrow(regions),
  enrichment_region_figure = "figures/wt_enriched_de_pseudotime_regions.png",
  session_workspace = workspace,
  matrix = normalized_matrix_id,
  de_operation = "de_wilcox",
  group_a = "WT", group_b = "KO", min_pct = 0.01,
  p_adjust_method = "BH", biological_replicates_per_group = 1L,
  inference_scope = "descriptive_cell_level_only",
  lineage_aggregate_score = list(
    source_url = "https://ars.els-cdn.com/content/image/1-s2.0-S2589004225011162-mmc2.xlsx",
    source_sheet = "Table S1B", source_sha256 = expected_lineage_sha,
    score = "per-cell arithmetic mean log1p-CP10K expression across exactly mapped unique genes",
    implementation = "explicitly_approved_report_only_calculation_outside_native_terminal_operation",
    cell_set = "set:qc_use_downstream", samples = sample_order,
    unique_source_genes = stats::setNames(lineage_coverage$unique_source_genes,
      lineage_coverage$lineage),
    mapped_genes = stats::setNames(lineage_coverage$mapped_genes,
      lineage_coverage$lineage),
    trajectory = list(pseudotime = "mapped_pt",
      query_support = c("supported", "borderline"), underlying_bins = 100L,
      smoothing = "centered_5_bin_cell_weighted_mean",
      minimum_cells_per_window = 20L,
      boundaries = "BB1_BB2_dashed_at_source_pt_medians")),
  parse_cytokine_aggregate_score = list(
    source_report = parse_signal_report,
    source_membership = parse_signal_membership_path,
    profiles = signal_panel_order,
    score = "per-cell arithmetic mean log1p-CP10K expression across exactly mapped signature genes",
    implementation = "same_explicitly_approved_report_only_calculation_outside_native_terminal_operation",
    cell_set = "set:qc_use_downstream", samples = sample_order,
    mapped_genes = stats::setNames(parse_signal_coverage$mapped_genes,
      paste(parse_signal_coverage$cytokine,
        toupper(parse_signal_coverage$direction))),
    trajectory_profiles = c("Flt3L UP", "IL-7 UP"),
    trajectory = list(pseudotime = "mapped_pt",
      query_support = c("supported", "borderline"), underlying_bins = 100L,
      smoothing = "centered_5_bin_cell_weighted_mean",
      minimum_cells_per_window = 20L,
      boundaries = "BB1_BB2_dashed_at_source_pt_medians")),
  violin_selection = list(p_val_adj_below = 0.05,
    absolute_avg_log1p_fc_strictly_above = 0.2,
    maximum_genes_per_region_and_direction = 10L,
    directions_displayed_separately = TRUE,
    ranking_when_capped = "within_direction_ascending_BH_adjusted_P_then_descending_absolute_effect",
    qualifying_genes_by_region_direction = stats::setNames(
      violin_totals$qualifying_genes_in_direction,
      paste(violin_totals$region_id, violin_totals$direction, sep = ":")),
    plotted_genes_by_region_direction = {
      plotted <- violin_genes[, .N, by = .(region_id, direction)]
      stats::setNames(plotted$N,
        paste(plotted$region_id, plotted$direction, sep = ":"))
    }),
  promotions = 0L, checks = checks)
write_json(completion, file.path(stage_dir, "completion.json"))

if (!file.rename(stage_dir, output_dir))
  stop("Failed to atomically finalize report directory")
log_event("Completed report", output_dir)
