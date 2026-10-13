#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: 52_flt3l_il7_ko_vs_wt_hallmark_gsea.R",
  "<flt3l_de.tsv> <il7_de.tsv> <output_dir> <seed>"))
flt3l_path <- normalizePath(args[[1L]], mustWork = TRUE)
il7_path <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
seed <- as.integer(args[[4L]])
stopifnot(is.finite(seed))
if (dir.exists(output_dir)) stop("Output exists: ", output_dir)
stage_dir <- paste0(output_dir, ".incomplete")
if (dir.exists(stage_dir)) stop("Incomplete output exists: ", stage_dir)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(fgsea)
  library(msigdbr)
  library(jsonlite)
})

write_json <- function(x, path) jsonlite::write_json(x, path,
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null",
  na = "null")
file_md5 <- function(path) unname(tools::md5sum(path))

hallmark <- as.data.table(msigdbr::msigdbr(db_species = "MM",
  species = "Mus musculus", collection = "MH"))
required_hallmark <- c("gene_symbol", "gs_name", "gs_description",
  "gs_collection", "gs_collection_name", "db_version")
stopifnot(all(required_hallmark %in% names(hallmark)),
  uniqueN(hallmark$gs_name) == 50L,
  unique(hallmark$db_version) == "2026.1.Mm",
  unique(hallmark$gs_collection) == "MH")
hallmark <- unique(hallmark[, ..required_hallmark])
setorder(hallmark, gs_name, gene_symbol)
fwrite(hallmark, file.path(stage_dir,
  "mouse_msigdb_2026.1.Mm_hallmark_memberships.tsv"), sep = "\t")
gmt <- hallmark[, .(description = first(gs_description),
  genes = paste(sort(unique(gene_symbol)), collapse = "\t")), by = gs_name]
writeLines(paste(gmt$gs_name, gmt$description, gmt$genes, sep = "\t"),
  file.path(stage_dir, "mh.all.v2026.1.Mm.symbols.resolved.gmt"))
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)
pathways <- lapply(pathways, unique)

break_primary_ties <- function(d) {
  setorder(d, ko_vs_wt_log2fc, ko_vs_wt_avg_log1p_fc, feature_id)
  unique_scores <- sort(unique(d$ko_vs_wt_log2fc))
  positive_gaps <- diff(unique_scores)
  minimum_gap <- if (length(positive_gaps)) min(positive_gaps) else 1
  d[, rank_score := ko_vs_wt_log2fc]
  d[, tie_size := .N, by = ko_vs_wt_log2fc]
  tie_values <- d[tie_size > 1L, unique(ko_vs_wt_log2fc)]
  for (value in tie_values) {
    idx <- which(d$ko_vs_wt_log2fc == value)
    lower_gap <- value - max(unique_scores[unique_scores < value],
      na.rm = TRUE)
    upper_gap <- min(unique_scores[unique_scores > value],
      na.rm = TRUE) - value
    local_gaps <- c(lower_gap, upper_gap)
    local_gaps <- local_gaps[is.finite(local_gaps) & local_gaps > 0]
    margin <- if (length(local_gaps)) min(local_gaps) / 4 else
      minimum_gap / 4
    d$rank_score[idx] <- value + seq(-margin, margin,
      length.out = length(idx))
  }
  stopifnot(!anyDuplicated(d$rank_score))
  d
}

prepare_ranking <- function(path, comparison_id, comparison_label) {
  d <- fread(path)
  required <- c("feature_id", "mean_cp10k_wt", "mean_cp10k_ko",
    "standard_log2fc", "avg_log1p_fc", "p_val", "p_val_adj",
    "pt_left", "pt_right", "n_wt_cells", "n_ko_cells")
  stopifnot(all(required %in% names(d)), !anyDuplicated(d$feature_id))
  all_gene_count <- nrow(d)
  d <- d[is.finite(mean_cp10k_wt) & is.finite(mean_cp10k_ko) &
    mean_cp10k_wt + mean_cp10k_ko > 0 & is.finite(standard_log2fc)]
  d[, `:=`(ko_vs_wt_log2fc = -standard_log2fc,
    ko_vs_wt_avg_log1p_fc = -avg_log1p_fc,
    comparison_id = comparison_id, comparison = comparison_label)]
  d <- break_primary_ties(d)
  setorder(d, -rank_score)
  d[, rank := seq_len(.N)]
  stopifnot(!anyDuplicated(d$feature_id), !anyDuplicated(d$rank_score),
    all(diff(d$rank_score) < 0))
  list(data = d, all_gene_count = all_gene_count,
    region_id = unique(d$region_id), pt_range = unique(d[, .(pt_left, pt_right)]),
    cell_counts = unique(d[, .(n_wt_cells, n_ko_cells)]))
}

comparisons <- list(
  Flt3L = prepare_ranking(flt3l_path, "Flt3L_KO_vs_WT",
    "FLT3LKO vs FLT3L WT"),
  IL7 = prepare_ranking(il7_path, "IL7_KO_vs_WT",
    "IL7KO vs IL-7 WT"))

run_gsea <- function(prepared, name, seed_offset) {
  d <- prepared$data
  stats <- d$rank_score
  names(stats) <- d$feature_id
  set.seed(seed + seed_offset)
  result <- fgsea::fgseaMultilevel(pathways = pathways, stats = stats,
    minSize = 15L, maxSize = 500L, eps = 0, scoreType = "std",
    nproc = 1L)
  result <- as.data.table(result)
  setnames(result, "pathway", "hallmark_set")
  result[, `:=`(comparison_id = d$comparison_id[[1L]],
    comparison = d$comparison[[1L]],
    enrichment_direction = fifelse(NES > 0, "KO higher", "WT higher"),
    leading_edge_genes = vapply(leadingEdge, paste, collapse = ",",
      FUN.VALUE = character(1)))]
  result[, leadingEdge := NULL]
  setcolorder(result, c("comparison_id", "comparison", "hallmark_set",
    "enrichment_direction", "size", "ES", "NES", "pval", "padj",
    "log2err", "leading_edge_genes"))
  result[, absolute_NES_for_order := abs(NES)]
  setorder(result, padj, -absolute_NES_for_order, hallmark_set)
  result[, absolute_NES_for_order := NULL]
  fwrite(d, file.path(stage_dir, paste0(tolower(name),
    "_ko_vs_wt_full_detected_gene_ranking.tsv.gz")), sep = "\t",
    compress = "gzip", na = "NA")
  fwrite(result, file.path(stage_dir, paste0(tolower(name),
    "_ko_vs_wt_hallmark_gsea.tsv")), sep = "\t", na = "NA")
  list(result = result, stats = stats)
}

runs <- list(Flt3L = run_gsea(comparisons$Flt3L, "flt3l", 1L),
  IL7 = run_gsea(comparisons$IL7, "il7", 2L))
all_results <- rbindlist(lapply(runs, `[[`, "result"), fill = TRUE)
fwrite(all_results, file.path(stage_dir,
  "all_comparisons_hallmark_gsea.tsv"), sep = "\t", na = "NA")

pretty_set <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  tools::toTitleCase(tolower(gsub("_", " ", x)))
}
direction_colours <- c("KO higher" = "#2166AC", "WT higher" = "#B2182B")
make_nes_plot <- function(result, file_stub) {
  plot_data <- copy(result)
  plot_data[, label := pretty_set(hallmark_set)]
  setorder(plot_data, NES)
  plot_data[, label := factor(label, levels = label)]
  plot_data[, significant := padj < 0.05]
  p <- ggplot(plot_data, aes(NES, label, colour = enrichment_direction,
      alpha = significant)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = NES, yend = label), linewidth = 0.55) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = direction_colours) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35),
      labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR ≥ 0.05")) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 6.7),
      legend.position = "bottom") +
    labs(title = result$comparison[[1L]],
      subtitle = "Mouse MSigDB Hallmark preranked GSEA",
      x = "Normalized enrichment score (positive = KO higher)",
      y = NULL, colour = NULL, alpha = NULL)
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".png")),
    p, width = 8.5, height = 11, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".pdf")),
    p, width = 8.5, height = 11, device = cairo_pdf, bg = "white")
}
make_nes_plot(runs$Flt3L$result, "flt3l_ko_vs_wt_hallmark_gsea_nes")
make_nes_plot(runs$IL7$result, "il7_ko_vs_wt_hallmark_gsea_nes")

heat <- copy(all_results)
heat[, label := pretty_set(hallmark_set)]
heat[, comparison_short := factor(comparison_id,
  levels = c("Flt3L_KO_vs_WT", "IL7_KO_vs_WT"),
  labels = c("FLT3LKO / WT", "IL7KO / WT"))]
heat[, max_abs := max(abs(NES)), by = hallmark_set]
heat[, mean_nes := mean(NES), by = hallmark_set]
setorder(heat, max_abs, mean_nes, hallmark_set)
heat[, label := factor(label, levels = unique(label))]
heat[, fdr_mark := fifelse(padj < 0.001, "***",
  fifelse(padj < 0.01, "**", fifelse(padj < 0.05, "*", "")))]
heat_limit <- max(abs(heat$NES))
p_heat <- ggplot(heat, aes(comparison_short, label, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = fdr_mark), size = 2.4) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, limits = c(-heat_limit, heat_limit),
    name = "NES") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6.8),
    axis.title.x = element_blank()) +
  labs(title = "Hallmark GSEA across cytokine-knockout comparisons",
    subtitle = "Positive/blue: KO higher; negative/red: WT higher; * FDR < 0.05, ** < 0.01, *** < 0.001",
    y = NULL)
ggsave(file.path(stage_dir, "figures",
  "flt3l_il7_ko_vs_wt_hallmark_gsea_heatmap.png"), p_heat,
  width = 8, height = 11, dpi = 300, bg = "white")
ggsave(file.path(stage_dir, "figures",
  "flt3l_il7_ko_vs_wt_hallmark_gsea_heatmap.pdf"), p_heat,
  width = 8, height = 11, device = cairo_pdf, bg = "white")

make_curve_figure <- function(run, file_stub) {
  result <- copy(run$result)
  selected <- rbind(result[NES > 0][order(padj, -NES)][seq_len(min(4L, .N))],
    result[NES < 0][order(padj, NES)][seq_len(min(4L, .N))])
  plots <- lapply(seq_len(nrow(selected)), function(i) {
    row <- selected[i]
    fgsea::plotEnrichment(pathways[[row$hallmark_set]], run$stats) +
      theme_minimal(base_size = 8) +
      labs(title = pretty_set(row$hallmark_set),
        subtitle = sprintf("NES %+.2f; FDR %.3g", row$NES, row$padj),
        x = "KO-higher ← ranked genes → WT-higher",
        y = "Running ES") +
      theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7.5))
  })
  figure <- wrap_plots(plots, ncol = 2) +
    plot_annotation(title = paste0(run$result$comparison[[1L]],
      ": leading Hallmark enrichment curves"))
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".png")),
    figure, width = 12, height = 13, dpi = 300, bg = "white")
  ggsave(file.path(stage_dir, "figures", paste0(file_stub, ".pdf")),
    figure, width = 12, height = 13, device = cairo_pdf, bg = "white")
}
make_curve_figure(runs$Flt3L, "flt3l_ko_vs_wt_top_hallmark_enrichment_curves")
make_curve_figure(runs$IL7, "il7_ko_vs_wt_top_hallmark_enrichment_curves")

summary_rows <- rbindlist(lapply(names(runs), function(name) {
  r <- runs[[name]]$result
  data.table(comparison_id = r$comparison_id[[1L]],
    ranked_detected_genes = length(runs[[name]]$stats),
    hallmark_sets_tested = nrow(r),
    fdr_below_0.05 = sum(r$padj < 0.05),
    ko_higher_fdr_below_0.05 = sum(r$padj < 0.05 & r$NES > 0),
    wt_higher_fdr_below_0.05 = sum(r$padj < 0.05 & r$NES < 0))
}))
fwrite(summary_rows, file.path(stage_dir, "gsea_summary.tsv"), sep = "\t")

diagnostics <- list(
  status = "complete_preliminary_unregistered_outside_scworkbench",
  rationale = paste("scWorkbench pathway_enrich is hypergeometric",
    "over-representation analysis and does not implement preranked GSEA"),
  user_approval_for_outside_analysis = TRUE,
  inputs = list(flt3l_de = flt3l_path, il7_de = il7_path,
    flt3l_md5 = file_md5(flt3l_path), il7_md5 = file_md5(il7_path)),
  ranking = list(primary = "negative stored standard_log2fc = KO versus WT",
    universe = "all genes with nonzero mean CP10K in WT or KO",
    ties = paste("preserve primary ordering; resolve ties by KO-minus-WT",
      "avg_log1p_fc then feature_id within one-quarter of local score gap"),
    positive_direction = "KO higher", negative_direction = "WT higher"),
  gene_sets = list(source = "msigdbr mouse-native MSigDB",
    db_species = "MM", species = "Mus musculus", collection = "MH",
    db_version = unique(hallmark$db_version), sets = uniqueN(hallmark$gs_name),
    memberships = nrow(hallmark)),
  gsea = list(method = "fgseaMultilevel", min_size = 15L,
    max_size = 500L, eps = 0, score_type = "std", nproc = 1L,
    seed = seed, p_adjustment = "BH within each comparison"),
  comparison_summary = summary_rows,
  biological_replicates_per_genotype = 1L,
  inference = paste("descriptive; cells are not biological replicates and",
    "GSEA FDR does not provide replicate-level inference"),
  bundle_mutated = FALSE, capsule_registered = FALSE,
  package_versions = list(R = R.version.string,
    fgsea = as.character(packageVersion("fgsea")),
    msigdbr = as.character(packageVersion("msigdbr"))))
write_json(diagnostics, file.path(stage_dir, "diagnostics.json"))
writeLines(c("# FLT3L and IL-7 KO-versus-WT Hallmark GSEA", "",
  "This approved outside-scWorkbench report performs preranked GSEA because scWorkbench's available `pathway_enrich` operation is hypergeometric over-representation analysis rather than GSEA.",
  "", "## Ranking and direction", "",
  "Each comparison uses every gene with nonzero mean CP10K in WT or KO in the shared-coverage pseudotime region. The primary statistic is KO-versus-WT standard log2 fold-change, obtained by reversing the stored WT-versus-KO sign. Therefore positive NES denotes KO-higher enrichment and negative NES denotes WT-higher enrichment.",
  "", "## Gene sets and algorithm", "",
  paste0("Mouse-native MSigDB ", unique(hallmark$db_version),
    " Hallmark collection (", uniqueN(hallmark$gs_name),
    " sets) was retrieved with msigdbr ", packageVersion("msigdbr"),
    " and frozen in this report. GSEA used fgseaMultilevel ",
    packageVersion("fgsea"), " with minSize=15, maxSize=500, eps=0, seed=",
    seed, ", and one process."),
  "", "## Interpretation limitation", "",
  "There is one biological sample per genotype. Cell-level DE ranks and pathway FDR values are descriptive and do not supply biological-replicate-level inference."),
  file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
if (!file.rename(stage_dir, output_dir)) stop("Could not finalize output")
cat("Completed Hallmark GSEA at", output_dir, "\n")
