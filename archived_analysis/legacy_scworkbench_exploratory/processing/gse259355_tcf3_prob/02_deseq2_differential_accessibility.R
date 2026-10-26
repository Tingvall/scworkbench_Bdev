#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(
  "Usage: 02_deseq2_differential_accessibility.R FEATURECOUNTS_TSV SAMPLE_MANIFEST OUTPUT_DIR SCRIPT_VERSION",
  call. = FALSE
)

count_path <- normalizePath(args[[1L]], mustWork = TRUE)
manifest_path <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
script_version <- args[[4L]]
if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

manifest <- fread(manifest_path)
stopifnot(nrow(manifest) == 6L, sum(manifest$condition == "control") == 3L,
  sum(manifest$condition == "TCF3_depleted") == 3L,
  !anyDuplicated(manifest$sample_id),
  identical(sort(manifest$treatment_time_h[manifest$condition == "TCF3_depleted"]), c(6L, 6L, 8L)))

fc <- fread(count_path, skip = "Geneid")
stopifnot(all(c("Geneid", "Chr", "Start", "End", "Strand", "Length") %in% names(fc)))
count_cols <- setdiff(names(fc), c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
count_samples <- sub("\\.filtered\\.bam$", "", basename(count_cols))
stopifnot(setequal(count_samples, manifest$sample_id))
count_cols <- count_cols[match(manifest$sample_id, count_samples)]
count_samples <- sub("\\.filtered\\.bam$", "", basename(count_cols))
stopifnot(identical(count_samples, manifest$sample_id))

counts <- as.matrix(fc[, ..count_cols])
storage.mode(counts) <- "integer"
rownames(counts) <- fc$Geneid
colnames(counts) <- count_samples
stopifnot(!anyNA(counts), all(counts >= 0), !anyDuplicated(rownames(counts)))

tested <- rowSums(counts) >= 10L
if (sum(tested) < 100L) stop("Too few consensus peaks pass row-sum >=10: ", sum(tested))
condition <- factor(manifest$condition, levels = c("control", "TCF3_depleted"))
coldata <- data.frame(row.names = manifest$sample_id, condition = condition,
  treatment_time_h = manifest$treatment_time_h)
dds <- DESeqDataSetFromMatrix(countData = counts[tested, , drop = FALSE],
  colData = coldata, design = ~ condition)
dds <- DESeq(dds, quiet = FALSE)
res <- results(dds, contrast = c("condition", "TCF3_depleted", "control"),
  alpha = 0.05, independentFiltering = TRUE)

res_dt <- as.data.table(as.data.frame(res), keep.rownames = "feature_id")
setnames(res_dt, "padj", "padj_bh")
coords <- fc[, .(feature_id = Geneid, chromosome = Chr,
  start_1based = Start, end_1based = End, peak_length = Length)]
res_dt <- merge(coords, res_dt, by = "feature_id", all.y = TRUE, sort = FALSE)
res_dt[, `:=`(
  log2fc_definition = "TCF3_depleted_over_control",
  significant_padj_0p05 = is.finite(padj_bh) & padj_bh < 0.05,
  direction = fifelse(is.finite(padj_bh) & padj_bh < 0.05 & log2FoldChange < 0,
    "TCF3_dependent_lower_after_depletion",
    fifelse(is.finite(padj_bh) & padj_bh < 0.05 & log2FoldChange > 0,
      "higher_after_TCF3_depletion", "not_significant")))]
res_dt[, abs_log2fc_sort := abs(log2FoldChange)]
setorder(res_dt, padj_bh, pvalue, -abs_log2fc_sort, feature_id, na.last = TRUE)
res_dt[, abs_log2fc_sort := NULL]
fwrite(res_dt, file.path(stage_dir, "tcf3_depleted_vs_control_differential_accessibility.tsv"),
  sep = "\t", na = "NA")

norm <- counts(dds, normalized = TRUE)
norm_dt <- data.table(feature_id = rownames(norm), as.data.frame(norm, check.names = FALSE))
fwrite(norm_dt, file.path(stage_dir, "deseq2_normalized_peak_counts.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip")
fwrite(manifest, file.path(stage_dir, "sample_manifest.tsv"), sep = "\t", na = "NA")
fwrite(data.table(sample_id = names(sizeFactors(dds)), size_factor = unname(sizeFactors(dds))),
  file.path(stage_dir, "deseq2_size_factors.tsv"), sep = "\t")

plot_dt <- copy(res_dt[is.finite(log2FoldChange) & is.finite(padj_bh)])
if (!nrow(plot_dt)) stop("No peaks have finite log2FC and adjusted P value")
min_positive <- min(plot_dt$padj_bh[plot_dt$padj_bh > 0], na.rm = TRUE)
if (!is.finite(min_positive)) min_positive <- .Machine$double.xmin
plot_dt[, volcano_y := -log10(pmax(padj_bh, min_positive / 10))]
plot_dt[, plot_class := factor(direction,
  levels = c("not_significant", "TCF3_dependent_lower_after_depletion",
    "higher_after_TCF3_depletion"))]
counts_by_class <- plot_dt[, .N, by = plot_class]
subtitle <- paste(sprintf("%s: %s", counts_by_class$plot_class, counts_by_class$N),
  collapse = " | ")
p <- ggplot(plot_dt, aes(log2FoldChange, volcano_y, colour = plot_class)) +
  geom_point(size = 0.7, alpha = 0.55) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey35") +
  geom_vline(xintercept = c(-1, 1), linetype = 3, colour = "grey55") +
  scale_colour_manual(values = c(
    not_significant = "grey72",
    TCF3_dependent_lower_after_depletion = "#0072B2",
    higher_after_TCF3_depletion = "#D55E00"), drop = FALSE) +
  labs(
    title = "Global chromatin accessibility after acute TCF3/E2A depletion in pro-B cells",
    subtitle = subtitle,
    x = "log2 fold change (TCF3-depleted / control)",
    y = expression(-log[10](BH-adjusted~P)), colour = NULL,
    caption = "GSE259355; three controls vs three depleted samples (two 6 h, one 8 h). DESeq2 design: ~ condition.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))
ggsave(file.path(stage_dir, "figures", "tcf3_depleted_vs_control_volcano.png"),
  p, width = 9, height = 6.5, dpi = 220)
ggsave(file.path(stage_dir, "figures", "tcf3_depleted_vs_control_volcano.pdf"),
  p, width = 9, height = 6.5)

vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca, "percentVar"))
p_pca <- ggplot(pca, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 3) + geom_text(vjust = -0.8, size = 3) +
  scale_colour_manual(values = c(control = "#555555", TCF3_depleted = "#0072B2")) +
  labs(title = "Variance-stabilized consensus-peak counts",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"), colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(stage_dir, "figures", "sample_pca.png"), p_pca,
  width = 7, height = 5.5, dpi = 200)

summary_counts <- list(
  consensus_peaks = nrow(counts),
  tested_peaks_row_sum_at_least_10 = nrow(dds),
  finite_padj = sum(is.finite(res_dt$padj_bh)),
  significant_padj_below_0.05 = sum(res_dt$significant_padj_0p05),
  tcf3_dependent_lower_after_depletion = sum(res_dt$direction == "TCF3_dependent_lower_after_depletion"),
  higher_after_tcf3_depletion = sum(res_dt$direction == "higher_after_TCF3_depletion")
)
write_json(list(
  schema_version = "gse259355_tcf3_prob_da.v1",
  status = "complete",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script_version = script_version,
  source = list(geo = "GSE259355", bioproject = "PRJNA1081365",
    genome = "GRCm38/mm10", cell_type = "bone marrow pro-B"),
  comparison = list(group_a = "TCF3_depleted", group_b = "control",
    log2fc = "log2(TCF3_depleted/control)", design = "~ condition",
    controls = 3, depleted = 3, depleted_times_h = c(6, 6, 8)),
  filtering = list(minimum_total_fragment_count = 10),
  testing = list(method = "DESeq2 negative-binomial GLM and Wald test",
    multiple_testing = "Benjamini-Hochberg", significance = "padj < 0.05"),
  counts = summary_counts,
  input_sha256 = list(featurecounts = digest(count_path, algo = "sha256", file = TRUE),
    sample_manifest = digest(manifest_path, algo = "sha256", file = TRUE))
), file.path(stage_dir, "completion.json"), pretty = TRUE, auto_unbox = TRUE,
  digits = NA, null = "null", na = "null")

report <- c(
  "# GSE259355 pro-B TCF3/E2A differential accessibility", "",
  "This analysis compares three untreated controls with three acute E2A/TCF3-depleted pro-B ATAC-seq samples.",
  "The source perturbation is auxin-inducible acute protein degradation, not a constitutive genetic knockout.",
  "Two depleted samples are labelled 6 h and one is labelled 8 h in GEO; they are combined in the user-approved 3-versus-3 primary comparison.", "",
  "## Model", "",
  "- DESeq2 design: `~ condition`.",
  "- Reported log2 fold change: TCF3-depleted / control.",
  "- Negative log2 fold changes indicate peaks whose accessibility is maintained by TCF3.",
  "- Peaks with fewer than 10 fragments across all six libraries are excluded before fitting.",
  "- Significance is BH-adjusted P < 0.05; no fold-change cutoff is imposed.", "",
  "## Results", "",
  paste0("- Consensus peaks: ", summary_counts$consensus_peaks),
  paste0("- Peaks tested: ", summary_counts$tested_peaks_row_sum_at_least_10),
  paste0("- Peaks with finite adjusted P values: ", summary_counts$finite_padj),
  paste0("- Significant peaks: ", summary_counts$significant_padj_below_0.05),
  paste0("- TCF3-dependent (lower after depletion): ", summary_counts$tcf3_dependent_lower_after_depletion),
  paste0("- Higher after depletion: ", summary_counts$higher_after_tcf3_depletion), "",
  "## Primary outputs", "",
  "- `tcf3_depleted_vs_control_differential_accessibility.tsv`", 
  "- `figures/tcf3_depleted_vs_control_volcano.png` and `.pdf`",
  "- `deseq2_normalized_peak_counts.tsv.gz`",
  "- `figures/sample_pca.png`",
  "- `completion.json`"
)
writeLines(report, file.path(stage_dir, "REPORT.md"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
if (!file.rename(stage_dir, output_dir)) stop("Could not publish output directory")
cat("Completed GSE259355 TCF3 differential accessibility: ", output_dir, "\n", sep = "")
