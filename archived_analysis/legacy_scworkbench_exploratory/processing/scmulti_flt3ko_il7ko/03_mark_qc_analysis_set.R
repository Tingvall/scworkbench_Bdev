#!/usr/bin/env Rscript

# Mark the downstream QC analysis universe without deleting any cells.
# Usage: Rscript --vanilla 03_mark_qc_analysis_set.R UNFILTERED_QC_DIR OUTPUT_DIR

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 03_mark_qc_analysis_set.R UNFILTERED_QC_DIR OUTPUT_DIR")
qc_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- args[[2L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
completion_path <- file.path(output_dir, "completion.json")
if (file.exists(completion_path)) stop("Refusing to overwrite completed selection output: ", output_dir)

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
qc_path <- file.path(qc_dir, "cell_qc_metrics.tsv.gz")
if (!file.exists(qc_path)) stop("Missing cell-level QC table: ", qc_path)
qc <- fread(cmd = paste("gzip -cd", shQuote(qc_path)))

required <- c(
  "cell_id", "library_id", "knockout_target", "genotype", "tissue",
  "gex_umis_count", "gex_genes_count", "matrix_mito_percent",
  "atac_fragments", "tss_enrichment", "nucleosome_signal"
)
stopifnot(all(required %in% names(qc)), !anyDuplicated(qc$cell_id), nrow(qc) == 28541L)

selection <- qc[, .(
  cell_id, library_id, knockout_target, genotype, tissue,
  rna_umi = gex_umis_count,
  detected_genes = gex_genes_count,
  mitochondrial_percent = matrix_mito_percent,
  atac_fragments,
  tss_enrichment,
  nucleosome_signal
)]
selection[, pass_rna_umi_gt_1000 := is.finite(rna_umi) & rna_umi > 1000]
selection[, pass_genes_gt_500 := is.finite(detected_genes) & detected_genes > 500]
selection[, pass_mito_lt_20 := is.finite(mitochondrial_percent) & mitochondrial_percent < 20]
selection[, pass_atac_fragments_gt_1000 := is.finite(atac_fragments) & atac_fragments > 1000]
selection[, pass_tss_gt_2 := is.finite(tss_enrichment) & tss_enrichment > 2]
selection[, pass_nucleosome_lt_2 := is.finite(nucleosome_signal) & nucleosome_signal < 2]

pass_columns <- c(
  "pass_rna_umi_gt_1000", "pass_genes_gt_500", "pass_mito_lt_20",
  "pass_atac_fragments_gt_1000", "pass_tss_gt_2", "pass_nucleosome_lt_2"
)
selection[, qc_use_downstream := Reduce(`&`, .SD), .SDcols = pass_columns]
reason_labels <- c(
  pass_rna_umi_gt_1000 = "rna_umi_le_1000",
  pass_genes_gt_500 = "detected_genes_le_500",
  pass_mito_lt_20 = "mitochondrial_percent_ge_20",
  pass_atac_fragments_gt_1000 = "atac_fragments_le_1000",
  pass_tss_gt_2 = "tss_enrichment_le_2",
  pass_nucleosome_lt_2 = "nucleosome_signal_ge_2"
)
pass_matrix <- as.matrix(selection[, ..pass_columns])
selection[, qc_failure_reasons := apply(pass_matrix, 1L, function(x) {
  failed <- reason_labels[!x]
  if (!length(failed)) "" else paste(unname(failed), collapse = ";")
})]
selection[, qc_status := fifelse(qc_use_downstream, "use_downstream", "retained_not_selected")]

sample_summary <- selection[, .(
  cells_retained = .N,
  cells_use_downstream = sum(qc_use_downstream),
  cells_retained_not_selected = sum(!qc_use_downstream),
  percent_use_downstream = 100 * mean(qc_use_downstream)
), by = .(library_id, knockout_target, genotype, tissue)]

criterion_summary <- rbindlist(lapply(pass_columns, function(column) {
  selection[, .(
    criterion = column,
    cells_retained = .N,
    cells_pass = sum(get(column)),
    cells_fail = sum(!get(column)),
    percent_pass = 100 * mean(get(column))
  ), by = .(library_id, knockout_target, genotype, tissue)]
}))
criterion_total <- rbindlist(lapply(pass_columns, function(column) {
  selection[, .(
    criterion = column,
    cells_retained = .N,
    cells_pass = sum(get(column)),
    cells_fail = sum(!get(column)),
    percent_pass = 100 * mean(get(column))
  )][, `:=`(library_id = "ALL", knockout_target = "ALL", genotype = "ALL", tissue = "BM")]
}))
criterion_summary <- rbind(criterion_summary, criterion_total, use.names = TRUE)

fwrite(selection, file.path(output_dir, "cell_qc_selection.tsv.gz"), sep = "\t", compress = "gzip", na = "NA")
fwrite(sample_summary, file.path(output_dir, "sample_selection_summary.tsv"), sep = "\t")
fwrite(criterion_summary, file.path(output_dir, "criterion_selection_summary.tsv"), sep = "\t")
writeLines(selection[qc_use_downstream == TRUE, cell_id], file.path(output_dir, "downstream_cell_ids.txt"))
writeLines(selection[qc_use_downstream == FALSE, cell_id], file.path(output_dir, "retained_not_selected_cell_ids.txt"))

thresholds <- list(
  combination = "all criteria must pass",
  inequalities_are_strict = TRUE,
  rna_umi = list(operator = ">", value = 1000),
  detected_genes = list(operator = ">", value = 500),
  mitochondrial_percent = list(operator = "<", value = 20),
  atac_fragments = list(operator = ">", value = 1000),
  tss_enrichment = list(operator = ">", value = 2),
  nucleosome_signal = list(operator = "<", value = 2)
)
write_json(thresholds, file.path(output_dir, "thresholds.json"), pretty = TRUE, auto_unbox = TRUE)
validation <- list(
  status = "ok",
  physical_cell_filtering_performed = FALSE,
  input_cells = nrow(qc),
  retained_rows = nrow(selection),
  unique_cell_ids = uniqueN(selection$cell_id),
  cells_use_downstream = sum(selection$qc_use_downstream),
  cells_retained_not_selected = sum(!selection$qc_use_downstream),
  partition_is_complete = sum(selection$qc_use_downstream) + sum(!selection$qc_use_downstream) == nrow(selection),
  cell_order_preserved = identical(as.character(selection$cell_id), as.character(qc$cell_id)),
  missing_selection_flags = sum(is.na(selection$qc_use_downstream))
)
write_json(validation, file.path(output_dir, "validation.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
write_json(
  list(
    status = "complete", completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    physical_cell_filtering_performed = FALSE,
    outputs = c(
      "cell_qc_selection.tsv.gz", "sample_selection_summary.tsv",
      "criterion_selection_summary.tsv", "downstream_cell_ids.txt",
      "retained_not_selected_cell_ids.txt", "thresholds.json", "validation.json", "session_info.txt"
    )
  ), completion_path, pretty = TRUE, auto_unbox = TRUE
)
