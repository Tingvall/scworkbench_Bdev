#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--config")
  stop("Usage: 01b_audit_hub_concentration.R --config <config.yaml>")

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(digest)
  library(ggplot2)
})
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")

cfg_path <- normalizePath(args[[2L]], mustWork = TRUE)
cfg <- yaml::read_yaml(cfg_path)
.libPaths(c(cfg$runtime$run_library, .libPaths()))
suppressPackageStartupMessages(library(scworkbench))
root <- cfg$storage$results_root
tables <- file.path(root, "01_tables")
figures <- file.path(root, "02_figures")
reports <- file.path(root, "03_reports")
logs <- file.path(root, "04_logs")
sink(file.path(logs, "02_hub_concentration_audit.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

source_manifest <- file.path(cfg$inputs$source_mapping_bundle$path, "manifest.json")
bdev_manifest <- file.path(cfg$inputs$approved_atac_umap$bundle_path, "manifest.json")
source_before <- sha256(source_manifest); bdev_before <- sha256(bdev_manifest)
stopifnot(source_before == cfg$inputs$source_mapping_bundle$manifest_sha256,
  bdev_before == cfg$inputs$approved_atac_umap$manifest_sha256)

ledger <- read.delim(file.path(tables, "06_pseudocell_neighbor_ledger.tsv"),
  check.names = FALSE, quote = "\"", na.strings = "NA")
bdev <- wb_open_bundle(cfg$inputs$approved_atac_umap$bundle_path,
  validate = FALSE, read_only = TRUE)
meta <- wb_meta(bdev, cells = unique(ledger$cell_id), cols = list(
  wb_ref_metadata_field("core:cell_id"), "orig.ident", "peaks_snn_res.1.2",
  "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"))
wb_close_bundle(bdev)
joined <- merge(ledger, meta, by = "cell_id", sort = FALSE)
stopifnot(nrow(joined) == nrow(ledger))

cluster_weight <- aggregate(normalized_weight ~ gsm + query_cell_id +
  peaks_snn_res.1.2, joined, sum)
sample_weight <- aggregate(normalized_weight ~ gsm + query_cell_id + orig.ident,
  joined, sum)
libraries <- unique(ledger$gsm)
composition_rows <- lapply(libraries, function(gsm) {
  c6 <- cluster_weight[cluster_weight$gsm == gsm &
    as.character(cluster_weight$peaks_snn_res.1.2) == "6", ]
  sw <- sample_weight[sample_weight$gsm == gsm, ]
  sample_fraction <- function(sample_id) {
    values <- sw$normalized_weight[as.character(sw$orig.ident) == sample_id]
    if (!length(values)) 0 else median(values)
  }
  data.frame(gsm = gsm,
    median_cluster6_neighbor_weight = if (!nrow(c6)) 0 else
      median(c6$normalized_weight),
    median_reference_sample1_weight = sample_fraction("1"),
    median_reference_sample2_weight = sample_fraction("2"),
    median_reference_sample3_weight = sample_fraction("3"),
    stringsAsFactors = FALSE)
})
composition <- do.call(rbind, composition_rows); rownames(composition) <- NULL
write_tsv(composition,
  file.path(tables, "13_library_neighbor_technical_composition.tsv"))

occurrence <- table(ledger$cell_id)
library_count <- tapply(ledger$gsm, ledger$cell_id,
  function(x) length(unique(x)))
total_weight <- tapply(ledger$normalized_weight, ledger$cell_id, sum)
hub <- meta[match(names(occurrence), meta$cell_id), ]
hub$pseudocell_occurrence <- as.integer(occurrence[hub$cell_id])
hub$pseudocell_occurrence_fraction <- hub$pseudocell_occurrence /
  length(unique(ledger$query_cell_id))
hub$libraries_selecting_cell <- as.integer(library_count[hub$cell_id])
hub$total_normalized_weight <- as.numeric(total_weight[hub$cell_id])
hub <- hub[order(-hub$pseudocell_occurrence, -hub$total_normalized_weight,
  hub$cell_id), ]
write_tsv(hub, file.path(tables, "14_recurrent_reference_hubs.tsv"))

top <- head(hub, 20L)
top$cell_label <- factor(top$cell_id, levels = rev(top$cell_id))
p <- ggplot(top, aes(pseudocell_occurrence_fraction, cell_label,
  fill = factor(peaks_snn_res.1.2))) +
  geom_col(color = "white", linewidth = .25) +
  geom_text(aes(label = paste0("sample ", orig.ident)), hjust = 1.03,
    color = "white", size = 3) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)) +
  labs(title = "Rank-sparse pseudo-cells repeatedly select the same reference hubs",
    subtitle = "Top 20 cells by occurrence across 2,800 pseudo-cell primary-30 neighborhoods",
    x = "Fraction of pseudo-cells selecting reference cell", y = NULL,
    fill = "Reference cluster") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank())
ggsave(file.path(figures, "05_reference_hub_concentration.png"), p,
  width = 10, height = 8, dpi = 220, bg = "white")
ggsave(file.path(figures, "05_reference_hub_concentration.pdf"), p,
  width = 10, height = 8, device = cairo_pdf, bg = "white")

source_after <- sha256(source_manifest); bdev_after <- sha256(bdev_manifest)
stopifnot(source_before == source_after, bdev_before == bdev_after)
summary <- list(
  schema_version = "bdev.immgen_label_blind_sparsity_adapter_hub_audit.v1",
  status = "complete_post_projection_audit", run_identity = cfg$run_identity,
  unique_reference_neighbors = nrow(hub),
  maximum_pseudocell_occurrence = max(hub$pseudocell_occurrence),
  maximum_pseudocell_occurrence_fraction = max(
    hub$pseudocell_occurrence_fraction),
  top_hub_cell = hub$cell_id[[1L]], top_hub_sample = hub$orig.ident[[1L]],
  top_hub_cluster = as.character(hub$peaks_snn_res.1.2[[1L]]),
  median_library_cluster6_weight = median(
    composition$median_cluster6_neighbor_weight),
  median_library_sample1_weight = median(
    composition$median_reference_sample1_weight),
  labels_used_for_adapter = FALSE, audit_ran_after_projection = TRUE,
  source_manifests_unchanged = TRUE, mapping_accepted = FALSE)
write_json(summary, file.path(reports, "02_hub_concentration_summary.json"))
cat("Completed post-projection hub audit for", nrow(hub),
  "unique reference neighbors\n")
