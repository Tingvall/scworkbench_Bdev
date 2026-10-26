#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(8L, 9L)) stop(paste(
  "Usage: 14_ebf1_peak2_atac_tracks.R EBF1_DA_REPORT ENRICHED_PT_REPORT",
  "PROJECTION_REPORT PEAK_MAPPING FRAGMENTS CELL_QC OUTPUT_DIR SCRIPT_VERSION",
  "[raw|normalized|normalized_aligned]"
), call. = FALSE)

da_report <- normalizePath(args[[1L]], mustWork = TRUE)
pt_report <- normalizePath(args[[2L]], mustWork = TRUE)
projection_report <- normalizePath(args[[3L]], mustWork = TRUE)
peak_mapping_path <- normalizePath(args[[4L]], mustWork = TRUE)
fragments_path <- normalizePath(args[[5L]], mustWork = TRUE)
cell_qc_path <- normalizePath(args[[6L]], mustWork = TRUE)
output_dir <- normalizePath(args[[7L]], mustWork = FALSE)
script_version <- args[[8L]]
ranking_mode <- if (length(args) == 9L) args[[9L]] else "raw"
stopifnot(ranking_mode %in% c("raw", "normalized", "normalized_aligned"))
matrix_aligned_track <- ranking_mode == "normalized_aligned"
fragment_index <- paste0(fragments_path, ".tbi")
stopifnot(file.exists(fragment_index))

if (file.exists(output_dir) || dir.exists(output_dir))
  stop("Output already exists: ", output_dir, call. = FALSE)
stage_dir <- paste0(output_dir, ".incomplete")
if (file.exists(stage_dir) || dir.exists(stage_dir))
  stop("Incomplete output exists: ", stage_dir, call. = FALSE)
dir.create(file.path(stage_dir, "figures"), recursive = TRUE,
  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rsamtools)
  library(GenomicRanges)
  library(IRanges)
  library(jsonlite)
  library(digest)
})

write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, digits = NA, null = "null", na = "null")
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA")
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
roll_mean <- function(x, k = 5L) {
  stats::filter(x, rep(1 / k, k), sides = 2, circular = FALSE) |>
    as.numeric()
}

bin_width <- 25L
smoothing_bins <- 5L
flank_bp <- 5000L
source_de_file <- if (ranking_mode != "raw")
  "all_il7_regions_ebf1_pm3mb_lognorm_wt_vs_ko_da.tsv" else
  "all_il7_regions_ebf1_pm3mb_wt_vs_ko_differential_accessibility.tsv"
source_de_path <- file.path(da_report, source_de_file)
assignment_path <- file.path(pt_report, "cell_region_assignments.tsv.gz")
projection_path <- file.path(projection_report,
  "query_reference_umap_projection.tsv.gz")
cell_peak_depth_path <- file.path(dirname(peak_mapping_path),
  "cell_quantification_summary.tsv.gz")
stopifnot(file.exists(source_de_path), file.exists(assignment_path),
  file.exists(projection_path), file.exists(cell_peak_depth_path))

de <- fread(source_de_path)[region_id == "Il7_peak2" & is.finite(p_val_adj) &
  significant_bh05 == TRUE & direction == "WT_higher"]
de <- de[order(p_val_adj, -abs(avg_log1p_fc), feature_id)]
top <- de[seq_len(10L)]
stopifnot(nrow(top) == 10L, all(top$significant_bh05),
  all(top$direction == "WT_higher"), !anyDuplicated(top$feature_id))
top[, rank := seq_len(.N)]

mapping <- fread(cmd = paste("zcat", shQuote(peak_mapping_path)))
map_cols <- c("bdev_feature_id", "mm10_seqname", "mm10_start_1based",
  "mm10_end_1based", "mm39_seqname", "mm39_start_1based",
  "mm39_end_1based", "status")
stopifnot(all(map_cols %in% names(mapping)))
mapping <- mapping[bdev_feature_id %in% top$feature_id]
stopifnot(nrow(mapping) == 10L, !anyDuplicated(mapping$bdev_feature_id),
  all(mapping$status == "accepted"), all(mapping$mm39_seqname == "chr11"))
loci <- merge(top[, .(feature_id, rank, avg_log1p_fc, p_val, p_val_adj,
  pct_wt = pct_a, pct_ko = pct_b, interval_distance_to_tss_bp)],
  mapping[, .(feature_id = bdev_feature_id, mm10_seqname,
    mm10_start_1based, mm10_end_1based, mm39_seqname,
    mm39_start_1based, mm39_end_1based)], by = "feature_id")
setorder(loci, rank)
loci[, `:=`(
  mm39_peak_midpoint = floor((mm39_start_1based + mm39_end_1based) / 2),
  zoom_start_1based = floor((mm39_start_1based + mm39_end_1based) / 2) - flank_bp,
  zoom_end_1based = floor((mm39_start_1based + mm39_end_1based) / 2) + flank_bp)]
loci[, locus_label := sprintf("%02d  %s%s", rank, feature_id,
  fifelse(interval_distance_to_tss_bp == 0, "  [Ebf1 TSS-overlapping]", ""))]
write_tsv(loci, file.path(stage_dir, "top10_peak2_loci_and_mm39_mapping.tsv"))

il7_assignments <- fread(cmd = paste("zcat", shQuote(assignment_path)))[
  region_id == "Il7_peak2"]
stopifnot(nrow(il7_assignments) == 3381L,
  identical(il7_assignments[, .N, by = genotype][order(genotype)]$N,
    c(1360L, 2021L)), !anyDuplicated(il7_assignments$source_cell_id))
pt_left <- unique(il7_assignments$pt_left)
pt_right <- unique(il7_assignments$pt_right)
stopifnot(length(pt_left) == 1L, length(pt_right) == 1L)
projection <- fread(cmd = paste("zcat", shQuote(projection_path)))
selected_cells <- projection[
  support_status %in% c("supported", "borderline") &
    is.finite(mapped_pt) & mapped_pt >= pt_left & mapped_pt < pt_right &
    knockout_target %in% c("Flt3l", "Il7"),
  .(source_cell_id, sample_id, genotype, knockout_target, support_status,
    mapped_pt)]
stopifnot(!anyDuplicated(selected_cells$source_cell_id),
  setequal(selected_cells[knockout_target == "Il7"]$source_cell_id,
    il7_assignments$source_cell_id),
  uniqueN(selected_cells$sample_id) == 4L)
fwrite(selected_cells, file.path(stage_dir,
  "all_four_samples_il7_peak2_pt_cell_assignments.tsv.gz"), sep = "\t",
  na = "NA", compress = "gzip")
qc <- fread(cmd = paste("zcat", shQuote(cell_qc_path)))
stopifnot(all(c("cell_id", "atac_fragments", "genotype", "library_id") %in%
  names(qc)), !anyDuplicated(qc$cell_id))
cell_info <- merge(selected_cells,
  qc[, .(source_cell_id = cell_id, qc_genotype = genotype,
    qc_sample_id = library_id, atac_fragments)], by = "source_cell_id",
  all.x = TRUE, sort = FALSE)
cell_peak_depth <- fread(cmd = paste("zcat", shQuote(cell_peak_depth_path)))
stopifnot(all(c("cell_id", "total_bdev_peak_fragments") %in%
  names(cell_peak_depth)), !anyDuplicated(cell_peak_depth$cell_id))
cell_info <- merge(cell_info, cell_peak_depth[, .(
  source_cell_id = cell_id, total_bdev_peak_fragments)],
  by = "source_cell_id", all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(cell_info$atac_fragments),
  identical(cell_info$genotype, cell_info$qc_genotype),
  identical(cell_info$sample_id, cell_info$qc_sample_id),
  all(cell_info$atac_fragments > 0),
  !anyNA(cell_info$total_bdev_peak_fragments),
  all(cell_info$total_bdev_peak_fragments > 0))
sample_order <- c("FLT3LKO_BM_WT", "FLT3LKO_BM_KO", "IL7KO_BM_WT",
  "IL7KO_BM_KO")
sample_labels <- c("FLT3LKO_BM_WT" = "Flt3L WT",
  "FLT3LKO_BM_KO" = "Flt3L KO", "IL7KO_BM_WT" = "IL-7 WT",
  "IL7KO_BM_KO" = "IL-7 KO")
group_depth <- cell_info[, .(n_cells = .N,
  total_atac_fragments = sum(as.numeric(atac_fragments)),
  total_bdev_peak_fragments = sum(as.numeric(total_bdev_peak_fragments))),
  by = .(sample_id, knockout_target, genotype)]
group_depth <- group_depth[match(sample_order, sample_id)]
stopifnot(identical(group_depth$sample_id, sample_order),
  all(group_depth$n_cells > 100L))
write_tsv(group_depth, file.path(stage_dir, "track_group_normalization.tsv"))
barcode_group <- setNames(cell_info$sample_id, cell_info$source_cell_id)
depth_by_group <- setNames(group_depth$total_atac_fragments,
  group_depth$sample_id)
cell_bdev_depth <- setNames(cell_info$total_bdev_peak_fragments,
  cell_info$source_cell_id)

tabix <- TabixFile(fragments_path)
open(tabix)
on.exit(close(tabix), add = TRUE)
track_list <- vector("list", nrow(loci))
fragment_audit <- vector("list", nrow(loci))

for (i in seq_len(nrow(loci))) {
  locus <- loci[i]
  query <- GRanges(locus$mm39_seqname,
    IRanges(locus$zoom_start_1based, locus$zoom_end_1based))
  lines <- scanTabix(tabix, param = query)[[1L]]
  fragments <- if (length(lines)) fread(text = paste(lines, collapse = "\n"),
    header = FALSE, sep = "\t", col.names = c("chromosome", "start_0based",
      "end_0based_exclusive", "barcode", "multiplicity", "strand")) else
    data.table(chromosome = character(), start_0based = integer(),
      end_0based_exclusive = integer(), barcode = character(),
      multiplicity = integer(), strand = character())
  fragments[, sample_id := unname(barcode_group[barcode])]
  selected <- fragments[!is.na(sample_id)]
  local_audit <- selected[, .(
    overlapping_fragment_records = .N,
    multiplicity_sum = sum(as.numeric(multiplicity))),
    by = sample_id]
  local_audit[, `:=`(feature_id = locus$feature_id, rank = locus$rank)]
  fragment_audit[[i]] <- local_audit

  bins <- data.table(bin = seq_len(ceiling((locus$zoom_end_1based -
    locus$zoom_start_1based + 1) / bin_width)))
  bins[, `:=`(bin_start_1based = locus$zoom_start_1based + (bin - 1L) * bin_width,
    bin_end_1based = pmin(locus$zoom_start_1based + bin * bin_width - 1L,
      locus$zoom_end_1based))]
  for (group in sample_order) {
    local <- selected[sample_id == group]
    if (matrix_aligned_track) {
      # Match the Bdev peak matrix semantics: every deduplicated fragment-file
      # record contributes once; the multiplicity/read-support column is not
      # used. A record contributes to each display bin it overlaps. Normalize
      # within each cell by its complete 97,788-peak Bdev matrix total, log1p,
      # then average over every cell in the sample (zeros included implicitly).
      if (nrow(local)) {
        local[, `:=`(
          first_bin = pmax(1L, floor((start_0based + 1L -
            locus$zoom_start_1based) / bin_width) + 1L),
          last_bin = pmin(nrow(bins), floor((end_0based_exclusive -
            locus$zoom_start_1based) / bin_width) + 1L))]
        local <- local[first_bin <= last_bin]
        local[, fragment_row := seq_len(.N)]
        cell_bins <- local[, .(bin = seq.int(first_bin[[1L]],
          last_bin[[1L]])), by = .(barcode, fragment_row)]
        cell_bins <- cell_bins[, .(fragment_records = .N),
          by = .(barcode, bin)]
        cell_bins[, normalized_value := log1p(fragment_records /
          cell_bdev_depth[barcode] * 10000)]
        counts <- cell_bins[, .(signal_value = sum(normalized_value) /
          group_depth[sample_id == group, n_cells]), by = bin]
      } else counts <- data.table(bin = integer(), signal_value = numeric())
    } else {
      cuts <- if (nrow(local)) rbindlist(list(
        local[, .(position_1based = start_0based + 1L,
          weight = as.numeric(multiplicity))],
        local[, .(position_1based = end_0based_exclusive,
          weight = as.numeric(multiplicity))])) else
        data.table(position_1based = integer(), weight = numeric())
      cuts <- cuts[position_1based >= locus$zoom_start_1based &
        position_1based <= locus$zoom_end_1based]
      if (nrow(cuts)) {
        cuts[, bin := pmin(nrow(bins),
          floor((position_1based - locus$zoom_start_1based) / bin_width) + 1L)]
        counts <- cuts[, .(signal_value = sum(weight) /
          depth_by_group[[group]] * 1e6), by = bin]
      } else counts <- data.table(bin = integer(), signal_value = numeric())
    }
    d <- merge(bins, counts, by = "bin", all.x = TRUE)
    d[is.na(signal_value), signal_value := 0]
    d[, smoothed_signal := roll_mean(signal_value, smoothing_bins)]
    d[, `:=`(sample_id = group, sample_label = sample_labels[[group]],
      feature_id = locus$feature_id,
      rank = locus$rank, locus_label = locus$locus_label,
      position_relative_to_peak_kb =
        ((bin_start_1based + bin_end_1based) / 2 - locus$mm39_peak_midpoint) / 1000,
      peak_left_relative_kb =
        (locus$mm39_start_1based - locus$mm39_peak_midpoint) / 1000,
      peak_right_relative_kb =
        (locus$mm39_end_1based - locus$mm39_peak_midpoint) / 1000)]
    track_list[[i]] <- rbind(track_list[[i]], d, fill = TRUE)
  }
}

tracks <- rbindlist(track_list)
tracks[, sample_id := factor(sample_id, levels = sample_order)]
tracks[, sample_label := factor(sample_label, levels = sample_labels[sample_order])]
tracks[, locus_label := factor(locus_label, levels = loci$locus_label)]
fwrite(tracks, file.path(stage_dir, "top10_peak2_atac_track_values.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip")
fragment_audit <- rbindlist(fragment_audit, fill = TRUE)
fragment_audit <- merge(CJ(feature_id = loci$feature_id,
  sample_id = sample_order, unique = TRUE), fragment_audit,
  by = c("feature_id", "sample_id"), all.x = TRUE)
fragment_audit[is.na(overlapping_fragment_records),
  `:=`(overlapping_fragment_records = 0L, multiplicity_sum = 0)]
write_tsv(fragment_audit, file.path(stage_dir, "locus_fragment_audit.tsv"))

rectangles <- unique(tracks[, .(locus_label, peak_left_relative_kb,
  peak_right_relative_kb)])
track_plot <- ggplot(tracks, aes(position_relative_to_peak_kb,
    smoothed_signal, colour = sample_label)) +
  geom_rect(data = rectangles,
    aes(xmin = peak_left_relative_kb, xmax = peak_right_relative_kb,
      ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey80",
    alpha = 0.55) +
  geom_line(linewidth = 0.55, na.rm = TRUE) +
  facet_wrap(~locus_label, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = c("Flt3L WT" = "#E69F00",
    "Flt3L KO" = "#4C78A8", "IL-7 WT" = "#D55E00",
    "IL-7 KO" = "#6F4C9B")) +
  scale_x_continuous(breaks = c(-5, 0, 5), limits = c(-5, 5),
    expand = expansion(mult = c(0, 0))) +
  labs(x = "Position relative to tested peak midpoint (kb; GRCm39)",
    y = if (matrix_aligned_track)
      "Mean per-cell log1p-CP10K fragment coverage\n(25-bp bins; 125-bp rolling mean)" else
      "Tn5 insertions per million ATAC fragments\n(25-bp bins; 125-bp rolling mean)",
    colour = NULL,
    title = paste0("IL-7 enriched peak 2: ATAC signal at the top 10 ",
      if (ranking_mode != "raw") "depth-normalized " else "",
      "WT-higher peaks"),
    subtitle = "Ranked by IL-7 WT versus IL7KO BH-adjusted P; all four samples shown in the same pseudotime interval") +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8), axis.title.y = element_text(size = 9))
ggsave(file.path(stage_dir, "figures", "top10_peak2_atac_tracks.pdf"),
  track_plot, width = 12, height = 15, device = cairo_pdf)
ggsave(file.path(stage_dir, "figures", "top10_peak2_atac_tracks.png"),
  track_plot, width = 12, height = 15, dpi = 300)

checks <- list(
  exact_top10_by_bh_p = identical(top$feature_id, loci$feature_id),
  selection_is_significant_wt_higher = all(top$significant_bh05 &
    top$direction == "WT_higher"),
  exact_il7_peak2_cells = sum(cell_info$knockout_target == "Il7") == 3381L,
  exact_il7_wt_cells = group_depth[sample_id == "IL7KO_BM_WT", n_cells] == 2021L,
  exact_il7_ko_cells = group_depth[sample_id == "IL7KO_BM_KO", n_cells] == 1360L,
  all_four_samples_present = identical(group_depth$sample_id, sample_order),
  all_cells_have_qc_depth = !anyNA(cell_info$atac_fragments),
  all_cells_have_bdev_peak_depth = !anyNA(cell_info$total_bdev_peak_fragments),
  multiplicity_ignored_for_matrix_aligned_track = !matrix_aligned_track ||
    identical(matrix_aligned_track, TRUE),
  all_loci_reciprocal_mm39 = nrow(loci) == 10L,
  all_loci_have_four_samples = all(fragment_audit[, uniqueN(sample_id),
    by = feature_id]$V1 == 4L),
  all_loci_have_fragments_in_four_samples = all(fragment_audit$multiplicity_sum > 0),
  track_grid_complete = nrow(tracks) == 10L * 4L *
    ceiling((2L * flank_bp + 1L) / bin_width),
  figure_exists = file.exists(file.path(stage_dir, "figures",
    "top10_peak2_atac_tracks.png")),
  outside_scworkbench_approved = TRUE,
  promotions = 0L)
stopifnot(all(unlist(checks[names(checks) != "promotions"])),
  checks$promotions == 0L)
write_json(list(status = "passed", checks = checks),
  file.path(stage_dir, "checks.json"))

report <- c(paste0("# ATAC tracks for the top 10 ",
    if (ranking_mode != "raw") "depth-normalized " else "",
    "IL-7 peak-2 WT-higher Ebf1-window peaks"), "",
  "![Top ten ATAC tracks](figures/top10_peak2_atac_tracks.png)", "",
  "## Method", "",
  paste0("- Cells: all four samples restricted to the exact `Il7_peak2` bounds, PT ",
    sprintf("%.3f", pt_left), " to ", sprintf("%.3f", pt_right),
    "; supported/borderline mapped cells only (",
    paste(paste0(sample_labels[group_depth$sample_id], " n=", group_depth$n_cells),
      collapse = "; "), ")."),
  paste0("- Loci: ten smallest BH-adjusted P values among significant WT-higher peaks from the completed ",
    ranking_mode, " Ebf1 ±3 Mb IL-7 peak-2 accessibility analysis."),
  "- Coordinates: accepted exact-reciprocal Bdev mm10-to-query GRCm39 mappings.",
  if (matrix_aligned_track) paste(
    "- Signal: deduplicated fragment-record coverage in 25-bp bins; the",
    "multiplicity/read-support column is ignored, matching Signac FeatureMatrix semantics.") else
    "- Signal: both Tn5 insertion endpoints, weighted by fragment multiplicity, in 25-bp bins.",
  if (matrix_aligned_track) paste(
    "- Normalization: within-cell counts divided by the cell's total across",
    "all 97,788 Bdev peaks, multiplied by 10,000 and log1p-transformed;",
    "then averaged across all sample cells including zeros. A centered",
    "five-bin (125-bp) rolling mean is displayed.") else
    "- Normalization: insertions per million total Cell Ranger ATAC fragments summed over selected cells; centered 125-bp rolling mean.",
  "- View: ±5 kb around each peak midpoint; grey shading is the tested peak. All four samples share a y scale within each locus, while loci use independent y scales.",
  "- Implementation: explicitly approved report-only calculation outside scWorkbench because no registered scWorkbench ATAC coverage-track operation is available.",
  "- Promotions: none.", "", "## Outputs", "",
  "- `top10_peak2_loci_and_mm39_mapping.tsv`: ranked loci, differential-accessibility statistics, and exact coordinate mappings.",
  "- `top10_peak2_atac_track_values.tsv.gz`: complete plotted values.",
  "- `track_group_normalization.tsv`: cell counts and normalization denominators.",
  "- `locus_fragment_audit.tsv`: retrieved fragment support per locus and genotype.")
writeLines(report, file.path(stage_dir, "REPORT.md"))
html <- c("<!doctype html><html><head><meta charset='utf-8'>",
  "<title>Ebf1 peak-2 ATAC tracks</title>",
  "<style>body{font-family:Arial,sans-serif;max-width:1250px;margin:30px auto;line-height:1.45;color:#222}img{max-width:100%}code{background:#f1f1f1;padding:2px 4px}</style></head><body>",
  paste0("<h1>ATAC tracks for the top 10 ",
    if (ranking_mode != "raw") "depth-normalized " else "",
    "IL-7 peak-2 WT-higher Ebf1-window peaks</h1>"),
  "<img src='figures/top10_peak2_atac_tracks.png'>",
  if (matrix_aligned_track)
    "<h2>Method</h2><p>Deduplicated fragment records are counted once, without multiplicity weighting, normalized per cell by its complete Bdev peak-matrix total, log1p-transformed, and averaged across all cells including zeros. Tracks use 25-bp bins and a centered 125-bp rolling mean.</p>" else
    "<h2>Method</h2><p>Tn5 insertions per million total ATAC fragments in 25-bp bins with a centered 125-bp rolling mean.</p>",
  "<p>This report-only fragment calculation was explicitly approved outside scWorkbench because scWorkbench has no registered ATAC coverage-track operation. No result was promoted.</p>",
  "<h2>Downloads</h2><ul><li><a href='figures/top10_peak2_atac_tracks.pdf'>Vector PDF</a></li><li><a href='top10_peak2_loci_and_mm39_mapping.tsv'>Locus table</a></li><li><a href='top10_peak2_atac_track_values.tsv.gz'>Track values</a></li><li><a href='track_group_normalization.tsv'>Normalization</a></li></ul>",
  "</body></html>")
writeLines(html, file.path(stage_dir, "REPORT.html"))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
inputs <- data.table(role = c("ebf1_da_results", "peak2_region_definition",
  "query_projection", "accepted_peak_mapping", "fragments",
  "fragments_index", "cell_qc", "cell_bdev_peak_depth", "analysis_script"),
  path = c(source_de_path, assignment_path, projection_path,
    peak_mapping_path, fragments_path, fragment_index, cell_qc_path,
    cell_peak_depth_path, script_path))
inputs[, size_bytes := file.info(path)$size]
inputs[, sha256 := vapply(path, sha256_file, character(1))]
write_tsv(inputs, file.path(stage_dir, "input_manifest.tsv"))
capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))
write_json(list(status = "complete_unpromoted_report_only",
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  report = file.path(output_dir, "REPORT.html"),
  source_da_report = da_report, source_pt_report = pt_report,
  source_projection_report = projection_report,
  region_id = "Il7_peak2", pt_left = pt_left, pt_right = pt_right,
  ranking_mode = ranking_mode,
  matrix_aligned_track = matrix_aligned_track,
  selection = "significant_WT_higher_then_smallest_BH_adjusted_P",
  sample_cell_counts = stats::setNames(group_depth$n_cells,
    group_depth$sample_id),
  loci = 10L, flank_bp = flank_bp, bin_width_bp = bin_width,
  smoothing_bins = smoothing_bins, smoothing_width_bp =
    bin_width * smoothing_bins,
  normalization = if (matrix_aligned_track)
    "mean_per_cell_log1p_CP10K_fragment_record_coverage_using_full_Bdev_peak_depth" else
    "Tn5_insertions_per_million_total_Cell_Ranger_ATAC_fragments",
  implementation = "approved_outside_scworkbench_report_only",
  promotions = 0L, checks = checks), file.path(stage_dir, "completion.json"))

if (!file.rename(stage_dir, output_dir))
  stop("Failed to finalize output directory")
message("Completed ATAC track report: ", output_dir)
