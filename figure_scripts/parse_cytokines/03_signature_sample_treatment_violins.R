#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 03_signature_sample_treatment_violins.R <bundle_dir> <signature_report_dir>")
}

bundle_dir <- normalizePath(args[[1L]], mustWork = TRUE)
report_dir <- normalizePath(args[[2L]], mustWork = TRUE)
status_path <- file.path(report_dir, "all_signature_status.tsv")
scores_path <- file.path(report_dir, "ucell_sample_medians.tsv")
artifacts_path <- file.path(report_dir, "committed_ucell_artifacts.tsv")
for (path in c(status_path, scores_path, artifacts_path)) {
  if (!file.exists(path)) stop("Required input is missing: ", path)
}

suppressPackageStartupMessages({
  library(scworkbench)
  library(DBI)
})

pdf_path <- file.path(report_dir, "ucell_sample_treatment_violins_all_candidates.pdf")
figure_dir <- file.path(report_dir, "ucell_sample_treatment_violins")
manifest_path <- file.path(report_dir, "ucell_sample_treatment_violin_manifest.tsv")
exclusions_path <- file.path(report_dir, "signature_exclusion_reasons.tsv")
supplement_path <- file.path(report_dir, "SIGNATURE_EXCLUSIONS_AND_VIOLINS.md")
outputs <- c(pdf_path, figure_dir, manifest_path, exclusions_path, supplement_path)
if (any(file.exists(outputs))) {
  stop("One or more output paths already exist; refusing to overwrite: ",
       paste(outputs[file.exists(outputs)], collapse = ", "))
}

status <- read.delim(status_path, check.names = FALSE, na.strings = c("NA", ""))
scores <- read.delim(scores_path, check.names = FALSE, na.strings = c("NA", ""))
artifacts <- read.delim(artifacts_path, check.names = FALSE, na.strings = c("NA", ""))

required_status <- c(
  "signature_id", "cytokine", "direction", "status", "full_eligible_genes",
  "selected_genes", "selected_size", "discovery_cv_response_auc",
  "discovery_cv_specificity_auc", "discovery_cv_balanced_auc",
  "discovery_fold_gene_jaccard", "ucell_response_auc",
  "ucell_specificity_auc", "ucell_balanced_auc", "final_status"
)
required_scores <- c(
  "signature_id", "cytokine", "direction", "sample_id", "condition",
  "replicate", "treatment_category", "sample_median_ucell"
)
missing_status <- setdiff(required_status, names(status))
missing_scores <- setdiff(required_scores, names(scores))
if (length(missing_status)) stop("Status table is missing: ", paste(missing_status, collapse = ", "))
if (length(missing_scores)) stop("Score table is missing: ", paste(missing_scores, collapse = ", "))
if (anyDuplicated(status$signature_id)) stop("Status table has duplicate signature IDs")
if (anyDuplicated(scores[c("signature_id", "sample_id")])) {
  stop("Score table has duplicate signature/sample rows")
}

candidates <- status[status$status == "validated", , drop = FALSE]
if (nrow(status) != 52L || nrow(candidates) != 31L) {
  stop("Expected 52 cytokine/direction rows and 31 scored discovery candidates")
}
if (!setequal(unique(scores$signature_id), candidates$signature_id)) {
  stop("Scored signature IDs do not match the 31 discovery candidates")
}
score_counts <- table(scores$signature_id)
if (any(score_counts != 81L)) stop("Every scored candidate must have exactly 81 sample medians")
if (any(!is.finite(scores$sample_median_ucell)) || any(scores$sample_median_ucell < 0)) {
  stop("UCell sample medians must be finite and non-negative")
}
sample_counts <- aggregate(sample_id ~ signature_id + condition, scores, length)
if (any(sample_counts$sample_id != 3L) || any(table(sample_counts$signature_id) != 27L)) {
  stop("Every candidate must have three samples in each of 27 conditions")
}

dir.create(figure_dir, recursive = FALSE)

safe_name <- function(x) {
  value <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("^_+|_+$", "", value)
}

fmt <- function(x, digits = 3L) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

plot_signature <- function(info, metadata, values, sample_dat) {
  target <- info$cytokine[[1L]]
  direction <- info$direction[[1L]]
  final_status <- info$final_status[[1L]]
  conditions <- unique(metadata$condition)
  if (!("control" %in% conditions) || !(target %in% conditions)) {
    stop("Control or target condition is absent for ", info$signature_id[[1L]])
  }
  condition_order <- c("control", target, sort(setdiff(conditions, c("control", target))))
  sample_metadata <- unique(metadata[, c("sample_id", "condition", "replicate"), drop = FALSE])
  sample_metadata <- sample_metadata[
    order(match(sample_metadata$condition, condition_order),
          suppressWarnings(as.numeric(sample_metadata$replicate)), sample_metadata$sample_id),
    , drop = FALSE
  ]
  if (nrow(sample_metadata) != 81L) stop("Expected 81 samples")

  ymax <- max(0.05, max(values) * 1.05)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 0.05
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(10.3, 4.8, 5.0, 1.3), mgp = c(2.8, 0.8, 0), xpd = FALSE)
  plot(
    NA_real_, NA_real_, xlim = c(0.35, length(condition_order) + 0.65),
    ylim = c(0, ymax), xaxt = "n", xlab = "", ylab = "Per-cell UCell score",
    las = 1, bty = "l", main = info$signature_id[[1L]],
    sub = paste0(
      target, " ", direction, " | ", final_status,
      " | ", info$selected_genes[[1L]], " genes | sample-level UCell AUC: response ",
      fmt(info$ucell_response_auc[[1L]], 2L), ", specificity ",
      fmt(info$ucell_specificity_auc[[1L]], 2L), ", balanced ",
      fmt(info$ucell_balanced_auc[[1L]], 2L),
      " | each narrow violin is one biological sample"
    ), cex.main = 1.05, cex.sub = 0.84
  )
  grid(nx = NA, ny = NULL, col = "#E6E6E6", lty = 1)
  control_median <- median(sample_dat$sample_median_ucell[sample_dat$condition == "control"])
  abline(h = control_median, col = "#555555", lty = 3, lwd = 1.2)

  fill_for <- function(condition) {
    if (condition == "control") return("#BDBDBD")
    if (condition == target) {
      if (final_status == "accepted") return("#D73027")
      return("#F46D43")
    }
    "#91BFDB"
  }
  point_for <- function(condition) {
    if (condition == "control") return("#4D4D4D")
    if (condition == target) return("#9E0142")
    "#2166AC"
  }

  sample_offsets <- c(-0.22, 0, 0.22)
  for (i in seq_along(condition_order)) {
    local_samples <- sample_metadata[sample_metadata$condition == condition_order[[i]], , drop = FALSE]
    if (nrow(local_samples) != 3L) stop("Expected three samples for ", condition_order[[i]])
    fill <- fill_for(condition_order[[i]])
    edge <- point_for(condition_order[[i]])
    for (j in seq_len(nrow(local_samples))) {
      sample_id <- local_samples$sample_id[[j]]
      sample_values <- values[metadata$sample_id == sample_id]
      x <- i + sample_offsets[[j]]
      if (!length(sample_values)) stop("No cells for sample ", sample_id)
      if (diff(range(sample_values)) > .Machine$double.eps^0.5) {
        den <- density(sample_values, n = 256, bw = "nrd0")
        keep <- which(
          den$x >= 0 & den$x <= ymax & is.finite(den$y) &
            den$y >= max(den$y, na.rm = TRUE) * 0.002
        )
        if (length(keep) >= 2L) {
          keep <- seq.int(min(keep), max(keep))
          half_width <- 0.09 * den$y[keep] / max(den$y[keep])
          polygon(
            c(x - half_width, rev(x + half_width)), c(den$x[keep], rev(den$x[keep])),
            col = adjustcolor(fill, alpha.f = 0.62), border = adjustcolor(edge, alpha.f = 0.88),
            lwd = 0.65
          )
        }
      } else {
        segments(x - 0.07, sample_values[[1L]], x + 0.07, sample_values[[1L]],
                 col = adjustcolor(edge, alpha.f = 0.9), lwd = 3)
      }
      points(x, median(sample_values), pch = 21, bg = "white", col = "#111111", cex = 0.62, lwd = 0.75)
    }
  }

  axis(1, at = seq_along(condition_order), labels = FALSE, tck = -0.012)
  usr <- par("usr")
  label_y <- usr[[3L]] - 0.047 * diff(usr[3:4])
  label_col <- ifelse(
    condition_order == target, "#9E0142",
    ifelse(condition_order == "control", "#333333", "#333333")
  )
  text(
    x = seq_along(condition_order), y = label_y, labels = condition_order,
    srt = 55, adj = 1, cex = 0.70, xpd = NA, col = label_col
  )
  mtext("Treatment condition (three biological samples per condition)", side = 1,
        line = 8.9, cex = 0.88)
  legend(
    "topright", inset = 0.005, bty = "n", cex = 0.78,
    legend = c("control samples", "target-treatment samples", "other single-cytokine samples",
               "sample median", "median of control sample medians"),
    pch = c(15, 15, 15, 21, NA), pt.bg = c("#BDBDBD", "#D73027", "#91BFDB", "white", NA),
    col = c("#4D4D4D", "#9E0142", "#2166AC", "#111111", "#555555"),
    lty = c(NA, NA, NA, NA, 3), lwd = c(NA, NA, NA, NA, 1.2)
  )
}

score_artifacts <- artifacts[artifacts$type == "cell_vector_family", , drop = FALSE]
if (nrow(score_artifacts) != 1L) stop("Expected exactly one committed UCell score-family artifact")
score_artifact_id <- score_artifacts$artifact_id[[1L]]
bundle <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
on.exit(if (!is.null(bundle)) try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
analysis_cells <- as.character(wb_cells(bundle, cell_set = "set:parse_qc_control_single"))
metadata <- wb_meta(
  bundle, cells = analysis_cells,
  cols = list(wb_ref_metadata_field("core:cell_id"), "sample_id", "condition",
              "replicate", "treatment_category")
)
if (length(analysis_cells) != 154209L || !identical(as.character(metadata$cell_id), analysis_cells) || anyNA(metadata)) {
  stop("Analysis cell set or metadata did not match the completed signature analysis")
}

read_committed_score <- function(signature_id) {
  score <- DBI::dbGetQuery(
    bundle$con,
    paste(
      "SELECT cell_id, pseudotime FROM artifact_pseudotime",
      "WHERE artifact_id = ? AND lineage_id = ?;"
    ),
    params = list(score_artifact_id, signature_id)
  )
  if (nrow(score) != length(analysis_cells) || anyDuplicated(score$cell_id)) {
    stop("Committed UCell payload has invalid coverage for ", signature_id)
  }
  index <- match(analysis_cells, as.character(score$cell_id))
  if (anyNA(index)) stop("Committed UCell payload is missing cells for ", signature_id)
  values <- as.numeric(score$pseudotime[index])
  if (anyNA(values) || any(!is.finite(values))) stop("Non-finite UCell values for ", signature_id)
  values
}

manifest <- vector("list", nrow(candidates))
pdf(pdf_path, width = 14, height = 8.5, onefile = TRUE, useDingbats = FALSE)
for (i in seq_len(nrow(candidates))) {
  info <- candidates[i, , drop = FALSE]
  dat <- scores[scores$signature_id == info$signature_id[[1L]], , drop = FALSE]
  values <- read_committed_score(info$signature_id[[1L]])
  calculated_medians <- tapply(values, metadata$sample_id, median)
  expected_medians <- dat$sample_median_ucell[match(names(calculated_medians), dat$sample_id)]
  if (anyNA(expected_medians) || max(abs(as.numeric(calculated_medians) - expected_medians)) > 1e-12) {
    stop("Cell-level scores do not reproduce stored sample medians for ", info$signature_id[[1L]])
  }
  plot_signature(info, metadata, values, dat)
  png_file <- paste0(safe_name(info$signature_id[[1L]]), ".png")
  png_path <- file.path(figure_dir, png_file)
  png(png_path, width = 2800, height = 1700, res = 200, type = "cairo")
  plot_signature(info, metadata, values, dat)
  dev.off()
  manifest[[i]] <- data.frame(
    page = i,
    signature_id = info$signature_id[[1L]],
    cytokine = info$cytokine[[1L]],
    direction = info$direction[[1L]],
    genes = info$selected_genes[[1L]],
    final_status = info$final_status[[1L]],
    cells = length(values),
    samples = nrow(dat),
    treatments = length(unique(dat$condition)),
    target_samples = sum(dat$condition == info$cytokine[[1L]]),
    zero_sample_medians = sum(dat$sample_median_ucell == 0),
    nonzero_sample_medians = sum(dat$sample_median_ucell > 0),
    png = file.path("ucell_sample_treatment_violins", png_file),
    stringsAsFactors = FALSE
  )
}
dev.off()
wb_close_bundle(bundle)
bundle <- NULL
manifest <- do.call(rbind, manifest)
write.table(manifest, manifest_path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")

score_summary <- do.call(rbind, lapply(split(scores, scores$signature_id), function(dat) {
  target <- unique(dat$cytokine)
  if (length(target) != 1L) stop("Inconsistent cytokine labels in score table")
  control <- dat$sample_median_ucell[dat$condition == "control"]
  target_values <- dat$sample_median_ucell[dat$condition == target]
  other_values <- dat$sample_median_ucell[dat$treatment_category == "single" & dat$condition != target]
  data.frame(
    signature_id = unique(dat$signature_id),
    scored_samples = nrow(dat),
    zero_sample_medians = sum(dat$sample_median_ucell == 0),
    nonzero_sample_medians = sum(dat$sample_median_ucell > 0),
    minimum_sample_median = min(dat$sample_median_ucell),
    maximum_sample_median = max(dat$sample_median_ucell),
    control_group_median = median(control),
    target_group_median = median(target_values),
    other_single_group_median = median(other_values),
    stringsAsFactors = FALSE
  )
}))

exclusions <- merge(status, score_summary, by = "signature_id", all.x = TRUE, sort = FALSE)
exclusions <- exclusions[match(status$signature_id, exclusions$signature_id), , drop = FALSE]
exclusions$gene_set_available <- exclusions$status == "validated"
exclusions$violin_available <- exclusions$signature_id %in% manifest$signature_id
exclusions$exclusion_stage <- ifelse(
  exclusions$final_status == "accepted", "none",
  ifelse(exclusions$final_status == "not_advanced_from_discovery", "discovery", "UCell transfer")
)
exclusions$exclusion_reason <- NA_character_
for (i in seq_len(nrow(exclusions))) {
  if (exclusions$final_status[[i]] == "accepted") {
    exclusions$exclusion_reason[[i]] <- "Accepted; passed discovery and sample-median UCell transfer criteria."
  } else if (exclusions$final_status[[i]] == "not_advanced_from_discovery") {
    n <- exclusions$full_eligible_genes[[i]]
    exclusions$exclusion_reason[[i]] <- paste0(
      "Only ", n, " gene", if (!is.na(n) && n == 1) "" else "s",
      " passed both contrasts, effect/FDR criteria, and three-fold direction stability; fewer than the prespecified minimum of 5 genes, so no gene set or UCell score was created."
    )
  } else {
    exclusions$exclusion_reason[[i]] <- paste0(
      "A discovery gene set was scored, but sample-median UCell transfer did not reach all prespecified cutoffs (response AUC >=0.60, specificity AUC >=0.60, balanced AUC >=0.65); observed ",
      fmt(exclusions$ucell_response_auc[[i]], 2L), "/",
      fmt(exclusions$ucell_specificity_auc[[i]], 2L), "/",
      fmt(exclusions$ucell_balanced_auc[[i]], 2L), "."
    )
  }
}

exclusion_columns <- c(
  "signature_id", "cytokine", "direction", "status", "full_eligible_genes",
  "selected_genes", "discovery_cv_response_auc", "discovery_cv_specificity_auc",
  "discovery_cv_balanced_auc", "discovery_fold_gene_jaccard",
  "ucell_response_auc", "ucell_specificity_auc", "ucell_balanced_auc",
  "scored_samples", "zero_sample_medians", "nonzero_sample_medians",
  "minimum_sample_median", "maximum_sample_median", "control_group_median",
  "target_group_median", "other_single_group_median", "final_status",
  "gene_set_available", "violin_available", "exclusion_stage", "exclusion_reason"
)
write.table(
  exclusions[exclusion_columns], exclusions_path, sep = "\t", quote = TRUE,
  row.names = FALSE, na = "NA"
)

md_escape <- function(x) {
  x <- ifelse(is.na(x), "NA", as.character(x))
  gsub("\\|", "\\\\|", x)
}

md_table <- function(dat) {
  header <- paste0("| ", paste(md_escape(names(dat)), collapse = " | "), " |")
  divider <- paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  rows <- apply(dat, 1L, function(x) paste0("| ", paste(md_escape(x), collapse = " | "), " |"))
  c(header, divider, rows)
}

discovery_excluded <- exclusions[exclusions$final_status == "not_advanced_from_discovery", , drop = FALSE]
transfer_excluded <- exclusions[exclusions$final_status == "ucell_transfer_not_confirmed", , drop = FALSE]
discovery_table <- data.frame(
  cytokine = discovery_excluded$cytokine,
  direction = discovery_excluded$direction,
  eligible_genes = discovery_excluded$full_eligible_genes,
  reason = "fewer than 5 eligible genes",
  stringsAsFactors = FALSE
)
transfer_table <- data.frame(
  cytokine = transfer_excluded$cytokine,
  direction = transfer_excluded$direction,
  genes = transfer_excluded$selected_genes,
  response_AUC = fmt(transfer_excluded$ucell_response_auc, 2L),
  specificity_AUC = fmt(transfer_excluded$ucell_specificity_auc, 2L),
  balanced_AUC = fmt(transfer_excluded$ucell_balanced_auc, 2L),
  zero_sample_medians = paste0(transfer_excluded$zero_sample_medians, "/81"),
  stringsAsFactors = FALSE
)

supplement <- c(
  "# Signature exclusions and per-treatment UCell distributions",
  "",
  paste0("Generated: `", format(Sys.time(), tz = "UTC", usetz = TRUE), "`"),
  "",
  "## What was excluded, and why",
  "",
  "The workflow evaluated 52 possible cytokine/direction combinations (26 single cytokines x induced/repressed). Exclusion occurred at one of two prespecified gates; there was no post hoc biological curation.",
  "",
  paste0("**Discovery gate (21 directions).** These directions produced fewer than the required five genes after requiring the target-versus-control and target-versus-other-single-cytokines contrasts, the prespecified effect-size and FDR limits, and stable direction in all three leave-one-replicate-out fits. Because no valid gene set existed, UCell scoring and a violin plot would be undefined."),
  "",
  md_table(discovery_table),
  "",
  paste0("**UCell transfer gate (8 gene sets).** These sets passed pseudobulk discovery but failed confirmation using biological-sample median UCell scores. Acceptance required response AUC >=0.60, specificity AUC >=0.60, and balanced AUC >=0.65. All eight missed those cutoffs; seven had a median score of zero in all 81 samples, and IL-7 up had zero in ", transfer_excluded$zero_sample_medians[transfer_excluded$signature_id == "parse_cytokine.il_7.up"], " of 81 samples."),
  "",
  md_table(transfer_table),
  "",
  "This does not establish that these genes or cytokine responses are biologically absent. Pseudobulk tests total expression across all cells in a sample and can detect a consistent response contributed by a responsive subset. The UCell gate asks a stricter, different question: is within-cell rank enrichment sufficiently broad and reproducible to shift the typical cell, summarized by each sample median? Sparse or tail-only enrichment can pass pseudobulk discovery while leaving the sample median at zero. These eight sets were therefore excluded conservatively as signatures that did not transfer to robust typical-cell scoring, not labeled biologically invalid.",
  "",
  "## Violin plots for every scored gene set",
  "",
  "The multi-page PDF contains all 31 discovery gene sets, including the 23 accepted and 8 UCell-transfer failures. Each treatment is a group of three narrow violins: one cell-level UCell distribution for each biological sample. The white-centered dot marks that sample's median, so biological replication remains explicit while the violin shows within-sample cellular heterogeneity. Control samples are gray, the gene set's intended target-treatment samples are red, and other single-cytokine samples are blue. The horizontal dotted line is the median of the three control-sample medians. Cell-level violin shapes are descriptive; the sample-level medians were the confirmation unit.",
  "",
  "- [All 31 gene sets, multi-page PDF](ucell_sample_treatment_violins_all_candidates.pdf)",
  "- [Individual PNG directory](ucell_sample_treatment_violins/)",
  "- [Figure manifest](ucell_sample_treatment_violin_manifest.tsv)",
  "- [Complete 52-direction exclusion table](signature_exclusion_reasons.tsv)",
  "",
  "The 21 discovery-failed directions are listed above and in the complete table, but are not plotted because they never produced a gene set and therefore have no UCell score. No Bdev, combined-stimulation cells, or newly chosen thresholds were used for this supplement. It reuses the Parse-only scores and statuses in this report directory."
)
writeLines(supplement, supplement_path, useBytes = TRUE)

message("Wrote 31 candidate gene-set violin pages and exclusion supplement to: ", report_dir)
