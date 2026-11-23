#!/usr/bin/env Rscript

# Reciprocally lift the exact frozen HVP-depth20k axis to GRCm39 and quantify
# only the bounded four-library KO fixture. Archived inputs are opened read-only.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(paste(
    "Usage: 02_liftover_and_quantify_fixture.R AXIS_TSV FIXTURE_TSV",
    "MM10_TO_MM39_CHAIN MM39_TO_MM10_CHAIN FRAGMENTS OUTPUT_DIR PROCESS_N"
  ), call. = FALSE)
}

axis_path <- normalizePath(args[[1L]], mustWork = TRUE)
fixture_path <- normalizePath(args[[2L]], mustWork = TRUE)
forward_chain_path <- normalizePath(args[[3L]], mustWork = TRUE)
reverse_chain_path <- normalizePath(args[[4L]], mustWork = TRUE)
fragments_path <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
process_n <- as.integer(args[[7L]])
if (is.na(process_n) || process_n < 1L) stop("PROCESS_N must be positive", call. = FALSE)
if (!file.exists(paste0(fragments_path, ".tbi"))) stop("Missing fragment index", call. = FALSE)
if (file.exists(output_dir) || dir.exists(output_dir)) {
  stop("Refusing to overwrite output: ", output_dir, call. = FALSE)
}

options(width = 160, warn = 1)
suppressPackageStartupMessages({
  library(BPCells)
  library(Signac)
  library(GenomicRanges)
  library(rtracklayer)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(scworkbench)
})

write_json <- function(value, path) {
  jsonlite::write_json(value, path, pretty = TRUE, auto_unbox = TRUE,
                       null = "null", na = "null", digits = NA)
}

log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}

axis <- fread(axis_path)
fixture <- fread(fixture_path)
stopifnot(
  nrow(axis) == 20000L,
  identical(axis$feature_order, seq_len(20000L)),
  !anyNA(axis$feature_id), !anyDuplicated(axis$feature_id),
  nrow(fixture) == 400L, !anyNA(fixture$cell_id),
  !anyDuplicated(fixture$cell_id),
  identical(sort(unique(fixture$library_id)),
            sort(c("FLT3LKO_BM_WT", "FLT3LKO_BM_KO", "IL7KO_BM_WT", "IL7KO_BM_KO"))),
  all(table(fixture$library_id) == 100L),
  setequal(unique(fixture$knockout_target), c("Flt3l", "Il7")),
  setequal(unique(fixture$genotype), c("WT", "KO"))
)

parent <- dirname(output_dir)
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
stage_dir <- file.path(parent, paste0(".", basename(output_dir), ".tmp-", Sys.getpid()))
if (dir.exists(stage_dir)) unlink(stage_dir, recursive = TRUE)
dir.create(stage_dir, recursive = FALSE)
on.exit(if (dir.exists(stage_dir)) unlink(stage_dir, recursive = TRUE), add = TRUE)
tables_dir <- file.path(stage_dir, "01_tables")
reports_dir <- file.path(stage_dir, "03_reports")
objects_dir <- file.path(stage_dir, "05_objects")
manifests_dir <- file.path(stage_dir, "06_manifests")
dir.create(tables_dir)
dir.create(reports_dir)
dir.create(objects_dir)
dir.create(manifests_dir)

feature_ids <- axis$feature_id
mm10 <- Signac::StringToGRanges(feature_ids, sep = c("-", "-"))
names(mm10) <- feature_ids
stopifnot(length(mm10) == 20000L, all(width(mm10) > 0L))

log_event("Importing reciprocal liftOver chains")
forward_chain <- rtracklayer::import.chain(forward_chain_path)
reverse_chain <- rtracklayer::import.chain(reverse_chain_path)
forward <- rtracklayer::liftOver(mm10, forward_chain)
forward_n <- S4Vectors::elementNROWS(forward)
single_forward <- which(forward_n == 1L)
mm39 <- unlist(forward[single_forward], use.names = FALSE)
names(mm39) <- feature_ids[single_forward]
standard <- grepl("^chr([1-9]|1[0-9]|X|Y)$", as.character(seqnames(mm39)))
mm39_key <- paste0(seqnames(mm39), ":", start(mm39), "-", end(mm39))
duplicated_mm39 <- duplicated(mm39_key) | duplicated(mm39_key, fromLast = TRUE)

reverse_n <- integer(length(mm39))
reciprocal_exact <- rep(FALSE, length(mm39))
eligible <- which(standard & !duplicated_mm39)
if (length(eligible)) {
  reverse <- rtracklayer::liftOver(mm39[eligible], reverse_chain)
  reverse_n[eligible] <- S4Vectors::elementNROWS(reverse)
  reverse_single <- which(S4Vectors::elementNROWS(reverse) == 1L)
  if (length(reverse_single)) {
    reverse_gr <- unlist(reverse[reverse_single], use.names = FALSE)
    candidate_position <- eligible[reverse_single]
    original_position <- single_forward[candidate_position]
    reciprocal_exact[candidate_position] <-
      as.character(seqnames(reverse_gr)) == as.character(seqnames(mm10[original_position])) &
      start(reverse_gr) == start(mm10[original_position]) &
      end(reverse_gr) == end(mm10[original_position])
  }
}

accepted_local <- standard & !duplicated_mm39 & reciprocal_exact
accepted_original <- single_forward[accepted_local]
status <- rep("forward_unmapped", length(feature_ids))
status[forward_n > 1L] <- "forward_multimapped"
status[single_forward] <- "candidate"
status[single_forward[!standard]] <- "nonstandard_mm39_chromosome"
status[single_forward[standard & duplicated_mm39]] <- "duplicate_mm39_interval"
status[single_forward[standard & !duplicated_mm39 & reverse_n != 1L]] <- "reverse_not_unique"
status[single_forward[standard & !duplicated_mm39 & reverse_n == 1L & !reciprocal_exact]] <- "reverse_not_exact"
status[accepted_original] <- "accepted"

mapping <- data.table(
  feature_order = axis$feature_order,
  feature_id = feature_ids,
  mm10_seqname = as.character(seqnames(mm10)),
  mm10_start_1based = start(mm10),
  mm10_end_1based = end(mm10),
  forward_mapping_count = forward_n,
  mm39_seqname = NA_character_,
  mm39_start_1based = NA_integer_,
  mm39_end_1based = NA_integer_,
  reciprocal_mapping_count = NA_integer_,
  status = status
)
mapping$mm39_seqname[single_forward] <- as.character(seqnames(mm39))
mapping$mm39_start_1based[single_forward] <- start(mm39)
mapping$mm39_end_1based[single_forward] <- end(mm39)
mapping$reciprocal_mapping_count[single_forward] <- reverse_n
fwrite(mapping, file.path(tables_dir, "03_hvp_depth20k_mm10_to_mm39.tsv.gz"),
       sep = "\t", quote = TRUE, na = "NA")

status_counts <- as.list(table(status))
liftover_report <- list(
  status = if (mean(status == "accepted") >= 0.9975) "passed_with_structural_zero_mask" else "failed_unavailable_fraction_gate",
  input_features = length(feature_ids),
  accepted_features = sum(status == "accepted"),
  unavailable_features = sum(status != "accepted"),
  unavailable_fraction = mean(status != "accepted"),
  status_counts = status_counts,
  exact_model_axis_required = TRUE,
  unavailable_feature_policy = "explicit_structural_zero_with_mask",
  maximum_unavailable_fraction = 0.0025
)
write_json(liftover_report, file.path(reports_dir, "02_liftover_audit.json"))

if (mean(status == "accepted") < 0.9975) {
  capture.output(sessionInfo(), file = file.path(manifests_dir, "02_session_info.txt"))
  if (!file.rename(stage_dir, output_dir)) stop("Could not publish failed liftOver audit", call. = FALSE)
  stop("Unavailable feature fraction gate failed: ", sum(status == "accepted"),
       " of ", length(status), " selected peaks accepted", call. = FALSE)
}

accepted <- mapping[status == "accepted"]
lookup <- GRanges(
  seqnames = accepted$mm39_seqname,
  ranges = IRanges(start = accepted$mm39_start_1based,
                   end = accepted$mm39_end_1based)
)
names(lookup) <- accepted$feature_id

query_cells <- fixture$cell_id
log_event("Creating read-only fragment object for", length(query_cells), "fixture cells")
fragment_object <- Signac::CreateFragmentObject(
  path = fragments_path,
  cells = setNames(query_cells, query_cells),
  validate.fragments = FALSE,
  verbose = FALSE
)
log_event("Quantifying exact HVP-depth20k axis")
counts <- Signac::FeatureMatrix(
  fragments = fragment_object,
  features = lookup,
  cells = query_cells,
  process_n = process_n,
  keep_all_features = FALSE,
  verbose = TRUE
)
if (!setequal(colnames(counts), query_cells)) stop("Query cell axis mismatch", call. = FALSE)
counts <- counts[, query_cells, drop = FALSE]
if (nrow(counts) != nrow(accepted)) stop("FeatureMatrix dropped exact lifted features", call. = FALSE)
rownames(counts) <- accepted$feature_id
stopifnot(
  identical(as.character(rownames(counts)), as.character(accepted$feature_id)),
  identical(as.character(colnames(counts)), as.character(query_cells)),
  all(counts@x >= 0), all(counts@x == floor(counts@x))
)

triplet <- summary(counts)
full_counts <- Matrix::sparseMatrix(
  i = accepted$feature_order[triplet$i],
  j = triplet$j,
  x = triplet$x,
  dims = c(length(feature_ids), length(query_cells)),
  dimnames = list(feature_ids, query_cells)
)
counts <- as(full_counts, "dgCMatrix")
stopifnot(
  identical(as.character(rownames(counts)), as.character(feature_ids)),
  identical(as.character(colnames(counts)), as.character(query_cells)),
  all(Matrix::rowSums(counts[status != "accepted", , drop = FALSE]) == 0)
)

matrix_dir <- file.path(objects_dir, "01_ATAC_hvp_depth20k_counts")
BPCells::write_matrix_dir(BPCells::convert_matrix_type(counts, "uint32_t"), matrix_dir)
persisted <- BPCells::open_matrix_dir(matrix_dir)
stopifnot(
  identical(dim(persisted), c(20000L, 400L)),
  identical(as.character(rownames(persisted)), feature_ids),
  identical(as.character(colnames(persisted)), as.character(query_cells))
)

cell_totals <- as.numeric(Matrix::colSums(counts))
cell_detected <- as.integer(Matrix::colSums(counts != 0))
feature_detected <- as.integer(Matrix::rowSums(counts != 0))
cell_summary <- copy(fixture[, .(cell_id, library_id, knockout_target, genotype, tissue)])
cell_summary[, total_hvp_depth20k_fragments := cell_totals]
cell_summary[, detected_hvp_depth20k_peaks := cell_detected]
fwrite(cell_summary, file.path(tables_dir, "04_fixture_quantification_summary.tsv.gz"),
       sep = "\t", quote = TRUE, na = "NA")

quantification <- list(
  status = "passed",
  dimensions = unname(dim(counts)),
  nonzero_entries = length(counts@x),
  features_with_signal = sum(feature_detected > 0L),
  unavailable_features_structural_zero = sum(status != "accepted"),
  unavailable_fraction = mean(status != "accepted"),
  unavailable_feature_policy = "explicit_structural_zero_with_mask",
  cells_with_zero_counts = sum(cell_totals == 0),
  cell_total_summary = as.list(summary(cell_totals)),
  cell_detected_feature_summary = as.list(summary(cell_detected)),
  exact_model_feature_order = TRUE,
  fixture_cells_exact = TRUE,
  source_archives_modified = FALSE,
  paper_analysis = FALSE
)
write_json(quantification, file.path(reports_dir, "03_quantification_audit.json"))
capture.output(sessionInfo(), file = file.path(manifests_dir, "02_session_info.txt"))

if (!file.rename(stage_dir, output_dir)) stop("Could not publish quantified fixture", call. = FALSE)
log_event("Published", output_dir)
