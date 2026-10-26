#!/usr/bin/env Rscript

# Translate the canonical Bdev mm10 peak universe to temporary mm39 lookup
# intervals. The original Bdev feature IDs remain canonical. Only unique,
# standard-chromosome, exactly reciprocal liftOver mappings are accepted.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "Usage: 06_prepare_bdev_peak_universe.R BDEV_PEAKS_BPCELLS MM10_TO_MM39_CHAIN MM39_TO_MM10_CHAIN OUTPUT_DIR MIN_RETAINED_FRACTION",
    call. = FALSE
  )
}

bdev_matrix_dir <- normalizePath(args[[1L]], mustWork = TRUE)
forward_chain_path <- normalizePath(args[[2L]], mustWork = TRUE)
reverse_chain_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
minimum_retained_fraction <- as.numeric(args[[5L]])
if (!is.finite(minimum_retained_fraction) || minimum_retained_fraction <= 0 ||
    minimum_retained_fraction > 1) {
  stop("MIN_RETAINED_FRACTION must be in (0, 1]", call. = FALSE)
}
if (file.exists(output_dir) || dir.exists(output_dir)) {
  stop("Output directory already exists; refusing to overwrite: ", output_dir, call. = FALSE)
}

options(width = 160, warn = 1)
suppressPackageStartupMessages({
  library(BPCells)
  library(Signac)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(jsonlite)
  library(scworkbench)
})

log_event <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}

write_json <- function(value, path) {
  jsonlite::write_json(value, path, pretty = TRUE, auto_unbox = TRUE,
                       null = "null", na = "null", digits = NA)
}

sha256 <- function(path) scworkbench::wb_source_sha256(path)

main <- function() {
  parent <- dirname(output_dir)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  stage_dir <- file.path(parent, paste0(".", basename(output_dir), ".tmp-", Sys.getpid()))
  if (dir.exists(stage_dir)) unlink(stage_dir, recursive = TRUE)
  dir.create(stage_dir, recursive = FALSE)
  on.exit(if (dir.exists(stage_dir)) unlink(stage_dir, recursive = TRUE), add = TRUE)

  log_event("Opening Bdev peaks_counts BPCells matrix read-only")
  bdev <- BPCells::open_matrix_dir(bdev_matrix_dir)
  peak_ids <- rownames(bdev)
  stopifnot(nrow(bdev) == 97846L, length(peak_ids) == nrow(bdev),
            !anyNA(peak_ids), !anyDuplicated(peak_ids))

  mm10 <- Signac::StringToGRanges(peak_ids, sep = c("-", "-"))
  names(mm10) <- peak_ids
  stopifnot(length(mm10) == length(peak_ids), all(width(mm10) > 0L))

  log_event("Importing official UCSC liftOver chains")
  forward_chain <- rtracklayer::import.chain(forward_chain_path)
  reverse_chain <- rtracklayer::import.chain(reverse_chain_path)
  log_event("Lifting", length(mm10), "Bdev peaks from mm10 to mm39")
  forward <- rtracklayer::liftOver(mm10, forward_chain)
  forward_n <- S4Vectors::elementNROWS(forward)
  single_forward <- which(forward_n == 1L)
  mm39 <- unlist(forward[single_forward], use.names = FALSE)
  names(mm39) <- peak_ids[single_forward]

  standard_mm39 <- grepl("^chr([1-9]|1[0-9]|X|Y)$", as.character(seqnames(mm39)))
  mm39_key <- paste0(seqnames(mm39), ":", start(mm39), "-", end(mm39))
  duplicated_mm39 <- duplicated(mm39_key) | duplicated(mm39_key, fromLast = TRUE)

  reciprocal_exact <- rep(FALSE, length(mm39))
  reverse_n <- integer(length(mm39))
  eligible_reverse <- which(standard_mm39 & !duplicated_mm39)
  if (length(eligible_reverse)) {
    log_event("Running reciprocal mm39-to-mm10 validation on", length(eligible_reverse), "candidate peaks")
    reverse <- rtracklayer::liftOver(mm39[eligible_reverse], reverse_chain)
    reverse_n[eligible_reverse] <- S4Vectors::elementNROWS(reverse)
    reverse_single_local <- which(S4Vectors::elementNROWS(reverse) == 1L)
    if (length(reverse_single_local)) {
      reverse_gr <- unlist(reverse[reverse_single_local], use.names = FALSE)
      candidate_position <- eligible_reverse[reverse_single_local]
      original_position <- single_forward[candidate_position]
      reciprocal_exact[candidate_position] <-
        as.character(seqnames(reverse_gr)) == as.character(seqnames(mm10[original_position])) &
        start(reverse_gr) == start(mm10[original_position]) &
        end(reverse_gr) == end(mm10[original_position])
    }
  }

  accepted_local <- standard_mm39 & !duplicated_mm39 & reciprocal_exact
  accepted_original <- single_forward[accepted_local]
  accepted_mm39 <- mm39[accepted_local]
  accepted_ids <- peak_ids[accepted_original]
  names(accepted_mm39) <- accepted_ids

  status <- rep("forward_unmapped", length(peak_ids))
  status[forward_n > 1L] <- "forward_multimapped"
  status[single_forward] <- "candidate"
  status[single_forward[!standard_mm39]] <- "nonstandard_mm39_chromosome"
  status[single_forward[standard_mm39 & duplicated_mm39]] <- "duplicate_mm39_interval"
  status[single_forward[standard_mm39 & !duplicated_mm39 & reverse_n != 1L]] <-
    "reverse_not_unique"
  status[single_forward[standard_mm39 & !duplicated_mm39 & reverse_n == 1L & !reciprocal_exact]] <-
    "reverse_not_exact"
  status[accepted_original] <- "accepted"

  mapping <- data.table(
    bdev_row = seq_along(peak_ids),
    bdev_feature_id = peak_ids,
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

  accepted <- mapping[status == "accepted"]
  stopifnot(identical(accepted$bdev_feature_id, accepted_ids),
            !anyDuplicated(accepted$bdev_feature_id),
            !anyDuplicated(paste(accepted$mm39_seqname,
                                 accepted$mm39_start_1based,
                                 accepted$mm39_end_1based, sep = ":")))
  retained_fraction <- nrow(accepted) / nrow(mapping)
  if (retained_fraction < minimum_retained_fraction) {
    stop(sprintf("Retained fraction %.4f is below gate %.4f", retained_fraction,
                 minimum_retained_fraction), call. = FALSE)
  }

  fwrite(mapping, file.path(stage_dir, "bdev_mm10_to_mm39_peak_mapping.tsv.gz"),
         sep = "\t", quote = TRUE, na = "NA")
  fwrite(accepted, file.path(stage_dir, "accepted_bdev_peak_mapping.tsv.gz"),
         sep = "\t", quote = TRUE, na = "NA")
  bed <- accepted[, .(
    chrom = mm39_seqname,
    start_0based = mm39_start_1based - 1L,
    end_0based_half_open = mm39_end_1based,
    bdev_feature_id
  )]
  fwrite(bed, file.path(stage_dir, "accepted_bdev_peaks_mm39_lookup.bed"),
         sep = "\t", col.names = FALSE, quote = FALSE)

  status_counts <- as.data.frame(table(mapping$status), stringsAsFactors = FALSE)
  names(status_counts) <- c("status", "peaks")
  completion <- list(
    status = "complete",
    stage = "bdev_peak_universe_liftover",
    canonical_feature_identity = "Bdev peaks_counts mm10 feature IDs",
    lookup_coordinate_system = "mm39/GRCm39",
    bdev_matrix_dir = bdev_matrix_dir,
    bdev_matrix_dimensions = unname(dim(bdev)),
    input_peaks = nrow(mapping),
    accepted_peaks = nrow(accepted),
    retained_fraction = retained_fraction,
    minimum_retained_fraction = minimum_retained_fraction,
    status_counts = split(status_counts$peaks, status_counts$status),
    chain_files = list(
      mm10_to_mm39 = list(path = forward_chain_path, sha256 = sha256(forward_chain_path)),
      mm39_to_mm10 = list(path = reverse_chain_path, sha256 = sha256(reverse_chain_path))
    ),
    acceptance = paste(
      "unique forward mapping; standard mouse chromosome; unique mapped interval;",
      "unique reverse mapping; exact reciprocal recovery of the original mm10 interval"
    ),
    bdev_parent_modified = FALSE
  )
  write_json(completion, file.path(stage_dir, "completion.json"))
  capture.output(sessionInfo(), file = file.path(stage_dir, "session_info.txt"))

  if (!file.rename(stage_dir, output_dir)) stop("Could not publish output directory", call. = FALSE)
  log_event("Accepted", nrow(accepted), "of", nrow(mapping), "Bdev peaks (",
            sprintf("%.2f%%", 100 * retained_fraction), ")")
}

main()

