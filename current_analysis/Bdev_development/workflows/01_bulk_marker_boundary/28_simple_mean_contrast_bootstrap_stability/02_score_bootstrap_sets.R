#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, width = 180)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config")) stop("Usage: 02_score_bootstrap_sets.R --config <config.yaml>")
suppressPackageStartupMessages({ library(scworkbench); library(yaml); library(jsonlite); library(digest) })

cfg <- yaml::read_yaml(normalizePath(args[[2L]], mustWork = TRUE))
root <- cfg$storage$results_root
intermediate <- cfg$storage$intermediate_root
dirs <- setNames(file.path(root, c("01_tables", "03_reports", "04_logs", "06_manifests")), c("tables", "reports", "logs", "manifests"))
session_root <- file.path(intermediate, "05_objects", "01_sessions")
invisible(lapply(c(dirs, sessions = session_root), dir.create, recursive = TRUE, showWarnings = FALSE))
sink(file.path(dirs[["logs"]], "01_bootstrap_scoring.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
run_ref <- function(run, output) wb_ref_run(run$run_id, output)

membership <- read.delim(gzfile(file.path(dirs[["tables"]], "01_bootstrap_set_membership.tsv.gz")), check.names = FALSE)
fixture <- read.delim(cfg$inputs$fixture_cells$path, check.names = FALSE)
baseline <- read.delim(cfg$inputs$baseline_decisions$path, check.names = FALSE, na.strings = character())
stopifnot(nrow(fixture) == as.integer(cfg$fixture$expected_cells), !anyDuplicated(fixture$cell_id))
stopifnot(nrow(baseline) == nrow(fixture), identical(sort(baseline$cell_id), sort(fixture$cell_id)))
rounds <- sort(unique(as.integer(membership$bootstrap_round)))
stopifnot(length(rounds) == as.integer(cfg$bootstrap$rounds))

wb_operation_pack_load("scwbGenePrograms")
pack_record <- wb_operation_packs()
stopifnot(nrow(pack_record) == 1L, identical(as.character(pack_record$manifest_sha256[[1L]]), cfg$runtime$operation_pack_manifest_sha256))
bundle <- wb_open_bundle(cfg$inputs$execution_bundle$path, validate = FALSE, read_only = FALSE)
on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
manifest_before <- wb_manifest(bundle)
stopifnot(identical(manifest_before$current_state$fingerprint, cfg$inputs$execution_bundle$fingerprint))

fixture_md <- wb_meta(bundle, cells = fixture$cell_id)
stopifnot(identical(as.character(fixture_md$cell_id), fixture$cell_id))
metadata <- fixture_md[, c("cell_id", cfg$fixture$sample_field, cfg$fixture$pseudotime_field), drop = FALSE]
names(metadata)[names(metadata) == cfg$fixture$sample_field] <- cfg$fixture$reported_sample_field

all_sets <- list(); round_keys <- vector("list", length(rounds))
for (i in seq_along(rounds)) {
  round_id <- rounds[[i]]
  rows <- membership[membership$bootstrap_round == round_id, , drop = FALSE]
  rows <- rows[order(rows$population, rows$rank), , drop = FALSE]
  sets <- split(as.character(rows$gene_symbol), as.character(rows$population))
  stopifnot(length(sets) == 64L, all(lengths(sets) == 20L))
  keys <- paste0(sprintf("B%03d__", round_id), names(sets))
  names(sets) <- keys
  all_sets <- c(all_sets, sets)
  round_keys[[i]] <- keys
}
stopifnot(length(all_sets) == length(rounds) * 64L, !anyDuplicated(names(all_sets)))
batched_session_root <- file.path(session_root, "simple_mean_contrast_bootstrap_stability_batched_v1")
score_paths <- Sys.glob(file.path(batched_session_root, "stage", "operation_pack_v2", "transition_*", "outputs", "signature_scores", "payload.rds"))
score_primitive <- NULL; score_path <- NULL; score_run_id <- NULL
for (candidate in score_paths) {
  item <- readRDS(candidate)
  if (inherits(item, "wb_primitive") && identical(item$kind, "cell_vector_family") &&
      length(item$key) == nrow(fixture) * length(all_sets) && length(unique(item$key)) == length(all_sets)) {
    score_primitive <- item; score_path <- candidate; break
  }
}
if (is.null(score_primitive)) {
  session <- wb_session_create(bundle, workspace_root = session_root, session_id = "simple_mean_contrast_bootstrap_stability_batched_v1")
  score_result <- wb_session_run_op(session, bundle, "scwbgeneprograms::ucell", inputs = list(
    matrix = cfg$scoring$matrix_id, cell_ids = fixture$cell_id, min_cells = 10L,
    features = all_sets, name = "SIMPLE_MEAN_BOOTSTRAP_BATCH", max_rank = as.integer(cfg$scoring$ucell_max_rank),
    ncores = as.integer(cfg$scoring$ncores), parallel_backend = cfg$scoring$parallel_backend,
    chunk_size = as.integer(cfg$scoring$chunk_size), missing_features_policy = cfg$scoring$missing_features_policy),
    notes = "One native UCell pass for all 30 traceable ImmGen bootstrap collections; round membership is encoded in signal keys.")
  score_primitive <- score_result$run$outputs$signature_scores
  score_run_id <- score_result$run$run_id
  message("completed batched native UCell scoring")
} else {
  receipt_path <- file.path(dirname(dirname(dirname(score_path))), "receipt.json")
  if (file.exists(receipt_path)) score_run_id <- jsonlite::fromJSON(receipt_path)$run_id
  message("reusing completed batched native UCell score artifact")
}
stopifnot(length(score_primitive$key) == nrow(fixture) * length(all_sets),
  setequal(unique(as.character(score_primitive$key)), names(all_sets)))

reduce_exact <- function(signal, keys, round_id) {
  rows <- as.character(signal$key) %in% keys
  x <- data.frame(cell_id = as.character(signal$cell_id[rows]), key = as.character(signal$key[rows]),
    value = as.numeric(signal$value[rows]), stringsAsFactors = FALSE)
  cells <- unique(x$cell_id); key_order <- keys[keys %in% unique(x$key)]
  mat <- matrix(NA_real_, length(cells), length(key_order), dimnames = list(cells, key_order))
  mat[cbind(match(x$cell_id, cells), match(x$key, key_order))] <- x$value
  unknown <- cfg$scoring$unknown_label; tol <- as.numeric(cfg$scoring$tie_tolerance)
  label <- rep(unknown, length(cells)); top <- second <- margin <- rep(NA_real_, length(cells)); reason <- rep("none", length(cells))
  for (j in seq_along(cells)) {
    vals <- mat[j, ]; ok <- is.finite(vals)
    if (!any(ok)) next
    ordered <- order(vals[ok], decreasing = TRUE, na.last = TRUE)
    vv <- vals[ok][ordered]; kk <- names(vals[ok])[ordered]
    top[[j]] <- vv[[1]]; second[[j]] <- if (length(vv) > 1L) vv[[2]] else NA_real_
    margin[[j]] <- if (length(vv) > 1L) vv[[1]] - vv[[2]] else NA_real_
    tied <- abs(vv - vv[[1]]) <= tol
    if (sum(tied) > 1L && identical(cfg$scoring$ties, "unknown")) { reason[[j]] <- "tie"; next }
    winner <- if (sum(tied) > 1L) sort(kk[tied], method = "radix")[[1]] else kk[[1]]
    if (vv[[1]] < as.numeric(cfg$scoring$minimum_raw_score)) { reason[[j]] <- "min_score"; next }
    if (is.finite(margin[[j]]) && margin[[j]] < as.numeric(cfg$scoring$minimum_raw_margin)) { reason[[j]] <- "min_margin"; next }
    label[[j]] <- winner
  }
  assigned <- label != unknown; label[assigned] <- sub("^B[0-9]{3}__", "", label[assigned])
  data.frame(cell_id = cells, label = label, top_score = top, runner_up_score = second,
    score_margin = margin, abstain_reason = reason, bootstrap_round = round_id, stringsAsFactors = FALSE)
}

decision_rows <- vector("list", length(rounds))
for (i in seq_along(rounds)) {
  decision_rows[[i]] <- reduce_exact(score_primitive, round_keys[[i]], rounds[[i]])
  message(sprintf("completed in-memory bootstrap reduction %d/%d", i, length(rounds)))
}
decisions <- do.call(rbind, decision_rows)
decisions <- merge(decisions, metadata, by = "cell_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(decisions) == nrow(fixture) * length(rounds), !anyNA(decisions[[cfg$fixture$reported_sample_field]]))
output_path <- file.path(dirs[["tables"]], "05_bootstrap_bdev_decisions.tsv.gz")
con <- gzfile(output_path, "wt")
write.table(decisions, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
close(con)

manifest_after <- wb_manifest(bundle)
stopifnot(identical(manifest_after$current_state$fingerprint, manifest_before$current_state$fingerprint))
native_equivalence <- FALSE
native_decision_paths <- Sys.glob(file.path(batched_session_root, "stage", "operation_pack_v2", "transition_*", "outputs", "decision_table", "payload.rds"))
for (candidate in native_decision_paths) {
  item <- readRDS(candidate)
  if (!inherits(item, "wb_primitive") || is.null(item$table) || nrow(item$table) != nrow(fixture)) next
  native <- item$table
  if (!any(startsWith(as.character(native$label), "B001__"))) next
  assigned <- native$label != cfg$scoring$unknown_label
  native$label[assigned] <- sub("^B[0-9]{3}__", "", native$label[assigned])
  calculated <- decision_rows[[1L]]
  native <- native[match(calculated$cell_id, native$cell_id), names(native), drop = FALSE]
  native_equivalence <- identical(as.character(native$label), calculated$label) &&
    identical(as.character(native$abstain_reason), calculated$abstain_reason) &&
    isTRUE(all.equal(as.numeric(native$top_score), calculated$top_score, tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(native$runner_up_score), calculated$runner_up_score, tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(native$score_margin), calculated$score_margin, tolerance = 1e-12))
  if (native_equivalence) break
}
stopifnot(native_equivalence)
receipt <- list(schema_version = "bdev.simple_mean_contrast_bootstrap_scoring_receipt.v1",
  status = "complete_without_promotion", run_identity = cfg$run_identity,
  rounds = length(rounds), execution_strategy = cfg$scoring$execution_strategy,
  score_run_id = score_run_id, score_payload_path = score_path,
  score_payload_sha256 = if (is.null(score_path)) NULL else digest::digest(file = score_path, algo = "sha256"),
  reducer_equivalence_check = "exact_against_completed_native_round_1", reducer_equivalence_passed = native_equivalence,
  bundle_fingerprint_unchanged = TRUE, artifacts_promoted = FALSE)
write_json(receipt, file.path(dirs[["manifests"]], "01_scoring_session_receipt.json"))
print(receipt[c("status", "rounds", "bundle_fingerprint_unchanged", "artifacts_promoted")])
