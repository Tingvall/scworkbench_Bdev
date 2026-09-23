#!/usr/bin/env Rscript
# Restore missing commit files from the authenticated, retained checkpoint.
# This performs no analysis and writes no registry rows. Every restored file
# must match the original registry's raw-file SHA-256 before publication.
restore_parse_materializations <- function(bundle_dir, report_dir, session_dir) {
  suppressPackageStartupMessages(library(scworkbench))
  bundle_dir <- normalizePath(bundle_dir, mustWork = TRUE)
  report_dir <- normalizePath(report_dir, mustWork = TRUE)
  session_dir <- normalizePath(session_dir, mustWork = TRUE)
  stopifnot(grepl("/bundle_recovered_[0-9]{8}T[0-9]{6}Z$", bundle_dir))
  sha <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
  pointer_path <- file.path(session_dir, "CURRENT.json")
  pointer <- jsonlite::read_json(pointer_path)
  stopifnot(grepl("^[A-Za-z0-9_.-]+$", pointer$generation))
  checkpoint_dir <- file.path(session_dir, "checkpoints", pointer$generation)
  metadata_path <- file.path(checkpoint_dir, "session.json")
  checkpoint_path <- file.path(checkpoint_dir, "session.rds")
  stopifnot(identical(sha(metadata_path), pointer$metadata_sha256))
  metadata <- jsonlite::read_json(metadata_path)
  stopifnot(identical(sha(checkpoint_path), metadata$checkpoint_sha256))
  # This is our own checkpoint, authenticated above. Opening its session would
  # mutate/resume the historical workspace; reading the saved outputs does not.
  session <- readRDS(checkpoint_path)
  stopifnot(identical(session$session_id, "parse_signature_ucell_20260921T152017"))
  commit_id <- "commit_20260921_155948_s32blwqt"
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(bundle_dir, "db/app.sqlite"),
                        flags = RSQLite::SQLITE_RO)
  expected <- tryCatch(DBI::dbGetQuery(con,
    "SELECT * FROM prov_materializations WHERE commit_id = ? ORDER BY run_id, output_name;",
    params = list(commit_id)), finally = DBI::dbDisconnect(con))
  stopifnot(nrow(expected) == 5L,
            all(expected$content_identity_scheme == "scworkbench.content.raw-file.sha256.v1"),
            all(expected$output_hash_scheme == "scworkbench.materialization.r-serialize-v3.sha256.v2"),
            setequal(expected$output_name, c("signature_scores", "signature_summary", "group_summary", "key_summary", "pairwise_summary")))
  candidate_dir <- tempfile("materialization-restore-")
  dir.create(candidate_dir)
  on.exit(unlink(candidate_dir, recursive = TRUE), add = TRUE)
  candidates <- character(nrow(expected))
  destinations <- file.path(bundle_dir, expected$path)
  for (i in seq_len(nrow(expected))) {
    row <- expected[i, ]
    stopifnot(identical(row$path, file.path("provenance", "materializations", commit_id,
                                         row$run_id, paste0(row$output_name, ".rds"))))
    prim <- session$runs[[row$run_id]]$outputs[[row$output_name]]
    stopifnot(!is.null(prim), identical(prim$kind, row$primitive_kind),
              row$primitive_kind %in% c("obs_table", "cell_vector_family"))
    # The engine's promotion_transitive.R leaves these primitive kinds intact
    # and uses exactly saveRDS(version = 3), followed by both identities below.
    stopifnot(identical(getFromNamespace(".wb_prov_materialization_hash", "scworkbench")(prim), row$output_hash))
    candidates[i] <- file.path(candidate_dir, paste0(i, ".rds"))
    saveRDS(prim, candidates[i], version = 3)
    stopifnot(identical(sha(candidates[i]), row$content_sha256))
    if (file.exists(destinations[i])) stopifnot(identical(sha(destinations[i]), row$content_sha256))
    cat("VERIFIED", row$run_id, row$output_name, row$content_sha256, "\n")
  }
  missing <- !file.exists(destinations)
  # Acquire the normal scWorkbench writer lease before restoring contained
  # files, and keep all registry identities unchanged across the operation.
  bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
  closed <- FALSE
  on.exit(if (!closed) wb_close_bundle(bundle), add = TRUE)
  before <- wb_manifest(bundle)$current_state$fingerprint
  for (i in which(missing)) {
    destination <- destinations[i]
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    stopifnot(!file.exists(destination),
              startsWith(normalizePath(dirname(destination)), paste0(bundle_dir, "/")))
    temporary <- tempfile(".restore-", tmpdir = dirname(destination))
    stopifnot(file.copy(candidates[i], temporary),
              identical(sha(temporary), expected$content_sha256[i]),
              !file.exists(destination), file.rename(temporary, destination))
  }
  stopifnot(all(vapply(seq_len(nrow(expected)), function(i) {
    identical(sha(destinations[i]), expected$content_sha256[i])
  }, logical(1))), identical(wb_manifest(bundle)$current_state$fingerprint, before))
  wb_close_bundle(bundle)
  closed <- TRUE
  stopifnot(identical(sha(metadata_path), pointer$metadata_sha256),
            identical(sha(checkpoint_path), metadata$checkpoint_sha256))
  jsonlite::write_json(list(status = "complete", commit_id = commit_id,
    checkpoint = checkpoint_path, checkpoint_sha256 = metadata$checkpoint_sha256,
    checkpoint_metadata_sha256 = pointer$metadata_sha256,
    restored_files = sum(missing), all_five_raw_file_sha256_exact = TRUE,
    all_five_materialization_identities_exact = TRUE,
    registry_state_unchanged = TRUE, checkpoint_unchanged = TRUE, files = expected),
    file.path(report_dir, "materialization_restoration.json"),
    pretty = TRUE, auto_unbox = TRUE, null = "null")
  cat("PASS: restored", sum(missing), "missing materialization files with exact original identities\n")
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 3L)
  restore_parse_materializations(args[1], args[2], args[3])
}
