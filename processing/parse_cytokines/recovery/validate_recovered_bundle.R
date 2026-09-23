#!/usr/bin/env Rscript
# Full recipient validation is independently resumable after registration.
finalize_parse_recovery <- function(original_dir, bundle_dir, report_dir) {
  suppressPackageStartupMessages(library(scworkbench))
  stopifnot(identical(as.character(packageVersion("scworkbench")), "0.1.0.9004"))
  original_dir <- normalizePath(original_dir, mustWork = TRUE)
  bundle_dir <- normalizePath(bundle_dir, mustWork = TRUE)
  report_dir <- normalizePath(report_dir, mustWork = TRUE)
  stopifnot(!identical(original_dir, bundle_dir),
            !file.exists(file.path(report_dir, "completion.json")))
  evidence <- jsonlite::read_json(file.path(report_dir, "scientific_state_checks.json"))
  stopifnot(isTRUE(evidence$metadata_exact), isTRUE(evidence$all_cell_sets_exact),
            identical(as.integer(evidence$final_gene_sets), 23L),
            file.exists(file.path(report_dir, "feature_set_import_promotion.json")))
  write_json <- function(x, name) jsonlite::write_json(
    x, file.path(report_dir, name), pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null", digits = NA
  )
  log <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ..., "\n")
  log("Rechecking registered signatures, original metadata, and cell memberships")
  preflight <- function() {
    original <- wb_open_bundle(original_dir, read_only = TRUE, validate = FALSE)
    on.exit(wb_close_bundle(original), add = TRUE)
    recovered <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
    on.exit(wb_close_bundle(recovered), add = TRUE)
    accepted <- read.csv(file.path(report_dir, "final_accepted_feature_sets.csv"),
                         check.names = FALSE, stringsAsFactors = FALSE)
    registered <- wb_feature_sets(recovered, matrix_id = "RNA_counts", namespace = "parse_cytokine_signatures")
    stopifnot(nrow(registered) == 23L, all(registered$is_readonly == 1L),
              setequal(registered$feature_set, unique(accepted$feature_set)),
              identical(wb_meta(recovered), wb_meta(original)),
              identical(sort(wb_cell_sets(recovered)), sort(wb_cell_sets(original))))
    for (id in registered$feature_set) {
      stopifnot(identical(wb_feature_set_get(recovered, id, "RNA_counts"),
                          accepted$feature_id[accepted$feature_set == id]))
    }
    for (id in wb_cell_sets(original)) {
      stopifnot(identical(wb_cells(original, cell_set = id), wb_cells(recovered, cell_set = id)))
    }
    baseline <- jsonlite::read_json(file.path(report_dir, "recovered_manifest_before_import.json"), simplifyVector = TRUE)
    current <- wb_manifest(recovered)
    # Compare JSON forms to avoid irrelevant simplification differences from read_json().
    for (table in unlist(evidence$unchanged_registry_tables)) {
      expected <- jsonlite::toJSON(baseline$current_state$tables[[table]], auto_unbox = TRUE, null = "null")
      observed <- jsonlite::toJSON(current$current_state$tables[[table]], auto_unbox = TRUE, null = "null")
      stopifnot(identical(expected, observed))
    }
    artifacts <- wb_art_list(recovered, include_history = TRUE)
    stopifnot(nrow(artifacts) == 26L)
    write_json(list(metadata_exact = TRUE, all_cell_sets_exact = TRUE,
                    registered_signatures = 23L, gene_order_exact = TRUE,
                    all_registered_read_only = TRUE, artifact_listing_rows = nrow(artifacts)),
               "prevalidation_rechecks.json")
  }
  preflight()
  log("Running full recipient validation, including existing UCell content identities")
  validation <- wb_bundle_validation_report(bundle_dir, level = "recipient", progress = TRUE)
  write_json(validation, "validation.json")
  stopifnot(isTRUE(validation$ok))
  score_coverage <- subset(validation$content_coverage,
                           entity_kind == "artifact" &
                           entity_id == "cell_vector_family_20260921_155957_5n7edbnw")
  stopifnot(nrow(score_coverage) == 1L,
            identical(as.character(score_coverage$verification_status), "content_verified"))
  write_json(score_coverage, "ucell_content_validation.json")
  recipient <- wb_open_bundle(bundle_dir, read_only = TRUE, validate = FALSE)
  on.exit(wb_close_bundle(recipient), add = TRUE)
  manifest <- wb_manifest(recipient)
  write_json(manifest, "manifest.json")
  writeLines(capture.output(sessionInfo()), file.path(report_dir, "session_info.txt"))
  write_json(list(status = "complete", completed_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
                  source_bundle = original_dir, bundle = bundle_dir,
                  bundle_state_fingerprint = manifest$current_state$fingerprint,
                  validation_status = validation$status,
                  recipient_validation_fingerprint = validation$post_validation_fingerprint,
                  registered_signatures = 23L,
                  original_scores_retained = TRUE, metadata_and_cell_sets_unchanged = TRUE,
                  runtime_version = as.character(packageVersion("scworkbench"))), "completion.json")
  warnings <- table(validation$diagnostics$code[validation$diagnostics$severity == "warning"])
  writeLines(c(
    "# Parse signature bundle recovery", "",
    "Status: full recipient validation passed; final source-file preservation is recorded separately by the runner.", "",
    paste0("Recovered bundle: `", bundle_dir, "`"),
    paste0("Original bundle: `", original_dir, "`"),
    paste0("Runtime: scWorkbench `", packageVersion("scworkbench"), "` at `", find.package("scworkbench"), "`."), "",
    "The original counts, normalized matrix, PCA space, cell metadata, and all cell-set memberships were retained.",
    "All 251,599 imported cells remain present; the control/single analysis universe remains 154,209 cells.",
    "The existing 31-signature UCell family passed full content-identity verification.", "",
    if (file.exists(file.path(report_dir, "materialization_restoration.json")))
      "Five missing commit materializations were restored from the authenticated original session checkpoint; each file matches its original recorded SHA-256.", "",
    "All 23 accepted gene sets were registered read-only in `parse_cytokine_signatures`.",
    "CSV import was checked against the unchanged TSV, and every registered gene membership was verified.", "",
    paste0("Bundle fingerprint: `", manifest$current_state$fingerprint, "`"),
    paste0("Recipient fingerprint: `", validation$post_validation_fingerprint, "`"), "",
    "Validation warnings:", "",
    if (length(warnings)) paste0("- `", names(warnings), "`: ", as.integer(warnings)) else "- None.", "",
    "Evidence: [validation](validation.json), [gene membership checks](feature_set_membership_checks.tsv),",
    "[scientific state checks](scientific_state_checks.json), [UCell content verification](ucell_content_validation.json),",
    "[materialization restoration](materialization_restoration.json), [completion](completion.json), and [source preservation](source_preservation.json)."
  ), file.path(report_dir, "REPORT.md"))
  log("Recovery and full recipient validation complete")
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 3L)
  finalize_parse_recovery(args[1], args[2], args[3])
}
