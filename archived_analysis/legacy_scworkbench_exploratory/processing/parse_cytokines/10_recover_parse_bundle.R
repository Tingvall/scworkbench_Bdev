#!/usr/bin/env Rscript
# Reconcile the existing manifest and import already accepted gene sets using
# public scWorkbench APIs, with the isolated serialization-fixed engine.
suppressPackageStartupMessages(library(scworkbench))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
original_dir <- normalizePath(args[1], mustWork = TRUE)
bundle_dir <- normalizePath(args[2], mustWork = TRUE)
report_dir <- normalizePath(args[3], mustWork = TRUE)
accepted_path <- normalizePath(args[4], mustWork = TRUE)
stopifnot(!identical(original_dir, bundle_dir),
          identical(as.character(packageVersion("scworkbench")), "0.1.0.9004"))
write_json <- function(x, name) jsonlite::write_json(
  x, file.path(report_dir, name), pretty = TRUE, auto_unbox = TRUE,
  null = "null", na = "null", digits = NA
)
log <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), ..., "\n")
main <- function() {
  source <- wb_open_bundle(original_dir, read_only = TRUE, validate = FALSE)
  on.exit(wb_close_bundle(source), add = TRUE)
  write_json(wb_manifest(source), "source_durable_manifest.json")
  log("Reconstructing original registry identity with bounded JSON serialization")
  baseline <- wb_manifest(source, view = "reconstructed")
  write_json(baseline, "source_reconstructed_manifest.json")
  sets <- sort(as.character(wb_cell_sets(source)))
  memberships <- setNames(lapply(sets, function(id) as.character(wb_cells(source, cell_set = id))), sets)
  stopifnot(length(memberships[["set:parse_cytokines"]]) == 251599L,
            length(memberships[["set:parse_qc_control_single"]]) == 154209L)
  source_meta <- wb_meta(source)

  log("Opening copied bundle and publishing its authoritative manifest")
  bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
  closed <- FALSE
  on.exit(if (!closed) try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  before <- wb_manifest(bundle)
  stopifnot(identical(before$bundle_id, "parse_cytokines"))
  unchanged_tables <- c("cells", "metadata_fields", "metadata_field_identities", "datasets", "features",
                         "matrices", "matrix_cells", "matrix_fragments", "matrix_fragment_cells",
                         "matrix_fragment_features", "spaces", "space_payloads", "space_models",
                         "cell_set_membership", "cell_set_recipes", "cell_set_styles", "cell_vectors")
  for (table in unchanged_tables) {
    stopifnot(identical(before$current_state$tables[[table]], baseline$current_state$tables[[table]]))
  }
  write_json(before, "recovered_manifest_before_import.json")
  interface <- wb_op_interface("feature_set_import", bundle = bundle)
  write_json(interface, "feature_set_import_interface.json")
  accepted <- read.delim(accepted_path, check.names = FALSE, stringsAsFactors = FALSE)
  expected_sets <- unique(as.character(accepted$feature_set))
  stopifnot(length(expected_sets) == 23L,
            !anyDuplicated(accepted[c("feature_set", "feature_id")]),
            nrow(wb_feature_sets(bundle, matrix_id = "RNA_counts",
                                namespace = "parse_cytokine_signatures")) == 0L)
  # The documented table file parser reads CSV. Preserve the reviewed TSV and
  # supply an exactly equivalent CSV, with an explicit round-trip check.
  accepted_csv <- file.path(report_dir, "final_accepted_feature_sets.csv")
  write.csv(accepted, accepted_csv, row.names = FALSE)
  stopifnot(identical(read.csv(accepted_csv, check.names = FALSE, stringsAsFactors = FALSE), accepted))
  log("Registering the 23 previously accepted signature gene sets")
  # This bounded import is cheaply repeatable. Its full evidence becomes
  # durable through promotion; no expensive scientific computation is rerun.
  session <- wb_session_new()
  imported <- wb_session_run_op(session, bundle, "feature_set_import", inputs = list(
    matrix = "RNA_counts", cell_set = "set:parse_qc_control_single", min_cells = 10L,
    source_path = accepted_csv, source_format = "table",
    namespace = "parse_cytokine_signatures", collection = "individual_cytokine",
    species = "Mus musculus", organism = "Mus musculus", gene_id_type = "gene_symbol",
    set_id_col = "feature_set", set_name_col = "set_name", gene_col = "feature_id",
    description_col = "description", split_by_collection = FALSE, overwrite = FALSE,
    strict = TRUE, created_by = "jtg", source_id = "parse_cytokine_signature_workflow_v1",
    source_version = "2026-09-21", origin = "derived", is_readonly = TRUE
  ), notes = paste(
    "Register the unchanged 23 accepted gene sets from signature_report_20260921T221910Z.",
    "Recovery preserves the original UCell artifact and all existing scientific state."
  ))
  staged_manifest <- imported$run$outputs$import_manifest$table
  stopifnot(nrow(staged_manifest) == 23L,
            setequal(as.character(staged_manifest$feature_set), expected_sets),
            all(staged_manifest$status == "imported"),
            all(staged_manifest$n_genes_raw == staged_manifest$n_genes_mapped_to_bundle))
  promotion <- wb_session_promote_current(bundle, imported$session, outputs = list(
    wb_ref_run(imported$run$run_id, "import_summary"),
    wb_ref_run(imported$run$run_id, "import_manifest"),
    wb_ref_run(imported$run$run_id, "rejected_sets")
  ), name_prefix = "parse_cytokine_final_signatures",
  description = "Final accepted individual-cytokine transcriptional signature feature sets.",
  created_by = "jtg", fail_fast = TRUE)
  write_json(list(promoted = promotion$promoted,
                  snapshot_ids = unique(as.character(promotion$promoted$snapshot_id)),
                  run_id = imported$run$run_id, session_type = "ephemeral_promoted"),
             "feature_set_import_promotion.json")
  registered <- wb_feature_sets(bundle, matrix_id = "RNA_counts",
                                namespace = "parse_cytokine_signatures")
  stopifnot(setequal(as.character(registered$feature_set), expected_sets),
            all(registered$is_readonly == 1L))
  membership_checks <- lapply(expected_sets, function(id) {
    expected <- as.character(accepted$feature_id[accepted$feature_set == id])
    observed <- wb_feature_set_get(bundle, id, "RNA_counts")
    stopifnot(!anyDuplicated(observed), setequal(observed, expected))
    data.frame(feature_set = id, genes = length(observed), members_exact = TRUE,
               order_exact = identical(observed, expected))
  })
  write.table(do.call(rbind, membership_checks), file.path(report_dir, "feature_set_membership_checks.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  stopifnot(identical(sort(as.character(wb_cell_sets(bundle))), sets),
            identical(wb_meta(bundle), source_meta))
  for (id in sets) stopifnot(identical(as.character(wb_cells(bundle, cell_set = id)), memberships[[id]]))
  after <- wb_manifest(bundle)
  for (table in unchanged_tables) {
    stopifnot(identical(before$current_state$tables[[table]], after$current_state$tables[[table]]))
  }
  write_json(list(metadata_exact = TRUE, all_cell_sets_exact = TRUE,
                  unchanged_registry_tables = unchanged_tables, final_gene_sets = 23L,
                  gene_set_source_sha256 = digest::digest(accepted_path, algo = "sha256", file = TRUE),
                  csv_roundtrip_exact = TRUE,
                  csv_source_sha256 = digest::digest(accepted_csv, algo = "sha256", file = TRUE)),
             "scientific_state_checks.json")
  readiness <- wb_finalize_recipient_bundle(bundle)
  write_json(readiness, "recipient_readiness.json")
  stopifnot(isTRUE(readiness$ready))
  wb_close_bundle(bundle)
  closed <- TRUE
  script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  source(file.path(dirname(script_path), "recovery", "restore_materializations.R"))
  restore_parse_materializations(bundle_dir, report_dir,
    file.path(dirname(original_dir), "sessions", "parse_signature_ucell_20260921T152017"))
  source(file.path(dirname(script_path), "recovery", "validate_recovered_bundle.R"))
  finalize_parse_recovery(original_dir, bundle_dir, report_dir)
}
main()
