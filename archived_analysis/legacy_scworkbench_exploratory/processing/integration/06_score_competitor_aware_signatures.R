#!/usr/bin/env Rscript

# Transfer the 33 reviewed competitor-aware full signatures into an isolated
# Bdev integration session, report Bdev coverage, and calculate native UCell
# scores for all 33 sets. No output is promoted.

source("src/competitor_aware_integration_workflow.R")

run_or_reuse <- function(session, bundle, label, op, inputs, notes) {
  marked <- paste0("[", label, "] ")
  existing <- names(session$runs)[vapply(session$runs, function(x) {
    identical(x$op_name, op) && startsWith(as.character(x$notes %||% ""), marked)
  }, logical(1))]
  stopifnot(length(existing) <= 1L)
  if (length(existing)) {
    log_event("Reusing completed checkpoint", label)
    list(session = session, run = session$runs[[existing[[1L]]]])
  } else {
    log_event("Running", label, op)
    wb_session_run_op(session, bundle, op, inputs = inputs,
      notes = paste0(marked, notes))
  }
}

main <- function() {
  completion_path <- file.path(report_dir, "scoring.json")
  if (file.exists(completion_path)) {
    stopifnot(identical(read_json("scoring.json")$status, "complete"))
    log_event("Completed scoring report already exists; leaving it unchanged")
    return(invisible(NULL))
  }

  parents <- check_competitor_aware_parents()
  source <- read_source_memberships()
  sets <- source$sets
  passing <- source$summary

  ref <- open_parent("bdev_reference")
  on.exit(try(wb_close_bundle(ref), silent = TRUE), add = TRUE)
  parse <- open_parent("parse_cytokines_rest_only")
  on.exit(try(wb_close_bundle(parse), silent = TRUE), add = TRUE)
  bundle <- wb_open_bundle(bundle_dir, validate = FALSE)
  on.exit(try(wb_close_bundle(bundle), silent = TRUE), add = TRUE)

  manifest_before <- wb_manifest(bundle)
  cells <- wb_cells(bundle, "all")
  ref_cells <- wb_cells(ref, "all")
  metadata <- wb_meta(bundle, cells = cells)
  ref_metadata <- wb_meta(ref, cells = ref_cells)
  stopifnot(identical(cells, ref_cells), length(cells) == 26723L,
    identical(as.numeric(metadata$pt), as.numeric(ref_metadata$pt)),
    identical(is.na(metadata$pt), is.na(ref_metadata$pt)),
    sum(is.finite(metadata$pt)) == 12642L, sum(is.na(metadata$pt)) == 14081L)

  # Confirm that the exported source table is exactly the registered accepted
  # membership before using it as the cross-bundle transfer input.
  for (id in unique(sets$feature_set)) {
    expected <- sets$feature_id[sets$feature_set == id]
    observed <- as.character(wb_feature_set_get(parse, id, "RNA_counts"))
    stopifnot(identical(observed, expected))
  }

  features <- wb_feature_ids(bundle, "RNA_counts")
  coverage <- do.call(rbind, lapply(unique(sets$feature_set), function(id) {
    genes <- sets$feature_id[sets$feature_set == id]
    row <- passing[match(id, passing$signature_id), , drop = FALSE]
    mapped <- genes[genes %in% features]
    data.frame(feature_set = id, cytokine = row$cytokine,
      direction = row$direction, selected_size = row$signature_size,
      original_genes = length(genes), mapped_genes = length(mapped),
      coverage = length(mapped) / length(genes),
      coverage_gate_pass = length(mapped) >= competitor_cfg$minimum_genes &&
        length(mapped) / length(genes) >= competitor_cfg$minimum_feature_coverage,
      included = TRUE,
      missing_genes = paste(setdiff(genes, features), collapse = ";"),
      stringsAsFactors = FALSE)
  }))
  stopifnot(nrow(coverage) == 33L, !anyDuplicated(coverage$feature_set))
  write_tsv(coverage, "feature_coverage.tsv")
  write_tsv(sets, "accepted_source_gene_sets.tsv")
  write_tsv(passing, "source_signature_summary.tsv")
  saveRDS(pt_family(bundle), file.path(report_dir, "exact_pt_family.rds"))

  for (op in c("feature_set_import", "feature_set_overlap_qc", "ucell",
    "cell_vector_family_qc")) {
    write_json(wb_op_interface(op, bundle = bundle),
      paste0("interface_", op, ".json"))
  }

  input_csv <- file.path(report_dir, "accepted_source_gene_sets.csv")
  write.csv(sets, input_csv, row.names = FALSE)
  stopifnot(identical(read.csv(input_csv, check.names = FALSE), sets))

  workspace <- file.path(session_root, "competitor_aware_scoring")
  session <- if (dir.exists(workspace)) {
    wb_session_resume(bundle, workspace)
  } else {
    wb_session_create(bundle, workspace_root = session_root,
      session_id = "competitor_aware_scoring")
  }
  if (!is.null(session$durable$promotion_receipt)) {
    stop("The review session was unexpectedly promoted; refusing to continue")
  }

  imported <- run_or_reuse(session, bundle, "source_import", "feature_set_import",
    list(matrix = "RNA_counts", cell_set = "all", min_cells = 10L,
      source_path = input_csv, source_format = "table",
      namespace = competitor_cfg$source_namespace,
      collection = "competitor_aware_full", species = "Mus musculus",
      organism = "Mus musculus", gene_id_type = "bundle_feature_id",
      set_id_col = "feature_set", set_name_col = "set_name",
      gene_col = "feature_id", description_col = "description",
      split_by_collection = FALSE, overwrite = FALSE, strict = FALSE,
      created_by = "scworkbench_Bdev integration review",
      source_id = "parse_cytokines_rest_only_competitor_aware",
      source_version = "2026-09-25", origin = "derived", is_readonly = TRUE),
    paste("Exact accepted memberships from parse_cytokines_rest_only;",
      "genes absent from Bdev RNA are retained in the source table and dropped from the staged import."))
  session <- imported$session
  import_manifest <- imported$run$outputs$import_manifest$table
  write_tsv(import_manifest, "signature_import_manifest.tsv")
  stopifnot(nrow(import_manifest) == 33L,
    setequal(import_manifest$feature_set, coverage$feature_set))

  qc <- run_or_reuse(session, bundle, "overlap_qc", "feature_set_overlap_qc",
    list(matrix = "RNA_counts", cell_set = "all",
      namespace = competitor_cfg$source_namespace,
      min_set_size = as.integer(competitor_cfg$minimum_genes),
      min_detected_fraction = competitor_cfg$minimum_feature_coverage,
      deduplicate_sets = FALSE),
    "Apply the prespecified minimum five mapped genes and 80 percent source-membership coverage gate.")
  session <- qc$session
  qc_table <- qc$run$outputs$overlap_qc$table
  write_tsv(qc_table, "signature_overlap_qc.tsv")
  write_tsv(qc$run$outputs$namespace_summary$table,
    "signature_overlap_namespace_summary.tsv")
  stopifnot(setequal(coverage$feature_set[coverage$coverage_gate_pass],
    qc_table$feature_set[qc_table$screening_status == "include"]))
  selected <- coverage$feature_set[coverage$included]

  signatures <- setNames(lapply(selected, function(id) {
    sets$feature_id[sets$feature_set == id]
  }), selected)
  scored <- run_or_reuse(session, bundle, "ucell_scores_all_33", "ucell",
    list(matrix = "RNA_counts", cell_set = "all", min_cells = 10L,
      features = signatures, name = "parse_cytokine_rest_only",
      max_rank = as.integer(competitor_cfg$max_rank), ncores = 4L,
      parallel_backend = "auto", chunk_size = 400L,
      missing_features_policy = "drop"),
    paste("Raw-count UCell scoring of all 33 passing competitor-aware signatures on all Bdev cells;",
      "missing Bdev features are dropped and coverage is reported, but is not an exclusion gate per user instruction;",
      "max_rank 1500, no graph smoothing, and no cross-signature calibration."))
  session <- scored$session
  scores <- scored$run$outputs$signature_scores
  score_summary <- scored$run$outputs$signature_summary$table
  stopifnot(length(scores$value) == length(cells) * length(selected),
    setequal(unique(scores$key), selected), all(is.finite(scores$value)),
    all(scores$value >= 0 & scores$value <= 1))
  for (id in selected) {
    stopifnot(identical(scores$cell_id[scores$key == id], cells))
  }
  saveRDS(scores, file.path(report_dir, "ucell_scores.rds"), compress = "gzip")
  write_tsv(score_summary, "ucell_summary.tsv")
  write_json(wb_session_dag(session, bundle = bundle, include_artifacts = TRUE,
    include_active = TRUE), "native_session_dag.json")

  # Writable session setup may republish manifest bookkeeping, while the
  # registry-state fingerprint must remain exactly unchanged without promotion.
  manifest_after <- wb_manifest(bundle)
  stopifnot(identical(manifest_after$current_state$fingerprint,
    manifest_before$current_state$fingerprint))
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
  write_json(list(status = "complete", promoted = FALSE,
    source_bundle = "parse_cytokines_rest_only",
    source_bundle_fingerprint = parents$parse_cytokines_rest_only$bundle_state_fingerprint,
    reference_bundle = "bdev_reference",
    reference_bundle_fingerprint = parents$bdev_reference$bundle_state_fingerprint,
    integration_bundle_fingerprint = manifest_before$current_state$fingerprint,
    source_signatures = nrow(coverage), scored_signatures = length(selected),
    below_prior_coverage_gate = coverage$feature_set[!coverage$coverage_gate_pass],
    coverage_policy = "diagnostic_only_all_33_passing_sets_retained",
    cells_scored = length(cells), trajectory_cells = sum(is.finite(metadata$pt)),
    missing_pt = sum(is.na(metadata$pt)), pseudotime_exact = TRUE,
    score_run_id = scored$run$run_id, session_workspace = workspace,
    source_membership_sha256 = wb_source_sha256(source_membership_tsv),
    source_signature_summary_sha256 = wb_source_sha256(source_signature_summary_tsv),
    bundle_manifest_unchanged = TRUE), "scoring.json")
  log_event("Competitor-aware scoring complete: all", length(selected),
    "passing Parse signatures retained")
}

main()
