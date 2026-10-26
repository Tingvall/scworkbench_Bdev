#!/usr/bin/env Rscript
source("src/parse_rest_only_workflow.R")
main <- function() {
  if (file.exists(file.path(report, "preparation.json"))) {
    stopifnot(read_json(file.path(report, "preparation.json"), simplifyVector = TRUE)$status == "complete")
    log_event("Verified existing preparation receipt"); return(invisible(NULL))
  }
  parent <- parent_check()
  b <- wb_open_bundle(parent$bundle_dir, read_only = TRUE, validate = FALSE)
  on.exit(wb_close_bundle(b), add = TRUE)
  cells <- wb_cells(b, cfg$source_cell_set)
  fields <- c("sample_id", "condition", "replicate", "treatment_category")
  md <- wb_meta(b, cells = cells, cols = c(list(wb_ref_metadata_field("core:cell_id")), as.list(fields)))
  stopifnot(length(cells) == 154209L, identical(md$cell_id, cells), !anyNA(md), !anyDuplicated(cells))
  metadata <- md
  metadata$source_bundle <- "parse_cytokines"
  features <- wb_feature_ids(b, "RNA_counts")
  raw <- wb_matrix(b, "RNA_counts")
  raw <- raw[, match(cells, colnames(raw)), drop = FALSE]
  stopifnot(identical(colnames(raw), cells), identical(rownames(raw), features), length(features) == 78298L)
  source_path <- file.path(report, "preparation_source.json")
  if (!file.exists(source_path)) write_json(list(parent = parent, source_cell_set = cfg$source_cell_set,
    matrix = "RNA_counts", cells = length(cells), features = length(features)), "preparation_source.json")
  source_identity <- list(dataset_name = cfg$bundle_id, parent_bundle = parent,
    extraction = "Public scWorkbench getters; exact raw RNA counts and sample metadata for the complete accepted control/single analysis population",
    inputs = list(wb_source_input(path = source_path, source_id = "accepted_parse_parent",
      accession = cfg$run_id, role = "parent_bundle_identity", format = "json", scope = "external_file",
      sha256 = wb_source_sha256(source_path))))
  wb_source_preflight(source_identity)
  derived <- if (dir.exists(cfg$bundle_dir)) wb_open_bundle(cfg$bundle_dir, validate = FALSE) else
    wb_init_bundle(cfg$bundle_dir, bundle_id = cfg$bundle_id,
      bundle_name = "Parse cytokine signatures relative to other individual stimuli: pending review", overwrite = FALSE)
  on.exit(wb_close_bundle(derived), add = TRUE)
  stopifnot(wb_manifest(derived)$bundle_id == cfg$bundle_id)
  if (!length(wb_cells(derived, "all"))) {
    # The public ingest adapter requires a Matrix payload; its atomic importer
    # cannot coerce the backed BPCells S4 object directly. Keep the same counts.
    log_event("Materializing a sparse Matrix payload for the scWorkbench importer")
    sparse <- methods::as(raw, "dgCMatrix")
    stopifnot(identical(dimnames(sparse), dimnames(raw)))
    log_event("Ingesting exact RNA counts; sparse payload GiB", round(as.numeric(object.size(sparse)) / 1024^3, 2))
    ingest <- wb_ingest_apply(derived, set_name = "parse_qc_control_single", cell_ids = cells,
      metadata = metadata, matrices = list(RNA_counts = list(mat = sparse, layer = "counts", feature_type = "gene",
        description = "Exact accepted Parse control/single RNA counts; source parent retained read-only")),
      source = "accepted parse_cytokines", source_identity = source_identity,
      created_by = "jtg", apply_mode = "atomic", emit_report_artifact = TRUE)
    write_json(attr(ingest, "ingest_report"), "ingest_report.json")
    rm(sparse); gc()
  }
  log_event("Verifying sample metadata, matrix identifiers and count aggregates")
  observed_cells <- wb_cells(derived, "all")
  observed_md <- wb_meta(derived, cells = observed_cells, cols = c(list(wb_ref_metadata_field("core:cell_id")), as.list(fields)))
  stopifnot(identical(cells, observed_cells), identical(features, wb_feature_ids(derived, "RNA_counts")))
  for (field in c("cell_id", fields)) stopifnot(identical(as.character(md[[field]]), as.character(observed_md[[field]])))
  observed <- wb_matrix(derived, "RNA_counts")
  stopifnot(identical(dim(raw), dim(observed)), identical(rownames(raw), rownames(observed)),
    identical(colnames(raw), colnames(observed)), identical(as.numeric(colSums(raw)), as.numeric(colSums(observed))))
  source_pb_path <- file.path(parent_cfg$signature_analysis$discovery$report, "pseudobulk_counts.rds")
  pb <- readRDS(source_pb_path)
  sums <- BPCells::pseudobulk_matrix(observed, cell_groups = factor(md$sample_id, levels = pb$sample_metadata$sample_id),
    method = "sum", threads = 4L)
  sums <- as.matrix(sums)
  stopifnot(identical(dim(sums), dim(pb$counts)), identical(rownames(sums), rownames(pb$counts)),
    identical(colnames(sums), colnames(pb$counts)), identical(as.numeric(sums), as.numeric(pb$counts)))
  stopifnot(nrow(wb_feature_sets(derived, namespace = cfg$namespace)) == 0L)
  write_tsv(data.frame(cell_id = cells, sample_id = md$sample_id, condition = md$condition,
    replicate = md$replicate, treatment_category = md$treatment_category), "analysis_cell_metadata.tsv", work)
  save_provenance("preparation", c(source_path, source_pb_path, "processing/parse_cytokines/12_prepare_rest_only_bundle.R"),
    list(bundle_state_fingerprint = wb_manifest(derived)$current_state$fingerprint,
      cells = length(cells), features = length(features), samples = length(unique(md$sample_id)),
      ordered_cell_and_feature_ids_exact = TRUE, sample_metadata_exact = TRUE,
      every_cell_total_exact = TRUE, all_gene_by_sample_pseudobulk_counts_exact = TRUE,
      individual_nonzero_entries_compared = FALSE, signal_namespace_registered = FALSE))
  log_event("Separate review bundle prepared")
}
main()
