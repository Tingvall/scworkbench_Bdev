# Validation of imported source values through scWorkbench public getters.
reference_values_equal <- function(a, b) {
  if (is.factor(a)) a <- as.character(a)
  if (is.factor(b)) b <- as.character(b)
  if (is.numeric(a) && is.numeric(b)) {
    return(identical(as.numeric(a), as.numeric(b)))
  }
  identical(unname(a), unname(b))
}

validate_reference_import <- function(bundle, seurat, mapping, reduction_names) {
  cells <- colnames(seurat)
  stopifnot(setequal(scworkbench::wb_cells(bundle, cell_set = "all"), cells))
  matrix_checks <- list()
  for (spec in mapping) {
    original <- SeuratObject::LayerData(seurat[[spec$assay]], layer = spec$layer)
    imported <- scworkbench::wb_matrix(bundle, spec$matrix_id)
    stopifnot(identical(dim(imported), dim(original)))
    stopifnot(identical(rownames(imported), rownames(original)))
    stopifnot(identical(colnames(imported), colnames(original)))
    # Compare distributed features across every cell without densifying a full assay.
    feature_idx <- unique(as.integer(round(seq(1L, nrow(original), length.out = min(128L, nrow(original))))))
    features <- rownames(original)[feature_idx]
    actual <- as.matrix(scworkbench::wb_matrix_read_features(
      bundle, spec$matrix_id, feature_ids = features, cell_ids = cells, missing = "NA"
    ))
    expected <- as.matrix(original[feature_idx, , drop = FALSE])
    stopifnot(identical(dim(actual), dim(expected)))
    stopifnot(identical(as.numeric(actual), as.numeric(expected)))
    matrix_checks[[spec$matrix_id]] <- list(
      dimensions = dim(original), ordered_cell_ids_equal = TRUE,
      ordered_feature_ids_equal = TRUE, sampled_features = length(features),
      sampled_cells = length(cells), sampled_values_exact = TRUE
    )
  }
  metadata <- scworkbench::wb_meta(bundle, cells = cells)
  original_metadata <- seurat[[]]
  stopifnot(identical(metadata$cell_id, cells))
  stopifnot(all(names(original_metadata) %in% names(metadata)))
  metadata_checks <- vapply(names(original_metadata), function(nm) {
    reference_values_equal(original_metadata[[nm]], metadata[[nm]])
  }, logical(1))
  stopifnot(all(metadata_checks))
  reduction_checks <- list()
  for (nm in reduction_names) {
    original <- SeuratObject::Embeddings(seurat[[nm]])
    imported <- scworkbench::wb_space(bundle, nm, cells = cells)
    if (is.data.frame(imported)) {
      stopifnot(identical(imported$cell_id, cells))
      imported <- as.matrix(imported[, c("x", "y"), drop = FALSE])
    } else {
      stopifnot(identical(rownames(imported), cells))
    }
    stopifnot(identical(dim(imported), dim(original)))
    stopifnot(identical(as.numeric(imported), as.numeric(original)))
    reduction_checks[[nm]] <- list(dimensions = dim(original), all_coordinates_exact = TRUE)
  }
  list(matrices = matrix_checks, metadata = as.list(metadata_checks), reductions = reduction_checks)
}

validate_reference_pseudotime <- function(bundle, artifact_id, source_metadata) {
  cells <- rownames(source_metadata)
  expected <- as.numeric(source_metadata$pt)
  metadata <- scworkbench::wb_meta(bundle, cells = cells)
  stopifnot(identical(as.numeric(metadata$pt), expected))
  primitive <- scworkbench::wb_art_load(bundle, artifact_id, expected_kind = "cell_vector")
  stopifnot(identical(primitive$value_type, "continuous"))
  stopifnot(identical(primitive$cell_id, cells[!is.na(expected)]))
  retained <- rep(NA_real_, length(cells))
  retained[match(primitive$cell_id, cells)] <- as.numeric(primitive$value)
  reader <- scworkbench::wb_read_cell_vector(bundle, paste0("artifact:", artifact_id), cells)
  stopifnot(identical(names(reader), cells))
  stopifnot(identical(as.numeric(reader), retained))
  # This image's wb_prim_cell_vector() uses as.character(x), which rounds
  # numeric values to R's character representation. Require that exact known
  # conversion, and retain the original doubles in metadata as authoritative.
  stopifnot(identical(retained, as.numeric(as.character(expected))))
  stopifnot(identical(is.na(retained), is.na(expected)))
  stopifnot(identical(rank(retained, na.last = "keep"), rank(expected, na.last = "keep")))
  list(source_field = "pt", artifact_id = artifact_id, n_cells = length(cells),
       n_defined = sum(!is.na(retained)), n_missing = sum(is.na(retained)),
       metadata_values_exact = TRUE, missingness_equal = TRUE, pseudotime_ranks_equal = TRUE,
       artifact_values_exact = identical(retained, expected),
       artifact_matches_scworkbench_numeric_character_conversion = TRUE,
       artifact_max_abs_difference = max(abs(retained - expected), na.rm = TRUE),
       precision_note = "Original pt doubles remain exact in metadata. scWorkbench cell_vector primitives use as.character(numeric), introducing the quantified rounding in this derived artifact.")
}

finish_reference_import <- function(bundle, source_metadata, mapping, reduction_names,
                                    source_path, source_sha, source_before, report_dir,
                                    artifact_id, promotion = NULL, run_id = NULL) {
  bundle_dir <- bundle$dir
  write_report <- function(value, name) {
    jsonlite::write_json(value, file.path(report_dir, name), pretty = TRUE,
                         auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  }
  pt_check <- validate_reference_pseudotime(bundle, artifact_id, source_metadata)
  write_report(c(pt_check, list(run_id = run_id, promotion = promotion)), "pseudotime.json")
  certificates <- list()
  for (mid in names(mapping)) certificates[[paste0("matrix:", mid)]] <- scworkbench::wb_entity_certificate(bundle, "matrix", mid)
  for (sid in reduction_names) certificates[[paste0("space:", sid)]] <- scworkbench::wb_entity_certificate(bundle, "space", sid)
  certificates$pseudotime <- scworkbench::wb_entity_certificate(bundle, "artifact", artifact_id)
  write_report(certificates, "entity_certificates.json")
  write_report(scworkbench::wb_manifest(bundle), "manifest.json")
  capture.output(sessionInfo(), file = file.path(report_dir, "session_info.txt"))
  cat(format(Sys.time()), "Finalizing bundle for read-only integrity validation\n")
  readiness <- scworkbench::wb_finalize_recipient_bundle(bundle)
  write_report(readiness, "recipient_readiness.json")
  scworkbench::wb_close_bundle(bundle)
  validation <- scworkbench::wb_bundle_validation_report(bundle_dir, level = "recipient", progress = TRUE)
  write_report(validation, "validation.json")
  stopifnot(isTRUE(validation$ok))
  stopifnot(identical(source_before, file.info(source_path)[, c("size", "mtime", "ctime")]))
  write_report(list(status = "complete", bundle_name = "bdev_reference", bundle_dir = bundle_dir,
                    source_sha256 = source_sha, source_stat_unchanged = TRUE,
                    matrices = names(mapping), reductions = reduction_names,
                    metadata_columns = names(source_metadata), n_cells = nrow(source_metadata),
                    pseudotime_artifact_id = artifact_id, pseudotime_precision = pt_check,
                    validation_status = validation$status), "completion.json")
  cat(format(Sys.time()), "Complete:", bundle_dir, "\n")
  invisible(validation)
}
