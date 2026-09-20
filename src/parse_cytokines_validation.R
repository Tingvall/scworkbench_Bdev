# Source-to-bundle validation for the Parse cytokine import.

parse_values_equal <- function(source, imported) {
  if (is.factor(source)) source <- as.character(source)
  if (is.factor(imported)) imported <- as.character(imported)
  identical(unname(source), unname(imported))
}

parse_distributed_indices <- function(n, target = 128L) {
  unique(as.integer(round(seq.int(1L, n, length.out = min(as.integer(target), n)))))
}

validate_parse_cytokines_import <- function(bundle, seurat,
                                             required_metadata = c("sample_id", "replicate", "condition")) {
  cells <- colnames(seurat)
  features <- rownames(seurat)
  source_metadata <- seurat[[]]
  source_counts <- SeuratObject::LayerData(seurat[["RNA"]], layer = "counts")

  stopifnot(
    identical(rownames(source_metadata), cells),
    identical(colnames(source_counts), cells),
    identical(rownames(source_counts), features),
    all(required_metadata %in% names(source_metadata))
  )

  bundle_cells <- scworkbench::wb_cells(bundle, cell_set = "all")
  stopifnot(identical(bundle_cells, cells))

  matrix_registry <- scworkbench::wb_matrices(bundle)
  stopifnot(
    nrow(matrix_registry) == 1L,
    identical(as.character(matrix_registry$matrix_id), "RNA_counts"),
    identical(as.character(matrix_registry$layer), "counts"),
    identical(as.character(matrix_registry$feature_type), "gene"),
    identical(as.integer(matrix_registry$n_features), nrow(source_counts)),
    identical(as.integer(matrix_registry$n_cells), ncol(source_counts))
  )

  imported_counts <- scworkbench::wb_matrix(bundle, "RNA_counts")
  stopifnot(
    identical(dim(imported_counts), dim(source_counts)),
    identical(rownames(imported_counts), features),
    identical(colnames(imported_counts), cells)
  )

  # Exercise every cell and every feature without densifying the complete
  # 78,298 x 251,599 matrix. One comparison spans distributed features over
  # all cells; the other spans all features over distributed cells.
  feature_indices <- parse_distributed_indices(nrow(source_counts))
  sampled_features <- features[feature_indices]
  actual_by_feature <- as.matrix(scworkbench::wb_matrix_read_features(
    bundle, "RNA_counts", feature_ids = sampled_features,
    cell_ids = cells, missing = "NA"
  ))
  expected_by_feature <- as.matrix(source_counts[feature_indices, , drop = FALSE])
  stopifnot(
    identical(dim(actual_by_feature), dim(expected_by_feature)),
    identical(rownames(actual_by_feature), rownames(expected_by_feature)),
    identical(colnames(actual_by_feature), colnames(expected_by_feature)),
    identical(as.numeric(actual_by_feature), as.numeric(expected_by_feature))
  )

  cell_indices <- parse_distributed_indices(ncol(source_counts))
  sampled_cells <- cells[cell_indices]
  actual_by_cell <- as.matrix(scworkbench::wb_matrix_read_features(
    bundle, "RNA_counts", feature_ids = features,
    cell_ids = sampled_cells, missing = "NA"
  ))
  expected_by_cell <- as.matrix(source_counts[, cell_indices, drop = FALSE])
  stopifnot(
    identical(dim(actual_by_cell), dim(expected_by_cell)),
    identical(rownames(actual_by_cell), rownames(expected_by_cell)),
    identical(colnames(actual_by_cell), colnames(expected_by_cell)),
    identical(as.numeric(actual_by_cell), as.numeric(expected_by_cell))
  )

  imported_metadata <- scworkbench::wb_meta(bundle, cells = cells)
  stopifnot(
    identical(imported_metadata$cell_id, cells),
    identical(as.character(imported_metadata$dataset_id), rep("parse_cytokines", length(cells)))
  )
  metadata_checks <- vapply(required_metadata, function(column) {
    column %in% names(imported_metadata) &&
      parse_values_equal(source_metadata[[column]], imported_metadata[[column]])
  }, logical(1L))
  stopifnot(all(metadata_checks))

  metadata_fields <- scworkbench::wb_metadata_fields(bundle, include_core = FALSE)
  stopifnot(
    nrow(metadata_fields) == length(required_metadata),
    setequal(as.character(metadata_fields$source_name), required_metadata)
  )

  source_mapping <- unique(data.frame(
    sample_id = as.character(source_metadata$sample_id),
    condition = as.character(source_metadata$condition),
    replicate = as.character(source_metadata$replicate),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
  imported_mapping <- unique(imported_metadata[, c("sample_id", "condition", "replicate"), drop = FALSE])
  source_mapping <- source_mapping[order(source_mapping$sample_id), , drop = FALSE]
  imported_mapping <- imported_mapping[order(imported_mapping$sample_id), , drop = FALSE]
  rownames(source_mapping) <- NULL
  rownames(imported_mapping) <- NULL
  stopifnot(
    identical(imported_mapping, source_mapping),
    nrow(source_mapping) == 96L,
    !anyDuplicated(source_mapping$sample_id)
  )

  list(
    matrix = list(
      matrix_id = "RNA_counts",
      dimensions = dim(source_counts),
      ordered_cell_ids_equal = TRUE,
      ordered_feature_ids_equal = TRUE,
      sampled_features_across_all_cells = length(feature_indices),
      sampled_cells_across_all_features = length(cell_indices),
      sampled_values_exact = TRUE
    ),
    metadata = list(
      required_columns = required_metadata,
      values_exact = as.list(metadata_checks),
      only_required_noncore_fields_imported = TRUE,
      sample_mappings_exact = TRUE,
      samples = nrow(source_mapping),
      conditions = length(unique(source_mapping$condition)),
      replicates = sort(unique(source_mapping$replicate))
    )
  )
}

