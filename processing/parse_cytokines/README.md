# Parse cytokine workflow

This workflow is planned. Add numbered scripts only after the source inspection
has established the delivered matrices, sample/replicate/condition metadata,
cell identifiers, checksums, and supported scWorkbench operations.

Expected phases are:

1. `01_inspect_source.R`
2. `02_create_bundle.R`
3. `03_qc_preprocess.R`
4. `04_dataset_analysis.R`
5. `05_validate_bundle.R`

The independently validated output is the `parse_cytokines` bundle. Reviewed
cytokine gene sets may later be consumed by `processing/integration/`.
