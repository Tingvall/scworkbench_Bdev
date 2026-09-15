# FLT3KO/IL7KO multiome workflow

This workflow is planned. Add numbered scripts only after delivered files and
their RNA/ATAC, fragment, genome, genotype, sample, and replicate identities
have been inventoried.

Expected phases are:

1. `01_inspect_source.R`
2. `02_create_bundle.R`
3. `03_qc_preprocess.R`
4. `04_dataset_analysis.R`
5. `05_validate_bundle.R`

The independently validated output is the `scmulti_flt3ko_il7ko` bundle. It is
the query bundle for later mapping to the Bdev reference.
