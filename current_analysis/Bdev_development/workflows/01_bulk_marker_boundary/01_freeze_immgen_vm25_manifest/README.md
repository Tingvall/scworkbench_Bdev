# Freeze the ImmGen GENCODE vM25 library manifest

This workflow validates the checksummed GSE109125 GENCODE vM25 count matrix,
maps every matrix column to a GEO RNA sample, reconciles populations to paper
Table S1, attaches the population-matched Table S1 ATAC GSM accessions, and
applies the approved 64-population discovery allowlist.

It does not perform normalization, differential expression, marker selection,
or any scWorkbench operation.

## Dependencies

The accepted manifest was generated with Python 3.10.12 and `openpyxl` 3.0.9.
All other imported Python modules are from the standard library.

## Inputs

Large source files remain outside Git. Their URLs and required SHA-256 hashes
are frozen in
`../../../inputs/01_bulk_marker_boundary/01_immgen_vm25_manifest_policy.json`.

## Reproduction

From the repository root, with the five source files downloaded to the shown
temporary paths:

```bash
python current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/01_freeze_immgen_vm25_manifest/01_build_immgen_vm25_manifest.py \
  --policy current_analysis/Bdev_development/inputs/01_bulk_marker_boundary/01_immgen_vm25_manifest_policy.json \
  --vm25-counts /tmp/GSE109125_Gene_count_table_GENCODE_vM25.csv.gz \
  --vm25-readme /tmp/GSE109125_readme_GENCODE_vM25.txt \
  --table-s1 /tmp/S0092867418316507-mmc1.xlsx \
  --rna-soft /tmp/GSE109125_samples.soft \
  --atac-soft /tmp/GSE100738_samples.soft \
  --output-dir current_analysis/Bdev_development/inputs/01_bulk_marker_boundary
```

The script fails if a source checksum, expected dimension, mapping cardinality,
count-integrity check, cohort size, or matched-ATAC requirement changes.
