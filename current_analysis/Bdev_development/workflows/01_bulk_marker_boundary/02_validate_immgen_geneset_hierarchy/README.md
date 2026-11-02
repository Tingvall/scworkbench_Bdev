# Validate the ImmGen gene-set hierarchy

This workflow checks the approved hierarchy and population roles against the
frozen vM25 library manifest. It verifies population partitions, lineage
counts, retained unique branches, broad-validation roles, the B-lineage branch
partition, symmetric local-stage neighborhoods, and stable edge and contrast
identifiers.

It performs no expression analysis.

```bash
python current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/02_validate_immgen_geneset_hierarchy/01_validate_immgen_geneset_hierarchy.py \
  --manifest current_analysis/Bdev_development/inputs/01_bulk_marker_boundary/02_immgen_vm25_library_manifest.tsv \
  --hierarchy current_analysis/Bdev_development/config/01_bulk_marker_boundary/01_immgen_geneset_hierarchy.yaml
```

The accepted validation environment used Python 3.10.12 and PyYAML 6.0.2.
