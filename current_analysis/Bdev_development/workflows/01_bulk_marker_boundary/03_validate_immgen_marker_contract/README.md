# Validate the ImmGen marker statistical contract

This workflow validates the approved statistical contract against the frozen
vM25 manifest and approved hierarchy. It checks discovery and held-out counts,
contrast families, thresholds, stability settings, deterministic selection,
and the outside-scWorkbench boundary.

It performs no normalization, model fitting, differential expression, or
gene-set selection.

```bash
python current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/03_validate_immgen_marker_contract/01_validate_immgen_marker_contract.py \
  --manifest current_analysis/Bdev_development/inputs/01_bulk_marker_boundary/02_immgen_vm25_library_manifest.tsv \
  --hierarchy current_analysis/Bdev_development/config/01_bulk_marker_boundary/01_immgen_geneset_hierarchy.yaml \
  --contract current_analysis/Bdev_development/config/01_bulk_marker_boundary/02_immgen_marker_statistical_contract.yaml
```

The accepted structural validation environment used Python 3.10.12 and PyYAML
6.0.2.
