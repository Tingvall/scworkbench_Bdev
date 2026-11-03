# Refine ImmGen population programs for default percentile margins

This ImmGen-only exploration starts from the completed global-local fold
rankings. It greedily refines each complete and fold-specific program using the
same percentile scale consumed by the default flat scWorkbench reducer. It
then applies the unchanged score 0.80 and margin 0.10 gates.

```bash
python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/18_percentile_margin_population_genesets/01_refine_percentile_margin.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/19_percentile_margin_population_genesets.yaml
```

No Bdev or scWorkbench bundle is read. Outputs are diagnostic candidates only.
