# Positive competitor-refinement pilot

This diagnostic pilot adapts the strongest-competitor idea from the archived
Parse cytokine analysis to ImmGen while permitting positive genes only. It
does not modify the approved marker contract, refit the production models, or
create an importable gene-set collection.

Run from the repository root:

```bash
python current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/09_positive_competitor_refinement_pilot/01_run_positive_competitor_refinement.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/09_positive_competitor_refinement_pilot.yaml
```

All generated diagnostic candidates and reports are written to paired active
storage under the configured workflow root.

