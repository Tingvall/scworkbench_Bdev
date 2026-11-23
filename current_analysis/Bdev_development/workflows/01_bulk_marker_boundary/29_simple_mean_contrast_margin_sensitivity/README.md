# Simple mean-contrast raw-margin sensitivity

This diagnostic rerun applies a raw UCell winner–runner-up margin of 0.05 to
the same 64-program scores and 1,746-cell fixture used by the selected 0.01
development default. It does not rescore cells or change the selected default.

The workflow produces a 0.05 label map, a direct 0.01-versus-0.05 comparison,
an assignment-loss map, and cell-, label-, and sample-level comparison tables.
All outputs remain unpromoted development diagnostics.

```bash
python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/29_simple_mean_contrast_margin_sensitivity/01_compare_margin005_atac_umap.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/30_simple_mean_contrast_margin_sensitivity.yaml
```
