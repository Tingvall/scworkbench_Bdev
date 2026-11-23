# Simple mean-contrast margin comparison

This diagnostic compares three reductions of the same authenticated 64-program
raw UCell score family: pure argmax, winner margin at least 0.01, and winner
margin at least 0.05. Exact top-score ties remain unassigned under all rules.

The first script reads the authenticated score payload without changing its
session or source bundle, reproduces the recorded 0.01 decisions exactly, and
exports all three decisions. The second script creates the three-panel ATAC
UMAP, a confidence-stratum UMAP, summary plots, tables, and the final manifest.

```bash
Rscript current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/30_simple_mean_contrast_margin_comparison/01_export_margin_decisions.R \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/31_simple_mean_contrast_margin_comparison.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/30_simple_mean_contrast_margin_comparison/02_plot_margin_comparison.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/31_simple_mean_contrast_margin_comparison.yaml
```

Outputs are unpromoted development diagnostics and do not replace the selected
0.01 provisional default.
