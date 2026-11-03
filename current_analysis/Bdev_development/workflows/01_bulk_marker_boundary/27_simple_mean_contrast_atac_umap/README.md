# Plot simple mean-contrast labels on the Bdev ATAC UMAP

This diagnostic step joins the exploratory raw-UCell decisions from workflow 26
to the existing, checksummed `umap.atac` coordinates. It produces a combined
label map and a four-sample faceted map. The source bundle is not changed and
the labels and figures are not promoted. It also creates a highlight grid and
standalone plots for every assigned label with more than 20 cells; the target
label is colored and every other cell is grey.

```bash
python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/27_simple_mean_contrast_atac_umap/01_plot_atac_umap.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/28_simple_mean_contrast_atac_umap.yaml
```
