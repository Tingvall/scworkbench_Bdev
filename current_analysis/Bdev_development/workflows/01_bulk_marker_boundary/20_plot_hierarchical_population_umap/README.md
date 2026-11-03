# Plot hierarchical ImmGen calls on the Bdev ATAC UMAP

This diagnostic plotting step reads the existing `umap.atac` coordinates from
the qualified development bundle in read-only mode and joins them to the three
step-19 policy decisions. It creates a status comparison and a direct-labeled
population-calibrated UMAP. No labels or plots are registered in scWorkbench.

```bash
/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --bind /mnt/sata2:/mnt/sata2 \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif \
  Rscript current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/20_plot_hierarchical_population_umap/01_export_umap.R \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/21_hierarchical_population_umap.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/20_plot_hierarchical_population_umap/02_plot_umap.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/21_hierarchical_population_umap.yaml
```
