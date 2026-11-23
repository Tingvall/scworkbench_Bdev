# Combined supported up-signature local-GSEA UMAPs

This post-qualification visualization consumes the checksummed outputs from
step 02 and combines the ten `up` signatures having at least one supported
local-GSEA interval into a 2-by-5 ATAC-UMAP grid. All panels share one NES
scale. The finalized step-01 and step-02 result roots remain unchanged.

Run from the repository root using the same qualified development SIF and R
library chain recorded in the configuration:

```bash
<apptainer> exec --cleanenv --bind /mnt/sata2:/mnt/sata2 \
  --env 'R_LIBS_USER=<03c-run-library>:<dependency-library>:<shared-library>' \
  <development-sif> Rscript \
  current_analysis/Bdev_development/workflows/03c_parse_cytokine_depth_adjusted_trajectory/03_combine_supported_up_geneset_gsea_umaps/01_plot_supported_up_geneset_gsea_umap_grid.R \
  --config current_analysis/Bdev_development/config/03c_parse_cytokine_depth_adjusted_trajectory/03_supported_up_geneset_gsea_umap_grid.yaml

python3 current_analysis/Bdev_development/workflows/03c_parse_cytokine_depth_adjusted_trajectory/03_combine_supported_up_geneset_gsea_umaps/02_finalize_run.py \
  --config current_analysis/Bdev_development/config/03c_parse_cytokine_depth_adjusted_trajectory/03_supported_up_geneset_gsea_umap_grid.yaml
```

Outputs are bounded development visualizations and are not promoted to
`Bdev_implementation`.
