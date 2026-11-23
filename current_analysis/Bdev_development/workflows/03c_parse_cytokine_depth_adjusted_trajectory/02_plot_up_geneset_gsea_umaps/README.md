# Up-signature significant-GSEA UMAPs

This post-qualification visualization step leaves the finalized parent 03c run
unchanged. It creates one ATAC-UMAP for each of the 21 Parse38 `up` signatures.
All Bdev cells are gray context. Only authenticated L4 fixture cells within a
supported merged local-GSEA interval are colored. Each colored cell receives
the NES of the nearest constituent supported-window center. Every map uses the
same symmetric blue-white-red NES scale.

Signatures without supported local-GSEA intervals still receive a map that
explicitly states that no region was supported.

Run from the repository root:

```bash
/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --cleanenv --bind /mnt/sata2:/mnt/sata2 \
  --env 'R_LIBS_USER=<03c-run-library>:<dependency-library>:<shared-development-library>' \
  <development-sif> Rscript \
  current_analysis/Bdev_development/workflows/03c_parse_cytokine_depth_adjusted_trajectory/02_plot_up_geneset_gsea_umaps/01_plot_up_geneset_gsea_umaps.R \
  --config current_analysis/Bdev_development/config/03c_parse_cytokine_depth_adjusted_trajectory/02_up_geneset_gsea_umaps.yaml

python3 current_analysis/Bdev_development/workflows/03c_parse_cytokine_depth_adjusted_trajectory/02_plot_up_geneset_gsea_umaps/02_finalize_run.py \
  --config current_analysis/Bdev_development/config/03c_parse_cytokine_depth_adjusted_trajectory/02_up_geneset_gsea_umaps.yaml
```

Generated results mirror this workflow under paired development storage. They
are development visualizations, are not promoted, and do not unlock
`Bdev_implementation`.
