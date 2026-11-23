# Depth-adjusted Parse38 cytokine trajectory qualification

This bounded development workflow mirrors 03b while replacing raw UCell-based
score analyses with reusable cell-level scores constructed from an authenticated
`RNA_log1p_cp10k` matrix. Each gene is adjusted for `nCount_RNA` while
preserving fitted native L4 structure with `ns(L4, df=8)`, standardized, and
averaged with equal gene weight. Dynamics and interval detection consume the
same adjusted score. The interval operation must report `already_adjusted` so
that depth is not removed twice. GSEA remains based on `SCT_data` with
`adjust_depth=FALSE`.

The analysis reuses the exact checksummed fitted-L4 pseudotime, model, and
deterministic fixture produced by 03b; it does not refit Slingshot. It uses
matched scales: 20 dynamics bins and 5%-wide windows at 2.5% steps for
intervals and GSEA. Analysis uses native fitted L4. Figures reverse only the
display direction to LSK (0) -> CD19 (1).

The isolated run library receives the pinned `scworkbench` core before the
pinned operation packs because `scwbTrajectory` 0.1.8 uses the core's
authenticated parameter-derived key-axis grant.

Run from the repository root:

```bash
python3 current_analysis/Bdev_development/workflows/03c_parse_cytokine_depth_adjusted_trajectory/01_qualify_depth_adjusted_cytokine_trajectory/01_prepare_execution_clone.py \
  --config current_analysis/Bdev_development/config/03c_parse_cytokine_depth_adjusted_trajectory/01_bdev_cytokine_depth_adjusted_qualification.yaml

python3 current_analysis/Bdev_development/workflows/03c_parse_cytokine_depth_adjusted_trajectory/01_qualify_depth_adjusted_cytokine_trajectory/02_install_operation_packs.py \
  --config current_analysis/Bdev_development/config/03c_parse_cytokine_depth_adjusted_trajectory/01_bdev_cytokine_depth_adjusted_qualification.yaml

/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --cleanenv --bind /mnt/sata2:/mnt/sata2 \
  --env 'R_LIBS_USER=<run-library>:<dependency-library>:<shared-development-library>' \
  <development-sif> Rscript \
  current_analysis/Bdev_development/workflows/03c_parse_cytokine_depth_adjusted_trajectory/01_qualify_depth_adjusted_cytokine_trajectory/03_run_depth_adjusted_cytokine_trajectory.R \
  --config current_analysis/Bdev_development/config/03c_parse_cytokine_depth_adjusted_trajectory/01_bdev_cytokine_depth_adjusted_qualification.yaml

/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --cleanenv --bind /mnt/sata2:/mnt/sata2 \
  --env 'R_LIBS_USER=<run-library>:<dependency-library>:<shared-development-library>' \
  <development-sif> Rscript \
  current_analysis/Bdev_development/workflows/03c_parse_cytokine_depth_adjusted_trajectory/01_qualify_depth_adjusted_cytokine_trajectory/04_plot_depth_adjusted_methods.R \
  --config current_analysis/Bdev_development/config/03c_parse_cytokine_depth_adjusted_trajectory/01_bdev_cytokine_depth_adjusted_qualification.yaml

python3 current_analysis/Bdev_development/workflows/03c_parse_cytokine_depth_adjusted_trajectory/01_qualify_depth_adjusted_cytokine_trajectory/05_finalize_run.py \
  --config current_analysis/Bdev_development/config/03c_parse_cytokine_depth_adjusted_trajectory/01_bdev_cytokine_depth_adjusted_qualification.yaml
```

Generated state mirrors this workflow under paired development storage. This is
not a paper analysis and does not unlock `Bdev_implementation`.
