# Bounded Parse38 cytokine trajectory qualification

This development workflow qualifies all three selected cytokine trajectory
methods on a deterministic, 20-bin Bdev L4 fixture. It is not a paper analysis
and promotes no result artifacts.

The common scale is 20 non-overlapping dynamics bins and 5%-wide overlapping
interval/GSEA windows with a 2.5% step. Signal intervals and GSEA use unblocked
permutations. `orig.ident` is retained only in the bin-composition audit.

Run from the repository root:

```bash
python3 current_analysis/Bdev_development/workflows/03b_parse_cytokine_trajectory_analysis/01_qualify_cytokine_trajectory_methods/01_prepare_execution_clone.py \
  --config current_analysis/Bdev_development/config/03b_parse_cytokine_trajectory_analysis/02_bdev_cytokine_trajectory_qualification.yaml

/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --cleanenv --bind /mnt/sata2:/mnt/sata2 \
  --env 'R_LIBS=<run-library>:<shared-development-library>' \
  <development-sif> Rscript \
  current_analysis/Bdev_development/workflows/03b_parse_cytokine_trajectory_analysis/01_qualify_cytokine_trajectory_methods/02_run_cytokine_trajectory.R \
  --config current_analysis/Bdev_development/config/03b_parse_cytokine_trajectory_analysis/02_bdev_cytokine_trajectory_qualification.yaml

python3 current_analysis/Bdev_development/workflows/03b_parse_cytokine_trajectory_analysis/01_qualify_cytokine_trajectory_methods/03_finalize_manifest.py \
  --config current_analysis/Bdev_development/config/03b_parse_cytokine_trajectory_analysis/02_bdev_cytokine_trajectory_qualification.yaml

/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --cleanenv --bind /mnt/sata2:/mnt/sata2 \
  <development-sif> Rscript \
  current_analysis/Bdev_development/workflows/03b_parse_cytokine_trajectory_analysis/01_qualify_cytokine_trajectory_methods/04_plot_cytokine_trajectory_methods.R \
  --config current_analysis/Bdev_development/config/03b_parse_cytokine_trajectory_analysis/02_bdev_cytokine_trajectory_qualification.yaml

python3 current_analysis/Bdev_development/workflows/03b_parse_cytokine_trajectory_analysis/01_qualify_cytokine_trajectory_methods/05_finalize_figures.py \
  --config current_analysis/Bdev_development/config/03b_parse_cytokine_trajectory_analysis/02_bdev_cytokine_trajectory_qualification.yaml

/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --cleanenv --bind /mnt/sata2:/mnt/sata2 \
  <development-sif> Rscript \
  current_analysis/Bdev_development/workflows/03b_parse_cytokine_trajectory_analysis/01_qualify_cytokine_trajectory_methods/06_plot_flt3l_il7_four_methods.R \
  --config current_analysis/Bdev_development/config/03b_parse_cytokine_trajectory_analysis/02_bdev_cytokine_trajectory_qualification.yaml

python3 current_analysis/Bdev_development/workflows/03b_parse_cytokine_trajectory_analysis/01_qualify_cytokine_trajectory_methods/07_finalize_focus_figures.py \
  --config current_analysis/Bdev_development/config/03b_parse_cytokine_trajectory_analysis/02_bdev_cytokine_trajectory_qualification.yaml

# Direction-corrected biological display: preserve native fitted L4 for analysis,
# but plot 0 = LSK and 1 = CD19, including the fitted curve on ATAC-UMAP.
/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --cleanenv --bind /mnt/sata2:/mnt/sata2 \
  --env 'R_LIBS=<run-library>:<shared-development-library>' \
  <development-sif> Rscript \
  current_analysis/Bdev_development/workflows/03b_parse_cytokine_trajectory_analysis/01_qualify_cytokine_trajectory_methods/08_plot_lsk_to_cd19_focus_figures.R \
  --config current_analysis/Bdev_development/config/03b_parse_cytokine_trajectory_analysis/02_bdev_cytokine_trajectory_qualification.yaml

python3 current_analysis/Bdev_development/workflows/03b_parse_cytokine_trajectory_analysis/01_qualify_cytokine_trajectory_methods/09_finalize_lsk_to_cd19_figures.py \
  --config current_analysis/Bdev_development/config/03b_parse_cytokine_trajectory_analysis/02_bdev_cytokine_trajectory_qualification.yaml
```

Generated state mirrors this workflow under paired development storage.

The original focused figures `04` and `05` retain native L4 orientation for
provenance. Figures `06` and `07` are the biological display intended for use:
their one-dimensional axes use `display_fraction = 1 - native_fraction`, so
`0 = LSK` and `1 = CD19`. Direction-dependent summary labels are transformed
consistently, while local scores and all underlying result tables are unchanged.
