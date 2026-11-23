# Full gated 20D all-cell fixed-UMAP view

Run `gated_20d_mapping_all_cells_fixed_umap_20261123T000113Z` completed as a
separate follow-up to the accepted-only visualization. It preserves the prior
run and displays all 23,252 query cells on the unchanged Bdev ATAC UMAP.

The 20,278 cells failing the complete support-and-annotation acceptance gate
are drawn in dark grey (`#3F3F3F`). The 2,974 accepted cells are drawn on top
and colored by accepted native L4 pseudotime. The reference panel retains the
authenticated native-L4 colors and native Slingshot orientation.

The PNG and matching PDF are stored at:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/05_trajectory_group_contrasts/06_visualize_flt3lko_il7ko_gated_20d_mapping/gated_20d_mapping_all_cells_fixed_umap_20261123T000113Z/02_figures/01_all_cells_fixed_bdev_atac_umap.{png,pdf}
```

The six-file run inventory validates without mismatch. This view inherits the
checksummed 30-neighbor fixed-UMAP placements from the preceding visualization;
no model, UMAP, or Slingshot refit occurred, and no genotype contrast or paper
analysis was performed.
