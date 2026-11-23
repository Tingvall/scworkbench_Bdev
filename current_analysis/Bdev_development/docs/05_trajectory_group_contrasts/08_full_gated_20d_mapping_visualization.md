# Full gated 20D mapping visualization

Run `gated_20d_mapping_visualization_20261122T235003Z` completed on
2026-11-22. It visualizes the finalized Step 05 gated mapping without
changing the mapping result or performing a genotype contrast.

Only rows with `accepted == TRUE` and finite `accepted_native_l4` are colored.
This is the complete support-and-annotation acceptance gate, not a separate
numeric cutoff applied to the pseudotime value. The figures therefore color
2,974 of 23,252 query cells (12.79%). Rejected cells are retained as a neutral
background on the native query UMAPs and are not assigned a displayed
pseudotime.

## Views

1. The authenticated Bdev L4 reference and accepted query placements on the
   unchanged Bdev `umap.atac`. Query coordinates are 30-neighbor
   distance-kernel barycenters calculated from the finalized 20D checkpoint
   ledgers.
2. Accepted cells on the read-only source Cell Ranger ARC ATAC UMAP.
3. Accepted cells on the read-only source Cell Ranger ARC GEX UMAP.
4. Descriptive accepted native-L4 histograms for each library.

The source ATAC and GEX UMAP CSV paths and checksums are declared in
`config/05_trajectory_group_contrasts/07_flt3lko_il7ko_gated_20d_mapping_visualization.yaml`.
They remain in archived storage and were read without modification.

## Outputs

The results root is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/05_trajectory_group_contrasts/06_visualize_flt3lko_il7ko_gated_20d_mapping/gated_20d_mapping_visualization_20261122T235003Z/
```

The four numbered figure pairs are:

- `02_figures/01_accepted_mapping_fixed_bdev_atac_umap.{png,pdf}`
- `02_figures/02_accepted_native_l4_query_atac_umap.{png,pdf}`
- `02_figures/03_accepted_native_l4_query_gex_umap.{png,pdf}`
- `02_figures/04_accepted_native_l4_distributions.{png,pdf}`

The combined cell-level placement table is
`01_tables/01_gated_20d_umap_placements.tsv.gz`; the library summary is
`01_tables/02_accepted_native_l4_summary.tsv`. The final inventory contains
13 checksummed files and validates without mismatch. No model, UMAP, or
Slingshot refit occurred, no paper analysis was performed, and
`Bdev_implementation` remains locked.
