# Step 4 full FLT3LKO/IL7KO mapping completion

Step 4 completed on 2026-11-18 under the corrective run identity
`flt3lko_il7ko_full_mapping_fixed_umap_20261118T031500Z`.

## Scope and execution contract

The workflow mapped all 23,252 QC-selected cells from the four FLT3LKO and
IL7KO libraries. It reused the qualified reciprocal liftOver feature map,
the already published exact HVP-depth20k query matrix, the frozen PeakVI model,
the fixed 26,723-cell Bdev `umap.atac`, and the authenticated fitted Slingshot
lineage L4 in its native orientation.

The corrective configuration is
`config/05_trajectory_group_contrasts/05_flt3lko_il7ko_full_mapping_fixed_umap_frozen_bundle_clone.yaml`
with SHA-256
`7f9632fbfb89db5255e7c833229d4b4b959923dc46ecfb74af3d713ddf378537`.
The mapping used only the configured primary neighborhood size, `k = 30`, and
processed the query in 47 bounded, checksummed batches (46 batches of 500
cells and one batch of 252 cells).

No PeakVI, UMAP, or Slingshot model was refit. The workflow performed no
genotype contrast or paper-level inference, built no SIF, and did not unlock
`Bdev_implementation`.

## Validated outcome

| Check | Result |
|---|---:|
| Reference cells | 26,723 |
| Query cells | 23,252 |
| FLT3LKO BM WT | 6,079 |
| FLT3LKO BM KO | 5,648 |
| IL7KO BM WT | 7,619 |
| IL7KO BM KO | 3,906 |
| Supported | 23,251 |
| Borderline | 1 |
| Out of reference | 0 |
| Finite accepted native L4 | 15,171 |
| Checkpoint files | 47 |

The support thresholds were calibrated on all 26,723 reference cells. The
primary-`k` thresholds exactly matched the retained multi-`k` calibration
baseline: supported upper distance `2.48994143907288` and borderline upper
distance `2.93105332005195`.

The frozen model-bundle manifest remained
`1f9e399600a9ab74b6037301af8a6b521bf4669b8cfd4b76f80f2d07ba437b98`.
The fixed Bdev execution-bundle manifest remained
`c54f0ef36effc37db30bf40b1f0f0d171ab701041b44dea41b26f3e58d1992f2`.
Both hashes were independently recomputed after completion and matched their
pinned configuration values.

## Durable outputs

The results root is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/05_trajectory_group_contrasts/04_map_flt3lko_il7ko_full/flt3lko_il7ko_full_mapping_fixed_umap_20261118T031500Z/
```

Key outputs are:

- `01_tables/05_full_support_aware_native_l4.tsv.gz`
- `01_tables/08_full_fixed_umap_placements.tsv.gz`
- `01_tables/09_full_mapping_group_summary.tsv`
- `02_figures/01_all_flt3lko_il7ko_fixed_atac_umap.png`
- `02_figures/01_all_flt3lko_il7ko_fixed_atac_umap.pdf`
- `03_reports/03_full_mapping_fixed_umap_summary.json`
- `04_logs/02_map_and_plot_fixed_umap.log`
- `05_objects/01_bounded_transfer_checkpoints/`
- `06_manifests/03_full_mapping_execution_receipt.json`

The corresponding derived mapping state is kept separately under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/intermediate_objects/05_trajectory_group_contrasts/04_map_flt3lko_il7ko_full/flt3lko_il7ko_full_mapping_fixed_umap_20261118T031500Z/
```

The source bundles were not modified. This completed run supersedes the failed
mapping attempt `flt3lko_il7ko_full_mapping_fixed_umap_20261118T022951Z` for
full-cell mapping, while the earlier run's qualified query quantification and
failure diagnostics remain traceable.
