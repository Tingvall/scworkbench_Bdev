# 03b Parse cytokine trajectory methods with scwbTrajectory 0.1.5

**Status:** bounded Bdev development qualification passed; immutable SIF build and
acceptance still pending

## Qualified runtime source

The optimized scwbTrajectory development source remains package version `0.1.5`
at published commit `34c792d39a7a9950294bf3e29e1cb65ebf1aba6c`. The source
operation-pack manifest file has SHA-256
`576156623f04d8d34f07062eb17a1612ce1e4d20926516afbc06a74af141c8f9`.
The loaded registry identity was
`5afd786c35a40e88972f800d34c8b96fc7b56f586104980554c13718c18c0a0d`.

The three qualified operations are:

1. `scwbtrajectory::trajectory_signal_dynamics_explicit`;
2. `scwbtrajectory::trajectory_signal_intervals`; and
3. `scwbtrajectory::trajectory_preranked_gsea`.

The run used `fgsea 1.34.2`, `BiocParallel 1.42.2`, and the optimized bounded
GSEA permutation implementation. All 12 focused interval/GSEA contracts passed.

## Bdev qualification

Run `bdev_cytokine_trajectory_20261105T153001Z` used the authenticated fitted
Slingshot lineage `L4` from the approved ATAC-UMAP recipe in native orientation.
The preserved scalar `pt` was unchanged and was not used as the trajectory.

The deterministic fixture contained 1,810 cells across 20 bins and all four
`orig.ident` samples, selected from 13,451 available L4 cells. All 38 frozen
Parse cytokine sets passed overlap checks for both `RNA_counts` and `SCT_data`.

All methods used the same nominal trajectory scale:

| Method | Scale |
|---|---|
| Dynamics | 20 non-overlapping bins, median summary |
| Signal intervals | width `0.05`, step `0.025` |
| Local GSEA | width `0.05`, step `0.025` |

Signal intervals and GSEA used unblocked permutations. `orig.ident` was retained
only for the bin-by-sample coverage audit and was not passed as a permutation
block. The bounded development run used 99 permutations per inferential method;
the frozen paper analysis may use the higher prespecified count after the new SIF
is built and accepted.

The run produced 53 supported signal intervals, 38 global GSEA results, 1,406
local GSEA results, and 50 GSEA enrichment intervals. fgsea reported tied-rank
and small-p-value precision warnings for some rankings; these are retained in
the log and were non-fatal.

## Method figures

The development run includes one comparison-aligned figure for each method:

1. a row-standardized 20-bin UCell dynamics heatmap;
2. an elevated/reduced signal-interval map with global-support markers; and
3. a local NES heatmap with an aligned global NES panel.

All figures retain the same dynamics-score signature order and native L4
orientation. PNG and vector PDF versions are stored under the run's
`02_figures/` directory. Their source-table and file checksums are recorded in
`06_manifests/04_figure_manifest.json`.

Focused FLT3L-up and IL-7-up summaries add a fourth aligned view: the qualified
L4 fixture's UCell scores over the approved `umap.atac` embedding, with all
other Bdev cells shown as a gray background. Each focused figure also contains
the 20-bin dynamics curve, unblocked signal intervals, and local/global GSEA.
Their provenance is recorded in `06_manifests/05_focus_figure_manifest.json`.

The retained preliminary-analysis figures use the biologically interpretable
display direction `LSK (0) -> CD19 (1)`. This is a visualization transform,
`display_fraction = 1 - native_L4_fraction`; it does not alter the authenticated
native L4 fit or any result table. Direction-dependent labels are transformed
consistently, and the fitted L4 curve is overlaid with an LSK-to-CD19 arrow on
the ATAC UMAP. The retained files are:

1. `02_figures/06_flt3l_up_lsk_to_cd19_four_method_summary.{png,pdf}`;
2. `02_figures/07_il7_up_lsk_to_cd19_four_method_summary.{png,pdf}`; and
3. `06_manifests/06_lsk_to_cd19_figure_manifest.json`.

The direction-figure manifest has SHA-256
`012fb7e08b652a4c77fef60355cccb1fd0c45b4ed9df4f71d7696797fdc6b329`.

## Preliminary-analysis handoff

This completed bounded run and its LSK-to-CD19 focused figures are retained as
the preliminary analysis to guide the later Bdev implementation. They define
the selected cytokine programs, the four-view presentation, the common
bin/window scale, and the display direction to reproduce.

They are not paper-analysis outputs and must not be copied or promoted as final
implementation results. After the new immutable SIF is built, qualified,
checksummed, and explicitly accepted, the implementation stage should rerun the
frozen analysis from its declared inputs and compare the new results with this
preliminary run as a reproducibility reference.

## Provenance boundary

The archived source bundle was not modified. Registry imports, fitted lineage,
signals, and method outputs were written only to the isolated development clone.
No outputs were promoted and `current_analysis/Bdev_implementation/` remains
locked pending a newly built, checksummed, qualified, and explicitly accepted
immutable SIF.

The final checksummed run manifest is stored at:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/03b_parse_cytokine_trajectory_analysis/01_qualify_cytokine_trajectory_methods/bdev_cytokine_trajectory_20261105T153001Z/06_manifests/03_final_run_manifest.json
```

Its SHA-256 is
`bd6ac82554a7cbacd63c6ef0a750f69b97b19295cecfa013562ea097c20f7085`.
