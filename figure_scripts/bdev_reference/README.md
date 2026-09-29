# Bdev reference reports

`01_bdev_reference_report.Rmd` renders the accepted Bdev reference dataset
report from the promoted bundle. Its UMAP panels use `umap.atac`, the peaks
assay embedding, and the accepted pseudotime artifact recorded in the bundle.

`02_il7_flt3_top25_core_atac_umap.R` is a read-only sanity check for the
Bdev-compatible Parse signature review. It plots log1p(CP10k) Bdev RNA
expression for every IL-7 and Flt3L top25 core gene on the accepted
`umap.atac` coordinates, with separate combined figures for up and down sets.
Run it with `run_02_il7_flt3_top25_core_atac_umap.sh`; outputs and provenance
are recorded under the report path in `config/bdev_reference.yml`.

`03_flt3_il7_refined_sct_umap.R` plots the final top-100 refinement-pilot
memberships as four separate figures: Flt3L up/down and IL-7 up/down. It reads
expression directly from Bdev `SCT_data` and uses the existing `umap.atac`
coordinates. Non-positive cells are gray; positive cells are overlaid as
larger, fully opaque points on a saturated blue-to-magenta scale. Run it with
`run_03_flt3_il7_refined_sct_umap.sh`. The final v2 output and provenance are
registered in `config/bdev_reference.yml`.

`04_flt3_il7_refined_aggregate_sct_umap.R` creates one four-panel UMAP for
the same final sets. A cell's aggregate expression score is the arithmetic
mean of its `SCT_data` values across all genes in the set. The raw per-cell
scores and panel summaries accompany the figure. Run it with
`run_04_flt3_il7_refined_aggregate_sct_umap.sh`; its report path is registered
in `config/bdev_reference.yml`.

`05_flt3_il7_refined_ucell_trajectory.R` reproduces the earlier native-UCell
trajectory method for the four refined sets. It scores raw Bdev `RNA_counts`
with `max_rank = 1500`, preserves the existing `pt` values and missingness, and
uses 5% sliding windows with 1% steps, median lines, and IQR bands. Each curve
point is placed at the midpoint of its actual window, including clipped endpoint
windows. Dotted vertical lines mark the median source `pt` for the `up1` and
`up2` groups in Bdev metadata field `bb_up`. Run it with
`run_05_flt3_il7_refined_ucell_trajectory.sh`; the unpromoted report and score
session are registered in `config/bdev_reference.yml`.

`06_flt3_il7_refined_ucell_umap.R` displays the saved native UCell scores from
the refined trajectory analysis on Bdev `umap.atac`. It creates a four-panel
overview and four separate figures. Raw UCell values are retained in the
output table; the bright color scale is normalized within each panel to the
positive-score 99th percentile. Run it with
`run_06_flt3_il7_refined_ucell_umap.sh`; its read-only report is registered in
`config/bdev_reference.yml`.

`07_flt3_il7_ucell_trajectory_all_windows_diagnostic.R` is a sensitivity
display of the refined Flt3L/IL-7 trajectories. It retains the six windows
below the primary 30-cell support threshold and marks their unchanged native
medians with hollow points. It does not recalculate UCell scores or window
summaries. Run it with
`run_07_flt3_il7_ucell_trajectory_all_windows_diagnostic.sh`; the diagnostic
report is registered in `config/bdev_reference.yml`.

# All-cytokine refined UCell report

After reviewing a completed all-cytokine Parse refinement report, run:

```bash
figure_scripts/bdev_reference/run_08_all_cytokine_refined_ucell_report.sh \
  /path/to/all_cytokine_competitor_refinement_report
```

The unpromoted report scores only final Parse sets passing the mean and/or
responder route. It uses native scWorkbench UCell on Bdev `RNA_counts`, the
unchanged `umap.atac` coordinates, and the unchanged `pt` metadata. Pseudotime
curves use 5% windows, 1% steps, a 20-cell minimum, median lines and type-7 IQR
bands. Missing `pt` values remain missing.

`09_all_cytokine_refined_driver_report.R` extends that completed 38-profile
review without recalculating UCell scores or rolling windows. For every passing
cytokine-direction profile it generates one combined figure containing the
UCell score on `umap.atac`, the existing 5% rolling-window median/IQR plot, and
five colored driver-gene lines with right-edge labels. The exact reviewed
member sets are staged read-only on Bdev `SCT_data`; native scWorkbench
`trajectory_geneset_dynamics` ranks member genes from unsmoothed mean-expression
curves by dynamic score (amplitude tie-breaker). A second native operation
calculates unsmoothed median-expression lines for those independently selected
five genes in 50 equal-width bins. For plotting, each gene's mean trajectory is
independently min-max scaled to 0-1 so sparse, low-expression trajectories
remain visible; raw mean and median curves remain in separate tables. No result is promoted. Run
`run_09_all_cytokine_refined_driver_report.sh`; the configured report contains
38 PNGs and a 38-page PDF.

`10_all_cytokine_balanced_score_comparison.R` retains the original UCell
results and adds the user-authorized expression-balanced sensitivity analysis.
For each signature gene, zero `SCT_data` expression contributes zero and
positive values are converted to within-gene percentiles; genes are then
averaged with exact weight `1 / signature_size`. Each of 38 comparison figures
shows old and new UMAPs, old and new 5% rolling median/IQR trajectories, and the
existing top-five driver curves on the same per-gene internal 0-1 scale. Driver
selection is unchanged. The script records cell scores, equal weights,
gene detection diagnostics, score correlations, window membership, checks, and
provenance. This custom scoring rule is outside a single built-in scWorkbench
score operation by explicit user approval; bundle access and rolling
distribution summaries remain native scWorkbench operations.

The same script now also implements the explicitly requested third,
trajectory-level consensus method. For every signature and existing 5% window,
it calculates each member gene's mean `SCT_data` expression, min-max scales
that gene's eligible-window trajectory to 0-1, and reports the median and IQR
across all signature genes. This robust center cannot be set by a single gene.
It has no per-cell value and therefore no UMAP. The three-method figures show
both per-cell UMAPs, all three trajectories, and the unchanged five driver
genes; labeled dashed BB1/BB2 lines are included in the driver panel.

`11_all_cytokine_neighborhood_coordination_report.R` adds a third per-cell
method designed to reward coordinated multi-gene activity. Native scWorkbench
`knn_neighbors` finds 50 neighbors in the existing 50-dimensional RNA-PCA
space; the nested first 15, 30, and 50 neighbors support sensitivity analysis.
For each neighborhood size, every signature gene's neighborhood-mean
`SCT_data` value is converted to a zero-aware within-gene percentile. The
primary score uses `k=30` and equals the equal-gene mean percentile (intensity)
times the fraction of genes at or above the 75th percentile (participation).
Sensitivity results also cover thresholds 0.60 and 0.90. Each figure compares
native UCell, the per-cell balanced score, and the primary coordination score
as UMAPs and rolling trajectories, followed by the unchanged driver panel with
BB1/BB2. The session is unpromoted and bundle state is fingerprint-checked.

`12_all_cytokine_neighborhood_rank_change_report.R` reuses the saved native
RNA-PCA neighbors and measures change gene by gene. It smooths every signature
gene over `k=30`, converts the result to a within-gene percentile, calculates
the gene's median rank in the existing 5% windows, and subtracts its value five
1%-step windows earlier. Only after differencing are genes combined by their
mean. Because every gene change is bounded from -1 to 1, one member contributes
at most `1 / signature_size`; a separate flag requires at least 20% of genes to
change by 0.05 in the score direction. The UMAP uses a cell-specific version of
the same comparison; the trajectory shows the mean and gene-wise IQR; and the final panel shows the
five unchanged driver identities as individual rank-change curves with BB1/BB2.
Run `run_12_all_cytokine_neighborhood_rank_change_report.sh`. The custom score
is unpromoted and preserves exact source `pt` values and missingness.

`13_all_cytokine_rankmean_comparison_report.R` runs the built-in scWorkbench
`trajectory_geneset_dynamics` operation with `score_method="rank_mean"` on the
same 38 signatures and exact source `pt`. It compares native UCell, rank-mean,
and the saved k=30 gene-wise mean-rank-change method as UMAPs and on the same
5% rolling-window memberships. Since the native operation exports binned
curves rather than per-cell values, the script reproduces its documented
within-cell rank formula for UMAPs and checks every reconstructed 50-bin median
against the native result. The durable session and report remain unpromoted.

`14_all_cytokine_breadth_intensity_report.R` implements the explicitly
authorized breadth/rank-intensity score on raw Bdev `RNA_counts`. Detection is
count greater than zero; all positive genes are ranked within each cell from
highest to lowest with average ranks for ties. With `R=1500`, Breadth is the
detected signature fraction, Intensity is the cutoff-linear rank-weight sum
divided by the full signature size, and Combined uses `alpha=0.5`. The script
reuses the validated native rank-mean cell scores, unchanged `umap.atac`, and
exact saved 5%-window memberships. Every page is ordered rank-mean, combined,
breadth, intensity in both the UMAP and trajectory rows. It verifies raw-count
detection against `nFeature_RNA` for every cell, preserves exact `pt` values
and missingness, and writes only an unpromoted report. Run
`run_14_all_cytokine_breadth_intensity_report.sh`.

For display consistency, the rank-mean UMAP exactly follows the preceding
native rank-mean report's purple raw-value p1/p99 scale. The Combined UMAP
follows the established score-map convention used by the preceding comparison
reports: zero is gray and positive values are divided by their within-panel
p99 on the blue-cyan-green-yellow-orange-magenta palette. Breadth and Intensity
retain their green and orange raw-value p1/p99 scales.

After review, `processing/bdev_reference/04_promote_flattened_rank_scores.R`
promotes Combined as the accepted default Bdev gene-set score and retains
Breadth and Intensity as separately addressable component families. The script
uses native `cell_vector_family_calibrate` with `method="zscore"`,
`center=FALSE`, and `scale=FALSE` as an identity materialization, verifies all
promoted values exactly, records one atomic promotion receipt, and runs full
recipient-level validation.
