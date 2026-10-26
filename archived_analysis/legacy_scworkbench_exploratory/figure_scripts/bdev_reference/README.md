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

`15_rna_atac_latent_progress.R` compares developmental movement in the imported
RNA `pca` and ATAC `lsi` spaces. It retains only cells with non-missing source
`pt`, orders those cells without changing `pt`, and calculates centroids in
full 200-cell sliding windows with a 20-cell stride. Consecutive-centroid
Euclidean distances are accumulated separately for each modality and divided
by the modality endpoint, producing relative distances from 0 to 1. The figure
plots ATAC relative distance against RNA relative distance, includes the
`y = x` reference, and highlights the windows nearest the existing `bb_up`
`up1`/`up2` median pseudotimes as BB1/BB2. The exact window membership, raw and
standardized distances, ordered trajectory cells, and diagnostics accompany
the PNG/PDF. `run_15_rna_atac_latent_progress.sh` executes the terminal custom
analysis and registers its script and principal outputs as an scWorkbench
analysis capsule; no source coordinates, metadata, or missing values are
changed.

`16_rna_atac_latent_progress_w30_s5.R` repeats that analysis with 30-cell
windows and a 5-cell stride. In addition to the relative-distance scatter, it
calculates the local centered finite-difference slope of RNA relative distance
with respect to ATAC relative distance. The slope figure uses a log10 y-axis,
marks slope 1 as matched local progress, shades the full source pseudotime
ranges of `bb_up` `up1`/`up2`, and marks their medians as BB1/BB2. Exact raw
slopes and BB ranges are retained in tables. Run and register the result with
`run_16_rna_atac_latent_progress_w30_s5.sh`.

`17_preliminary_smoothed_slope_gene_expression.R` is an explicitly preliminary,
read-only sensitivity display built on the 30-cell/5-cell windows. It overlays
a robust LOESS curve on log10 local slope and adds same-figure panels for Ebf1,
Pax5, and Ezh2 mean `SCT_data` expression with the same smoothing span. Raw
window values remain in exported tables and BB1/BB2 source ranges are shaded.
The runner intentionally skips capsule registration and full-bundle validation;
use it only for preliminary review.

`18_preliminary_slope_turning_regions_umap.R` finds sign changes in the centered
finite-difference derivative of that smoothed log10 slope. It reports all
candidate extrema, retains those with at least 1.10-fold prominence relative to
both adjacent extrema, and defines each highlighted region as the exact source
30-cell window at the selected extremum. Matching colors mark those windows on
the slope, Ebf1/Pax5/Ezh2 curves, and the existing Bdev `umap.atac`. The runner
is read-only and intentionally omits capsule registration and full validation.

`19_preliminary_turning_summary_ebf1_ezh2.R` presents the same selected T1–T4
windows with the smoothed slope, Ebf1, and Ezh2 vertically aligned on one
pseudotime axis. A companion RNA-versus-ATAC relative-distance panel marks the
turning windows on the trajectory, and `umap.atac` highlights their exact cells.
Pax5 is intentionally omitted. This remains a read-only preliminary display.

`20_preliminary_modality_ratio_umap_five_genes.R` extends that display with
Nipbl, Mau2, and Wapl. For each trajectory cell, the UMAP color is the mean
`log2(smoothed delta-RNA/delta-ATAC)` over every overlapping window containing
the cell. Red indicates RNA-faster, blue ATAC-faster, and white equal local
progress; non-trajectory cells are gray. T1–T4 cells are outlined and labeled.
The current expanded version also includes Ctcf, Eed, Tcf7, and Lef1. Tcf7 uses
a recorded Gaussian-LOESS fallback because robust LOESS fails numerically on
its low, sparse window-mean trajectory; all other genes retain robust LOESS.

`21_preliminary_slope_gene_set_signals.R` places the same smoothed RNA/ATAC
slope above four aligned aggregate trajectories: B-lineage, T-lineage,
Flt3L-UP, and IL-7-UP. Scores are arithmetic means of Bdev `SCT_data` across
the accepted memberships, summarized in the same 30-cell/5-cell windows and
smoothed with robust LOESS. Explicit Bdev mapping retains 236/246 B-lineage,
214/233 T-lineage, 45/45 Flt3L-UP, and 20/20 IL-7-UP genes; missing IDs are
exported rather than silently discarded.

`22_preliminary_lsi1_sensitivity.R` repeats the ATAC path calculation with
LSI1 excluded while holding RNA coordinates, trajectory cells, window
membership, pseudotime, stride, and smoothing fixed. It overlays the full-LSI
and LSI2–50 RNA-versus-ATAC paths and slopes, exports every paired window, and
uses the existing T1–T4 windows only as fixed visual references rather than
re-calling them.

`23_preliminary_tsix_expression.R` makes a separate two-panel Tsix figure from
the read-only Bdev bundle. One panel maps positive-cell `SCT_data` expression
to the existing `umap.atac`, with positive values clipped at their 99th
percentile for display. The other shows a robust-LOESS curve of 30-cell means
at 5-cell stride along the preserved Bdev `pt`. Both panels mark the fixed
T1–T4 windows from the RNA–ATAC distance analysis; the windows are not
re-called. The report also includes a standalone UMAP with no T1–T4 overlay;
it uses a conventional linear white-to-crimson expression scale from zero to
the positive-cell p99 so higher-expressing cells are easy to locate. Run it
with `run_23_preliminary_tsix_expression.sh`. A second standalone pseudotime
figure replaces T1–T4 with the exact source BB1/BB2 pseudotime ranges and
median lines, and fixes the expression y-axis minimum at zero.

`24_preliminary_log2_modality_ratio.R` displays the saved smoothed local
progression slope explicitly as `log2(delta RNA / delta ATAC)`, centered at
zero for matched progression. A secondary axis shows the exact inverse
`log2(delta ATAC / delta RNA)`, which is the negative of the primary curve and
not a separately fitted signal. T1–T4 remain fixed visual references. The
report also records a proposal for direct and inverse gene-pattern screens.

`25_preliminary_modality_gene_pattern_screen.R` implements the two proposed
screens across sufficiently detected, variable Bdev `SCT_data` genes. The
state approach correlates 30-cell window-mean expression with the saved log2
RNA/ATAC target; the change approach correlates the centered expression
derivative with that target. Both use one fixed centered smoother spanning 8%
of ordered windows before scoring. Negative scores are the exact ATAC-faster inverse
interpretation. Circular shifts preserve trajectory autocorrelation, pooled
null tails provide empirical P values, and BH correction is applied separately
per approach. This remains descriptive because feature-fold cross-fitting of
the RNA latent trajectory is not part of this preliminary screen.

`26_preliminary_state_top20_separate_figures.R` makes separate faceted figures
for the top 20 RNA-faster state genes and the top 20 ATAC-faster inverse state
genes. Each gene panel shows its fixed 8%-window-smoothed mean plus the smoothed
within-window cell-level IQR and the matching RNA/ATAC or inverse ATAC/RNA
target. A companion figure shows Angpt1, Adgrl4, Gcnt2, and Mecom in two tiers:
expression with IQR above and the screened expression derivative against the
RNA-faster target below.

`27_preliminary_modality_turn_pattern_screen.R` replaces whole-trajectory
state similarity with an event-focused up/down screen. It requires a gene to
have p95-p05 smoothed `SCT_data` range of at least 0.50, weights agreement over
the strongest 30% of RNA/ATAC target changes, and separately checks the two
flanks of each fixed T1-T4 turn. Selection requires six of eight matching
flanks and at least two complete turns, so long flat portions cannot dominate
the ranking. Direct RNA/ATAC and inverse ATAC/RNA results, every per-turn
decision, amplitude-cutoff sensitivity, and separate top-20 figures with
within-window IQR are exported by
`run_27_preliminary_modality_turn_pattern_screen.sh`.

`28_preliminary_epigenetic_pattern_candidates.R` plots the chromatin and
cohesin-associated candidates highlighted from the event-focused screen.
Separate direct and inverse figures overlay each standardized gene trajectory
and its smoothed within-window IQR on the appropriate modality target. A third
figure audits Nipbl/Mau2, core cohesin-associated factors, and PRC2 genes.
Every facet explicitly states `PASS` or `NOT PASS` and gives the p95-p05 range,
matched T1-T4 flanks, and complete-turn count. Run it with
`run_28_preliminary_epigenetic_pattern_candidates.sh`.

`29_preliminary_inverse_cohesin_turn_audit.R` replots eleven
cohesin-associated genes against the inverse ATAC/RNA target. The trajectory
figure marks every T1-T4 pre/post flank with a green pass or red fail dot; its
companion heatmap shows all eight flank calls and four complete-turn calls
explicitly. Wapl is highlighted as the principal inverse near-hit, passing
both flanks at T2 and T3 plus the T4 pre-flank. Run it with
`run_29_preliminary_inverse_cohesin_turn_audit.sh`.

`30_preliminary_t1_t3_continuous_gene_screen.R` performs an unbiased
genome-wide continuous shape ranking after excluding T4. Its raw score is 60%
target-change-weighted derivative cosine across the complete T1-T3 interval
and 40% cosine of the six T1-T3 flank changes. No flank or complete-turn call
is used as a cutoff. Low-amplitude or heterogeneous profiles are instead
smoothly downweighted using the p95-p05 range and conservative cell-level
spread derived from exact within-window IQR and SD. It exports all 11,745 gene
scores, separate direct/inverse top-100 tables, top-20 IQR figures, prior-gene
scores, and a Nipbl rank diagnostic. Run it with
`run_30_preliminary_t1_t3_continuous_gene_screen.sh`.

`31_preliminary_turn_local_balanced_gene_screen.R` removes full-trajectory
similarity from the ranking and scores the six T1-T3 pre/post flanks locally.
Each flank combines target-change-weighted derivative shape and oriented
endpoint change, with smooth amplitude and exact within-window IQR/SD
precision weights. Pre/post flanks contribute equally within a turn, and T1,
T2, and T3 contribute equally with penalties for imbalance, preventing one
strong peak from dominating the result. The top-20 direct and inverse figures
display the complete pseudotime range while shading the scored +/-100-window
regions; T4 is shown only as unscored context. Run it with
`run_31_preliminary_turn_local_balanced_gene_screen.sh`.

`32_preliminary_inverse_hits_depth_diagnostic.R` compares the inverse top 20
with total RNA (`nCount_RNA`) and total ATAC (`nCount_peaks`) along the same
trajectory. It exports a two-panel total-signal plot with cell-level IQR and a
20-panel overlay of each gene, total RNA, and total ATAC. Both display the full
pseudotime range and shade the fixed +/-100-window T1-T3 regions. The companion
table reports full-path, scored-region, cell-level, and pseudotime-residual
correlations to distinguish shared developmental trends from residual depth
association. Run it with
`run_32_preliminary_inverse_hits_depth_diagnostic.sh`.

`33_preliminary_ratio_total_signal_explained.R` quantifies how much of the
log2 RNA/ATAC ratio is associated with smoothed log-total RNA and ATAC. It
reports in-sample level R2, contiguous-block prediction, derivative R2 in the
fixed T1-T3 neighborhoods, per-turn derivative fits, and leave-one-turn-out
shape prediction. Figures show the full pseudotime range, the fitted and
residual ratio, and the local residual derivative with the same shaded
+/-100-window regions. Run it with
`run_33_preliminary_ratio_total_signal_explained.sh`.

`34_preliminary_residual_ratio_bb_gene_stack.R` makes two full-pseudotime,
vertically aligned figures headed respectively by the total-signal-adjusted
RNA/ATAC ratio and adjusted ratio derivative. Ebf1, Ezh2, Nipbl, and Wapl are
stacked underneath with smoothed within-window cell-level IQR. Every panel
shades the exact source `bb_up` BB1/BB2 pseudotime ranges and marks their
medians. Run it with
`run_34_preliminary_residual_ratio_bb_gene_stack.sh`.

`35_preliminary_residual_ratio_umap.R` projects the full-path residual
log2 RNA/ATAC ratio onto the existing `umap.atac`. Window residuals are
linearly interpolated from window-center rank to each cell's exact trajectory
rank, avoiding the five-cell stripes produced by membership-mean projection;
the original membership mean remains in the exported comparison table. The
blue-white-red scale is centered at zero; both an absolute-p99-clipped figure
and an unclipped figure using the observed maximum absolute value are
exported. Cells outside the trajectory are gray. Run it with
`run_35_preliminary_residual_ratio_umap.sh`.

`36_preliminary_residual_ratio_weighted_gene_screen.R` screens all detected
`SCT_data` genes for direct and inverse up/down agreement with the corrected ratio
without using BB or T-region boundaries. It conditionally removes each gene's
RNA-depth component after controlling for preserved Bdev pseudotime. Every
pseudotime window is scored, with a nonzero weight floor and extra continuous
weight assigned from the corrected ratio's magnitude, slope, and curvature.
The ratio and gene curves receive the same 3% smoother. Slopes are normalized
by each curve's own scale and tanh-clipped before comparison, so magnitude is
allowed to differ. Eight equal-interest trajectory segments and a smooth
between-segment consistency penalty reduce one-feature hits. It exports
complete scores, top-100 tables, and separate top-20 plots with cell-level IQR. Run it with
`run_36_preliminary_residual_ratio_weighted_gene_screen.sh`.

`37_preliminary_residual_ratio_full_direction_matches.R` applies the literal
full-pattern audit to that screen and retains genes whose normalized slope has
the requested direction in every one of the eight balanced trajectory
segments. It exports the complete direct and inverse concordant sets and plots
them over the standardized target with adjusted-expression IQR. Run it with
`run_37_preliminary_residual_ratio_full_direction_matches.sh`.

`38_preliminary_adjusted_cohesin_prc2_profiles.R` plots the corrected ratio
above full-pseudotime, RNA-depth-adjusted expression for the prior cohesin and
PRC2 audit sets. Ebf1 appears in both figures. Ratio and expression use the
same 3% trajectory smoother, expression is standardized so magnitude may
differ, and ribbons show adjusted-expression cell-level IQR. Panel labels give
direct and inverse agreement counts across the eight balanced trajectory
segments. Run it with
`run_38_preliminary_adjusted_cohesin_prc2_profiles.sh`.

`39_preliminary_adjusted_selected_genes_bb.R` makes a compact two-column
landscape figure containing the corrected ratio and only Ebf1, Tcf3, Stat5a,
Stat5b, Ezh2, Nipbl, and Wapl. Expression uses the same conditional RNA-depth
adjustment and matched 3% smoothing as the direction screen, ribbons show
cell-level adjusted-expression IQR, and every panel shades the exact source
BB1/BB2 ranges and marks their medians. Run it with
`run_39_preliminary_adjusted_selected_genes_bb.sh`.

`40_preliminary_nipbl_wapl_ezh2_umap.R` maps Nipbl, Wapl, and Ezh2 `SCT_data`
onto the existing `umap.atac`. Each panel uses a white-to-crimson scale
normalized to that gene's positive-cell p99, clips higher values for display,
and draws higher-expression cells last. Absolute p99 values and detection
fractions remain in the subtitles and exported scale table. Run it with
`run_40_preliminary_nipbl_wapl_ezh2_umap.sh`.

`41_preliminary_cox6a1_cox6a2_ube2c_umap.R` maps Cox6a1, Cox6a2, and Ube2c
`SCT_data` onto the existing Bdev `umap.atac`. It uses the same gene-specific
white-to-crimson zero-to-positive-cell-p99 display and draws high-expression
cells last. The uncapped cell-level values, color thresholds, PNG, and PDF are
written to a separate read-only report. Run it with
`run_41_preliminary_cox6a1_cox6a2_ube2c_umap.sh`.
