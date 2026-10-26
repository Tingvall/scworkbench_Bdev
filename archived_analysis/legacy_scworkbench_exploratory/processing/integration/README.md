# Bdev-centered integration workflow

The latest review branch repeats the pooled Bdev trajectory analysis with all
33 full signatures accepted by the competitor-aware
`parse_cytokines_rest_only` workflow. Run ID
`competitor_aware_pt_all33_20260925` uses native raw-count UCell scores on all
26,723 Bdev cells and the exact 12,642 defined `Bdev@meta.data$pt` values. The
same 5%, 10%, and 15% sliding windows, 1% steps, 30-cell minimum, endpoint
clipping, medians, and IQR bands are used so the output is directly comparable
with the prior trajectory report.

All 33 signatures are retained. The former transfer rule of at least five
mapped genes and 80% membership is reported as a diagnostic and does not filter
this branch. Twenty-four signatures meet that rule; nine are scored with lower
Bdev mapping after missing genes are dropped. The report lists the exact mapped
membership for each set. Scores, imports, and native operation history remain
in the durable review session and have not been promoted. The reproducible
entry point is
`bash processing/integration/run_competitor_aware_trajectory.sh score|figures competitor_aware_pt_all33_20260925`.
The report is under
`reports/competitor_aware_pt_all33_20260925/trajectory_review/` in the
configured integration storage root.

The approved [cytokine pseudotime strategy](CYTOKINE_PSEUDOTIME_STRATEGY.md)
starts with the accepted Bdev and Parse parents. KO integration follows after
that parent is accepted. The registry distinguishes active from future parents.

The earlier `cytokine_pt_20260923T2231Z` run is also awaiting user review of its
figure report. Analysis artifacts
had already been promoted into the working derived bundle before the request
to review figures first. No additional promotion, finalization or acceptance is
authorized while this review is pending. The runner enforces
`review_status: awaiting_user_review` by allowing only the read-only `figures`
stage. Change that status only after explicit user approval to resume.

That earlier scope used pooled data with sliding windows. Its figures stage
selects `03_sliding_window_cytokine_report.R`, which computes native window
medians, quartiles and SD from retained UCell scores and writes a separate
`sliding_window_review/` report. The main width is 10% of the observed pt range,
with 1% steps; 5% and 15% widths are sensitivity checks. Windows require 30 cells
and are clipped at observed endpoints. Continuous median lines and IQR bands
connect retained centers; low-coverage centers remain marked. All source
samples and cell-cycle states are combined. The fixed-bin `pooled_review/`
report remains available for comparison.
The historical stratified workflow in `04_cytokine_pseudotime.R` refuses to run
under the current pooled configuration. Its already-retained outputs are
preserved but are not analyzed or displayed in this revision.

Integration creates the separate `bdev_integrated` bundle. It must never add
derived integration state to any parent bundle.

The initial cytokine workflow consisted of:

1. `01_preflight_inputs.R`: check accepted parent fingerprints and validation
   evidence, gene coverage, source samples, and exact pseudotime values.
2. `02_create_bdev_integrated.R`: import Bdev RNA counts, metadata and all five
   coordinate spaces through scWorkbench; verify every RNA value and coordinate;
   define pseudotime, sample, and cell-cycle subsets through registered sessions.
3. `03_score_parse_signatures.R`: import all 23 source sets, check the prespecified
   coverage gate, and retain UCell scores for the 22 eligible sets.
4. `04_cytokine_pseudotime.R`: run pooled, sample and phase dynamics, score QC,
   RNA/cell-cycle covariate curves, and secondary member-gene expression.
5. `figure_scripts/bdev_integrated/01_cytokine_pseudotime_report.R`: render the
   historical report. The current runner uses the pooled sliding-window `03`
   report and `src/sliding_window_signature_summary.R`, with all bundles read-only.
6. After review approval, `05_validate_cytokine_bundle.R`: verify scientific
   preservation, finalize recipient readiness, and perform complete recipient
   validation before registry acceptance.

Run stages with `bash processing/integration/run_stage.sh STAGE RUN_ID`, where
`STAGE` is `preflight`, `create`, `score`, `dynamics`, `figures`, or `validate`.
The initial run ID is `cytokine_pt_20260923T2231Z`. Do not rerun analysis stages
while report review is pending. Inputs are mounted read-only; the figures stage
makes only its report directory writable. Other authorized stages make the
integration storage root writable. The runner
uses the registered runtime, four CPU threads, and a 64-GiB address-space ceiling.
Completed operation receipts are reused and incomplete durable sessions resumed.
Trajectory comparisons share one durable session, checkpoint each operation,
and promoted their selected outputs together in this initial run. A future
analysis requiring new promotion must first expose its staged results for
figure review; the existing analysis scripts are not a pre-promotion preview.
Bundle creation refuses to replace retained data; rerunning creation after an
interruption verifies the existing import before continuing.

The primary pseudotime family is constructed from exact metadata doubles.
The older scalar artifact is recorded as provenance but is not used as the
numeric source because its character conversion rounds values. No values are
imputed, flipped, rescaled, or refitted. The 14,081 missing entries remain in
metadata. Gene-set membership is transferred by exact feature ID; missing genes
are reported and no synonym mapping is inferred.

In the initial branch, coverage eligibility is at least five genes and 80% of
the original membership.
IFN-gamma up has 35/50 genes and is excluded from primary scoring. The source set
and its mapped membership remain retained. Secondary member-gene curves use
mean log1p(CP10K) expression to retain sparse RNA support; primary UCell curves
use bin medians. All results are descriptive.

For the Parse parent, signature provenance refers to the accepted
`set:parse_qc_control_single` universe and its promoted
`parse_treatment_category` vector. Combined-stimulation cells remain retained
in the parent bundle but are outside the integration input scope.

Reports, interfaces, operation receipts, figures, and validation evidence are
under the configured integration root's `reports/RUN_ID/`; durable sessions are
under `sessions/RUN_ID/`. Source fingerprints and runtime build evidence are
recorded in preflight and completion reports. Actual operation interfaces were
inspected in the supplied image with its registered runtime before execution.

## KO mapping stage

KO mapping now uses its own derived bundle,
`bdev_integrated_ko_mapping_v1`, and does not alter the cytokine bundle or
either parent. `07_create_ko_mapping_bundle.R` combines the accepted Bdev cells
and accepted KO query cells over the 97,788 canonical Bdev mm10 peak IDs. KO
counts come from the exact-reciprocal mm39 lookup intervals prepared by the KO
workflow. Cell IDs are namespaced, and the exact Bdev metadata `pt` values and
missing entries are retained as `reference_pt`.

`08_peakvi_mapping_smoke.R` performs a runtime/contract test only: it trains a
one-epoch, five-dimensional CPU PeakVI model on 256 deterministic Bdev cells,
projects 256 stratified KO cells without refitting, row-binds coordinates with
the model-bound union operation, and fits a display UMAP. All coordinate,
ordering, projection-mode, and finite-value checks pass. Results remain in an
unpromoted durable session and are not scientific mapping evidence.

Run `bash processing/integration/run_ko_mapping.sh create|smoke`. Creation is
refuse-to-overwrite; the smoke stage resumes and verifies its durable session.
The mapping bundle passed recipient validation with only expected warnings for
external provenance inputs outside the recipient trust boundary. Full PeakVI
training and projection have not started and require a separate retained run.

Before the retained fit, `09_peakvi_full_reference_benchmark.R` runs the real
full-reference path for five GPU epochs in a separate unpromoted session.
`run_peakvi_gpu_benchmark.sh` records one-minute utilization, memory,
temperature, power, and process observations. Without requiring administrator
configuration, its courtesy monitor pauses this analysis process group when a
foreign GPU compute process appears and resumes it after the GPU is free.

`10_peakvi_full_mapping.R` runs the retained full-reference fit, frozen-model
projection of every accepted KO cell, model-bound union, and display UMAP in a
separate durable session. `run_peakvi_full_mapping.sh` samples GPU state every
10 seconds and pauses the analysis if a compute process owned by another user
appears. Full outputs remain unpromoted for mapping QC and review.
If another user's long-lived CUDA process would otherwise pause the mapping
indefinitely, `run_peakvi_low_duty_share.sh` adds a no-admin cooperative mode:
while foreign work is present, this analysis runs for 15 seconds and yields for
45 seconds; it runs continuously when no foreign process is present.

The full `peakvi_full_mapping_v1` session completed without promotion. It
contains the 26,723-cell Bdev reference latent space, frozen-model projection
of all 23,252 KO cells, their 49,975-cell model-bound union, and a joint UMAP.
All registered dimension, cell-order, finite-value, projection-without-refit,
and union-without-transformation checks passed.

`11_peakvi_mapping_qc.R` performs the approved read-only exception outside a
registered scWorkbench operation. The registered neighbor implementation is a
brute-force R calculation and cannot do the required Bdev-restricted searches
at this scale. The script therefore uses exact RANN Euclidean searches in the
20-dimensional staged PeakVI space, with Bdev as the only neighbor pool. It
uses k=30 as primary and k=15/50 as sensitivity checks, and calibrates label
transfer and median-30-distance out-of-reference gates with five deterministic,
peak-cluster-stratified Bdev folds. Pseudotime uses only the 12,642 Bdev cells
with defined `reference_pt`, with a weighted median and weighted interquartile
range; values are withheld from the primary output for out-of-reference cells.

Run `bash processing/integration/run_peakvi_mapping_qc.sh`. The completed,
unpromoted report is under `reports/ko_mapping_qc_v1_20260930/`. At primary
k=30, held-out Bdev cluster accuracy is 0.800 (balanced 0.787), phase accuracy
is 0.756 (balanced 0.536), and pseudotime Spearman correlation is 0.990.
Among 23,252 KO cells, 20,595 are supported, 1,876 borderline, and 781
out-of-reference. Exact neighbor tables, per-cell confidence and sensitivity,
calibration predictions, sample summaries, and figures are retained as report
artifacts. Nothing was promoted, and the exact Bdev `reference_pt` vector and
its missing entries remain unchanged.

`figure_scripts/bdev_integrated/06_peakvi_reference_umap_projection.R` provides
the reference-anchored display requested after review of the joint UMAP. It
keeps all native `bdev_reference` `umap.atac` coordinates unchanged and assigns
each query cell the distance-weighted barycenter of its exact Bdev-restricted
PeakVI neighbors. The primary display uses k=30, with k=15/50 sensitivity and
weighted neighbor-coordinate spread retained per cell. The 781 calibrated
out-of-reference cells retain raw audit coordinates but are withheld from the
primary display. Run it with
`bash processing/integration/run_peakvi_reference_umap_projection.sh`. The
completed report-only output is under
`reports/ko_mapping_reference_umap_v1_20261001/`; no space was registered or
promoted and the Bdev bundle fingerprint remained unchanged.

### Current default Bdev mapping and enrichment regions

At present, the default query-to-Bdev mapping uses the Bdev-only PeakVI model
frozen before query projection; query cells are transformed without refitting.
The existing general-support gate is applied first. Trajectory membership then
requires a score of at least 0.25, where the score is the distance-weighted
fraction of the 30 Bdev neighbors that have defined `pt`. This score was
calibrated by leaving out each original Bdev sample in turn. Query pseudotime is
the distance-weighted median of defined neighbor `pt` values. The exact
`Bdev@meta.data$pt` values and their missing entries remain unchanged.

Mapped enrichment is evaluated in 100 Bdev pseudotime bins. Within each target
and genotype, the retained-query conditional cell fraction is smoothed with a
centered three-bin rolling mean and divided by the corresponding Bdev rolling
fraction. A genotype-specific background envelope is the bin-wise 95th
percentile from 1,000 multinomial draws from the Bdev pseudotime distribution,
matching that genotype's retained cell count (seed 1729). Consecutive WT bins
above the WT-specific envelope define the primary enrichment peaks. KO
intervals and individual pseudotime bins above the separately size-matched KO
envelope are retained for possible later use, but do not redefine the primary
WT peaks. Significant Bdev residual-UCell cytokine regions are display-only
overlays and are not used in peak calling.

Run `bash processing/integration/run_84_bdev_cutoff025_mapping_report.sh`. The
current report and machine-readable WT/KO region artifacts are under
`reports/bdev_cutoff025_mapping_report_v4_20261012/`. This is a report-only
default and does not mutate or promote any source bundle.

`85_wt_enrichment_rolling_window_cell_support.R` applies the reviewed support
rule for the motif-analysis display domain. Within each WT-defined enrichment
peak it finds maximal consecutive runs in which every source bin has at least
10 WT cells. Runs shorter than three bins are audited and excluded because a
centered three-bin window cannot be formed without extending outside the run.
Displayed counts sum the left, middle, and right source bins and are located at
the middle-bin pseudotime. KO support is assessed from the unsmoothed middle
source bin. Regions without any middle bin containing at least 10 KO cells are
excluded; in each retained region, one orange vertical line separates the last
supported middle bin from the subsequent unsupported tail. Run with
`bash processing/integration/run_85_wt_enrichment_rolling_window_cell_support.sh`.
The completed read-only report is under
`reports/wt_enrichment_rolling_3bin_cell_support_v2_20261013/`.

`86_rolling_window_motif_trajectory_heatmaps.R` applies those exact windows to
the motif analysis. Within each retained region, Bdev alone defines accessible
peaks, chromVAR background and expectation, temporal-variability ranking,
top-30 motif selection, centering, scaling, and row order. WT and KO enter only
afterward. Each figure aligns Bdev, WT, KO, and descriptive WT-minus-KO panels;
the agreed KO middle-bin support boundary is drawn in every panel. The display
uses centered three-bin pseudobulks, and trajectories for all 746 eligible
motifs are exported. Run with
`bash processing/integration/run_86_rolling_window_motif_trajectory_heatmaps.sh`.
The completed report is under
`reports/rolling_3bin_motif_trajectory_heatmaps_v1g_20261013/`.

`run_87_single_bin_motif_trajectory_heatmaps.sh` runs the same Bdev-defined
motif workflow without rolling aggregation. Each displayed pseudobulk contains
exactly one of the retained source bins: 18 positions for Flt3L WT peak 1 and
13 for IL-7 WT peak 1. The top-30 motifs are reselected from Bdev single-bin
temporal variability, while calibration, scaling, and row order remain
Bdev-only. KO bins with zero cells, and their WT-minus-KO comparisons, are
explicitly gray rather than imputed. Filled tracks above the heatmaps show each
group's within-group accepted-cell fraction divided by the Bdev defined-pt
cell fraction in the same bin. A value of one is Bdev-like occupancy; the
comparison panel overlays WT and KO. The completed report is under
`reports/single_bin_motif_trajectory_heatmaps_v2_20261013/`.

`run_88_fine_rolling_motif_trajectory_heatmaps.sh` provides the narrowed
rolling-window version. It halves the original source-bin width and aggregates
three fine bins per centered pseudobulk, so each window spans 1.5 original
bins and advances by 0.5 original bin. This yields 31 Flt3L and 21 IL-7 focal
positions. WT support remains at least 22 and 43 cells per window,
respectively. The Bdev-only motif workflow, density tracks, KO boundary, and
explicit gray display for zero-cell KO windows are otherwise unchanged. KO
and WT-minus-KO motif tiles after the support line are masked, while their
computed values remain exported with a display-support flag. The
completed report is under
`reports/fine_rolling_motif_trajectory_heatmaps_v2_20261013/`. In v2, the KO
line is aligned to the fine grid: it falls halfway between the last half-width
middle bin with at least 10 KO cells and the first subsequent unsupported tail
bin.

`89_single_cell_axis_15nn_motif_heatmaps.R` removes the fixed plotting grid:
every Bdev cell in each retained region supplies a focal position at its
preserved `pt`. At every position, Bdev, WT, and KO activity is the mean of
that group's 15 nearest pseudotime-cell chromVAR deviations. The top-30 motif
identities are fixed from the Bdev-only fine rolling analysis for a direct
resolution comparison; per-cell chromVAR calibration and trajectory scaling
remain Bdev-only. A group is masked if its 15th neighbor is more than 0.25 pt
away, and KO plus WT-minus-KO are masked after the fine-bin KO boundary.
Filled tracks show 15-NN local density relative to Bdev. Run with
`bash processing/integration/run_89_single_cell_axis_15nn_motif_heatmaps.sh`.
The completed report is under
`reports/single_cell_axis_15nn_motif_heatmaps_v2_20261013/`.

`run_90_single_cell_axis_30nn_motif_heatmaps.sh` repeats that exact workflow
with 30 nearest pseudotime cells per group. It is smoother but has a wider
effective bandwidth. Neighbor radius is exported for auditing but does not
mask values; the 0.25 threshold applies only to trajectory-support selection
of mapped query cells. KO and WT-minus-KO remain masked after the agreed KO
boundary. The completed comparison report is under
`reports/single_cell_axis_30nn_motif_heatmaps_v2_20261013/`.

`91_ebf1_pm2mb_shared_support_da.R` tests differential accessibility for all
251 Bdev peaks overlapping Ebf1 TSS +/-2 Mb (mm10 chr11:42,617,317-46,617,317).
For each cytokine comparison, cells are limited to the maximal consecutive
original-bin interval in which every bin contains at least 10 accepted mapped
WT cells and 10 accepted mapped KO cells: bins 14-25 for Flt3L and 22-30 for
IL-7. The complete peak matrix is normalized with native scWorkbench
log1p-CP10K and each target is tested with native `de_wilcox`; BH adjustment is
within the 251-peak target test. Because there is one biological sample per
genotype, results are descriptive cell-level differences. Run with
`bash processing/integration/run_91_ebf1_pm2mb_shared_support_da.sh`. The
read-only, unpromoted report is under
`reports/ebf1_pm2mb_shared_support_da_v1_20261013/`.

`92_add_wt_ko_da_to_tcf3_placseq_tracks.R` adds the shared-support Flt3L and
IL-7 differential-accessibility results to the existing E2A/TCF3 degradation
coverage and Ebf1-promoter PLAC-seq display. The two aligned DA tracks show
accessible-cell fraction in WT minus KO for all 251 tested peaks. Peaks passing
BH < 0.05 within their cytokine comparison are highlighted red for WT-higher
and blue for KO-higher; nonsignificant peaks remain gray. Both the Ebf1 ±3 Mb
and chr11:44-46 Mb layouts are rendered. The source coverage and PLAC-seq
values are reused unchanged. Run with
`bash processing/integration/run_92_add_wt_ko_da_to_tcf3_placseq_tracks.sh`.
The completed unpromoted report is under
`reports/e2a_tcf3_placseq_with_shared_support_wt_ko_da_v2_20261013/`.

`93_ebf1_depth_complexity_matched_da.R` tests whether the strong directional
skew in the Ebf1-window DA is explained by sample-level ATAC depth and
complexity. Within every retained pseudotime bin, WT and KO cells are matched
without replacement on an oriented PC1 of log total counts and log detected
peaks across all 97,788 Bdev peaks. Matching requires common support and a
0.2-target-SD caliper. It retains 475 pairs for Flt3L and 255 for IL-7 while
reducing both covariate absolute standardized differences below 0.05. Native
scWorkbench log1p-CP10K and `de_wilcox` are then reapplied to the matched cells.
Only 5 Flt3L peaks remain BH-significant (1 WT-higher and 4 KO-higher), while
none remain for IL-7. The Ebf1 TSS-overlapping peak remains Flt3L WT-higher
(BH q=0.037) but not IL-7 significant (q=0.892). This is a sensitivity analysis,
not a substitute for biological replicates. Run with
`bash processing/integration/run_93_ebf1_depth_complexity_matched_da.sh`; the
completed report is under
`reports/ebf1_pm2mb_depth_complexity_matched_da_v1_20261013/`.

`94_ebf1_locus_matched_da_e2a_placseq.R` renders the final Ebf1 ±2 Mb locus
view from the matched sensitivity results. Separate control and E2A/TCF3-
depleted coverage tracks are followed by aligned matched Flt3L and IL-7
WT-minus-KO accessibility tracks and the Ebf1-promoter PLAC-seq arcs. The five
remaining Flt3L BH-significant peaks are colored and labeled at the Ebf1 TSS,
-0.43 Mb, +0.50 Mb, +1.19 Mb, and +1.51 Mb. The IL-7 track contains no
significant highlight. Original E2A/TCF3-dependent peaks retain their P1-P7
labels within this interval. Run with
`bash processing/integration/run_94_ebf1_locus_matched_da_e2a_placseq.sh`;
the completed report is under
`reports/ebf1_locus_matched_da_e2a_placseq_v1_20261013/`.

`95_bdev_five_region_ebf1_dynamic_peaks.R` partitions the preserved Bdev
trajectory into five consecutive local states around the current shared
WT/KO-supported cytokine domains: before Flt3L (bins 10-13), Flt3L only
(14-21), dual (22-25), IL-7 only (26-30), and after IL-7 (31-34). The two
four-bin flanks match the width of the dual interval. Across all 251 Bdev peaks
within Ebf1 TSS ±2 Mb, native scWorkbench `de_wilcox` tests each adjacent
state pair after full-matrix log1p-CP10K normalization. A peak is statistically
dynamic when any adjacent comparison has BH q < 0.05; the robust list also
requires an absolute accessible-cell-fraction change of at least 0.05 and a
concordant normalized-effect direction. This identifies 108 statistically
dynamic and 40 robust dynamic peaks. In the heatmap, Bdev accessibility below
5% is white and values from 5% upward use a sequential blue-to-purple scale.
Two aligned columns at right show the existing bulk E2A/TCF3 control and
depleted mean fragment-coverage CPM over each peak, using a separate orange
scale. Run with
`bash processing/integration/run_95_bdev_five_region_ebf1_dynamic_peaks.sh`;
the read-only, unpromoted report is under
`reports/bdev_five_region_ebf1_dynamic_peaks_v3_20261013/`.

`96_plot_two_dynamic_tcf3_peaks_bdev_umap.R` displays the two statistically
dynamic Bdev elements that overlap significant E2A/TCF3-depletion peaks:
`chr11-45089034-45090047` and `chr11-45847272-45848407`. All 26,723 Bdev
cells are shown on the unchanged native `umap.atac`. Accessibility is the
previously exported native scWorkbench log1p-CP10K value; gray denotes zero,
and both panels share a pooled nonzero-p99 display cap so color intensity is
comparable. Run with
`bash processing/integration/run_96_plot_two_dynamic_tcf3_peaks_bdev_umap.sh`;
the read-only report is under
`reports/two_dynamic_tcf3_peaks_bdev_umap_v1_20261013/`.

`97_tcf3_peaks_detection_normalized_wt_ko_ratio.R` evaluates all eight
significant E2A/TCF3-depletion peaks in the shared-support Flt3L bins 14-25
and IL-7 bins 22-30. For each cell and peak, its score is 10,000 times the
binary accessibility indicator divided by the number of detected peaks across
all 97,788 Bdev-aligned peaks. WT/KO is the ratio of genotype mean cell
scores. No pseudocount is used, so peaks absent from one genotype retain exact
WT-only or KO-only labels. The output also retains raw accessible-cell
percentages and cell-level scores. Run with
`bash processing/integration/run_97_tcf3_peaks_detection_normalized_wt_ko_ratio.sh`;
the read-only report is under
`reports/tcf3_peaks_detection_normalized_wt_ko_ratio_v2_20261013/`.

`98_ebf1_tcf3_bdev_dynamics_wtko_placseq_locus.R` integrates the eight
significant E2A/TCF3-depletion peaks in a focused Ebf1 locus view extending
from -0.25 to +2.10 Mb relative to the TSS, which retains P8 at +2.018 Mb.
Aligned panels show five-region Bdev accessible-cell percentages, Flt3L and
IL-7 detected-peak-normalized WT/KO ratios, original TCF3 depletion log2FC,
and 28 Ebf1-promoter H3K4me3 PLAC-seq arcs. P3 is marked as statistically
Bdev-dynamic and P7 as robustly Bdev-dynamic; P8 is explicitly marked as
outside the original ±2 Mb dynamic screen. Run with
`bash processing/integration/run_98_ebf1_tcf3_bdev_dynamics_wtko_placseq_locus.sh`;
the report is under
`reports/ebf1_tcf3_bdev_dynamics_wtko_placseq_locus_v1_20261013/`.

`figure_scripts/bdev_integrated/07_ko_cytokine_positive_signal_report.R`
combines the reviewed mapping with the accepted positive Flt3L and IL-7 score
state. It produces exactly two five-panel composite figures: accepted combined
flattened-rank score on native Bdev `umap.atac`, its Bdev trajectory, the
established five driver-gene trajectories, WT/KO cells on the reference-
anchored UMAP, and WT/KO query scores along mapped Bdev pseudotime. The same
raw-count score formula is applied to query RNA counts; both signatures have
complete query feature coverage. Query curves use the existing aligned 5%
windows and require 20 cells per genotype, leaving explicit gaps where support
is insufficient. BB1/BB2 are dashed in every trajectory panel. Run with
`bash processing/integration/run_ko_cytokine_positive_signal_report.sh`.
The completed unpromoted report is under
`reports/ko_cytokine_positive_signal_v1_20261001/`.

`figure_scripts/bdev_integrated/08_ko_cytokine_positive_distribution_report.R`
is the v2 display revision. It produces exactly two four-panel composites,
retains the accepted positive-score UMAP at a larger display size with its scale
on the right and its Bdev
trajectory, removes the driver panel, and keeps the reference-anchored WT/KO
UMAP with the v1 point styling and no redundant genotype legend. The two rows have
equal height. Panel D shares the exact horizontal span of panel B and shows a
centered three-bin rolling mean of genotype-normalized cell fractions over 100
underlying mapped-pseudotime bins rather than query scores. Deliberate whitespace
around A provides the alignment.
Supported and borderline query cells are included; out-of-reference cells are
withheld. BB1/BB2 are dashed in panels B and D. Run with
`bash processing/integration/run_ko_cytokine_positive_distribution_report.sh`.
The completed unpromoted report is under
`reports/ko_cytokine_positive_distribution_v2_20261001/`; v1 is preserved.

The same script supports the report-only `bdev_normalized` mode used for v3.
It enlarges panel C's projected points with v1-like styling and divides each
genotype's centered three-bin query fraction by the centered three-bin fraction
of the 12,642 Bdev cells with defined `pt`. Panel D is therefore a density ratio
where 1 denotes agreement with the Bdev reference distribution. Run with
`bash processing/integration/run_ko_cytokine_positive_distribution_report.sh <output_dir> bdev_normalized`.
The v3 output is under
`reports/ko_cytokine_positive_distribution_bdev_normalized_v3_20261001/`; all
earlier report versions remain preserved.

Adding the `gene_expression_panels` flag produces the v4 six-panel revision.
Panels A/B show `Flt3` for Flt3L or `Il7r` for IL-7 on native Bdev `umap.atac`
and as mean `SCT_data` expression in the saved aligned 5% windows. The v3
panels shift to C--F, and trajectory panels B, D, and F use the same horizontal
span and pseudotime limits. Run with
`bash processing/integration/run_ko_cytokine_positive_distribution_report.sh <output_dir> bdev_normalized gene_expression_panels`.
The output is under
`reports/ko_cytokine_gene_expression_bdev_normalized_v4_20261001/`.

`processing/integration/12_wt_enriched_peak_de.R` derives trajectory regions
directly from panel F of that v4 report: maximal contiguous WT bins whose
Bdev-normalized rolling fraction is strictly above 1. This yields one Flt3L
region and two separately labelled IL-7 regions. The same bounds mark WT and
KO supported/borderline query cells, and native scWorkbench `de_wilcox`
compares WT (group A) with KO (group B) on a staged log1p-CP10K RNA matrix.
Run with
`bash processing/integration/run_wt_enriched_peak_de.sh`.
The unpromoted report-only result is under
`reports/ko_wt_enriched_peak_rna_de_v1_20261001/`. Each genotype/cytokine has
only one sample, so the tests are explicitly descriptive cell-level contrasts,
not biological-replicate-level inference.
The report includes a two-panel visual audit of the exact shaded pseudotime
intervals used for DE, overlaid on the WT/KO Bdev-normalized rolling curves
with the ratio-1 reference and dashed BB1/BB2 boundaries.
The same report includes WT/KO expression violins for every significant gene
with `|avg_log1p_fc| > 0.2`. WT-higher and KO-higher genes are displayed in
separate panels, retaining up to the 10 smallest BH-adjusted P values in each
direction for every region.
The report additionally extracts the B- and T-lineage lists from `Table S1B`
of the cited supplementary workbook and shows all four KO-study samples in a
single two-panel violin figure. Exact mapping retains 246/251 unique B-lineage
genes and 233/234 T-lineage genes. With explicit user approval, the per-cell
score is calculated as the arithmetic mean of staged log1p-CP10K expression
because the native scWorkbench gene-set trajectory operation does not export
per-cell arithmetic-mean scores. Source rows, mapping coverage, cell scores,
and sample summaries are retained in the report.
The same approved calculation is also applied to the accepted Parse-derived
Flt3L-UP, Flt3L-DOWN, IL-7-UP, and IL-7-DOWN profiles. All 105 signature genes
map exactly to the KO RNA matrix. A single four-panel figure shows all four
samples in every profile; exact memberships, per-cell scores, coverage, and
sample summaries are exported.

The default GO trajectory analysis is the completed residual-UCell
dynamic-range permutation run produced by
`run_73_go_ucell_rna_feature_residual_dynamic_range.sh`. It residualizes each
native-compatible UCell score for `log1p(nFeature_RNA)` without pseudotime in
the adjustment model, summarizes the unchanged Bdev `pt` trajectory in 97
overlapping windows, and uses 1,000 permutations. The accepted downstream
pattern resolution is BP=25, MF=20, and CC=20 via
`figure_scripts/bdev_integrated/run_36_significant_region_residual_ucell_kmeans_k25_20_20.sh`.
Prior gene-set runs are retained for sensitivity and provenance, not selected
as the default.

The same selected protocol is complete for the 38 cytokine signatures in
`cytokine_ucell_rna_feature_residual_dynamic_range_perm1000_v1_20261012`.
For each signature, it fits
`UCell ~ ns(log1p(nFeature_RNA), df=4)` without pseudotime, carries the per-cell
residuals into the identical 5%-width/1%-step windows, and recomputes both the
1,000-permutation peak-to-trough test and trajectory-wide positive-region test.
The run finds 32 dynamically significant signatures, 56 supported positive
peaks, and 34 merged positive regions. It is produced by
`run_74_cytokine_ucell_rna_feature_residual_dynamic_range.sh`; native-score
reconstruction and bundle-fingerprint checks passed. The earlier residual
cytokine curves under
`cytokine_ucell_dynamic_range_perm1000_v1_20261012/diagnostics/` validate curve
shape and reduced RNA-complexity association only; their completion record
explicitly sets `permutation_calls_changed` to false and remains historical.

The predeclared WT-versus-KO GSEA region policy uses `Flt3l_shared_peak1` and
`Il7_shared_peak2` from the mapped shared-enrichment report. These intervals
retain thousands of mapped cells in both genotypes and will be matched by
pseudotime bin before constructing KO-minus-WT gene ranks. Bdev residual-UCell
cytokine regions are annotation context only: KO coverage falls to 35 versus
672 WT cells in the early FLT3L region and 20 versus 274 WT cells in the IL-7
region, consistent with failure to progress into those states. The late FLT3L
region contains only 17 query cells in total. It must not be used for primary
WT-versus-KO GSEA.

`processing/integration/75_mapped_enrichment_pseudotime_matched_gsea.R`
implements the completed analysis. It reads the retained scWorkbench
log1p-CP10K matrix from the existing session stage without mutation, calculates
WT and KO linear-CP10K means in each selected mapped-pseudotime bin, and ranks
genes by the Bdev-density-weighted mean bin-wise KO-versus-WT log2FC. The report
tests Hallmark, GO BP/MF/CC, and the 38 accepted cytokine signatures using
`fgseaMultilevel`; BH adjustment is separate by comparison and collection.
The runner is `run_75_mapped_enrichment_pseudotime_matched_gsea.sh`, and the
completed output is
`reports/flt3l_il7_ko_vs_wt_pseudotime_matched_gsea_v1_20261012/`.
Primary-versus-sparse-bin sensitivity correlations are 0.988 and 0.991 for
Hallmark NES in FLT3L and IL-7, respectively. Source bundle fingerprints are
unchanged and no results are promoted.

`processing/integration/76_mapped_enrichment_dynamic_change_gsea.R` extends
that matched-region analysis from an average KO-versus-WT effect to a
genotype-by-pseudotime change contrast. For every detected gene it calculates
`(WT late - WT early) - (KO late - KO early)` on bin-level
log2(mean linear CP10K + 1), then runs preranked GSEA on Hallmark, GO BP/MF/CC,
and cytokine sets. The primary endpoints are the first and last two bins with
at least 20 cells in each genotype; a separate sensitivity analysis includes
the sparse terminal bin. Gene-set curves center each genotype to its own early
baseline. Positive GSEA is only classified as failed induction when the
member-gene summary rises in WT and rises less in KO. The runner is
`run_76_mapped_enrichment_dynamic_change_gsea.sh`; completed output is under
`reports/flt3l_il7_ko_vs_wt_dynamic_change_gsea_v1_20261012/`.

`77_summarize_dynamic_go_top_terms.R` keeps the complete BP/MF/CC testing in
that parent run but produces a focused report with at most five significant
terms per GO group and comparison. Reporting is restricted to terms whose
member-gene summary rises in WT and rises less or declines in KO; each row also
shows the sparse-terminal-bin sensitivity result. Run with
`run_77_summarize_dynamic_go_top_terms.sh`. Output is under
`reports/flt3l_il7_ko_vs_wt_dynamic_change_top_go_v1_20261012/`.

`78_summarize_go_failure_to_silence.R` tests the reciprocal strict pattern in
the same complete GO results: a term must decrease in WT, change upward
relative to WT in KO, have negative differential-change NES, and pass FDR.
The report shows at most five terms per BP/MF/CC group and plots curves for the
reported terms that also pass sparse-terminal sensitivity. Run with
`run_78_summarize_go_failure_to_silence.sh`; output is under
`reports/flt3l_il7_ko_vs_wt_go_failure_to_silence_v1_20261012/`.

`79_individual_gene_divergent_trajectories.R` applies both strict trajectory
directions to individual genes. A gene must retain the sign definition under
both the reliable-bin and sparse-terminal endpoint policies; genes are ranked
by the smaller absolute divergence across the two. The report lists the top 20
genes and plots the top 12 per direction and comparison. Plot lines retain the
log2(mean CP10K + 1) trajectories used for selection, while ribbons show the
cell-level log2(CP10K + 1) IQR. The ribbons describe expression heterogeneity,
not confidence intervals. No gene-level P values are assigned because there is
one biological sample per genotype. Run with
`run_79_individual_gene_divergent_trajectories.sh`; output is under
`reports/flt3l_il7_ko_vs_wt_individual_gene_divergence_v1_20261012/`.

`processing/integration/13_ebf1_enhancer_il7_da.R` tests accessibility in all
259 quantified Bdev peaks overlapping the GRCm38/mm10 Ebf1 TSS ±3 Mb window
(`chr11:41,617,317-47,617,317`). It reuses the exact saved cell assignments for
the two IL-7 WT-enriched pseudotime regions above and keeps the regions as
separate WT-versus-KO contrasts. Native scWorkbench `de_wilcox` runs on the
`ATAC_bdev_mm10_peaks` count matrix with `min_pct=0`; BH adjustment is separate
within each region among peaks with finite tests. Run with
`bash processing/integration/run_ebf1_enhancer_il7_da.sh`. The completed,
unpromoted report-only output is under
`reports/ebf1_enhancer_il7_enriched_pt_da_v1_20261001/`. As above, one sample
per genotype makes this descriptive cell-level evidence rather than
biological-replicate-level genotype inference.

`processing/integration/14_ebf1_peak2_atac_tracks.R` renders one combined,
ten-locus ATAC track figure for the ten smallest BH-adjusted P values in the
IL-7 peak-2 Ebf1-window contrast. All four samples are restricted to the exact
IL-7 peak-2 pseudotime interval and to supported/borderline projected cells.
The report uses accepted reciprocal GRCm39 intervals, ±5 kb views, 25-bp Tn5
insertion bins, a centered 125-bp rolling mean, and normalization by the total
Cell Ranger ATAC-fragment depth of the selected cells in each sample. Run with
`bash processing/integration/run_ebf1_peak2_atac_tracks.sh`. The completed,
unpromoted report-only output is under
`reports/ebf1_enhancer_il7_peak2_atac_tracks_v1_20261001/`. This fragment-based
plotting step was explicitly approved outside scWorkbench because the supplied
version has no registered ATAC coverage-track operation.

`processing/integration/15_ebf1_enhancer_il7_lognorm_da.R` repeats the two
Ebf1-window IL-7 comparisons after native per-cell ATAC depth normalization.
scWorkbench `matrix_log_normalize` creates a staged log1p-CP10K matrix using
each cell's total across all 97,788 Bdev peaks; native `de_wilcox` then tests
the same 259 peaks in the same exact saved cells, with BH correction separately
within each pseudotime region. The report also compares every normalized result
with the earlier raw-count result. Run with
`bash processing/integration/run_ebf1_enhancer_il7_lognorm_da.sh`. The completed
unpromoted report is under
`reports/ebf1_enhancer_il7_enriched_pt_lognorm_da_v1_20261001/`. Normalization
addresses per-cell peak-matrix depth but does not overcome the one-sample-per-
genotype inference limitation.

`run_ebf1_peak2_lognorm_atac_tracks.sh` uses the normalized IL-7 peak-2 result
to select the ten significant WT-higher peaks with the smallest BH-adjusted P
values and renders the same four-sample ±5 kb track layout. Peak selection uses
only IL-7 WT versus IL7KO; Flt3L WT/KO are display-only context within the same
pseudotime interval. The v2 report is under
`reports/ebf1_enhancer_il7_peak2_lognorm_atac_tracks_v2_20261002/`; the earlier
raw-ranked v1 remains preserved.

`run_ebf1_peak2_matrix_aligned_atac_tracks.sh` is the corrected v3 track
workflow. It ignores the fragment file's multiplicity/read-support column, as
does Signac `FeatureMatrix`; normalizes each cell by its complete 97,788-peak
Bdev matrix total; applies log1p-CP10K; and averages fragment-record coverage
over all cells, including zeros. This makes track height directly comparable to
the normalized differential-accessibility analysis. The current interpretive
report is under
`reports/ebf1_enhancer_il7_peak2_matrix_aligned_atac_tracks_v3_20261002/`.
Earlier v1/v2 tracks remain preserved but are superseded for interpretation.

`run_il7_peak2_top_de_genes.sh` resumes the saved report-only RNA session used
for the IL-7 peak-2 `Il7r` follow-up. It reuses the complete 33,696-gene native
scWorkbench `de_wilcox` tables and the staged log1p-CP10K RNA matrix without
running a new operation. The report ranks genes separately by effect direction
for FLT3L WT versus FLT3LKO and IL-7 WT versus IL7KO, and reads only the forty
selected profiles for plotting. The completed unpromoted report is under
`reports/il7_peak2_top_de_genes_v1_20261002/`.

`run_il7_peak2_fc_filtered_de_genes.sh` reuses the same saved normalization and
complete Wilcoxon tables, calculates standard log2FC from mean linear CP10K,
and retains genes with genome-wide BH-adjusted P < 0.05 and
`|log2FC| >= log2(1.5)`. It creates paginated violin panels separately for WT-
higher and KO-higher genes in both IL-7 WT versus IL7KO and FLT3L WT versus
FLT3LKO, always within the exact IL-7 peak-2 cell interval. The completed
unpromoted report is under
`reports/il7_peak2_fc1p5_bh05_de_genes_v1_20261002/`.

`17_shared_enriched_peak_rna_de.R` defines new comparison regions from maximal
contiguous bins where both WT and KO Bdev-normalized rolling fractions exceed
1; IL-7 candidates are restricted to the original `Il7_peak2`. It reuses the
staged native log1p-CP10K matrix and runs new resumable scWorkbench Wilcoxon
tests across all 33,696 genes because the cell sets change. The report exports
complete normalized WT/KO mean CP10K, standard log2FC, P and genome-wide BH
tables, highlights the selected pseudotime regions, shows a new two-comparison
volcano, and plots up to 20 most significant genes in each direction passing
adjusted P < 0.05 and 1.5-fold change, with WT-higher and KO-higher genes in
separate figures. Run with
`bash processing/integration/run_shared_enriched_peak_rna_de.sh`; the completed
unpromoted output is under
`reports/ko_shared_enriched_peak_rna_de_v1_20261002/`.
The final figures use orange for WT, blue for FLT3LKO, and green for IL7KO.
Violin facet strips include BH-adjusted P values. Volcano plots label every
gene with BH-adjusted P < 0.05 and `|log2FC| >= 1`.
Targeted `Flt3` and `Il7r` violin figures are included for both shared-region
comparisons, with standard log2FC and genome-wide BH-adjusted P in each facet.
The two Parse-derived UP aggregate scores (Flt3L-UP and IL-7-UP) and the B- and
T-lineage aggregate scores are also plotted for both comparisons. Facets show
the WT-minus-KO mean difference, two-sided cell-level Wilcoxon P, and BH P
adjusted jointly across the eight requested score tests; these are descriptive
cell-level results because there is one biological sample per group.
Both aggregate-score families also have all-four-sample profiles along mapped
Bdev pseudotime. The lineage figure contains B and T panels; the cytokine
figure contains Flt3L-UP and IL-7-UP panels. Profiles use 100 underlying bins,
a centered five-bin cell-weighted mean with at least 20 cells, explicit gaps
for unsupported intervals, and dashed BB1/BB2 boundaries.

`52_flt3l_il7_ko_vs_wt_hallmark_gsea.R` uses the complete shared-region DE
tables for a user-approved analysis outside scWorkbench. This is necessary
because the available scWorkbench `pathway_enrich` operation is hypergeometric
over-representation analysis rather than preranked GSEA. Every gene with
nonzero mean CP10K in WT or KO is ranked by KO-over-WT standard log2FC;
positive NES therefore means KO-higher enrichment. `fgseaMultilevel` tests the
50 mouse-native MSigDB 2026.1.Mm Hallmark sets, which are frozen with the
report along with leading edges and full rankings. Run with
`bash processing/integration/run_52_flt3l_il7_ko_vs_wt_hallmark_gsea.sh`.
The completed unpromoted report is under
`reports/flt3l_il7_ko_vs_wt_hallmark_gsea_v1_20261010_run2/`. Results remain
descriptive because there is one biological sample per genotype.

`34_ebf1_putative_enhancers_three_way.R` evaluates all 259 Bdev peaks in the
mm10 Ebf1 TSS ±3 Mb window using three independent evidence routes. Route 1
maps GSE259355 peaks with `abs(log2FC) > 1` and DESeq2 BH-adjusted P < 0.05 to
the Bdev peak universe by direct interval overlap. Route 2 tests Ebf1
log1p-CP10K expression against log1p-CP10K peak accessibility in the exact 971
Bdev cells in `Il7_shared_peak2`; partial Spearman tests control for nonlinear
source `pt`, RNA depth, and ATAC depth, with BH correction over the 259 local
peaks. Route 3 selects GSE162858 230-238 H3K4me3 PLAC-seq interactions whose
5-kb anchor contains the Ebf1 TSS and maps Bdev peaks to the opposite anchor
after expanding each border by 20 kb. Run with
`bash processing/integration/run_ebf1_putative_enhancers_three_way.sh`.

The completed unpromoted report is under
`reports/ebf1_putative_enhancers_three_way_v2_20261005/`. It retains a complete
one-row-per-peak audit plus all one-to-many source mappings. Eight Bdev peaks
have TCF3 differential-accessibility support, 63 have PLAC-seq support, and no
peak has BH-significant partial correlation with Ebf1 expression in the
restricted IL-7 interval. The distal union contains 64 putative enhancers; four
are supported by both TCF3 loss and Ebf1-promoter contact, and none by all three
routes. No parent bundle or result was promoted.

`35_plot_ebf1_tcf3_peaks_bdev_umap.R` plots the eight route-1 peaks across all
26,723 Bdev cells on the unchanged native `umap.atac` embedding. It reads the
existing staged scWorkbench log1p-CP10K peak matrix and does not rerun
normalization or modify the reference. Small gray points denote zero
accessibility, while larger colored points show nonzero cells. Each peak's
positive values are capped at its nonzero 99th percentile and scaled
independently for the spatial display; the uncapped values are exported.
Run with
`bash processing/integration/run_plot_ebf1_tcf3_peaks_bdev_umap.sh`. The
completed unpromoted report is under
`reports/ebf1_tcf3_peaks_bdev_umap_v3_20261005/`. This v3 uses the revised
±20 kb PLAC-seq classification and retains the larger nonzero-cell display;
the earlier reports remain preserved.

`36_tcf3_ebf1_coverage_placseq_tracks.R` creates normalized bulk ATAC coverage
tracks for the six GSE259355 libraries. Each filtered paired-end BAM is
normalized independently to fragment-span coverage per million filtered
first-mate records, then the three replicate tracks are averaged with equal
weight within control and TCF3-depleted conditions. The depleted mean retains
the source design of two 6-hour and one 8-hour libraries. The full mm10 Ebf1
TSS ±3 Mb figure uses 1-kb bins, a centered 5-kb rolling mean, and all 28
230238 H3K4me3 PLAC-seq Ebf1-promoter links. A second figure shows ±10 kb
around each of the eight TCF3-positive Bdev peaks using 50-bp bins and a
centered 250-bp rolling mean. Condition-average 50-bp bedGraph tracks and all
replicate-level values are exported. This user-approved report-only coverage
calculation is outside scWorkbench because no registered coverage-track
operation is available; no bundle is modified or promoted. Run with
`bash processing/integration/run_tcf3_ebf1_coverage_placseq_tracks.sh`. The
completed report is under
`reports/tcf3_ebf1_coverage_placseq_tracks_v4_20261005/`. The v2 report added a
coverage-plus-PLAC-seq view restricted to mm10 chr11:44-46 Mb; v3 displays
the PLAC-seq arcs downward from their anchor baseline; v4 adds full-region and
44-46 Mb layouts with control and depleted means in separate stacked tracks
on a shared y-axis. The current in-place v4 display revision uses filled ATAC
areas, halves the ATAC panel-height allocation, and omits the PLAC distal-anchor
vertical bars. Earlier versions remain preserved.

`42_bdev_selected_motif_pt_heatmaps.R` creates separate motif-accessibility
heatmaps for `Flt3l_peak1` and `Il7_peak2`. For each region, the accessible
peak universe, chromVAR background and expectation, top-30 temporal-variability
selection, row ordering, centering, and scaling are all fixed from Bdev alone.
The mapped WT and knockout pseudobulks are evaluated only afterward and cannot
affect motif selection. Bins with fewer than 20 cells are blank, and the final
panel is descriptive WT-minus-KO in Bdev standard-deviation units. The approved
report-only workflow reads the exact pre-existing Bdev JASPAR match matrix
because the bundles contain no registered motif library/match set; it does not
modify or promote a bundle. Run with
`bash processing/integration/run_42_bdev_selected_motif_pt_heatmaps.sh`. The
completed output is under
`reports/bdev_selected_motif_pt_heatmaps_v1_20261010_run5/`. The displayed
heatmaps remain top-30 summaries, while run 5 also exports standardized
trajectories for all 746 eligible motifs for unbiased downstream screening.

`47_bdev_up_down_motif_trajectory_adjusted_wt_ko.R` screens all 746 eligible
motifs rather than imposing a top-30 limit. A strong Bdev directional pattern
requires an absolute last-third minus first-third change of at least 1.5 Bdev
temporal SD, absolute Spearman correlation with Bdev pseudotime of at least
0.6, and agreement between the two directions. WT-minus-KO is then calculated
within matched mapped-pseudotime bins and averaged using the Bdev cell
distribution, removing genotype occupancy differences along the interval.
All passing motifs are exported and shown in paginated up/down figures.
Leave-one-bin-out intervals measure trajectory-stratum stability, not
biological-replicate uncertainty. Run with
`bash processing/integration/run_47_bdev_up_down_motif_trajectory_adjusted_wt_ko.sh`.
The completed report is under
`reports/bdev_strong_directional_motif_trajectory_adjusted_wt_ko_v1_20261010/`.

`48_top_adjusted_wt_ko_motif_heatmaps.R` visualizes the 20 largest positive
and 20 largest negative trajectory-adjusted WT-minus-KO motif effects in each
cytokine region. Each figure aligns Bdev, WT, KO, and matched-bin WT-minus-KO
heatmaps across pseudotime. Values are in Bdev temporal-SD units, row labels
show the adjusted effect, and stars mark motifs passing the prespecified
effect and sign-consistency thresholds. Run with
`bash processing/integration/run_48_top_adjusted_wt_ko_motif_heatmaps.sh`.
The completed report is under
`reports/top_adjusted_wt_ko_motif_heatmaps_v1_20261010/`.

`51_single_cell_30nn_top_adjusted_motif_heatmaps.R` retains those same motif
rows and their order but replaces pseudotime bins with one displayed column
per Bdev cell, ordered by the preserved Bdev `pt`. At every Bdev focal
pseudotime, Bdev, WT, and KO values are independently averaged over that
group's 30 closest pseudotime cells; matched-position WT-minus-KO is then
calculated on the common focal axis. Query positions outside each genotype's
observed mapped-pseudotime support are audited. Figures retain only the largest
contiguous interval where both WT and KO 30-neighbor radii are at most 0.25
pseudotime units; full-axis values and support flags remain exported. Run with
`bash processing/integration/run_51_single_cell_30nn_top_adjusted_motif_heatmaps.sh`.
The completed report is under
`reports/single_cell_30nn_top_adjusted_motif_heatmaps_v1_20261010_run4/`.

`49_selected_motif_bdev_umap_query_violins.R` plots SNAI2 (MA0745.2), TCF3
(MA0522.3), RUNX1 (MA0002.2), SREBF2 (MA0596.1), and STAT5A::STAT5B
(MA0519.1). Bdev panels use native `chromvar_data` on the unchanged
`umap.atac`, with motif-specific symmetric blue-white-red scales capped at the
absolute p99 for display. Query-cell deviations are recomputed independently
for Flt3L peak 1 and IL-7 peak 2 using region-specific Bdev-defined accessible
peaks, background peaks, and expectation. Violin plots show raw per-cell WT
and KO distributions; facet labels add matched-bin WT-minus-KO effects weighted
to the Bdev pseudotime distribution. Run with
`bash processing/integration/run_49_selected_motif_bdev_umap_query_violins.sh`.
The completed report is under
`reports/selected_motif_bdev_umap_query_violins_v2_20261010/`.

`run_50_tcf3_ebf1_zeb_stat5_runx1_composite.sh` reuses the same read-only
workflow for TCF3, EBF1, ZEB-family, STAT5A::STAT5B, and RUNX1. The single tall
composite keeps every native Bdev UMAP in the left column and places the
corresponding Flt3L and IL-7 WT/KO violins directly opposite on the right.
Because the 746-motif set has no ZEB2 entry, ZEB1 (MA0103.3) is explicitly
labeled as a ZEB-family proxy. The completed report is under
`reports/tcf3_ebf1_zeb_stat5_runx1_composite_v1_20261010/`.

`43_selected_tf_motif_pt_curves.R` plots Bdev (black), WT (red), and knockout
(blue) motif-accessibility trajectories together for EBF1, PAX5, TCF3, CEBP,
HLF, BACH2, and LEF1. One fixed motif is selected per requested TF/family by
the highest mean Bdev temporal SD across `Flt3l_peak1` and `Il7_peak2`, so both
cytokine figures use exactly the same seven motif IDs. Selection and chromVAR
calibration remain Bdev-only, and bins with fewer than 20 cells are retained as
gaps. Run with
`bash processing/integration/run_43_selected_tf_motif_pt_curves.sh`. The
completed output is under
`reports/selected_tf_motif_pt_curves_v1_20261009/`.

`44_bdev_variable_gene_expression_pt_heatmaps.R` mirrors the motif heatmap
screen at the RNA level. It applies native scWorkbench full-matrix
log1p-CP10K normalization to Bdev and reuses the identically normalized query
matrix. For each of `Flt3l_peak1` and `Il7_peak2`, the top 30 genes are chosen
independently by temporal SD across Bdev bin means; WT and KO never affect
selection, clustering, centering, or scaling. Every nonempty bin is calculated,
including bins with fewer than 20 cells. Those low-cell bins have a black
border and cross, while zero-cell bins are gray. Run with
`bash processing/integration/run_44_bdev_variable_gene_expression_pt_heatmaps.sh`.
The completed output is under
`reports/bdev_variable_gene_expression_pt_heatmaps_v1_20261009_run4/`.

`45_query_vs_bdev_gene_expression_diff_heatmaps.R` reuses those fixed top-30
gene lists and Bdev-scaled bin means to calculate matched-bin WT-minus-Bdev and
KO-minus-Bdev heatmaps. Positive values indicate higher query expression and
negative values higher Bdev expression. Every calculable bin is included; a
black border and cross flags comparisons where either constituent bin has
fewer than 20 cells, while gray denotes a noncalculable comparison. Run with
`bash processing/integration/run_45_query_vs_bdev_gene_expression_diff_heatmaps.sh`.
The completed output is under
`reports/query_vs_bdev_gene_expression_diff_heatmaps_v1_20261009_run2/`.

`46_ko_specific_gene_expression_diff_heatmaps.R` starts from all 30,063 common
RNA features but restricts the candidate pool using Bdev alone. Within each
region, genes must have a Bdev maximum bin mean of at least 0.1, expression in
at least three Bdev bins, and temporal SD in the top quartile. The retained
genes are then screened for large KO-minus-Bdev but small WT-minus-Bdev
divergence in matched pseudotime bins. The selection score is the
reliability-weighted mean absolute KO-minus-Bdev difference minus its WT
counterpart. Primary hits require KO
mean absolute difference of at least 0.2 log1p-CP10K and a KO/WT divergence
ratio of at least 1.5. Every nonempty bin contributes; bins below 20 cells are
down-weighted, not excluded, and remain marked in the heatmaps. No ribosomal,
mitochondrial, or other technical-category filter is applied. Each figure puts
the Bdev mean-expression heatmap beside the WT-minus-Bdev and KO-minus-Bdev
heatmaps using the same genes, order, and pseudotime bins. Run with
`bash processing/integration/run_46_ko_specific_gene_expression_diff_heatmaps.sh`.
The completed output is under
`reports/ko_specific_bdev_variable_gene_expression_diff_heatmaps_v1_20261009/`.
