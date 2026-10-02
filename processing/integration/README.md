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
