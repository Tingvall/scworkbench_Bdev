# Integrated Bdev reports

`05_competitor_aware_sliding_window_report.R` is the current report source. It
uses the unpromoted native UCell session created by
`processing/integration/06_score_competitor_aware_signatures.R` and profiles all
33 accepted competitor-aware full signatures. Bdev gene coverage is retained
as a diagnostic: 24 pass the former five-gene/80% transfer rule and all 33 are
plotted. The report uses the same pooled 5%/10%/15% sliding-window design as the
last iteration and writes overview PNG/PDF files plus a multi-page PDF with one
page per cytokine. Run both stages through
`processing/integration/run_competitor_aware_trajectory.sh` with run ID
`competitor_aware_pt_all33_20260925`.

`03_sliding_window_cytokine_report.R` renders the earlier pooled sliding-window
review from the working `bdev_integrated` bundle. It reuses retained UCell scores
and runs native scWorkbench distribution summaries on overlapping cell subsets.
Run it through
`bash processing/integration/run_stage.sh figures RUN_ID`.

It writes an offline HTML report and PNG/PDF figures with continuous pooled
median curves and 25th–75th percentile bands (IQR). The primary window is 10%
of the observed pseudotime range, advanced in 1% steps; 5% and 15% widths are
included for sensitivity. Windows use at least 30 cells and are clipped at the
observed endpoints. Native `cell_vector_family_qc` computes spread; the reusable
implementation is `src/sliding_window_signature_summary.R`. Exact membership,
window bounds, cell counts, quartiles and SD are exported. These bands show
cell-to-cell variability, not confidence intervals; neighboring windows share
cells. Source `pt` and missingness are preserved. Unsupported centers are marked
gray, with connecting line segments used only for display.

The report also includes an all-gene-set 5% window median/IQR overview in PNG
and PDF. `04_add_five_percent_iqr.R` adds it from saved native summaries through
`bash processing/integration/run_stage.sh figures_5pct RUN_ID`, without rerunning
analysis. The standard report script also generates this overview on a full run.

All source samples and cell-cycle states are pooled. No sample or phase split
is analyzed or displayed in the current revision. The earlier
`01_cytokine_pseudotime_report.R` and fixed-bin `02_pooled_cytokine_report.R`,
with their generated reports, remain historical records. The fixed-bin report
is linked from the sliding-window report for comparison.

Only the report directory is writable during this stage; every bundle and
session is mounted read-only. The earlier report is
`intermediate_objects/bdev_integrated/reports/cytokine_pt_20260923T2231Z/sliding_window_review/REPORT.html`
under mirrored storage. The latest report is
`intermediate_objects/bdev_integrated/reports/competitor_aware_pt_all33_20260925/trajectory_review/REPORT.html`.
The integration remains working and unaccepted pending user review. The latest
scoring session and trajectory review perform no promotion.

`06_peakvi_reference_umap_projection.R` is the report-only KO mapping display.
It reads the native Bdev `umap.atac` through `wb_space()` and reuses the exact
Bdev-restricted PeakVI neighbor table from mapping QC. Query coordinates are
distance-weighted Bdev-neighbor means at k=30; k=15/50 projections, neighbor
UMAP spread, and cross-k displacement are retained as diagnostics. Native Bdev
coordinates are never changed, and calibrated out-of-reference query cells are
withheld from the primary display. The runner is
`processing/integration/run_peakvi_reference_umap_projection.sh`.

`07_ko_cytokine_positive_signal_report.R` generates the final two-image Flt3L
and IL-7 review. Each composite contains panels A-E: the accepted positive
combined score on Bdev UMAP, Bdev score trajectory, established top-five driver
trajectories, the reference-anchored WT/KO projection, and the WT/KO combined-
score trajectory along mapped Bdev pseudotime. All trajectory panels redraw
the saved BB1/BB2 medians as dashed lines. Unsupported query windows remain
gaps. The runner is
`processing/integration/run_ko_cytokine_positive_signal_report.sh`.

`08_ko_cytokine_positive_distribution_report.R` is the four-panel v2 display
revision. A enlarges the accepted Bdev positive-score UMAP and moves its scale
to the right;
B retains its score trajectory; C shows the reference-anchored WT/KO
projection with the v1 point styling and no genotype legend; and D, sharing B's
exact horizontal span, shows a centered three-bin rolling mean of supported/borderline
WT and KO cell fractions over 100 underlying mapped-pseudotime bins. The two rows
have equal height. Deliberate whitespace around A provides the alignment. The driver and
query-score panels are omitted. Panels B and D include dashed BB1/BB2 lines.
The runner is
`processing/integration/run_ko_cytokine_positive_distribution_report.sh`.

With the `bdev_normalized` flag, the same script generates v3 without changing
v2. Panel C uses enlarged v1-like projected points. Panel D divides the centered
three-bin query fraction by the centered three-bin Bdev reference fraction over
the same 100 bins; its dotted line at 1 marks equal query and Bdev density. The
12,642 Bdev cells with defined original `pt` form the denominator, and their
values and missing entries remain unchanged.

With both `bdev_normalized` and `gene_expression_panels`, the script creates the
six-panel v4 report. A/B display `Flt3` or `Il7r` `SCT_data` on the native Bdev
ATAC UMAP and as a mean trajectory using the exact saved 5%-wide windows. The
former A--D panels become C--F. B, D, and F share one horizontal span and exact
pseudotime limits, with BB1/BB2 dashed in every trajectory panel.

`09_il7_peak2_top_de_genes.R` resumes the saved IL-7 peak-2 RNA session and
renders four separate ten-gene violin panels: WT-higher and KO-higher for each
of FLT3L WT versus FLT3LKO and IL-7 WT versus IL7KO. Selection uses the already
completed 33,696-gene native scWorkbench Wilcoxon tables, ranked within
direction by genome-wide BH-adjusted P, raw P, absolute effect, and feature ID.
No normalization or differential test is repeated.

`10_il7_peak2_fc_filtered_de_genes.R` calculates standard log2FC from the saved
log1p-CP10K matrix by returning values to linear CP10K before group averaging.
It intersects a 1.5-fold threshold with genome-wide BH-adjusted P < 0.05 from
the saved native Wilcoxon tables, then plots every retained gene in paginated
violin panels for both IL-7 and FLT3L WT/KO comparisons within `Il7_peak2`.
