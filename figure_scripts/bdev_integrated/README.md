# Integrated Bdev reports

`03_sliding_window_cytokine_report.R` renders the current pooled sliding-window
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
session is mounted read-only. The current report is
`intermediate_objects/bdev_integrated/reports/cytokine_pt_20260923T2231Z/sliding_window_review/REPORT.html`
under mirrored storage. The integration remains working and unaccepted pending
user review. Its analysis artifacts were already promoted before the user
requested report review; the pooled revision performs no further promotions.
