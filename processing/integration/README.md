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

Future KO mapping must use its own accepted inputs and inspected interfaces.
It is not part of the cytokine-scoring run and cannot change its parent inputs.
