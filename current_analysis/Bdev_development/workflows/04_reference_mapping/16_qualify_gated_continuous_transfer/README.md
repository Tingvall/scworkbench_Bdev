# Qualify gated continuous reference transfer

This bounded development workflow qualifies the generic
`scwbReferenceMap` 0.2.0 continuous-acceptance policy before any new full KO
mapping. It clones the existing balanced 28-cell fixture bundle, fits a new
20-dimensional PeakVI model to the frozen depth-adjusted HVP20k Bdev matrix,
and projects the already quantified query fixture without refitting the
reference model.

Reference calibration uses five deterministic folds balanced within sample,
complete L4-domain annotation, and native-L4 trajectory strata. A reference
cell is `in_domain` only when native L4 is finite and its fitted L4 weight is
at least 0.5; every other reference cell is explicitly `out_domain`.
Continuous L4 is accepted only when general reference support passes, the
30-neighbor annotation vote is `in_domain`, its confidence is at least 0.51,
and the raw estimate is finite. Raw estimates and an explicit rejection reason
remain recorded for every query.

This is a qualification fixture only. It performs no genotype comparison,
does not map all KO cells, does not write to `Bdev_implementation`, and does
not build or accept a SIF.

Status: passed on run `gated_continuous_transfer_20261122T192735Z`.
Calibration covered all 26,723 reference cells across five folds, achieved
0.821 domain-label accuracy at 30 neighbors, and had native-L4 median absolute
error 0.568. Of the 28 balanced query cells, all had an explicit raw-estimate
row, seven passed the full gate, and 21 were withheld because their local vote
was `out_domain`. No query was out of general reference support and no value
was rejected solely for low confidence. The source bundle manifest remained
unchanged.

The first transfer invocation stopped before producing mapping values because
the scWorkbench retained-parameter contract requires calibrated support
thresholds to be passed explicitly. The clean-run script now does so. The run
resumed from the already promoted 20D model and calibration using action 3;
the final inventory records the resumed path.

Run in order:

1. `01_prepare_runtime_and_clone.py --config <config>`
2. `02_fit_20d_calibrate_and_transfer.R --config <config>` inside the configured SIF and run library
3. If an older scWorkbench invocation stops because explicit retained support
   thresholds were omitted, resume with `03_resume_gated_transfer.R --config
   <config>`; a clean execution of step 2 does not need this action.
4. `04_finalize.py --config <config>`
