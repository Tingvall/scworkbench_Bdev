# ImmGen-to-Bdev labeling acceptance contract

**Step:** `15_qualify_immgen_bdev_labeling`  
**Status:** technical scoring chain complete; flat 64-way labels not accepted  
**Collection:** `immgen_vm25_20261102_all_population_parse` (64 populations)

This step tests the selected diagnostic ImmGen collection on a deterministic,
bounded Bdev fixture. It does not execute the paper analysis and cannot accept
cell labels or unlock `Bdev_implementation`.

The archived `bdev_reference` bundle is read-only. Because the baseline runtime
requires a newer schema, execution uses an active-storage copy. Only that copy
is migrated and mutated. The workflow verifies the source manifest identity,
all 26,723 cell identifiers, and the exact authoritative `pt` values and missing
entries before scoring.

The fixture contains up to 50 deterministically selected cells from every
`peaks_snn_res.1.2` cluster by `pt`-defined stratum. This yields 1,746 cells and
retains both trajectory and non-trajectory portions of the reference.

The supplied scWorkbench labeling recipe defines the native chain:

1. import the pinned GMT into a read-only feature-set namespace;
2. screen matrix overlap, requiring at least five mapped genes and 50% source
   membership detected;
3. score included sets on `RNA_counts` with UCell (`max_rank = 1500`);
4. inspect raw score-family QC and cluster-grouped behavior;
5. percentile-calibrate each population score independently; and
6. reduce to provisional labels with `min_score = 0.8`, `min_margin = 0.1`,
   tie abstention, and the explicit label `unassigned`.

All raw scores, overlap exclusions, calibrated evidence, score margins, and
abstention reasons remain attached to the result. Singleton, two-fold, and
held-out-derived gene sets are reviewed separately. Winner cell sets are not
registered because these labels remain diagnostic.

## Result

The deterministic fixture contained 1,746 cells across 40 cluster-by-`pt`
strata. The workflow verified all 26,723 source cell identifiers and the exact
authoritative pseudotime digest
`2a9e9b9b8cd85daa00a85c3897cb2706844a306513efe70956bda7cd4745d6ee`
after clone-only schema migration. The archive manifest remained unchanged.

All 64 gene sets imported and passed overlap QC. Median detected membership was
100%, with a median of 44 detected genes per set. The native UCell, raw-family
QC, percentile calibration, and abstaining reduction operations all completed
and were promoted with their evidence into the development execution clone.

At `min_score = 0.8` and `min_margin = 0.1`, only 21 of 1,746 cells (1.2%)
received a provisional label. The other 1,725 were unassigned: 1,685 failed the
margin gate, 39 were ties, and one failed the score gate. Most of the 21 labels
were macrophage, T-cell, stromal, or other non-B programs; only one provisional
`proB.FrA.BM` call represented a B-development stage. The flat result therefore
does not provide biologically credible Bdev annotations.

Raw score families also showed substantial redundancy, including Spearman
correlations of 0.78 between paired NK programs, 0.75 between MPP3 and MPP4,
and 0.75 between CLP and pro-B FrA. Several expected developmental programs
were strongly correlated with authoritative `pt`; this supports retention of
the continuous scores but does not validate forced discrete winners.

## Decision

The complete 64-population collection remains the working scoring and lookup
default. Flat 64-way winner labeling is not accepted, thresholds must not be
relaxed, and no winner cell sets may be registered from this run.

The next labeling design should retain all 64 scores but introduce a reviewed
ontology-aware hierarchy: use broad non-B programs as exclusion/confounding
evidence, then reduce only among biologically eligible B/progenitor stages.
That scientific policy is a separate approval gate; it is not inferred from
this negative acceptance result.
