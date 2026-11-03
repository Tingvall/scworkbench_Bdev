# Broad-state-first ImmGen labeling exploration

**Step:** `16_qualify_immgen_bdev_hierarchy`  
**Run:** `immgen_bdev_hierarchical_exploration_20261102T1835Z`  
**Status:** exploratory comparison complete; no labels accepted

This bounded exploration responds to the flat 64-way labeling failure by
separating lineage eligibility from population competition. It uses the same
deterministic 1,746-cell Bdev fixture and preserves the exact authoritative
`pt` values and missing entries.

Nine existing, globally contrasted ImmGen broad-lineage panels are scored
first. A cell retains every broad state with a percentile-calibrated score of
at least 0.80 and within 0.10 of its highest broad score. Population
competition is then restricted to the 64 panels belonging to those eligible
states, with the unchanged 0.80 score and 0.10 margin gates.

The population step compares two variants:

1. the pinned original 64-set collection; and
2. a filtered variant removing a gene when it occurs in population panels
   from at least three distinct broad lineages.

The filter removed 158 genes. All 64 panels remained scoreable, retaining at
least 13 genes, so no population was lost through the overlap floor.

## Result

The broad gate retained one state for 741 cells and multiple borderline states
for 658 cells; 347 cells did not meet the broad score gate. Among all fixture
cells, 303 retained B and 312 retained Stem/Prog as eligible; these categories
did not co-occur under the declared 0.10 window.

The hierarchy substantially increased provisional assignment relative to the
flat run, from 21 cells to 357 cells (20.4%) with the original panels. Of these,
29 were assigned a B population and 86 a Stem/Prog population. The dominant
B call was pro-B FrBC (23 cells); the few later calls were B FrE, T1, T2, and
T3. Among cells with defined `pt`, stem/progenitor calls tended to occur at
earlier values and pro-B FrBC and transitional-B calls later, which is
directionally coherent but not independent validation.

Shared-gene filtering did not improve the primary exploratory outcome. It
assigned 332 cells (19.0%), 25 fewer than the original panels. The variants
agreed for 1,436 cells; filtering gained 139 assignments, lost 164, and changed
seven assigned labels. This instability and lower yield do not support making
the filtered panels the working default.

## Decision

Broad-state-first competition is a useful exploratory direction and is much
better behaved than flat 64-way reduction, but the present calls remain too
sparse and reference-dependent for acceptance. Keep the original 64 panels as
the scoring and lookup default. Do not adopt the simple cross-lineage
prevalence filter. A future refinement should target biologically reviewed
nuisance genes or reconstruct fine panels within the accepted hierarchy rather
than select a filter from this single-cell result.

No feature sets were registered, no scWorkbench artifacts were promoted, no
winner cell sets were created, and implementation remains locked. scWorkbench
requires a writable bundle handle for durable sessions, so the final run used
a copy-on-write session bundle; its registry fingerprint remained exactly
equal to the parent fingerprint
`16b3d4160164bce7ad152f6d821bf9bef17fa6fa8b957500dd890c5fa445985b`.

Two earlier exploratory attempts are retained as incomplete storage evidence:
the first was stopped during unnecessary import promotion, and the second
confirmed that durable sessions cannot be created from a read-only bundle.
Neither produced labels or modified the completed parent bundle.

The final checksum manifest is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/01_bulk_marker_boundary/16_qualify_immgen_bdev_hierarchy/immgen_bdev_hierarchical_exploration_20261102T1835Z/06_manifests/04_final_run_manifest.json
SHA-256: 4d11068462dd78868e397185e47b1b3219cce355525f8dbcb676443e757dd1d4
```
