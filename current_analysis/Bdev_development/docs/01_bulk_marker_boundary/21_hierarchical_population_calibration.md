# Hierarchical population-specific calibration test

**Step:** `19_hierarchical_population_calibration`  
**Run:** `immgen_global_local_hierarchy_calibration_20261102T201305Z`  
**Status:** diagnostic test complete; no labels accepted

This test combined the nine-state broad gate with the global-local population
programs. Fine competition was limited to eligible broad states. Population
thresholds were then estimated from held-out ImmGen predictions by requiring
at least two predictions with perfect observed precision.

## Reference result

The broad gate retained the correct lineage for 115/118 held-out libraries
(97.5%). Without a margin, 101 libraries received non-tied fine winners and 78
were correct (77.2% precision). The fixed 0.10 margin accepted 12/118, all
correct.

Only 32/64 populations had enough prediction evidence to estimate an
exploratory population-specific threshold. Applying those thresholds back to
the same held-out inventory accepted 67/118 calls, all correct. This is
resubstitution performance, not independent validation: the thresholds were
selected and evaluated on the same small inventory.

## Bdev fixture comparison

On the bounded 1,746-cell fixture:

- no margin produced 1,227 fine calls, 172 broad-only calls, and 347 fully
  unassigned cells;
- the fixed 0.10 margin produced 312 fine calls, 1,087 broad-only calls, and
  347 fully unassigned cells;
- population-specific thresholds produced 385 fine calls, 1,014 broad-only
  calls, and 347 fully unassigned cells.

The calibrated policy therefore recovers 73 fine calls relative to the fixed
margin while retaining broad-only fallback. It includes 104 B-lineage fine
calls versus 32 with the fixed margin. Bdev has no independent fine-label
truth in this test, so these counts measure coverage, not accuracy.

## Decision

The hierarchy is useful: broad-state eligibility is robust and provides a
safe fallback. Population-specific calibration is promising but cannot be
accepted from two or three observations per population or from resubstitution
precision. Retain broad states and scores as diagnostic evidence. Fine labels
require a larger independent calibration inventory or a separate nested
validation design.

No feature sets were registered, no scWorkbench artifacts were promoted, and
no accepted labels or winner cell sets were emitted. The exact Bdev `pt` field
was inherited unchanged from the checksum-verified parent bundle. The final
manifest is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/01_bulk_marker_boundary/19_hierarchical_population_calibration/immgen_global_local_hierarchy_calibration_20261102T201305Z/06_manifests/03_final_run_manifest.json
SHA-256: 88cbff73366ca656005ef9a7197bd0c249221f7eac5736880fc6f10061c281e1
```
