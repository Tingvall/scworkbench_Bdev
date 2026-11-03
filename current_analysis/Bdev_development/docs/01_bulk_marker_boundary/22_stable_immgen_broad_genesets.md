# Stable ImmGen broad-lineage gene-set exploration

**Run:** `immgen_stable_broad_lineages_20261102T205257Z`  
**Status:** complete, diagnostic only; not accepted or promoted

This test used only ImmGen data for discovery. Replicates were collapsed to
population medians so that every ImmGen population had equal weight. Candidate
genes had to be consistently rank-visible within the target lineage and exceed
the typical value in every competing lineage. Final sets prioritized genes
repeatedly selected when each target population was left out. Bdev did not
affect gene selection or internal evaluation.

## Internal result

Nested leave-one-population-out winners were correct for 90/126 ImmGen
libraries (71.4%). The default-compatible score >= 0.80 and margin >= 0.10 gate
accepted 50/126 libraries, of which 45 were correct (90.0% precision).

The main failures were structural rather than threshold-only: conventional
alpha/beta T was correct for 9/24 libraries (37.5%), and the combined myeloid
state was correct for 4/16 (25.0%). These broad labels combine heterogeneous
programs that cannot be represented reliably by one universal positive-marker
set. Even among the 20 final genes per lineage, only 4 DC genes, 4 monocyte
genes, 7 myeloid genes, and 6 stroma genes were selected in at least 80% of
their target-population holdouts.

## Bdev result by sampleId

The fixed gate accepted 564/1,746 fixture cells, left 936 borderline, and left
246 unassigned for low score.

| sampleId | cells | accepted | borderline | unassigned | accepted B |
|---|---:|---:|---:|---:|---:|
| 1 | 741 | 242 | 433 | 66 | 0 |
| 2 | 286 | 93 | 153 | 40 | 9 |
| 3 | 396 | 111 | 199 | 86 | 44 |
| 4 | 323 | 118 | 151 | 54 | 58 |

Only 111/564 accepted calls were B, and sampleId 1 had no accepted B calls.
The large sample-to-sample shift and many accepted non-B calls are incompatible
with adopting these sets for Bdev labeling.

## Interpretation

Removing broad or unstable genes is not sufficient. A single positive gene set
per broad ImmGen lineage loses subtype-specific programs and still confuses
biologically adjacent categories. The next useful ImmGen-only design would
retain multiple subtype modules within each broad lineage and combine their
scores at the broad-state reducer, with explicit conflicting-lineage evidence
and an unknown outcome. This run does not justify changing the current default.

The exact Bdev `pt` values and missing entries were read without modification.
No feature set, label, score artifact, or bundle state was promoted.
