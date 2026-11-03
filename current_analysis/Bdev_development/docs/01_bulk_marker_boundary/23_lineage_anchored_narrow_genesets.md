# Lineage-anchored narrow ImmGen gene-set exploration

**Run:** `immgen_lineage_anchored_narrow_20261102T213318Z`  
**Status:** complete, diagnostic only; not accepted or promoted

All 64 ImmGen populations remained direct labeling options. Each population
received a target-specific cross-lineage anchor module and a within-lineage
local module. Candidate scoring used twice the anchor score plus the local
score minus the strongest anchor from an incompatible lineage. No broad label
was assigned and no broad gate was used.

All gene selection, score-weight comparison, and target-membership holdout
evaluation used ImmGen only. Bdev was scored only after the modules and rule
were frozen.

## ImmGen internal result

The unthresholded winner was the exact population for 102/118 held-out
libraries (86.4%) and belonged to the correct lineage for 114/118 (96.6%). The
fixed rule requiring anchor score >= 0.80, non-negative anchor contrast, and
composite margin >= 0.10 accepted 78/118 libraries. All 78 accepted winners
had the correct lineage, and 75/78 had the exact population. The three accepted
errors were confined to neighbouring populations in the same lineage.

This result is exploratory rather than independent validation: competing
programs were derived from the complete atlas, and the same holdout inventory
was used to compare composite weights.

## Bdev result

| sampleId | cells | accepted | borderline | unknown | accepted B |
|---|---:|---:|---:|---:|---:|
| 1 | 741 | 219 | 165 | 357 | 6 |
| 2 | 286 | 56 | 82 | 148 | 16 |
| 3 | 396 | 99 | 95 | 202 | 39 |
| 4 | 323 | 80 | 79 | 164 | 28 |
| **Total** | **1,746** | **454** | **421** | **871** | **89** |

Only 89/454 accepted calls were B-lineage. SampleId 1 contained 219 accepted
calls but only six accepted B calls; its most frequent accepted lineages were
Stem/Prog, myeloid, stroma, and ILC. Therefore the strong within-ImmGen result
did not transfer to Bdev.

## Interpretation

The architecture successfully reduced cross-lineage errors inside ImmGen, but
the selected target-specific anchors still contain ImmGen-specific tissue,
laboratory, maturity, and population-isolation signals. Those signals can be
highly discriminative within the source atlas without being portable to Bdev.
Further margin adjustment is not justified. A subsequent ImmGen-only attempt
would need to explicitly remove genes associated with organ or laboratory and
require marker stability across the available contexts, while acknowledging
that many narrow populations are confounded with a single organ or laboratory
and therefore cannot support that identification from this reference alone.

The exact Bdev `pt` values and missing entries were preserved. No feature set,
label, score artifact, or bundle state was promoted.
