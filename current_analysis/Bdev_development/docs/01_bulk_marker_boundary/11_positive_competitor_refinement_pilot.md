# Positive-only competitor-refinement pilot

**Capability:** `01_bulk_marker_boundary`  
**Step:** `09_positive_competitor_refinement_pilot`  
**Run:** `immgen_vm25_20261031_positive_competitor_refinement_pilot`  
**Status:** diagnostic improvement; not qualified or accepted  
**Executed:** 2026-10-31

## Question

Could the positive multigene, strongest-competitor strategy used in the
archived Parse cytokine work perform better than requiring every ImmGen gene
to pass FDR against every competitor or neighbor?

The pilot was deliberately restricted to positive genes. It constructed no
down sets, used no negative features, and performed no subtraction of a
negative module.

The machine-readable definition is
[`09_positive_competitor_refinement_pilot.yaml`](../../config/01_bulk_marker_boundary/09_positive_competitor_refinement_pilot.yaml),
and the workflow is
[`09_positive_competitor_refinement_pilot`](../../workflows/01_bulk_marker_boundary/09_positive_competitor_refinement_pilot/).

## Parse-derived idea and ImmGen adaptation

The archived Parse workflow ranked positive response genes, scored a multigene
set, identified its strongest competing condition, and greedily removed genes
only when removal improved the standardized target-to-competitor score gap.
The most useful transferable principle is therefore set-level discrimination,
not the Parse up/down pairing.

This ImmGen diagnostic retained the approved expression filter, annotation,
target support, and every positive effect threshold. FDR was used as ranking
evidence instead of a hard candidate-pool gate. Each target contributed at
most 100 positive candidates. Profiles were scored by the mean within-profile
expression percentile rank of the positive genes. One gene at a time was
removed only when doing so strictly improved the standardized gap between the
target median score and the strongest declared competitor median score. The
approved family minimum size was the deletion floor.

Population-level families used counts summed within discovery populations
before rank scoring. Local stages used discovery RNA libraries directly.

## Results

The positive-only pilot produced full-data positive separation and at least
80% target-unit support for all 27 targets. It produced 14 diagnostic panels
where the production workflow emitted no size-valid set:

- broad `abT`, `gdT`, monocyte, myeloid, and DC;
- committed marrow B, splenic transitional B, and mature conventional B
  branches; and
- local `proB.FrA.BM`, `B.FrE.BM`, `B.T1.Sp`, `B.T2.Sp`, `B.T3.Sp`, and
  `B.Fo.Sp`.

For every production panel except `proB.CLP.BM`, the optimized pilot set had a
higher discovery-data gap under this pilot's rank-score metric. That result is
expected to be optimistic because the same discovery profiles ranked genes,
guided deletion, and measured the final gap.

Thirteen of the 27 target results remain explicitly low-replication. The Parse
workflow had three biological samples per condition and could intersect three
fold-specific top-100 rankings. ImmGen has only one or two target populations
or libraries for many B branches and local stages. A faithful three-fold Parse
construction is therefore impossible for those targets, and the apparent
full-data separation cannot substitute for replication.

## Held-out broad-lineage result

The three validation-only populations were excluded from candidate discovery
and refinement. The positive pilot improved their top-one broad-lineage calls
from 1/3 to 2/3:

| Population | Expected | Production | Positive pilot | Pilot correct |
|---|---|---|---|---|
| `B.Sp` | B | B | B | yes |
| `MF.PC` | myeloid | Stem&Prog | myeloid | yes |
| `Tgd.Sp` | gdT | ILC | myeloid | no |

This is real evidence in favor of the set-level approach, particularly for the
myeloid panel, but three held-out populations are too few for qualification.
The `Tgd.Sp` failure also shows that discovery separation does not guarantee a
portable lineage classifier.

## Biological review observations

The pilot broad B panel contains coherent positive B-cell genes such as
`Vpreb3`, `Cd79a`, `Pou2af1`, `Blnk`, `Cd79b`, and `Tnfrsf13c`. This supports
continued evaluation of the method.

Other pilot panels need stronger scrutiny. For example:

- the committed marrow B panel includes stromal or compartment-associated
  genes such as `Il33`, `Reln`, `Thbs2`, `Alpl`, and `Akap12`, alongside
  B-development genes such as `Vpreb1`, `Igll1`, and `Lrmp`;
- the early lymphoid progenitor panel includes genes associated with several
  non-B immune programs, consistent with the broad CLP input population; and
- the local `B.T2.Sp` panel contains numerous poorly characterized, sensory-
  receptor, or structural genes and cannot yet be treated as a coherent
  transitional-B identity program.

Thus, optimizing a rank-score gap can amplify compartment or source-population
peculiarities as well as genuine lineage biology. The method reduces the
single-gene FDR bottleneck but does not by itself solve biological specificity,
confounding, or sparse replication.

## Decision

The Parse-like positive approach performs better as a diagnostic classifier:
it recovers all missing discovery panels and improves the genuinely held-out
broad result from 1/3 to 2/3. It is not yet a qualified gene-set definition.

No production threshold or confidence status is changed. No set is accepted,
no GMT or other import artifact is created, and scWorkbench is not invoked.
The pilot remains separate from the approved production results.

## Generated outputs

Generated artifacts are outside Git under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/09_positive_competitor_refinement_pilot/
```

| Artifact | SHA-256 |
|---|---|
| `01_tables/01_candidate_pool_membership.tsv.gz` | `7d0f19bc2058f120eceda4475efd6193d943c0811a62b66c0ab765beec3502f4` |
| `01_tables/02_refinement_summary.tsv` | `103ffeed0246e529ce160cd0a6f771393607c23c2baa3fb49999e8db11e16344` |
| `01_tables/03_final_membership.tsv.gz` | `24d7370a0d01fd45ac9ec599c47df490590836be7a96e1f834f881d08593384f` |
| `01_tables/04_removal_trace.tsv` | `0896db49ce03748034b616ca1408fa0644e2df7f13c09f0bffa2233e5c395dfa` |
| `01_tables/05_profile_scores.tsv.gz` | `c552ac418e793d6849ddd98ddb410bccd0ca26f55677b8cc5b213394863281c7` |
| `01_tables/06_production_comparison.tsv` | `7aaaa21b01b56e1838e7a41c0e41c0b0cbd870fd63e1d0be15085eec9c49751a` |
| `01_tables/07_heldout_broad_predictions.tsv` | `effe4b61660ca9e1a764615960ed88a64dd0929c9808286dac9dc58a4e374502` |
| `03_reports/01_pilot_summary.json` | `ab0f32fa1d88c99ccbbea1f4791a3527689c5bcadf01cc3857d8d890dce11e8a` |
| `06_manifests/01_pilot_manifest.json` | `1324037605a9b43ad9e7cd49909ae0c9ae6875746157174181a923facea6c865` |

The closure audit checked 20 input, script, and output records with no checksum
mismatch.

## Proposed next step

Retain the positive-only constraint and qualify the approach in two tiers:

1. For targets with at least three independent populations or libraries,
   construct positive pools independently in held-out folds, refine only on
   training profiles, and measure target rank in the untouched unit.
2. Keep targets with fewer than three independent units as low-replication
   candidates regardless of full-data separation. They may be represented by
   broader positive stages, but should not be promoted through optimization.

The nested positive-only qualification is now complete and documented in
[`12_positive_refinement_nested_qualification.md`](12_positive_refinement_nested_qualification.md).
It supports the procedure for eight targets but leaves 13 targets
non-qualifiable for low replication and six evaluable targets below the method
screen. A frozen biological/confounding review is required before proposing a
candidate full-data membership rule.

