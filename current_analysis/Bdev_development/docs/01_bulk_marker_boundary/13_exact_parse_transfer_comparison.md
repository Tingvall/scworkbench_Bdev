# Exact applicable Parse-method transfer and comparison

**Capability:** `01_bulk_marker_boundary`  
**Step:** `11_exact_parse_transfer_comparison`  
**Run:** `immgen_vm25_20261031_exact_parse_transfer_comparison`  
**Status:** diagnostic comparison complete; no accepted gene sets  
**Executed:** 2026-10-31

## Question

Does the latest Parse construction perform better on ImmGen than the strict
production selector and the earlier positive-only ImmGen refinements?

The machine-readable definition is
[`11_exact_parse_transfer_comparison.yaml`](../../config/01_bulk_marker_boundary/11_exact_parse_transfer_comparison.yaml),
and the workflow is
[`11_exact_parse_transfer_comparison`](../../workflows/01_bulk_marker_boundary/11_exact_parse_transfer_comparison/).

## Method transferred from Parse

The transfer follows the latest archived Parse construction where it is
defined for ImmGen:

1. fit three cross-fitted target-versus-equal-weight-all-other specificity
   models;
2. retain positive genes with log2 fold change greater than 0.5, using FDR as
   ranking evidence rather than a hard inclusion gate;
3. take the intersection of the three fold-specific top-100 lists, requiring at
   least 20 genes;
4. score raw ImmGen libraries with the Parse UCell capped-rank formula using
   `max_rank = 1500`;
5. greedily delete one gene at a time only when deletion strictly improves the
   standardized target-to-strongest-competitor gap, with a 20-gene floor; and
6. pass the mean route only when the target group median is strictly greater
   than every other group median in the family.

No single-cell expression value, label, embedding, or pseudotime was read.
ImmGen alone determined all memberships and refinements.

Two structural adaptations were unavoidable and frozen before execution.
ImmGen does not have matched three-replicate blocks across every group, so the
independent target populations or libraries were deterministically partitioned
into three folds. ImmGen bulk also has no matched control-cell distribution, so
the Parse responder-fraction route has no valid analogue and was not run.

The exact archived method and implementation are declared as read-only
provenance inputs from:

```text
archived_analysis/legacy_scworkbench_exploratory/processing/parse_cytokines/ALL_CYTOKINE_COMPETITOR_REFINEMENT_AND_BDEV_REPORTING.md
archived_analysis/legacy_scworkbench_exploratory/processing/parse_cytokines/17_flt3_il7_competitor_refinement_pilot.R
```

## Construction results

Fourteen of the 27 approved targets had at least three independent ImmGen
units. All 42 fold-specific edgeR fits completed, and all 14 eligible targets
formed a three-fold top-100 intersection of at least 20 genes. The other 13
targets remain explicitly non-constructible for low replication.

All 14 refined sets ranked their target group first on the complete discovery
profiles. This is an optimistic full-data test because the same discovery
profiles guided greedy refinement, matching the limitation of the archived
Parse review method.

| Family | Target | Intersection | Refined set | Discovery UCell rank |
|---|---|---:|---:|---:|
| Superlineage | stem/progenitor | 79 | 75 | 1 |
| Superlineage | lymphoid | 69 | 55 | 1 |
| Superlineage | myeloid/DC | 56 | 52 | 1 |
| Superlineage | stroma | 40 | 39 | 1 |
| Broad lineage | Stem/Prog | 78 | 63 | 1 |
| Broad lineage | B | 72 | 71 | 1 |
| Broad lineage | abT | 64 | 42 | 1 |
| Broad lineage | gdT | 66 | 34 | 1 |
| Broad lineage | ILC | 83 | 71 | 1 |
| Broad lineage | myeloid | 45 | 33 | 1 |
| Broad lineage | DC | 26 | 22 | 1 |
| Broad lineage | stroma | 36 | 35 | 1 |
| B branch | splenic transitional B | 42 | 30 | 1 |
| Local B stage | `B.Fo.Sp` | 90 | 28 | 1 |

Broad monocyte has only two independent target populations and therefore has
no exact-transfer set. The other twelve low-replication exclusions are four B
branches and eight local B stages.

## Common-score comparison

All memberships were rescored with the same UCell-style capped-rank score.
This comparison therefore differs from the original reports, which used each
method's native score.

| Method | Sets available | Full-data mean passes | Broad panel sets | Validation-only broad calls |
|---|---:|---:|---:|---:|
| Strict production | 13/27 | 12 | 4 | 1/3 |
| Positive full-data pilot | 27/27 | 22 | 9 | 3/3 |
| Positive nested stable memberships | 8/27 | 8 | 3 | 2/3 |
| Exact applicable Parse transfer | 14/27 | 14 | 8 | 3/3 |

On the 14 targets eligible for the exact transfer, it passed 14/14 full-data
mean tests, compared with 13/14 for the positive full-data pilot, 8/14
available and passing nested-stable sets, and 7/14 available and passing strict
production sets. The positive pilot failure under the common score was local
`B.Fo.Sp`, which ranked behind `proB.CLP.BM` and `proB.FrBC.BM`.

The positive pilot previously called `Tgd.Sp` incorrectly under its original
mean-expression-percentile score. Its unchanged membership calls `Tgd.Sp`
correctly under the common capped-rank score. Thus, its validation result changes
from the previously reported 2/3 to 3/3 because the scoring rule changed, not
because its genes were reconstructed.

## Independent broad validation

The three validation-only populations remained completely excluded from gene
ranking and refinement.

| Population | Expected | Strict production | Positive pilot | Nested stable | Exact Parse transfer |
|---|---|---|---|---|---|
| `B.Sp` | B | B | B | B | B |
| `MF.PC` | myeloid | B | myeloid | B | myeloid |
| `Tgd.Sp` | gdT | B | gdT | gdT | gdT |

The exact-transfer winning-score margins over the runner-up were 0.342 for
`B.Sp`, 0.239 for `MF.PC`, and 0.076 for `Tgd.Sp`. This is the best result in
the comparison, but three validation populations are too few for acceptance.

## Interpretation and decision

The comparison supports the user's intended construction principle: require a
positive per-gene effect and cross-fold ranking stability, then optimize and
test discrimination at the gene-set level rather than requiring every gene to
pass every pairwise FDR threshold.

The exact applicable Parse transfer is the strongest candidate method in this
comparison. It reconstructs eight of nine broad-lineage panels, passes every
eligible discovery target, and correctly classifies all three untouched broad
populations under the common score. The incomplete monocyte panel, limited
independent validation, full-data refinement optimism, unresolved confounding,
and absent biological review prevent acceptance.

No method is selected using single-cell performance. No set is accepted, no
GMT or other import artifact is created, and scWorkbench is not invoked. A
candidate method and its biological/confounding review policy must be frozen
before any evaluation on single-cell expression data.

## Generated outputs

Generated artifacts are outside Git under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/11_exact_parse_transfer_comparison/
```

| Artifact | SHA-256 |
|---|---|
| `01_tables/01_fold_plan.tsv` | `a410301ad1ddfbbe55e75f0f50c24bd99cd46822d33d5d36e7c2859beca18f5b` |
| `01_tables/02_fold_top100_rankings.tsv.gz` | `a00fead06dff1d532f3a48f3f034c5d4ff969c3751745064ef7628ccc5e949a6` |
| `01_tables/03_candidate_membership.tsv.gz` | `f48cb9a7b9886b6e76e74cea9c3f24bdefe65d7a90f4870647bbca432d7a6c4c` |
| `01_tables/04_model_audit.tsv` | `7ae482f9ff8663f01191fd03beefb71b8e7b535057bda233169e54f742c70927` |
| `01_tables/05_final_membership.tsv.gz` | `fa6d8cd2a02579cc6c4160bc8ef5f0c736224dbff6509fd5140105c78e8ca024` |
| `01_tables/06_refinement_trace.tsv.gz` | `261d0d5602939cee7d69b6c15fbdc52be3b5eb7e7118c2827fbad16de01d7fbf` |
| `01_tables/07_method_target_comparison.tsv` | `7b44731a2310cd7f150eb7c63f5e857be52a00ad6c03c2d38de131b9affac49e` |
| `01_tables/08_heldout_broad_comparison.tsv` | `c88521956ff40a16365e0146e460cddc8004a92de2cef380d3e3d1edad399811` |
| `01_tables/09_method_summary.tsv` | `486232456d26f1dd9007a3996eead02e656b53581615002ee68b8f41c9f8c7b5` |
| `03_reports/01_model_summary.json` | `e402119ff22f20bfe03a5b37b157bd39aa076f9bd24497e1bdac4d92bacf492b` |
| `03_reports/02_comparison_summary.json` | `721a022ba9d552989e5375d5eb4649f39d32562c829b6f228caaca985f0fbd80` |
| `06_manifests/01_comparison_manifest.json` | `8c272c3b6b65d7c9f12511880f74a7f0d6c6e8e7f5107cb03999f530f2733443` |

The closure audit verified 14 input, two script, and 11 output checksum
records. Every candidate membership exactly equals its three fold-specific
top-100 intersection, and every refined set is a subset of that intersection
with at least 20 genes.
