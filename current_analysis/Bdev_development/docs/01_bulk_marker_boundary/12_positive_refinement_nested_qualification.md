# Nested qualification of positive-only competitor refinement

**Capability:** `01_bulk_marker_boundary`  
**Step:** `10_qualify_positive_refinement_nested`  
**Run:** `immgen_vm25_20261031_positive_refinement_nested_qualification`  
**Status:** method screen complete; no accepted gene sets  
**Executed:** 2026-10-31

## Scope

This step tests whether the positive-only set-level method from step 09
generalizes when gene ranking and greedy refinement cannot see the evaluated
target unit. It does not qualify sparse targets by extrapolation and does not
create an importable gene-set collection.

The machine-readable definition is
[`10_positive_refinement_nested_qualification.yaml`](../../config/01_bulk_marker_boundary/10_positive_refinement_nested_qualification.yaml),
and the workflow is
[`10_qualify_positive_refinement_nested`](../../workflows/01_bulk_marker_boundary/10_qualify_positive_refinement_nested/).

## Frozen qualification rules

A target was evaluable only with at least three independent discovery units:
populations for superlineage, broad-lineage, and B-branch targets, or RNA
libraries for local B stages.

For every held-out target unit, the workflow:

1. refit edgeR using training units only;
2. rebuilt the positive candidate ranking using the approved expression,
   identity, support, and positive-effect gates, with FDR used as ranking
   evidence;
3. greedily refined the positive set using training rank-score profiles only;
   and
4. passed the fold only when the untouched target-unit score was strictly
   greater than every training competitor-group median.

A target passed the method screen only if every fold constructed at least the
family minimum, at least 80% of held-out units passed, and at least the family
minimum number of genes appeared in 80% or more of the fold-final sets.

No negative features, down modules, Bdev data, or full-data pilot membership
entered fold selection.

## Execution closure

Fourteen targets were evaluable and 13 were excluded for low replication. All
127 planned edgeR folds completed under the immutable SIF and pinned edgeR
4.6.3 library. Every fold constructed at least its family minimum. The nested
set procedure passed 125 of 127 held-out target units.

Eight targets passed both discrimination and membership stability:

| Family | Target | Held-out passes | Stable genes at ≥80% | Minimum | Median fold-set Jaccard |
|---|---|---:|---:|---:|---:|
| Superlineage | stem/progenitor | 5/5 | 41 | 40 | 0.531 |
| Superlineage | lymphoid | 39/39 | 70 | 40 | 0.785 |
| Superlineage | myeloid/DC | 12/12 | 40 | 40 | 0.552 |
| Superlineage | stroma | 5/5 | 57 | 40 | 0.499 |
| Broad lineage | B | 10/10 | 33 | 30 | 0.597 |
| Broad lineage | gdT | 6/6 | 35 | 30 | 0.462 |
| Broad lineage | stroma | 5/5 | 47 | 30 | 0.487 |
| Local B stage | `B.Fo.Sp` | 4/4 | 33 | 15 | 0.353 |

The broad B result improves on the fixed production set under the same fold
score check: the production set passed 9/10 population folds, while the nested
positive procedure passed 10/10 and retained 33 stable genes.

## Failures

Five broad targets discriminated their held-out units often enough but lacked
the required stable membership:

| Target | Held-out pass fraction | Stable genes at ≥80% | Required |
|---|---:|---:|---:|
| Stem/Prog | 1.000 | 23 | 30 |
| abT | 0.923 | 7 | 30 |
| ILC | 1.000 | 20 | 30 |
| myeloid | 1.000 | 18 | 30 |
| DC | 1.000 | 5 | 30 |

The abT failure included `NKT.Sp`, whose held-out score did not exceed gdT.
The low DC overlap confirms the population-specific instability seen in the
production leave-one-population-out analysis even though each optimized fold
could discriminate its own holdout.

Splenic transitional B failed nested discrimination: 2/3 populations passed,
below the 80% rule. The failed unit was `B.T3.Sp`, which ranked below the early
lymphoid progenitor competitor. Only five genes were stable in at least 80% of
fold sets, also below the 20-gene branch minimum.

The other 13 targets remain `not_qualified_low_replication`. This includes
broad monocyte, four of five B branches, and eight of nine local B-stage
targets. Full-data separation from step 09 does not override this result.

## Interpretation

The Parse-like positive set-level approach is supported for broad positive
programs, especially broad B, and for the four superlineages. It does not, with
the available ImmGen replication, solve most close B-development stages.

The discovery gdT procedure passes the nested screen, but the independent
validation-only `Tgd.Sp` population was still misclassified by the step-09
full panel. Therefore gdT cannot be promoted based on discovery cross-fitting
alone. Likewise, the four follicular-B libraries support a reproducible local
score within `B.Fo.Sp`, but they do not provide independent populations or a
complete adjacent-stage validation panel.

Method-screen success is not gene-set acceptance. Confounding, external
validation, and biological coherence remain separate requirements.

## Generated outputs

Generated artifacts are outside Git under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/10_positive_refinement_nested_qualification/
```

| Artifact | SHA-256 |
|---|---|
| `01_tables/01_fold_plan.tsv` | `8a8fe5d0302c1263d8f070fd9d8a6682444ca76ada953a21ec3357bb204efe34` |
| `01_tables/02_fold_candidate_pools.tsv.gz` | `2d1cda90fa8ff4fa812ea8b6b2c242790f32e251391ae19e57666f9781f398fc` |
| `01_tables/03_fold_model_audit.tsv` | `173194e23811b6d0eb5554d7700e0b98f5dc28a8e9bbe93b70cc61ff93330a12` |
| `01_tables/04_fold_results.tsv` | `3871c73cc98ed85506691d932dbccb3ed60b92edd295b2c6461c41f8d0a56f79` |
| `01_tables/05_fold_final_membership.tsv.gz` | `78d7d85569eb2240e3aefb0b2b66082223834e934a9b56cd9bcaefbb5aa3c09f` |
| `01_tables/06_fold_removal_trace.tsv.gz` | `7edadc499d851ca79d1c16e05fff08ff1ba5efa9c2ff47c450bdafd6abb1296a` |
| `01_tables/07_target_summary.tsv` | `f90cc6abf90b51f28a58d9f2d164bfffb26fe9c4a5e9c19dc79cb33ca415c3b7` |
| `01_tables/08_stable_membership.tsv.gz` | `47fc60db58604e34b7102ad2e2d49b535dbd33431d408cf031702e27ce9ccace` |
| `03_reports/01_nested_model_summary.json` | `f695dade2ac4adacfcb046e30b049bb2e94b2697689da2182d0c95efc5f8400f` |
| `03_reports/02_qualification_summary.json` | `30d7b7a9032a6a058f17ae1fac328960cdcc28315e77b5293ec97193d895dc4e` |
| `06_manifests/01_qualification_manifest.json` | `349e247aa18be7df39fb31b271de8f4c9dce2027d056f47faae8e801ead381ba` |

The closure check covered 24 input, script, and output records with no
mismatch. The 6.6 GB immutable SIF was checksum-verified during execution.

## Decision and next step

No set is accepted or imported, and the approved production contract is not
changed.

The next step should freeze and apply a positive-only biological/confounding
review policy to the stable memberships of the eight method-screen targets.
The review should prioritize broad B, the only broad immune lineage with both
10/10 nested discovery performance and a correct validation-only population.
It should separately reject or downgrade genes associated with compartment,
state, contradictory lineage, or unresolved source separation. Only after
that review should a candidate full-data membership rule be proposed.

