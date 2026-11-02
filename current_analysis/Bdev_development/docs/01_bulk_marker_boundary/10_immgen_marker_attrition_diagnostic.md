# ImmGen marker attrition diagnostic

**Capability:** `01_bulk_marker_boundary`  
**Step:** `08_diagnose_immgen_marker_attrition`  
**Run:** `immgen_vm25_20261031_seed20261031_attrition_diagnostic`  
**Status:** complete; biological review required  
**Executed:** 2026-10-31

## Scope

This read-only diagnostic reconstructs every primary marker gate from the
complete statistics produced in step 07, then follows the preliminary genes
through balanced stability, leave-one-out stability, and the deterministic
size rule. It does not refit a model, modify a threshold, use Bdev single-cell
data, or create an importable gene set.

The machine-readable definition is
[`08_immgen_marker_attrition_diagnostic.yaml`](../../config/01_bulk_marker_boundary/08_immgen_marker_attrition_diagnostic.yaml),
and the workflow is
[`08_diagnose_immgen_marker_attrition`](../../workflows/01_bulk_marker_boundary/08_diagnose_immgen_marker_attrition/).

## Closure checks

The diagnostic streamed all 4,180,836 complete-statistics rows covering 27
targets, 76 contrasts, and 55,011 source genes per target. It exactly
reproduced:

- all 8,250 primary candidate rows, including gene identity for every target;
- every target's step-07 emitted candidate count; and
- the low-replication rule that retains otherwise eligible genes as candidates
  without treating them as accepted.

This establishes that the attrition results below are an audit of the frozen
production selection, not an alternate implementation.

## Main diagnostic result

The shared expression and annotation filters remove many source rows, as
expected, but they do not explain why particular targets fail while others in
the same family succeed. The discriminating losses occur later:

| Target group | Decisive result |
|---|---|
| Broad `abT`, `gdT`, monocyte, and myeloid | Each falls below the 30-gene minimum at the all-competitor FDR gate. The most restrictive comparisons are `abT` vs `gdT`, `gdT` vs `abT`, monocyte vs myeloid, and myeloid vs monocyte. |
| Broad DC | Has 130 primary genes and 119 balanced-stable genes, but only 13 pass population leave-one-out stability; the minimum-size failure is therefore a population-robustness failure rather than a primary-threshold failure. |
| Committed marrow B branch | The all-competitor FDR gate reduces 126 remaining genes to zero; committed marrow B vs splenic transitional B is the most restrictive declared contrast. |
| Splenic transitional and mature conventional B branches | The global FDR gate reduces both targets to zero. With only ten B populations across five branch groups and five residual degrees of freedom, these branch-level positive programs are not supported under the frozen design. |
| Local `B.T1.Sp`, `B.T2.Sp`, `B.T3.Sp`, and `B.Fo.Sp` | Each reaches the neighbor-effect gate but has zero genes after requiring FDR against every neighbor. The zero-pass pairs include T2 vs T3, T3 vs follicular, and follicular vs T3. |
| Local `B.FrE.BM` | Falls below the 15-gene minimum at the neighbor FDR gate: seven primary genes and three after leave-one-library-out stability. Its most restrictive comparison is against `B.T1.Sp`. |
| Local `proB.FrA.BM` | Has 22 primary and 20 balanced-stable genes, but only three pass leave-one-library-out stability. This is a robustness failure, with `proB.FrA.BM` vs `proB.CLP.BM` the most restrictive primary comparison. |

The two low-replication B-branch candidates remain unchanged: early lymphoid
progenitor retains 509 genes before truncation to 75, and B1b retains 65, but
their required population holdouts cannot be assessed. They are not accepted
sets.

The size-valid superlineage, broad B/ILC/Stem&Prog/stroma, and selected local
sets also remain candidates. Attrition success does not resolve the step-07
confounding, validation, compartment-gap, or biological-specificity findings.

## Interpretation

The failed close-stage targets are not explained by a single global filtering
mistake. Two distinct limitations are present:

1. Closely related targets often have genes with the required direction and
   effect but lack FDR support against every declared neighbor or competitor.
   This is especially clear for transitional/mature splenic stages and the
   paired broad T or myeloid lineages.
2. DC and `proB.FrA.BM` have enough primary genes but lose most of them under
   exact leave-one-out testing, indicating population- or library-specific
   rather than robust target identity signal.

It would therefore be unsafe to solve the missing panels by lowering one
universal cutoff. Any contract change must instead be target/family justified
and requalified. Following the user's subsequent decision, any alternative
must remain positive-only; positive multigene, competitor-aware scoring is
evaluated separately in step 09.

## Generated outputs

Generated state is stored outside Git under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/08_diagnose_immgen_marker_attrition/
```

| Artifact | Purpose | SHA-256 |
|---|---|---|
| `01_tables/01_attrition_funnel.tsv` | Sequential and independent gate counts for every target | `57aeffba56335058194bf53e0b8754b424d1cfad1dc375de2c684905bc89ff03` |
| `01_tables/02_first_failure_reasons.tsv` | Mutually exclusive first-failure counts | `4e7902f78154658def1800b8d410ed7dd7f6939b7188d12760a53554c9bce023` |
| `01_tables/03_contrast_bottlenecks.tsv` | Effect and FDR pass counts for all 76 contrasts | `2b14ab5124741d85a10e9f3661a039b553d38ac0b98c73af2b49ca301e4170ac` |
| `01_tables/04_post_primary_attrition.tsv` | Stability, holdout, minimum, and truncation losses | `282aa5ffe382077e9665c8ef4aca48be25331fe040a125a9dece25e21137158e` |
| `01_tables/05_target_diagnosis.tsv` | Concise target-level bottleneck classification | `b3ef90ae8766b33974dcb83697fa6ceb7238481004690969084111fbc05b15a0` |
| `03_reports/01_diagnostic_summary.json` | Machine-readable closure summary | `5a3f77652c4138a090637c1bd1c393860f0aa9fe0b4a14b0c6f1c5a774978306` |
| `06_manifests/01_diagnostic_manifest.json` | Input, script, and output checksums | `1f75984497dcbef982ad56f4d4e9aefa450827ab1ac26f077f81b578f9c83e6e` |

## Decision and next step

No thresholds are revised, no candidate is promoted, and nothing is imported
into scWorkbench.

The requested positive-only Parse-like comparison is now complete and
documented in
[`11_positive_competitor_refinement_pilot.md`](11_positive_competitor_refinement_pilot.md).
It supports further positive set-level qualification, while confirming that
biological/confounding review and independent selection validation are still
required before any narrowly scoped contract revision.
