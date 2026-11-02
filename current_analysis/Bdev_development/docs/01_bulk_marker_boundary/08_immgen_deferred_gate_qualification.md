# ImmGen deferred-gate qualification

**Capability:** `01_bulk_marker_boundary`  
**Step:** `06_qualify_immgen_deferred_gates`  
**Status:** deferred gate mechanics qualified; full ImmGen execution required  
**Qualified:** 2026-10-31

## Scope

This step qualified the remaining stability, confounding-audit, and internal
validation mechanics using only checksummed synthetic outputs from the prior
core qualification. It preserved the exact 126-library cohort structure and
did not read an ImmGen expression value.

The machine-readable definition is
[`06_immgen_deferred_gate_qualification.yaml`](../../config/01_bulk_marker_boundary/06_immgen_deferred_gate_qualification.yaml),
and the workflow is
[`06_qualify_immgen_deferred_gates`](../../workflows/01_bulk_marker_boundary/06_qualify_immgen_deferred_gates/).

## Stability qualification

The balanced schedule selected exactly one library from each of the 61
discovery populations in each of 200 iterations, producing 12,200 audited
selections. A seeded permutation and cycling schedule kept the number of times
each library was selected within one of every other library from its
population.

Each iteration recalculated TMM factors and every required contrast effect.
Effect stability used normalized expression contrasts because a
one-library-per-population local-stage resample has no residual degrees of
freedom for a quasi-likelihood refit. The separate leave-one-out gates used
full robust edgeR quasi-likelihood refits, including filtering, TMM,
`glmQLFit`, `glmTreat`, FDR, expression, and direction checks.

All 834 injected markers passed their family-specific 0.80 stability
requirement. The maximum stability of a negative control was zero.

Population leave-one-out produced 1,147 gene-level results: 955 were
assessable and 192 were correctly labeled low-replication candidates. The
non-assessable targets were monocytes and the four B branches represented by
fewer than three populations. The three-population splenic transitional branch
was assessable.

Local leave-one-library-out produced 144 results: 126 were assessable and 18
for singleton `B.T3.Sp` were labeled low-replication. `B.T1.Sp` had already
failed the minimum set size without padding in step 05 and therefore supplied
no selected set to this gate.

## Confounding qualification

Population-by-organ and population-by-laboratory contingency tables were
generated from the frozen metadata. Every one of the 76 declared contrasts
received explicit target/competitor separation flags. Twenty-five contrasts
were completely separated by organ and 25 by laboratory; these flags are
evidence for later downgrading or independent-support review, not a claim that
the confounding was statistically adjusted.

Five synthetic controls exercised and passed both required decisions:

- reproduced tissue bias in an unrelated lineage rejects a core-identity
  candidate, while absent reproduced bias does not; and
- a completely separated comparison without independent support becomes
  `candidate_confounding`, while independent support or lack of complete
  separation passes the confounding gate.

## Internal validation qualification

Training-standardized broad-lineage scores achieved 100% top-one accuracy for
all 61 discovery populations under leave-one-population-out validation. Every
lineage sensitivity was 1.0, exceeding the contract thresholds of 0.80 overall
accuracy and 0.67 per-lineage sensitivity.

All three validation-only populations received their declared parent lineage:

| Population | Expected and predicted lineage |
|---|---|
| `B.Sp` | B |
| `Tgd.Sp` | gdT |
| `MF.PC` | myeloid |

Nineteen local-neighbor AUROCs were evaluated. Twelve had sufficient
replication and all had AUROC 1.0. Seven exercised the required
insufficient-replication candidate path.

## Acceptance boundary

This qualification demonstrates workflow mechanics, not biological marker
validity. It emitted no `accepted_core` set and no GMT-like import artifact.
scWorkbench was not invoked, and no scWorkbench core or operation-pack change
is required.

## Reproducibility and outputs

Generated artifacts are outside Git under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/06_qualify_immgen_deferred_gates/
```

A complete repeat reproduced all 12 tables and both reports byte-for-byte. Key
checksums are:

| Artifact | SHA-256 |
|---|---|
| `01_tables/01_balanced_resample_schedule.tsv` | `1a7ea0a8c7909de5b3a6eafe47cc70d4985fced1c91f26cb74c6553de9a80c24` |
| `01_tables/02_balanced_stability.tsv.gz` | `bcc7076f98378aaf9c3b2b0d0b62e2e4ba6b8fcdc018395670749d342422e84c` |
| `01_tables/03_population_loo_stability.tsv` | `a79599f49a1bcd47d5aad04fafb6579a0e40a9e9546fb3a249b1bb4899a7a298` |
| `01_tables/04_library_loo_stability.tsv` | `cb0e622a1cb5a8885f46c8832b089020bd18d1b47e744ff1946748d627a6afcc` |
| `01_tables/11_contrast_confounding_audit.tsv` | `b179db681a29604bb5eed9767349c831e485f03a1674bb2ebeaeb113870ba67b` |
| `03_reports/02_deferred_gate_qualification_summary.json` | `2723f2bfdcf28b57fb8389a52bae4931a479601e5a2e2809943d0273f1e723d6` |
| `06_manifests/01_deferred_gate_qualification_manifest.json` | `0e13a1a41c3f47a11c23cde279fc6c106c9c19c81557edd6de867a2f609e4bb4` |

The manifest closes the checksums for all configurations, prior synthetic
inputs, both workflow scripts, and every generated output.

## Next step

The real execution subsequently completed and is documented in
[`09_immgen_vm25_marker_execution.md`](09_immgen_vm25_marker_execution.md).
Its incomplete broad panel, failed internal validation, sparse fine-stage
coverage, and confounding flags require diagnostic and biological review before
any gene set can be accepted for scWorkbench import.
