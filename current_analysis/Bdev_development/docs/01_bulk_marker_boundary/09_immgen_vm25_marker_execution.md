# ImmGen vM25 marker execution

**Capability:** `01_bulk_marker_boundary`  
**Step:** `07_execute_immgen_vm25_markers`  
**Run:** `immgen_vm25_20261031_seed20261031`  
**Status:** execution complete; candidates require scientific review  
**Executed:** 2026-10-31

## Scope

The synthetically qualified outside-scWorkbench workflow was executed against
the frozen GSE109125 GENCODE vM25 matrix. It used 126 approved steady-state RNA
libraries, including 121 discovery libraries from 61 populations and five
held-out libraries from `B.Sp`, `Tgd.Sp`, and `MF.PC`.

The machine-readable run definition is
[`07_immgen_vm25_marker_execution.yaml`](../../config/01_bulk_marker_boundary/07_immgen_vm25_marker_execution.yaml),
and the workflow is
[`07_execute_immgen_vm25_markers`](../../workflows/01_bulk_marker_boundary/07_execute_immgen_vm25_markers/).

## Model execution

All 76 declared contrasts completed. The four designs were full rank and
retained positive residual degrees of freedom:

| Family | Unit | Samples | Groups | Residual df | Genes passing `filterByExpr` |
|---|---|---:|---:|---:|---:|
| Superlineage | Population profile | 61 | 4 | 57 | 20,913 |
| Broad lineage | Population profile | 61 | 9 | 52 | 26,055 |
| B branch | Population profile | 10 | 5 | 5 | 19,496 |
| Local B stage | RNA library | 18 | 9 | 9 | 17,896 |

The complete-statistics output contains 4,180,836 rows: all 55,011 source genes
for every contrast, including filter failures, mapping outcomes, state flags,
identity exclusions, effects, FDR, support, and AUROC.

Primary thresholds retained 8,250 preliminary target-gene rows. After 200
balanced resamples and exact population/library holdouts, 1,301 rows remained
in size-valid review sets. Population support in both primary and holdout
models requires population CPM and expression in at least one underlying
library, as specified by the contract.

## Set outcomes

No set is accepted. The real classifications are:

| Family | Outcome |
|---|---|
| Superlineage | Stem/progenitor, lymphoid, myeloid/DC, and stroma each retain 150 review-ready candidates. |
| Broad lineage | B (120), ILC (94), Stem/Prog (120), and stroma (120) are `candidate_confounding`; abT, gdT, monocyte, myeloid, and DC fail minimum size. |
| B branch | Early lymphoid progenitor (75) and B1b (65) are `candidate_low_replication`; committed marrow B, splenic transitional B, and mature conventional B fail minimum size. |
| Local B stage | `proB.CLP.BM` (27) and `B.MZ.Sp` (40) are review-ready; `proB.FrBC.BM` (40) remains `candidate_gap_or_compartment`; the other six targets fail minimum size. |

In total, six sets are computationally review-ready, four are confounding
candidates, two are low-replication candidates, one is a gap/compartment
candidate, and 14 are rejected for insufficient eligible genes.

The broad B list begins with coherent B-cell genes including `Cd19`, `Cd79a`,
`Pax5`, `Vpreb3`, `Fcrla`, `Spib`, `Ebf1`, `Bank1`, `Blnk`, and `Cd22`.
However, the broad B comparisons are structurally separated by laboratory for
multiple required contrasts and therefore cannot be accepted from marker
plausibility alone.

The fine-stage results also require biological review. For example, the
`proB.CLP.BM` list includes `Hp`, `Tgfbi`, `Clec5a`, `Chil3`, and `Trem2`, which
may indicate myeloid or compartment-associated signal in the broad CLP gate
rather than a clean B-development identity program. This is one reason the
computationally passing label remains review-only.

## Internal validation

Only four of nine broad lineages produced size-valid sets. Scoring with this
incomplete panel gave overall leave-one-population-out top-one accuracy of
30/61, or 49.2%.

| Expected lineage | Sensitivity |
|---|---:|
| B | 1.00 |
| ILC | 1.00 |
| Stem/Prog | 1.00 |
| stroma | 1.00 |
| abT, gdT, monocyte, myeloid, DC | 0.00 |

The zero values primarily reflect absent size-valid gene sets for those
lineages, not evidence that the retained sets form a complete classifier.

Of the three validation-only populations, only `B.Sp` received its declared
parent lineage. `MF.PC` was scored as Stem/Prog instead of myeloid, and
`Tgd.Sp` as ILC instead of gdT. The broad panel therefore fails both its
overall and held-out validation requirements.

Five local neighbor comparisons were available for emitted local sets. All
had AUROC 1.0, but they cover only `proB.CLP.BM`, `proB.FrBC.BM`, and
`B.MZ.Sp`; they do not establish discrimination across the full B-development
sequence.

## Confounding audit

All 76 contrasts received explicit organ and laboratory separation audits.
Twenty-five were completely separated by organ and 25 by laboratory. The four
size-valid broad sets and B1b have at least one required contrast needing
confounding review. No claim of covariate adjustment is made.

## Acceptance decision

The results are insufficient for an accepted comprehensive ImmGen reference:

- the broad panel is incomplete and fails internal validation;
- most B branches and local developmental stages fail minimum set size;
- several retained sets are low-replication, gap-spanning, or structurally
  confounded; and
- some computationally passing fine-stage candidates have questionable
  biological specificity.

Accordingly, the run emitted no `accepted_core` status and no GMT-like import
artifact. scWorkbench was not invoked, and no scWorkbench core or operation
pack change is indicated.

## Reproducibility and outputs

Generated outputs are stored outside Git under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/07_execute_immgen_vm25_markers/
```

| Artifact | SHA-256 |
|---|---|
| `01_tables/04_complete_statistics.tsv.gz` | `893286b0b8a4932b96d477bfc259a75b2d2f41263ec744149a1879ea0625a4db` |
| `01_tables/19_candidate_genesets.tsv.gz` | `dd3e53660cec304c8f3ddf1ff90ba1f4760616ea01940e97815043fd0cd0440f` |
| `01_tables/20_set_summary.tsv` | `6236a623d1b78bc89bb0b58c72f964a66c6bc28c8f7b2d98e0642a3b010f8c6e` |
| `03_reports/03_execution_summary.json` | `93f9f209cc15f3d7db21c24d48b5f43268f1a20517d1220bc4a5ffe9ed7d8912` |
| `06_manifests/01_execution_manifest.json` | `bd334b9f93e47b18d0ffce24c729dc56c86a41c8bbb450e64a9b3c3dc58b3efb` |

The final manifest closes 38 input, script, and output checksum records with no
mismatches. The exact-support stability phase was rerun after enforcing the
underlying-library expression criterion and retained the same 1,301
validation-set rows.

## Next step

The requested attrition diagnostic is now complete and documented in
[`10_immgen_marker_attrition_diagnostic.md`](10_immgen_marker_attrition_diagnostic.md).
Do not import these candidates into scWorkbench. The remaining next step is a
structured biological and confounding review of retained and near-miss genes
before proposing any narrowly scoped contract revision.
