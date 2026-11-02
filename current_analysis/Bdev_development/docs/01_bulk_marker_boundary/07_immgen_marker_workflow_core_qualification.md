# ImmGen marker-workflow core qualification

**Capability:** `01_bulk_marker_boundary`  
**Step:** `05_qualify_immgen_marker_workflow`  
**Status:** model and deterministic-selection core qualified; deferred gates required  
**Qualified:** 2026-10-31

## Scope

This step qualified the outside-scWorkbench marker-analysis core using
synthetic integer counts and the exact frozen ImmGen cohort structure. It used
all 126 selected library identities, including the 121 discovery libraries
from 61 populations and the five held-out libraries. It did not read any
ImmGen expression value.

The machine-readable qualification definition is
[`05_immgen_marker_workflow_qualification.yaml`](../../config/01_bulk_marker_boundary/05_immgen_marker_workflow_qualification.yaml),
and the executable workflow is
[`05_qualify_immgen_marker_workflow`](../../workflows/01_bulk_marker_boundary/05_qualify_immgen_marker_workflow/).

## Qualified behavior

The workflow exercised:

- exact selected-column, discovery, and held-out-role handling;
- integer aggregation to one profile per population for superlineage, broad
  lineage, and B-branch models;
- library-level local B-stage models;
- full-rank designs with positive residual degrees of freedom;
- all 76 declared global, competitor, branch, and local-neighbor tests;
- `filterByExpr`, TMM normalization, robust `glmQLFit`, and `glmTreat` in the
  pinned R environment;
- complete output for filtered and tested synthetic genes;
- annotation eligibility, deterministic ranking, size limits, and no-padding;
  and
- exclusion of low-count, state-flagged, and ambiguous-mapping controls.

The synthetic fixture contained 890 genes and 126 libraries. All four model
families passed their design checks:

| Family | Unit | Samples | Groups | Rank | Residual df |
|---|---|---:|---:|---:|---:|
| Superlineage | Population profile | 61 | 4 | 4 | 57 |
| Broad lineage | Population profile | 61 | 9 | 9 | 52 |
| B branch | Population profile | 10 | 5 | 5 | 5 |
| Local B stage | RNA library | 18 | 9 | 9 | 9 |

The run produced 67,640 complete statistic rows. It recovered every injected
marker for every target that met its family minimum. The deliberate B.T1.Sp
control had 10 eligible markers against a minimum of 15; the workflow emitted
no set and did not pad it. The lymphoid and myeloid/DC superlineage controls
had 204 and 152 eligible genes, respectively, and were deterministically
truncated to the configured maximum of 150. No low-count, state, or ambiguous
decoy was selected.

## Acceptance boundary

This is an execution qualification, not a biological result. Synthetic
selection rows are labeled `synthetic_selected_not_importable`. The workflow
emitted no `accepted_core` set and no GMT-like import artifact. scWorkbench was
not invoked, and no scWorkbench core or operation-pack change is indicated.

The following gates remain mandatory:

- 200-iteration balanced stability resampling;
- leave-one-population-out and leave-one-library-out stability;
- tissue and laboratory confounding audits;
- held-out scoring thresholds; and
- full execution against the frozen ImmGen count matrix.

## Reproducibility and outputs

Generated artifacts are stored outside Git under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/05_qualify_immgen_marker_workflow/
```

A complete repeat reproduced every generated table and report byte-for-byte.
Key checksums are:

| Artifact | SHA-256 |
|---|---|
| `01_tables/01_synthetic_counts.tsv.gz` | `491736aa85ded20302fa05fb13e4d455421cb1c5d914640022eace42e97c0b12` |
| `01_tables/06_synthetic_complete_statistics.tsv.gz` | `d7f9ed0fa3b1fffa713bd9f08788b458949955754f50abd342500c1d2cfb8cfd` |
| `01_tables/08_synthetic_selected_genesets.tsv` | `59f4f459a46d601583f836fd3c0b53f0f71bd9ef0a4a8fd8a6e4acedfc946c08` |
| `03_reports/01_synthetic_qualification_summary.json` | `0828bd2764a744692abf35f86c8842fa0651713ccfea12d6767b996ddbc1c14c` |
| `06_manifests/01_synthetic_qualification_manifest.json` | `8697106f22d911db674e9a74d953e7c427e4d60ae8fbb926d0368b06e43108df` |

The final manifest records the qualification configuration, statistical
contract, hierarchy, software lock, all three workflow scripts, and every
generated output.

## Next step

The deferred stability, confounding, and held-out validation gates subsequently
passed the synthetic qualification documented in
[`08_immgen_deferred_gate_qualification.md`](08_immgen_deferred_gate_qualification.md).
The real execution subsequently completed and is documented in
[`09_immgen_vm25_marker_execution.md`](09_immgen_vm25_marker_execution.md).
No set was accepted because the real broad and fine-stage results require
diagnostic and biological review.
