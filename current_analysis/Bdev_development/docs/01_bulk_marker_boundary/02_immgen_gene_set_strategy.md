# ImmGen steady-state lineage gene-set strategy

**Capability:** `01_bulk_marker_boundary`  
**Bdev step:** derive ImmGen lineage gene sets outside scWorkbench  
**Status:** cohort, count source, hierarchy, statistical contract, and analysis
preflight frozen; marker models not run  
**Approved:** 2026-10-31

## Aim

Derive stable reference gene sets for steady-state immune lineages represented
in the Bdev analysis. The gene sets should reflect lineage identity rather than
sex, experimental treatment, infection, engineered genotype, or a single fine
population.

Bulk RNA data and statistical analysis remain outside scWorkbench. Only the
versioned gene sets and their complete provenance are imported into
scWorkbench for scoring and labeling the Bdev single-cell dataset.

## Approved cohort

Use only the 64 steady-state RNA populations with population-matched ATAC
profiles defined in
[`01_immgen_study_sample_inventory.md`](01_immgen_study_sample_inventory.md).
All 64 are present in the selected GENCODE vM25 count matrix, which contains
126 corresponding RNA library columns. The difference from the 127 libraries
listed for these populations in paper Table S1 reflects changed replicate
coverage in the updated matrix and must remain explicit in the input manifest.

The selected vM25 cohort contains untreated, non-reporter reference
populations from the following Table S1 lineages:

| Lineage | Populations | vM25 RNA libraries |
|---|---:|---:|
| Stem&Prog | 5 | 10 |
| B | 11 | 21 |
| abT | 13 | 24 |
| gdT | 7 | 13 |
| ILC | 10 | 18 |
| monocyte | 2 | 4 |
| myeloid | 8 | 16 |
| DC | 3 | 8 |
| stroma | 5 | 12 |
| **Total** | **64** | **126** |

Do not admit any vM25 column solely because it is present in the updated
matrix. Exclude the 22 Table S1 populations carrying explicit treatment,
infection, sex-specific sampling, or engineered reporter/transgenic status,
and exclude all later vM25 populations outside the 64-population allowlist.
Do not use any `T.act` population for discovery.

## Approved count source

Use the GEO `GSE109125` GENCODE vM25 matrix
`GSE109125_Gene_count_table_GENCODE_vM25.csv.gz` as the primary expression
input. The file retrieved on 2026-10-31 had SHA-256
`cc513819170cbd2e6cfbff6415f0533b114e66fa1c3588f438f232ba491bcb90`,
205 library columns, and 103 population labels before cohort restriction.

Select columns by an audited allowlist that maps each vM25 label to paper
Table S1. The updated matrix includes later populations and changed replicate
coverage, so selecting all columns or relying only on lineage-like name
patterns is prohibited. The original paper-era count matrix is not a discovery
input; retain it only as a provenance and optional sensitivity reference. Do
not combine counts across annotation versions.

The initial population audit requires these aliases for selected populations:

| GENCODE vM25 label | Paper Table S1 RNA population |
|---|---|
| `Ep.MECHi.Th` | `Ep.MEChi.Th` |
| `FRC.CD140a+.Madcam-.CD35-.SLN` | `FRC.SLN` |
| `MF.102+480+.PC` | `MF.ICAM+480hi.PC` |
| `MMP3.48+.BM` | `MPP3.48+.BM` |
| `MMP4.135+.BM` | `MPP4.135+.BM` |

All 64 selected populations are represented after applying these aliases.
Eleven have a different library count from paper Table S1: `BEC.SLN`,
`DC.4+.Sp`, `DC.8+.Sp`, `FRC.SLN`, and `ILC3.CCR6+.SI` each gain one vM25
column; `MF.PC`, `T.4.Nve.Sp`, `T.DN4.Th`, `T.ISP.Th`, `preT.DN2a.Th`, and
`preT.DN3.Th` each lose one. This produces the net change from 127 to 126
selected libraries and must be resolved at the library/GSM level in the frozen
manifest.

## Roles of excluded populations

Excluded samples are not discovery replicates. Use them only for declared
sensitivity checks:

- compare `B.Fem.Sp`, `T.4.Nve.Fem.Sp`, and `MF.Fem.PC` with their matched
  reference populations to detect sex-driven candidate genes;
- test terminal Prdm1-GFP B populations and reporter Tregs only for lineage
  generalization;
- use stimulated and infected populations only to identify gene sets that
  collapse under activation or become dominated by treatment-response genes.

Exclude Y-linked genes and `Xist` a priori. Remove other candidates with a
strong effect in the available sex comparisons. Report that only three
lineages have an explicit male/female comparison in this atlas, so complete
sex generalizability cannot be established from these data alone.

## Analysis principles

1. Use raw library-level RNA counts with an explicit sample-design table.
2. Preserve fine population, lineage, tissue, laboratory, sex/sex-specific
   status, treatment, genotype, and replicate identifiers.
3. Treat fine populations as biological strata, not interchangeable
   replicates. Prevent populations with more libraries from dominating a
   lineage contrast.
4. Compare each target lineage with reviewed competitor lineages while
   retaining complete statistics for every tested gene.
5. Require marker direction and effect to be consistent across multiple fine
   populations within the target lineage.
6. Test stability by leaving out one fine population at a time. A marker that
   depends on one population, tissue, or laboratory should not enter the core
   lineage set.
7. Apply reviewed expression, effect-size, false-discovery, specificity, and
   minimum-population-support thresholds. Freeze these thresholds in
   configuration before producing the accepted gene sets.
8. Preserve full results and exclusion reasons; do not export only selected
   genes.

The exact differential-expression engine and design formula will be selected
after the count matrix and metadata have been audited. No normalized
population-average matrix may substitute for replicate-level raw counts in
the statistical analysis.

## Outputs

The outside workflow must produce:

- a checksummed library and population manifest;
- library-level QC and inclusion/exclusion tables;
- complete long-form differential-expression and specificity statistics;
- leave-one-population-out stability results;
- sex, treatment, and engineered-background sensitivity results;
- a versioned long-form gene-set table with source and selection provenance;
- a deterministic GMT-like export for scWorkbench import; and
- a run manifest recording inputs, dependencies, configuration, random seeds
  where applicable, and output checksums.

## scWorkbench acceptance boundary

No bulk-analysis operation or scWorkbench core change is required. Acceptance
inside scWorkbench begins only after the outside gene-set version is frozen.
The subsequent Bdev step will import that gene-set artifact, run feature
overlap QC, score the Bdev single-cell expression layer, and evaluate
winner/abstention behavior without importing the bulk samples themselves.

## Immediate next step

The vM25 RNA input manifest is frozen and passed the checks documented in
[`03_immgen_vm25_manifest_freeze.md`](03_immgen_vm25_manifest_freeze.md).
The population hierarchy and contrast neighborhoods are approved and validated
in
[`04_immgen_geneset_hierarchy.md`](04_immgen_geneset_hierarchy.md). Next,
the replicate-aware statistical contract documented in
[`05_immgen_marker_statistical_contract.md`](05_immgen_marker_statistical_contract.md)
is approved, and its annotation, state-gene, and software preflight passed as
documented in
[`06_immgen_analysis_preflight.md`](06_immgen_analysis_preflight.md). The next
the model and deterministic-selection core has passed the synthetic
qualification documented in
[`07_immgen_marker_workflow_core_qualification.md`](07_immgen_marker_workflow_core_qualification.md).
The deferred stability, confounding, and held-out validation gates passed the
synthetic qualification documented in
[`08_immgen_deferred_gate_qualification.md`](08_immgen_deferred_gate_qualification.md).
The real execution subsequently completed and is documented in
[`09_immgen_vm25_marker_execution.md`](09_immgen_vm25_marker_execution.md).
No set was accepted; target-level attrition and biological specificity now
require review.
