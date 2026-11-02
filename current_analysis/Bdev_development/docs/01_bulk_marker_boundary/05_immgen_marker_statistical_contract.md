# Approved ImmGen marker statistical contract

**Capability:** `01_bulk_marker_boundary`  
**Bdev step:** freeze the marker-discovery method before execution  
**Status:** approved; analysis preflight subsequently passed  
**Proposed:** 2026-10-31

## Scope

This contract defines how the frozen GENCODE vM25 counts would be converted
into reproducible superlineage, broad-lineage, B-branch, and local B-stage gene
sets. The bulk analysis remains outside scWorkbench. A synthetic model using
the exact cohort structure has qualified the model and selection core; no
ImmGen expression model has been run and no biological genes have been
selected.

The machine-readable contract is
[`02_immgen_marker_statistical_contract.yaml`](../../config/01_bulk_marker_boundary/02_immgen_marker_statistical_contract.yaml).

## Cohort roles

- All 126 approved libraries remain available for source QC.
- Discovery models use 121 libraries from 61 populations.
- The five libraries from `B.Sp`, `Tgd.Sp`, and `MF.PC` remain completely held
  out from model fitting and are used only after gene selection to test broad
  lineage generalization.
- RNA and ATAC remain population-matched, never replicate-paired.

## Analysis layers

Use two complementary count representations:

1. **Library level:** retain biological libraries for local B-stage contrasts,
   within-population consistency, and resampling. Fit a population-factor model
   so every target is compared separately with each declared local neighbor.
2. **Population level:** sum integer counts within each fine population and
   recalculate normalization. Use the resulting 61 population profiles for
   superlineage, broad-lineage, and B-branch contrasts with equal population
   weighting.

## Statistical engine

Use the edgeR v4 robust quasi-likelihood pipeline:

- apply `filterByExpr` separately for each analysis layer and design;
- calculate TMM effective library sizes;
- fit robust `glmQLFit` models; and
- use `glmTreat` against a minimum absolute log2 fold change of 0.5.

The robust quasi-likelihood pipeline is chosen for count modeling and moderated
gene-wise variability. Exact R, Bioconductor, edgeR, limma, annotation, and
dependency versions must be locked before execution.

## Marker requirements

### Superlineage and broad lineage

Broad genes must be expressed across multiple target populations, not merely
one subtype. Broad-lineage genes require:

- support in at least 60% of target populations;
- global target-versus-other-lineage log2 fold change at least 1;
- log2 fold change at least 0.75 against every declared closest competitor;
- global FDR at most 0.05 and each closest-competitor FDR at most 0.10; and
- required direction and effect in at least 80% of stability resamples.

Accepted broad sets contain 30–120 eligible genes. Superlineage sets use a
40–150 range and analogous global and stability requirements.

### B branches

Branch genes require target support across the branch, log2 fold change at
least 0.75 globally and 0.5 against each sibling branch, FDR requirements, and
80% stability. Accepted sets contain 20–75 genes. The single-population B1b
branch can only produce a candidate set.

### Local B stages

For each of the nine local B targets, a positive marker must have log2 fold
change at least 0.75 and FDR at most 0.10 against every declared neighbor. Sets
contain 15–40 genes and are never padded with failing genes.

`B.T1.Sp` and `B.T3.Sp` each have only one selected vM25 RNA library, so their
sets can only be candidates. Sets depending on the missing `preB.FrD` interval
or the bone-marrow-to-spleen edge also remain candidates until orthogonally
validated.

## Stability and confounding

Run 200 deterministic balanced resamples with seed `20261031`, selecting one
library per population on every iteration. Broad and branch sets also undergo
leave-one-population-out testing when at least three target populations exist.

Do not claim that laboratory or organ effects are adjusted when they are
structurally confounded with population identity. Instead:

- report laboratory and organ contingency tables;
- flag completely confounded contrasts;
- screen tissue-spanning candidates for the same tissue bias in unrelated
  lineages; and
- downgrade unresolved sets to candidate status.

## Gene policy

The exact GENCODE vM25 GTF must be checksummed before execution. Full statistics
retain all 55,011 matrix rows and their filtering and mapping outcomes. Core
identity sets are positive, uniquely mapped protein-coding genes and exclude:

- Y-linked genes and `Xist`;
- mitochondrial, ribosomal-protein, and hemoglobin genes;
- immunoglobulin and TCR variable, diversity, and joining segments; and
- predefined cell-cycle, interferon, immediate-early, heat-shock, and stress
  genes.

State-related genes remain available in complete statistics and may be exported
only in separately named auxiliary state or transition sets. Their exact
versioned lists must also be frozen before execution.

## Deterministic selection and validation

Eligible genes are ranked by stability, target-population support, worst-case
effect, target-versus-competitor AUROC, worst required FDR, and finally gene
symbol. A set is not emitted as accepted when it contains fewer than its
minimum required number of genes.

Selection and threshold tuning use ImmGen only. Bdev expression, labels, and
the preserved `pt` field cannot influence gene selection. Internal validation
uses held-out population profiles and the three broad-validation populations.

## Structural validation

The contract validator passed and confirmed:

- 126 selected source libraries, 121 discovery-model libraries, and five
  held-out validation libraries;
- 61 discovery populations and all three declared held-out populations;
- four contrast families;
- nine local B targets, including singleton `B.T1.Sp` and `B.T3.Sp`;
- 200 balanced resamples with the frozen seed; and
- no requested scWorkbench core or operation-pack change.

## Approval and next action

The contract was approved on 2026-10-31. The subsequent preflight froze and
validated the exact GENCODE vM25 GTF, state-gene exclusion lists, and
outside-analysis software lockfile. The model, deterministic selection,
stability resampling, confounding audits, and held-out validation then passed
synthetic qualification. Full execution against the frozen ImmGen matrix then
completed as documented in
[`09_immgen_vm25_marker_execution.md`](09_immgen_vm25_marker_execution.md).
No set was accepted because the broad panel and fine-stage coverage failed the
frozen acceptance requirements.

## Method sources

- edgeR User's Guide:
  <https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf>.
- Chen et al., *From reads to genes to pathways: differential expression
  analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood
  pipeline*:
  <https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html>.
