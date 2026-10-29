# Proposed analysis strategy

**Status:** provisional foundation; not an approved execution plan  
**Scope:** B-cell development redo using the updated, module-oriented
scworkbench  
**Execution authority:** none; every analysis step requires separate user
review and approval before implementation or execution

## Purpose

This document records a possible organization for the new analysis. It maps
the requested scientific steps to relevant archived work, distinguishes
scworkbench operations from analysis performed outside the framework, and
identifies likely gaps that may require new operation packs.

The proposal is intentionally revisable. Inclusion of a step, parameter,
method, bundle, or operation pack below does not select or authorize it. Before
each step begins, its scientific question, inputs, method, parameters,
validation criteria, expected outputs, and scworkbench boundary must be
approved.

## Governing constraints

- Preserve `Bdev@meta.data$pt` as the selected Bdev pseudotime. Preserve its
  values and missing entries exactly when importing and using the reference.
- Any regenerated trajectory is a validation candidate and must not overwrite
  the authoritative `pt` field unless the user explicitly changes this project
  decision.
- Treat `archived_analysis/` and its paired storage as read-only.
- Put all new repository work under `current_analysis/` and all generated
  bundles and results under the paired active storage root.
- Use scworkbench by default. Small initial or unsupported analyses may run
  outside scworkbench only after their boundary is approved and with recorded
  inputs, dependencies, parameters, seeds, and outputs.
- Preserve source bundles during integration. Place joint or transferred state
  only in separate derived bundles.
- Use capsules, or the corresponding built-in reproducible reporting
  mechanism, for terminal figures and reports. Capsules must not hide
  undeclared statistical analysis or threshold selection.

## Current state and development baseline

The active area remains a scaffold. No biological workflow or custom operation
pack has been implemented for the redo.

The corrected image received on 2026-10-29 is the immutable baseline for
`Bdev_development` and the scWorkbench update:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif
sha256:2c9aa71f69d5775166f7dd0625176a8795d00caa498500eae7c39c42c2a7695a
```

This image supersedes the immutable `20261022` candidate for development because
the earlier image did not provide the required CPU PeakVI runtime. The corrected
image is **not eligible for the paper implementation**. Its standard profile
activates five core operations and seven operation packs, providing 62 operations
in total, including `scwblatent::peakvi_fit` and
`scwblatent::peakvi_project`. The supplied Python environment imports scvi-tools
and executes on CPU. A full CPU PeakVI lifecycle and the 17-check operation-pack
tutorial have not yet been rerun for this exact image and are not claimed as
qualification evidence.

Development uses this SIF as its toolchain baseline while missing reusable
capabilities are implemented and qualified in the owning scWorkbench
repositories. Changed source is tested against this image and installed into a
new derived SIF; the corrected baseline itself is not reconstructed. The paper
workflows begin only after development qualifies that new immutable SIF and its
accepted path, checksum, baseline identity, changed source commits, pack
identities, operation catalog, and qualification record are pinned in
`current_analysis/Bdev_implementation/config/`.

## Interpretation of status

- **Not started:** no corresponding project analysis was found.
- **Archived:** the exploratory generation contains relevant work; this does
  not mean the result is accepted for the redo.
- **Archived, unpromoted:** work or reports exist but were not accepted as
  durable state.
- **Proposed:** a possible new implementation, pending approval.

## Proposed analysis steps

### 1. Derive marker gene sets from ImmGen bulk RNA data

**Current status:** not started. No project-specific ImmGen bulk marker workflow
or accepted result was found in the archived analysis.

**Proposed boundary:** begin with a small analysis outside scworkbench because
the available framework operations are cell-oriented and do not provide a
replicate-aware bulk RNA differential-expression contract.

The possible workflow would:

1. declare checksummed raw and normalized count tables and complete sample
   metadata;
2. use raw counts and a reviewed replicate-aware statistical design for marker
   discovery;
3. use normalized values only for QC and display unless an approved method
   explicitly requires them;
4. compare each selected population with approved competitor populations;
5. retain complete statistics, including effect sizes, p-values, adjusted
   p-values, specificity criteria, and exclusion reasons; and
6. export versioned long-form and GMT-like feature-set definitions for import
   into scworkbench.

An operation pack such as `scwbBulkMarkers` should be considered only after the
bulk contrast and marker-selection contract is stable and expected to be
reused. A one-off initial analysis should remain a focused, reproducible outside
script rather than force bulk samples into an inappropriate single-cell
interface.

**Approval gate:** approve the ImmGen release, population selection, sample
design, statistical method, contrasts, gene identifiers, and marker filters.

### 2. Label Bdev cells with ImmGen marker sets

**Current status:** not started.

**Proposed boundary:** scworkbench-native after the marker sets from step 1 are
approved. Candidate operations are:

```text
feature_set_import
scwbgeneprograms::feature_set_overlap_qc
scwbgeneprograms::ucell
scwbgeneprograms::cell_vector_winners
scwbgeneprograms::winner_cell_sets
```

The output should retain all scores, not only the winning label. Winner calls
should have reviewed rules for minimum score, minimum winning margin, ties,
missing features, and an explicit unassigned state. Feature-set source
checksums and identifier conversions must be recorded.

**Approval gate:** approve the feature-set version, expression layer, scoring
parameters, overlap requirements, winner/abstention rules, and validation
plots.

### 3. Regenerate the original pseudotime trajectory and compare it with `pt`

**Current status:** the archived import preserved and validated the existing
`pt` metadata; it did not regenerate the trajectory. Relevant scripts are:

- [`01_inspect_bdev_reference.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/bdev_reference/01_inspect_bdev_reference.R)
- [`02_create_bdev_reference.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/bdev_reference/02_create_bdev_reference.R)
- [`03_validate_bdev_reference.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/bdev_reference/03_validate_bdev_reference.R)

The external source recipe is
[`scMulti_analysis.Rmd`](https://github.com/Tingvall/Bdev_epigenetics/blob/main/script/scMulti/scMulti_analysis.Rmd).
The relevant recipe excludes ATAC clusters 9, 15, 18, 19, 22, and 23; starts
Slingshot at cluster 5; uses `omega = TRUE` and `omega_scale = 5`; selects
lineage 4; and reverses its sign.

**Proposed boundary:** validation only. The current native Slingshot operations
do not expose all parameters needed for an exact reproduction. Two candidates
require approval:

1. run an exact, dependency-pinned outside reproduction and adopt the resulting
   model or vector into scworkbench; or
2. implement a small `scwbTrajectoryCompat` operation pack exposing the
   required omega parameters and lineage/orientation selection.

Comparison should cover the exact defined and missing cell IDs, selected
lineage membership, direction, rank correlation, numerical differences, and
the effects of package-version changes. The authoritative `pt` values remain
unchanged regardless of the comparison result unless a later explicit decision
states otherwise.

**Approval gate:** approve the purpose of the comparison, exact dependency
versions, outside-script versus pack implementation, and numerical equivalence
criteria.

### 4. Project selected bulk ATAC reference populations and align them to Bdev

**Current status:** not started. The archived PeakVI work concerns single-cell
knockout queries, not bulk ATAC population libraries.

**Proposed boundary:** scworkbench for peak-space preparation, reference-only
PeakVI fitting, typed latent spaces, and projection where the data contract is
valid. A limited outside proof of concept may be needed because projection of a
bulk library into a model trained on single cells is not automatically a
validated biological mapping.

The provisional design is to:

1. resolve genome builds and quantify or map all libraries to an approved Bdev
   peak universe;
2. fit PeakVI to Bdev only and freeze the reference model;
3. project bulk libraries without allowing them to affect the reference fit;
4. validate the behavior using Bdev pseudobulks and depth-matched pseudo-cells;
5. map supported query points to the authoritative Bdev pseudotime; and
6. report uncertainty and out-of-reference behavior instead of forcing every
   sample onto the trajectory.

A reusable bulk-query projection operation may be added to
`scwbReferenceMap` after the proof of concept defines valid normalization,
projection, and calibration contracts.

**Approval gate:** approve the bulk ATAC datasets and populations, genome and
peak-universe policy, projection semantics, validation controls, and support
thresholds.

### 5. Identify gene sets that vary along the Bdev trajectory

**Current status:** multiple archived implementations exist. They include
native trajectory scoring combined with bespoke permutation or standardized
outside analyses, as well as classical preranked GSEA:

- [`53_trajectory_geneset_permutation_group.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/53_trajectory_geneset_permutation_group.R)
- [`62_trajectory_geneset_standardized_min20.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/62_trajectory_geneset_standardized_min20.R)
- [`75_mapped_enrichment_pseudotime_matched_gsea.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/75_mapped_enrichment_pseudotime_matched_gsea.R)
- [`76_mapped_enrichment_dynamic_change_gsea.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/76_mapped_enrichment_dynamic_change_gsea.R)

**Proposed boundary:** first evaluate the native
`scwbtrajectory::trajectory_geneset_dynamics` operation. It can compute an
empirical dynamic-score p-value and an adjusted q-value by permuting cell-level
scores relative to fixed pseudotime/bin structure.

This test is not classical preranked GSEA. If the approved scientific output
requires normalized enrichment scores, leading-edge genes, or a ranked-feature
null model, implement a separate `scwbTrajectoryGSEA` operation pack around a
pinned GSEA engine. Do not label dynamic-score permutations as GSEA.

**Approval gate:** decide whether the scientific target is trajectory dynamics,
classical GSEA, or both; approve gene-set collections, statistic, null model,
number of permutations, multiple-testing families, support requirements, and
reporting thresholds.

### 6. Reprocess Parse data and derive cytokine signatures

**Current status:** extensively archived. Ingestion and replicate-aware
signature discovery included approved outside analysis; import, QC, state
management, clustering, and UCell scoring used scworkbench.

Representative scripts are:

- outside preprocessing: [`01_readparsebio_unfiltered.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/parse_cytokines/01_readparsebio_unfiltered.R) and [`02_prefilter_min500counts.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/parse_cytokines/02_prefilter_min500counts.R);
- scworkbench import and QC: [`03_create_parse_cytokines.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/parse_cytokines/03_create_parse_cytokines.R) through [`07_build_signature_states.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/parse_cytokines/07_build_signature_states.R);
- outside replicate-aware discovery: [`08_discover_specific_signatures.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/parse_cytokines/08_discover_specific_signatures.R); and
- scoring and later competitor-aware refinements: scripts 09 through 18 in the
  same archived directory.

**Proposed boundary:** re-run from checksummed raw inputs rather than promote or
copy the recovered legacy bundle. Keep ReadParseBio ingestion and the initial
replicate-aware pseudobulk discovery outside scworkbench unless an approved
bulk-marker pack exists. Perform bundle import, QC, feature-set import, UCell
scoring, and winner calls inside the new framework.

The archived original, rest-only, and all-cytokine competitor-aware signature
collections should be treated as alternatives for review, not silently merged
or assumed to be accepted for the redo.

**Approval gate:** approve input release, QC criteria, treatment map, comparator
design, replicate handling, signature-discovery method, feature-set selection,
and which archived refinements should be reproduced.

### 7. Examine cytokine signals along the Bdev trajectory

**Current status:** archived work exists, but several sessions and reports were
unpromoted or awaiting review. Representative processing scripts are:

- [`04_cytokine_pseudotime.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/04_cytokine_pseudotime.R)
- [`05_validate_cytokine_bundle.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/05_validate_cytokine_bundle.R)
- [`06_score_competitor_aware_signatures.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/06_score_competitor_aware_signatures.R)

Additional reports are under the archived `figure_scripts/bdev_reference/`
directory, particularly scripts 05 and 08 through 14.

**Proposed boundary:** use native UCell scoring and native trajectory signal or
gene-set dynamics. The explicit trajectory operations expect a pseudotime
family, while the authoritative Bdev field is a scalar metadata vector. A small
adapter operation should form a one-lineage family over the cells with defined
`pt`, without changing or imputing the original metadata.

**Approval gate:** approve the cytokine signature collection, score layer,
pseudotime domain, smoothing/binning method, significance test, display scaling,
and behavior for missing `pt`.

### 8. Map the knockout single-cell dataset and derive enrichment intervals

**Current status:** archived PeakVI fit and projection completed, followed by
outside mapping QC and exploratory enrichment-region analyses. Much of this
state remained unpromoted.

Relevant scripts include:

- [`07_create_ko_mapping_bundle.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/07_create_ko_mapping_bundle.R)
- [`09_peakvi_full_reference_benchmark.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/09_peakvi_full_reference_benchmark.R)
- [`10_peakvi_full_mapping.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/10_peakvi_full_mapping.R)
- [`11_peakvi_mapping_qc.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/11_peakvi_mapping_qc.R)
- [`84_bdev_cutoff025_mapping_report.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/84_bdev_cutoff025_mapping_report.R)

The archived source is described as **Flt3l ligand knockout** plus IL7 knockout,
not necessarily Flt3 receptor knockout. That identity must be resolved before
durable bundle names, labels, or biological conclusions are created.

**Proposed boundary:** use native CPU PeakVI and bounded k-nearest-neighbor
operations. Add a `scwbReferenceMap` pack for the currently missing cohesive
contract: reference-only cross-validation, calibrated similarity or distance
thresholds, weighted pseudotime transfer, uncertainty, support categories, and
out-of-reference calls. A trajectory-occupancy operation may then calculate
genotype enrichment intervals and their null envelopes from the accepted
mapped cells.

**Approval gate:** confirm dataset/genotype identity; approve source-cell QC,
peak universe, PeakVI configuration, candidate reference pool, neighbor count,
distance/similarity statistic, cross-validation plan, mapping threshold, and
enrichment-interval test.

### 9. Compare motifs through the FLT3L-enriched trajectory interval

**Current status:** extensive archived exploratory figures exist, but motif
registration, chromVAR preparation, and trajectory/genotype comparisons were
largely outside scworkbench. Representative scripts are:

- [`42_bdev_selected_motif_pt_heatmaps.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/42_bdev_selected_motif_pt_heatmaps.R)
- [`47_bdev_up_down_motif_trajectory_adjusted_wt_ko.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/47_bdev_up_down_motif_trajectory_adjusted_wt_ko.R)
- [`51_single_cell_30nn_top_adjusted_motif_heatmaps.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/51_single_cell_30nn_top_adjusted_motif_heatmaps.R)
- [`86_rolling_window_motif_trajectory_heatmaps.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/integration/86_rolling_window_motif_trajectory_heatmaps.R)

**Proposed boundary:** register the exact genome provider, motif-library
version, and peak-to-motif match set, then use `scwbatac::chromvar` for native
deviation scoring. Native trajectory signal dynamics may describe motif curves.
A `scwbTrajectoryContrast` operation is likely needed for support-aware WT/KO
comparisons within accepted enrichment intervals, including p-values and
multiple-testing correction.

**Approval gate:** approve the motif database and version, genome build,
background/matching policy, enrichment interval, support mask, WT/KO comparison
model, multiple-testing family, and display-selection rule.

## Provisional bundle organization

The following seven-bundle arrangement is a starting point, not an approved
bundle plan.

### Immutable source bundles

1. `01_bdev_reference` — Bdev assays, reductions, metadata, and exact `pt`.
2. `02_immgen_bulk_rna` — bulk counts, normalized table, sample design, and marker
   provenance.
3. `03_bulk_atac_reference` — selected bulk libraries, metadata, genome build, and
   peak mapping.
4. `04_parse_cytokines` — reprocessed Parse source and QC state.
5. `05_scmulti_flt3lko_il7ko` — knockout query source; final name awaits genotype
   confirmation.

### Derived bundles

6. `06_bdev_annotation_trajectory` — Bdev plus imported ImmGen and Parse feature
   libraries, UCell scores, winner labels, trajectory comparison, and approved
   trajectory gene-set results.
7. `07_bdev_atac_mapping` — Bdev reference model, separate bulk and single-cell
   query projections, mapping QC, transferred pseudotime, enrichment intervals,
   chromVAR state, and motif contrasts.

The source bundles must remain independent and unmodified. ImmGen and Parse
feature sets should be imported into a Bdev-derived bundle rather than merging
their observations with Bdev. PeakVI should be fit on Bdev only and then frozen
before query projection. Bulk and single-cell query union spaces should remain
separate even if they share a reference model.

If bulk and knockout mappings need separate review and release cycles, split
the seventh bundle into two derivatives, producing eight bundles. Whether the
framework can transfer one frozen model between those bundles without
duplication must be verified before choosing that structure.

## Provisional scWorkbench development plan

Reusable paper methods should be added to scWorkbench core or official
operation packs during `Bdev_development`, rather than retained as
implementation-local code. scWorkbench and PyBench are core infrastructure and
remain unchanged by default. New scientific operations go into the existing
official pack that fits their domain, or a new official pack when none fits.
Core changes require a pressing, documented missing primitive or infrastructure
capability and the user's explicit approval for the specific proposed change.
Before any scWorkbench or PyBench edit, development must state the evidence,
explain why an existing or new pack cannot solve the gap, motivate the minimal
core change, describe its qualification impact, and ask for approval. The
current likely ownership is:

| Capability | Preferred owner | Development note |
|---|---|---|
| Exact Slingshot omega behavior and scalar pseudotime-to-lineage adaptation | `scwbTrajectory` | Add separate qualified operations where adoption and refitting have different semantics. |
| Cross-validated query mapping, support thresholds, pseudotime transfer, uncertainty, and out-of-reference calls | Existing representation/latent packs where their contracts fit; otherwise a new official `scwbReferenceMap` pack | Keep reference fitting, projection, calibration, and transfer provenance explicit. |
| Classical preranked trajectory GSEA | `scwbTrajectory` | Keep it distinct from the existing cell-score permutation dynamics operation. |
| Trajectory occupancy and group-aware signal or motif contrasts | `scwbTrajectory` | Scientifically distinct alternatives, including alternative cytokine-signal analyses, should remain separate selectable operations rather than hidden modes. |
| Replicate-aware bulk marker discovery | `scwbStatistics` or a new official bulk-analysis pack | Add only after the ImmGen/Parse design and reusable sample-level contract stabilize. |

Core changes are appropriate only when a generic primitive, reader, writer,
axis, registry, provenance, or publication capability is missing, and even then
must not begin without specific user approval. Each approved operation should
have typed inputs and outputs, deterministic seeds where
applicable, complete parameter capture, dependency versions, machine-readable
QC and statistical tables, a small fixture, contract tests, and at least one
negative test. Source changes live on branches in the owning repositories and
enter implementation only through a qualified immutable SIF.

## Reporting and capsules

Figures should be generated only after the corresponding analytical state is
reviewed and promoted. A capsule should declare:

- exact bundle and artifact inputs;
- configuration and selection tables;
- figure/report source checksum;
- package and runtime identity; and
- output files and checksums.

Marker discovery, model fitting, mapping-threshold selection, enrichment-region
selection, and statistical testing should remain explicit workflow operations,
not unrecorded calculations inside a figure capsule.

Possible report families are Bdev annotation/trajectory validation, cytokine
trajectory signals, reference and knockout mapping, and motif contrasts. Their
number and content require approval and do not determine the bundle count.

## Step approval record

Before starting any step, add a dated decision record under
`current_analysis/Bdev_implementation/docs/` containing at least:

```text
step:
status: approved | revised | deferred | rejected
decision_date:
approved_by:
scientific_question:
input_manifests:
runtime_pin:
scworkbench_operations:
approved_outside_analysis:
parameters_and_thresholds:
validation_criteria:
expected_outputs:
bundle_targets:
report_or_capsule_targets:
notes:
```

Approval of one step does not imply approval of later steps. Parameters learned
from an exploratory or pilot run require explicit promotion before they become
defaults for a production run.

## Decisions still required

At minimum, resolve the following before implementation:

1. approve the immutable scworkbench runtime;
2. confirm the selected ImmGen and bulk ATAC datasets and sample identities;
3. confirm whether the knockout is Flt3, Flt3l, or another recorded condition;
4. select which archived Parse signature design, if any, should be reproduced;
5. define whether “trajectory GSEA” means dynamic feature-set testing,
   classical preranked GSEA, or both;
6. select mapping-support and enrichment-region methods;
7. select genome and motif-library versions; and
8. approve the final bundle graph and operation-pack boundaries.

Until those decisions are recorded, this document serves only as a foundation
for discussion and stepwise approval.
