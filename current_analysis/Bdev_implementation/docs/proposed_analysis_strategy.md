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

- Use authenticated fitted Slingshot lineage `L4` as the default Bdev
  pseudotime trajectory, in native Slingshot orientation. Use `-L4` only for an
  explicit comparison with the legacy scalar `pt` orientation.
- Preserve `Bdev@meta.data$pt` values and missing entries exactly. Retain the
  scalar field as validation evidence, a fallback when fitted L4 is
  unavailable, or an explicitly requested legacy source; never overwrite it.
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

**Current status:** selected as the proposed later-implementation strategy on
2026-11-02 after bounded development exploration; not active in the reserved
implementation stage. Use the checksummed `GSE109125` GENCODE vM25 count
matrix restricted to the 126 approved libraries representing 64 steady-state
ImmGen populations.

The proposed discovery rule is intentionally simple and population-balanced:

1. transform counts to `log2(CPM + 0.5)`;
2. take the median across libraries within each ImmGen population;
3. for every gene and target population, subtract the equally weighted mean of
   the other 63 population medians from the target median;
4. require target expression of at least 1 on the transformed scale and a
   positive contrast; and
5. retain the top 20 genes for each of all 64 populations.

Gene reuse across populations is allowed. Neither Bdev expression nor the
user-supplied comparison gene sets may be used to select or tune genes. This
keeps every ImmGen population eligible and prevents dataset-informed removal of
apparently implausible labels. The earlier broad-lineage gate, percentile
calibration, shared-gene prevalence filter, and replicate-aware edgeR proposal
are not part of this selected strategy.

Marker discovery remains reproducible outside scworkbench because the available
framework operations are cell-oriented. The versioned GMT, ranked candidates,
membership table, configuration, run manifest, and checksums form the import
boundary; bulk libraries must not be represented as single cells. The selected
development record is
[`15_default_immgen_geneset_collection.yaml`](../../Bdev_development/config/01_bulk_marker_boundary/15_default_immgen_geneset_collection.yaml),
and its method summary is
[`24_simple_mean_contrast_development_default.md`](../../Bdev_development/docs/01_bulk_marker_boundary/24_simple_mean_contrast_development_default.md).

**Activation gate:** before implementation execution, the selected artifacts
must be pinned in an implementation input manifest and rerun or verified with
the accepted immutable SIF. This proposal does not itself approve labels,
activate a workflow, or unlock `Bdev_implementation`.

### 2. Label Bdev cells with ImmGen marker sets

**Current status:** exact development strategy selected as the proposed later
implementation default; implementation not started.

**Proposed boundary:** scworkbench-native after the marker sets from step 1 are
approved. Candidate operations are:

```text
feature_set_import
scwbgeneprograms::feature_set_overlap_qc
scwbgeneprograms::ucell
scwbgeneprograms::cell_vector_winners
scwbgeneprograms::winner_cell_sets
```

After feature import and overlap QC, score `RNA_counts` with native UCell using
`maxRank = 1500`. Reduce the 64 raw scores by uncalibrated argmax, with no broad
or lineage gate and no percentile calibration. Require a winner-minus-runner-up
raw-score margin of at least `0.01`; ties and smaller margins are `unassigned`.
The minimum raw score is `0.0`.

Retain all 64 scores, the winner, runner-up, winner and runner-up scores, raw
margin, assignment or abstention reason, missing-feature information, and
source gene-set checksums. Preserve the existing Bdev `pt` field and its
missing entries exactly. Use the ATAC UMAP for label visualizations, but never
use its geometry to select genes or labels.

The 30-round within-population ImmGen bootstrap produced 71.2% overall
cell-round agreement with the baseline and a 76.7% median per-cell agreement.
The B-cell-heavy manifold was substantially more stable than adjacent
stem/progenitor labels. Therefore carry numeric bootstrap agreement as an
uncertainty field and avoid hard interpretation of unstable sublabels; do not
add a new confidence threshold unless it is separately justified. Eight
singleton ImmGen populations remain structurally unassessable for
target-replicate stability. The evidence is recorded in
[`25_simple_mean_contrast_bootstrap_stability.md`](../../Bdev_development/docs/01_bulk_marker_boundary/25_simple_mean_contrast_bootstrap_stability.md).

**Activation gate:** pin the qualified runtime, exact gene-set version and
checksums, expression layer, overlap QC, reducer parameters, bootstrap evidence,
and validation outputs before execution. No scWorkbench core or operation-pack
change is currently required because the existing import, UCell, and winner
reducer operations support this strategy.

### 3. Regenerate the original pseudotime trajectory and compare it with `pt`

**Current status:** recreated successfully as a bounded development validation
with scwbTrajectory 0.1.3; not activated in the reserved implementation stage.
The archived import preserved and validated the existing `pt` metadata.
Relevant archived scripts are:

- [`01_inspect_bdev_reference.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/bdev_reference/01_inspect_bdev_reference.R)
- [`02_create_bdev_reference.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/bdev_reference/02_create_bdev_reference.R)
- [`03_validate_bdev_reference.R`](../../../archived_analysis/legacy_scworkbench_exploratory/processing/bdev_reference/03_validate_bdev_reference.R)

The external source recipe is
[`scMulti_analysis.Rmd`](https://github.com/Tingvall/Bdev_epigenetics/blob/main/script/scMulti/scMulti_analysis.Rmd).
The relevant recipe excludes ATAC clusters 9, 15, 18, 19, 22, and 23; starts
Slingshot at cluster 5; uses `omega = TRUE` and `omega_scale = 5`; selects
lineage 4; and reverses its sign.

**Development result:** the audit completed on 2026-11-02 and the required
updates were implemented in core commit `da9477d` and scwbTrajectory 0.1.3
commit `080a192`. The full 26,723-cell reference fit retained 23,710 cells and
produced eight lineage families. Negated L4 shared 12,510 finite cells with
preserved `pt`, with Spearman correlation 0.999817 and Pearson correlation
0.999942. Membership and exact values differed slightly (941 recreated-only
cells, 132 preserved-only cells, mean absolute difference 0.5732), so this is a
near-exact ordering recreation rather than byte-identical recovery. Preserved
`pt` values and missingness were unchanged. The user accepted this as
sufficiently close and selected fitted L4 as the default pseudotime trajectory
for subsequent Bdev work. The audit, result, and decision are recorded
in
[`01_slingshot_semantics_gap_analysis.md`](../../Bdev_development/docs/02_slingshot_trajectory_semantics/01_slingshot_semantics_gap_analysis.md),
[`02_bdev_trajectory_recreation.md`](../../Bdev_development/docs/02_slingshot_trajectory_semantics/02_bdev_trajectory_recreation.md),
and
[`03_l4_default_trajectory.md`](../../Bdev_development/docs/02_slingshot_trajectory_semantics/03_l4_default_trajectory.md).

**Proposed boundary:** the reusable updates now exist in the official
`scwbTrajectory` pack, supported by the approved generic mixed-scalar core
change. Fitted lineage-family pseudotime is the default trajectory source, with
lineage keys selected explicitly. The preserved scalar `pt` is used only when
lineage output is unavailable, when specifically requested, or as the
validation reference for the published trajectory. PyBench is unchanged.

Comparison should cover the exact defined and missing cell IDs, selected
lineage membership, direction, rank correlation, numerical differences, and
the effects of package-version changes. The authoritative `pt` values remain
unchanged regardless of the comparison result unless a later explicit decision
states otherwise.

**Activation gate:** include the pinned core and pack commits in a new immutable
SIF, qualify and accept that SIF, then explicitly activate this workflow in the
implementation stage. The development session and figures are validation
evidence, not paper-analysis outputs.

### 4. Project selected bulk ATAC reference populations and align them to Bdev

**Current status:** proposed later-implementation method selected on 2026-11-22
after successful development bulk validation; implementation and paper execution
remain locked pending the qualified SIF.

Use the exact frozen 20-dimensional depth-adjusted HVP20k PeakVI reference
qualified in run `gated_continuous_transfer_20261122T192735Z`. The selected
model is `model_ref_20261122_194532_hka11fvs`, and the selected calibration is
`model_ref_20261122_195229_ewshintf`. The same model and calibration are used
for both single-cell multiome queries and bulk ATAC libraries; do not fit a
separate bulk-specific reference model by default.

Bulk libraries must be quantified on the exact 20,000-feature order and
projected without refitting. Apply the same primary-30 reference limits:
supported through `2.02620492264996`, borderline through
`2.56846898217638`, and out of reference above that value. Publish native L4
only when the local domain vote is `in_domain`, confidence is at least `0.51`,
and the estimate is finite. Retain the raw estimate and rejection reason for
every library. Fixed-UMAP figures use neighbor-based placement on the unchanged
Bdev `umap.atac`; neither UMAP nor Slingshot is refit.

The proposed implementation method is `gated_20d_shared_reference_mapping`.
Every valid query receives a distance-support call and fixed-UMAP placement.
All raw native-L4 estimates and rejection reasons are retained, but an accepted
native-L4 value is published only when the complete support, local-domain,
confidence, and finite-estimate gate passes. In the fixed UMAP, queries that do
not pass the complete L4 gate are shown in dark grey (`#343434`) and must not be
interpreted as accepted trajectory positions. This reporting contract applies
to both single-cell multiome and bulk ATAC queries using the same exact frozen
model and calibration.

The complete selection is pinned in
[`32_qualified_20d_reference_mapping_default.yaml`](../../Bdev_development/config/04_reference_mapping/32_qualified_20d_reference_mapping_default.yaml)
and mirrored as an inactive implementation default in
[`reference_mapping.yaml`](../config/reference_mapping.yaml).

Development validation run `immgen_bulk_gated_20d_20261123T001821Z` projected
all 28 frozen ImmGen ATAC libraries through that exact model. All were within
the calibrated distance-support limit; four passed the complete native-L4
domain gate and 24 retained raw estimates but had accepted values withheld.
This validates execution of the shared model and gate, not a paper result.

**Activation gate:** pin and accept the new immutable SIF, verify the exact
model and calibration payloads, and approve the bulk input manifest before
implementation execution. Development figures remain qualification evidence,
not paper outputs.

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

**Selection recorded 2026-11-02:** use the 38 passing signatures from the
all-cytokine top-100 competitor refinement. A development-stage standalone
reconstruction now reproduces its candidate membership, final membership, and
selected-set identity exactly; see
[`01_standalone_parse38_reconstruction.md`](../../Bdev_development/docs/06_parse_cytokine_signatures/01_standalone_parse38_reconstruction.md).
This records the collection for later implementation but does not activate the
implementation stage or authorize reuse of the archived bundle.

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
gene-set dynamics. Use authenticated fitted L4 pseudotime in native Slingshot
orientation by default. If fitted L4 is unavailable, or scalar `pt` is
specifically requested, the fallback adapter should form a one-lineage family
over cells with defined `pt` without changing or imputing the original
metadata.

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

**Identity decision (2026-11-17):** the user confirmed **Flt3l ligand
knockout (FLT3LKO)** plus **Il7 knockout (IL7KO)**, with corresponding WT
libraries. This is not an *Flt3* receptor knockout. The permanent dataset ID
`05_scmulti_flt3lko_il7ko` is therefore retained.

**Development evidence (2026-11-17):** the bounded current-axis fixture passed
native CPU PeakVI projection, exact k-nearest-neighbor mapping, Bdev-only
support calibration, and native-L4 transfer through `scwbReferenceMap`. This
qualifies the mapping contract, not the paper analysis. A trajectory-occupancy
operation is still required for fixed-bin group summaries and declared
reference-relative envelopes.

The available source has one library for each target/genotype combination and
no independent mouse/donor field. Cells cannot serve as biological replicates.
Unless independent sample provenance is supplied, the final workflow may show
explicitly descriptive/reference-relative occupancy only; it must not report
cell-level genotype p-values or call WT/KO differences significant.

**Approval gate:** confirm any independent biological-replicate provenance;
approve source-cell QC, peak universe, PeakVI configuration, candidate
reference pool, neighbor count, distance/similarity statistic,
cross-validation plan, mapping threshold, and whether occupancy remains
descriptive or uses a replicate-aware inferential model.

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
5. `05_scmulti_flt3lko_il7ko` — confirmed FLT3LKO/IL7KO knockout query source.

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
| Replicate-aware bulk marker discovery | Approved outside preprocessing | Keep bulk data and statistical analysis outside scworkbench; import only versioned generated gene sets for use on single-cell datasets. No pack or core change is required. |

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
3. preserve the confirmed Flt3l-ligand/Il7 knockout identities in every input
   and result manifest;
4. implement the selected 38-set Parse design from its frozen standalone
   development reconstruction after the qualified runtime is accepted;
5. define whether “trajectory GSEA” means dynamic feature-set testing,
   classical preranked GSEA, or both;
6. select mapping-support and enrichment-region methods;
7. select genome and motif-library versions; and
8. approve the final bundle graph and operation-pack boundaries.

Until those decisions are recorded, this document serves only as a foundation
for discussion and stepwise approval.
