# Proposed Parse-to-Bdev cytokine integration

Status: approved by the user, implementation started 2026-09-23. This document
records the design; execution status and observed coverage are recorded in
`config/bdev_integrated.yml` and the run reports described in `README.md`.

## Scientific question and first phase

Locate transcriptional programs resembling the accepted Parse cytokine
responses along Bdev's existing `pt` coordinate. Transfer the gene sets into a
separate `bdev_integrated` scWorkbench bundle and score Bdev RNA. This first
phase needs the accepted `bdev_reference` and `parse_cytokines` parents.
The planned `scmulti_flt3ko_il7ko` parent enters a later integration phase after
its acceptance. Update the integration configuration and parent gates to this
staged dependency structure when implementing the workflow.

The first phase uses Bdev cells and Parse-derived feature sets. Joint cell
embedding, Parse cell projection, and a newly fitted trajectory are unnecessary
for this question. Signature similarity is evidence of a response-like
transcriptional program; it does not establish cytokine exposure or causality.

## Accepted inputs

- Bdev: 26,723 cells; `RNA_counts` has 32,285 gene-symbol features. Preserve all
  metadata and existing embedding coordinates, including `umap.atac`.
- Pseudotime: authoritative `Bdev@meta.data$pt`, represented by reference
  artifact `cell_vector_20260914_162517_7q4y8xbb`. There are 12,642 defined
  values and 14,081 missing entries. Preserve values, direction, cell IDs, and
  the missing-value mask exactly. No fitted Slingshot model was found.
- Parse: the accepted `parse_cytokine_signatures` namespace contains 23 sets
  across 14 cytokines: nine up-regulated and fourteen down-regulated sets.
  Nine cytokines have both directions; five have only an accepted down set.
  Use the final accepted membership, not all 31 discovery candidates. Discovery
  and acceptance used control/single-stimulation cells and did not use Bdev.
- Source paths, accepted fingerprints, signature exports, and runtime are
  resolved from `config/bundles.yml` and the parent configurations. Use the
  accepted scWorkbench 0.1.0.9004 runtime declared for Parse and record its
  build provenance in the derived analysis.

## Analysis sequence

1. **Preflight and create the derived context.** Verify accepted parent
   identities, retained validation evidence, source/sample IDs, RNA feature
   identifiers, and exact pseudotime alignment. Create `bdev_integrated` through
   the documented bundle/ingest interfaces, with Bdev as the reference parent.
   Keep the complete source signature membership and its checksums alongside
   the target feature mapping. All parent bundles remain read-only.

2. **Check signature coverage before interpreting scores.** Import the accepted
   feature sets, retaining cytokine, direction, Parse acceptance evidence, and
   source namespace. Run `feature_set_overlap_qc` against Bdev `RNA_counts`.
   Proposed eligibility is at least five mapped genes and at least 80% of the
   original set represented in the matrix. Record missing genes and all failed
   sets; do not select replacement genes based on Bdev trends. Inspect identical
   mapped sets and overlapping membership as specificity diagnostics. The
   engine's overlap QC measures membership in the feature universe, not actual
   nonzero RNA expression; inspect expression support separately.

3. **Score with the same method used in Parse.** Use registered `ucell` on Bdev
   `RNA_counts`, `max_rank=1500`, without graph smoothing. Score all retained
   Bdev cells for spatial context; restrict trajectory summaries to the 12,642
   cells with defined `pt`. Report scoring coverage and missing features.
   Keep the 23 directional scores separate. Higher up-set expression together
   with lower down-set expression supports a concordant response-like pattern.
   A down-only signature supplies one-sided evidence. Do not manufacture an
   up signature, negate a down score into a calibrated activity estimate, or
   select a winning cytokine from raw scores across unrelated sets.

4. **Summarize along the existing coordinate.** Use
   `trajectory_signal_dynamics` with the UCell score family and explicit
   pseudotime, `pseudotime_direction=as_is`. Its typed pseudotime input is a
   `cell_vector_family`; represent the existing scalar vector under one
   descriptive key through the public primitive constructor and verify exact
   cell/value equality. The key identifies the supplied coordinate, not an
   inferred lineage. Retain the source artifact identity and full metadata
   missingness. Confirm this adapter against the live operation interface
   before execution.

   Start with 30 equal-width bins in original `pt` units, median scores, a
   minimum of 30 cells per retained pooled bin, and no smoothing or bootstrap.
   Report bin counts and unsupported intervals explicitly. Repeat at 20 and
   50 bins to identify unstable features. Record early/late changes, amplitude,
   and peak/trough locations descriptively. Optional loess overlays must retain
   the original binned values and must not imply support across empty regions.
   Do not label pseudotime as elapsed time or assign biological stage boundaries
   until the source annotations support them.

5. **Check sample consistency and competing explanations.** Bdev contains four
   `orig.ident` values (1–4). Inspect each sample's pseudotime coverage and repeat
   summaries by sample, initially requiring 20 cells per bin. Native binning
   uses each selected population's range: compare actual pseudotime positions
   and shared support, never equate bin numbers across sample-specific runs.
   Confirm whether these IDs represent independent biological replicates
   before making replicate-level claims. Inspect `nCount_RNA`, `nFeature_RNA`,
   `S.Score`, `G2M.Score`, and `Phase` alongside score trends; stratify where
   needed and label unresolved confounding. Inspect correlated signatures and
   member-gene expression to distinguish shared programs from cytokine-specific
   evidence. Use `cell_vector_family_qc` for score relationships and, as a
   secondary expression view, `trajectory_geneset_dynamics` with explicit `pt`
   on a recorded RNA log-normalized matrix for member-gene curves. Its gene-set
   expression summaries are distinct from the primary UCell scores.

6. **Retain and validate.** Save per-cell scores, feature mapping/QC, pooled and
   sample-specific curves, descriptive dynamics, source identities, parameters,
   runtime identity, and acceptance decisions in the separate integration
   context. Use durable sessions and review artifact sizes before promotion.
   Validate cell alignment, exact reference `pt` and missingness, unchanged
   coordinates, finite score coverage, retained memberships, parent
   fingerprints, and recipient readability before accepting the derived bundle.

## Figures and interpretation

- A 23-signature pseudotime heatmap, with directions labeled and within-signature
  scaling used only for display. Preserve raw-score curves beside it; heatmap
  colors do not compare absolute cytokine activity.
- Fourteen cytokine panels pairing available up/down curves, sample-specific
  curves, and cell-coverage information.
- Selected score overlays on the existing `umap.atac` embedding, with a separate
  display of defined versus missing pseudotime.
- Member-gene panels and a coverage/robustness table documenting unsupported
  sets, sample disagreement, and sensitivity to binning.

The first report is descriptive. scWorkbench's trajectory bootstrap resamples
cells within bins; those bands would not estimate between-animal uncertainty.
Do not use them as biological-replicate confidence intervals. Sample-aware
formal inference can be designed after confirming the experimental units. If
the required model is unsupported by scWorkbench, explain that specific gap
and obtain approval for outside analysis before running it.

## Capability evidence

Reviewed the documentation and implementation extracted from `/opt/scworkbench`
inside the supplied SIF during the accepted Parse recovery. The extraction is
under the configured runtime's sibling `source/` directory; its origin is
recorded in `build_provenance.json`. Relevant references are
`docs/op_catalog_reference.md`, `docs/bdev_external_user_workflows.md`,
`man/wb_op_trajectory_signal_dynamics.Rd`,
`man/wb_prim_cell_vector_family.Rd`, and `R/op_trajectory_smooth.R`.
This is a project-specific design, not a replacement for engine documentation.

The model-dependent `trajectory_signature_smooth` is unsuitable for the
retained Bdev state. The explicit-pseudotime mode of
`trajectory_signal_dynamics` supports the proposed primary analysis without
fitting or inventing a trajectory model. Feature coverage, sample-by-pseudotime
coverage, and runtime execution of the proposed workflow remain to be checked.
# Current scope amendment

The user subsequently selected pooled analysis only, with continuous lines and
a measure of spread along each trajectory. The current implementation uses
pooled medians with native scWorkbench 25th–75th percentile bands (IQR), exports
within-bin SD, and connects retained bins while marking low cell support.
Sample and cell-cycle stratification below is part of the original plan and is
superseded for the current analysis. No further promotion is authorized before
review of the pooled report.
