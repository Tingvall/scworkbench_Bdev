# scWorkbench B-cell development use case

This repository contains the B-cell development workflows for the scWorkbench
paper. It has three independently processed datasets and one derived integration
bundle:

1. `bdev_reference`: the reference B-cell development Seurat/Signac state and
   its accepted `pt` pseudotime.
2. `parse_cytokines`: Parse RNA cytokine-stimulation data used to derive
   cytokine response gene sets.
3. `scmulti_flt3ko_il7ko`: IL7KO/FLT3KO RNA and ATAC multiome data used as a
   query and validation dataset.
4. `bdev_integrated`: a derived integration bundle with `bdev_reference` as
   the reference. Its current cytokine phase uses the accepted Bdev and Parse
   parents; the KO parent is reserved for a later phase after acceptance.

The three parent bundles remain independent and read-only during integration.
All joint state is written to `bdev_integrated`.

## Repository organization

```text
config/                       Bundle registry and project-specific settings
processing/
  bdev_reference/             Reference import and validation
  parse_cytokines/            Parse inspection, processing, and validation
  scmulti_flt3ko_il7ko/       KO multiome inspection, processing, and validation
  integration/                Bdev-centered derived integration
figure_scripts/<bundle_id>/   Reproducible reports and final figures
figure_outputs/<bundle_id>/   Generated, Git-ignored outputs
src/                          Reusable validation and configuration helpers
```

Numeric script prefixes define execution order within each workflow directory.
`config/bundles.yml` is the central registry for bundle IDs, roles, paths,
accepted fingerprints, and parent relationships. The per-bundle configuration
files record source boundaries and workflow state without duplicating the
scWorkbench documentation supplied inside the project image.

Large datasets, bundles, run workspaces, and validation evidence remain under
`/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/`. New run-specific
evidence should use
`intermediate_objects/<bundle_id>/reports/<run_id>/`.

Mirrored storage is dataset-centered directly under
`intermediate_objects/<bundle_id>/`. Each dataset directory contains
`prepared/`, `reports/`, and `work/`, plus `bundle/` after its scWorkbench
bundle is created. Initialize missing dataset directories with
`processing/00_initialize_storage_layout.sh`.

## Bdev reference bundle

The reference source is the `Bdev` Seurat/Signac object in
`/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/datasets/scmulti_bdev/Bdev.Rdata`.
Generated data remain in mirrored storage under
`/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/`:

- `bdev_reference/bundle/`: the scWorkbench bundle.
- `bdev_reference/prepared/inspection/complete/`: the read-only source inventory.
- `bdev_reference/prepared/import/`: source identity, import mappings, and validation reports.

`processing/bdev_reference/01_inspect_bdev_reference.R` inspects the original
object. `processing/bdev_reference/02_create_bdev_reference.R` creates a new
bundle using scWorkbench
inside the supplied image; its arguments are the source file, bundle directory,
and report directory. Mount the source read-only and the bundle's parent
directory writable. The script checks the source checksum and refuses to
overwrite an existing bundle. General scWorkbench instructions remain inside
the image; see `AGENTS.md` for its location.

The import retains the nine available assay/layer matrices, all 68 metadata
columns, and `pca`, `lsi`, `umap.atac`, `umap.rna`, and `umap.wnn` coordinates.
Original assay/layer names, reduction names, feature and cell identifiers, and
the source checksum are recorded. Metadata factor labels are imported as text.
The existing `pt` metadata defines pseudotime: its 12,642 defined values form
the trajectory vector, while all 14,081 missing values remain missing. The
metadata retains the exact original numeric values. This image converts
cell-vector values to text internally, causing rounding of at most
`5e-14` in the derived pseudotime artifact; its ordering and missingness are
checked against the original. Use the metadata `pt` for exact source values.

`src/reference_validation.R` checks all metadata and reduction coordinates,
matrix dimensions and identifier ordering, and 128 distributed features per
matrix across every cell. The creation script also verifies the complete `pt`
vector and runs scWorkbench recipient validation. `completion.json` in the
import report directory records successful completion.
`processing/bdev_reference/03_validate_bdev_reference.R` accepts the same three
arguments to validate and finalize an existing import without recreating its
data or vector.

The reproducible dataset report source is
`figure_scripts/bdev_reference/01_bdev_reference_report.Rmd`. Generated HTML
and PNG outputs are in `figure_outputs/bdev_reference/`.

The IL-7/Flt3L top25-core sanity check maps log1p(CP10k) RNA expression onto
the accepted ATAC-derived UMAP. Its separate up and down figures, expression
detection table, and provenance are under
`intermediate_objects/bdev_reference/reports/il7_flt3_top25_core_atac_umap_20260925/`.

## Parse cytokine bundle

The accepted `parse_cytokines` bundle is built from the checksummed
pre-filtered Seurat RDS at
`intermediate_objects/parse_cytokines/prepared/parse_cytokines_prefiltered_min500counts_seurat.rds`.
It contains the raw `RNA_counts` matrix for 251,599 cells and 78,298 features.
Ingestion retained only the required `sample_id`, `replicate`, and `condition`
metadata for all 96 samples; the reviewed treatment step later added the
derived `treatment_category` and `cytokine_count` fields. At ingestion it
contained no normalized matrix or dimensional reduction. The subsequent
signature-state workflow added `RNA_log1p_cp10k` and `parse.signature_pca`.
The imported-cell universe remains the unchanged `set:parse_cytokines` set.

`processing/parse_cytokines/03_create_parse_cytokines.R` performs atomic
scWorkbench ingestion, compares source and imported identifiers, metadata, and
distributed count values, finalizes the bundle, and runs recipient validation.
The original ingestion bundle and evidence remain under
`intermediate_objects/parse_cytokines/bundle/` and
`intermediate_objects/parse_cytokines/reports/bundle_import_20260919T233232Z/`.
The current accepted copy is
`intermediate_objects/parse_cytokines/bundle_recovered_20260923T200125Z/`,
as recorded in `config/bundles.yml`.

`processing/parse_cytokines/04_calculate_rna_qc.R` uses the registered
scWorkbench `rna_qc_metrics` operation to calculate total counts, detected
features, and explicitly matched mitochondrial and ribosomal-prefix counts and
percentages for every imported cell. The promoted metrics, exact feature
membership, and summary passed recipient validation. The distribution report
is under `figure_outputs/parse_cytokines/rna_qc_20260920T210025Z/`.

After distribution review, the user selected `total_counts >= 1000`,
`detected_features >= 500`, and `mitochondrial_percent <= 10` as inclusive
acceptance criteria. `processing/parse_cytokines/05_apply_rna_qc_filter.R`
wrote the combined exclusion preview and sample/condition effects before
deriving any accepted set. It then promoted `set:parse_qc_excluded` with 57,384
cells and `set:parse_qc_accepted` with 194,215 cells. These are a validated,
disjoint, exhaustive partition of the 251,599 imported cells; all imported
cells remain retained in the unchanged `set:parse_cytokines` universe and raw
matrix. All 96 source samples retain accepted cells. The run evidence is under
`intermediate_objects/parse_cytokines/reports/rna_qc_filter_20260921T035031Z/`.

The reviewed treatment map classifies all 96 samples and 251,599 imported
cells as `control` (zero cytokines), `single` (one cytokine), or `combined`
(more than one cytokine). The classification is retained as the imported
`treatment_category` and `cytokine_count` metadata fields and the promoted
`parse_treatment_category` categorical vector. It includes 3 control, 78
single-cytokine, and 15 combined samples.

`processing/parse_cytokines/06_classify_treatments_and_select.R` wrote the
complete classification and treatment-scope exclusion preview before deriving
`set:parse_qc_control_single`. That accepted downstream universe contains
154,209 QC-accepted control/single cells across 81 samples and 27 conditions.
The 40,006 QC-accepted combined-stimulation cells remain in the bundle but are
outside that downstream set. The imported and both QC-derived cell sets remain
unchanged. Run evidence is under
`intermediate_objects/parse_cytokines/reports/treatment_selection_20260921T145825Z/`.

The individual-cytokine signature workflow used only that control/single
population and did not use Bdev. Replicate-aware response and specificity
testing produced 31 discovery candidates; sample-level transfer of raw-count
scWorkbench UCell scores retained 23 final signatures (9 induced and 14
repressed) across 14 cytokines. The report, heatmap, distributions, exact
sample scores, state-composition diagnostics, and final TSV/GMT gene sets are
under
`intermediate_objects/parse_cytokines/reports/signature_report_20260921T221910Z/`.

The same directory now also contains a 31-page treatment-by-sample UCell
violin collection covering every scored discovery gene set (23 accepted plus
8 UCell-transfer failures) and a complete 52-direction exclusion audit. The
21 remaining directions did not yield the prespecified minimum five eligible
genes, so they have no gene set or UCell score to plot. See
`SIGNATURE_EXCLUSIONS_AND_VIOLINS.md` in that report directory.

The Parse bundle is finalized. Recovery repaired scWorkbench's handling of
the existing 1.398-GB provenance record in a separate runtime built from the
supplied image's source. All 23 accepted gene sets are registered read-only in
`parse_cytokine_signatures`, and full recipient validation verified the existing
31-signature UCell family. Metadata, matrices, PCA, and cell-set memberships
were preserved; every original bundle file retained its SHA-256. Five missing
commit files were restored byte-for-byte from the authenticated saved session.

The accepted copy and its required scWorkbench `0.1.0.9004` library are pinned
in `config/bundles.yml`. Processing and figure runners select that runtime
automatically. The SIF itself is unchanged and still contains the original
limitation. Recovery and validation evidence is under
`intermediate_objects/parse_cytokines/reports/signature_recovery_20260923T200125Z/finalization/`.
Earlier signature reports retain their historical blocked status; this recovery
report records its resolution.

The separate `parse_cytokines_rest_only` accepted bundle rebuilds cytokine
signatures without treatment-versus-control selection. It retains the
target-versus-other-single-stimuli contrast, performs fresh membership and size
selection, and uses native UCell for confirmation. The promoted competitor-aware
collection contains 33 passing full signatures and 33 corresponding stable core
sets in separate read-only namespaces. Its configuration is
`config/parse_cytokines_rest_only.yml`; the source workflow is documented in
`processing/parse_cytokines/README.md`. Reports and GMT/TSV exports are under
`intermediate_objects/parse_cytokines_rest_only/reports/rest_only_20260924/`.
These directions describe differences from other individual stimuli, not
changes relative to unstimulated control.

A separate, unpromoted version tests Bdev-compatible 25--50 gene `full` sets
and fixed `top25` sets while preserving the accepted collection above. It
constructed both variants for all 52 cytokine directions; 24 full variants
and 27 top25 variants passed full-data competitor-aware discrimination. Exact
memberships, three-fold cores, holdout results, and the complete table are in
`intermediate_objects/parse_cytokines_rest_only/reports/rest_only_20260924/competitor_aware_bdev25_dual_v1/`.
The implementation and method are in `processing/parse_cytokines/` under
`16_competitor_aware_bdev25_dual.R` and
`COMPETITOR_AWARE_BDEV25_DUAL_WORKFLOW.md`.

A further unpromoted pilot refines Flt3L and IL-7 up/down candidates beginning
from the intersection of three fold-specific top-100 lists. It caches exact
sample-level UCell rank contributions for the 651-gene union, removes genes
only when deletion improves the full-data target-to-strongest-competitor gap,
and verifies initial and final sets with native scWorkbench UCell. All four
final sets pass the internal mean route; the report, exact memberships,
metrics, and Bdev `SCT_data` ATAC-UMAP figures are in
`intermediate_objects/parse_cytokines_rest_only/reports/rest_only_20260924/flt3_il7_competitor_refinement_top100_v1/`.

The fixed pilot method is implemented for all 26 cytokines and both directions
as a two-stage, unpromoted review workflow. First run
`processing/parse_cytokines/run_18_all_cytokine_competitor_refinement.sh` to
produce the complete 52-direction Parse inventory. After reviewing that result,
pass its directory to
`figure_scripts/bdev_reference/run_08_all_cytokine_refined_ucell_report.sh` to
score only passing signatures in Bdev and generate UCell UMAP and 5%
sliding-window pseudotime reports. The full workflow completed as an unpromoted
review on 2026-09-25.
Of 52 attempted directions, 43 constructed a stable top-100 intersection and
38 passed the final Parse mean and/or responder tests. Native Bdev UCell scored
those 38 sets across all 26,723 cells; the 5% pseudotime summaries preserve all
12,642 defined source `pt` values and 14,081 missing entries. Parse results are
under
`intermediate_objects/parse_cytokines_rest_only/reports/all_cytokine_competitor_refinement_top100_v1_20260925/`;
Bdev UMAP and trajectory results are under
`intermediate_objects/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925/`.

The follow-up driver report uses those same 38 passing profiles, saved UCell
scores, `umap.atac` coordinates, and 5% rolling-window summaries. Within each
reviewed signature, native scWorkbench mean-expression trajectory dynamics
select five drivers; a separate native median summary supplies the plotted
gene-expression lines. Every combined figure contains the score UMAP, rolling
median/IQR curve, and five colored, directly labeled driver curves. The 38 PNGs,
38-page PDF, driver tables, checks, and provenance are under
`intermediate_objects/bdev_reference/reports/all_cytokine_refined_drivers_v5_mean_internal_scale_bb_20260928/`.
The plotted mean trajectory for each driver is independently min-max scaled to
0-1 so sparse, low-expression genes are visible; raw mean and median curves
remain in separate tables and driver identities and ranks are unchanged. The
driver panel includes labeled dashed lines at the saved BB1/BB2 median
pseudotimes.

The expression-balanced sensitivity analysis retains the original UCell result
beside a score that gives every member gene exactly `1 / signature_size` weight.
For each gene, zero expression contributes zero and positive `SCT_data` values
are converted to within-gene percentiles before averaging. Each of the 38
figures compares old and new UMAPs and 5% rolling median/IQR trajectories and
retains the same five driver-expression curves. The report, cell scores, weight
table, score correlations, checks, and provenance are under
`intermediate_objects/bdev_reference/reports/all_cytokine_balanced_vs_ucell_v3_mean_internal_scale_20260928/`.

The subsequent three-method report adds a robust trajectory-level consensus.
Within each existing 5% pseudotime window, every signature gene's mean
`SCT_data` trajectory is independently scaled to 0-1; the plotted center and
band are the median and IQR across all member genes. This method has no
per-cell score or UMAP. Each figure compares its trajectory with native UCell
and the per-cell balanced score, and the driver panel marks the saved BB1/BB2
median boundaries with labeled dashed lines. Results are under
`intermediate_objects/bdev_reference/reports/all_cytokine_three_method_comparison_v1_20260928/`.

The neighborhood-coordination sensitivity analysis produces a third per-cell
score. It uses native scWorkbench neighbors in the existing 50-dimensional
RNA-PCA space, smooths each gene over local cells, converts genes independently
to zero-aware percentiles, and multiplies equal-gene mean intensity by the
fraction of participating genes. The primary setting is `k=30` with a 0.75
activity-percentile threshold; `k=15/50` and thresholds `0.60/0.90` are retained
as sensitivity variants. Its UMAPs, rolling trajectories, driver panels,
per-cell values, and sensitivity tables are under
`intermediate_objects/bdev_reference/reports/all_cytokine_neighborhood_coordination_v1_20260928/`.

The neighborhood mean-rank-change analysis retains gene-level changes as the unit
of evidence. For each signature gene, `SCT_data` expression is averaged over
the cell and its 29 nearest RNA-PCA neighbors and converted to a within-gene
percentile. Each gene's median percentile in a 5% pseudotime window is compared
with that gene five 1%-step windows earlier; the signature change is the mean
of those bounded gene-level differences. Per-cell UMAP values compare each cell's gene
percentiles with the corresponding lagged-window gene medians before taking
the signature mean. A separate coordination flag requires at least 20% of genes
to change by at least 0.05 in the score direction and is not multiplied into
the continuous score. Each of the 38 figures contains the blue-white-red change
UMAP, the change trajectory, and the unchanged top-five drivers plotted as
individual gene-rank changes with BB1/BB2 boundaries. Results are under
`intermediate_objects/bdev_reference/reports/all_cytokine_neighborhood_mean_rank_change_v2_20260928/`.

The native rank-mean comparison runs scWorkbench
`trajectory_geneset_dynamics(score_method="rank_mean")` on the same 38 sets and
compares it with saved UCell and k=30 gene-wise mean-rank-change scores. The
built-in operation ranks the union of selected genes within each cell. Because
it exports trajectory bins rather than per-cell scores, the report reconstructs
its documented cell-level formula for UMAP display and validates those values
exactly against all native 50-bin medians. The main figures use the same saved
5%-wide rolling-window memberships for all methods. Results are under
`intermediate_objects/bdev_reference/reports/all_cytokine_ucell_rankmean_genechange_v1_20260928/`.

The breadth/rank-intensity comparison applies the separately authorized custom
score to raw `RNA_counts` for the same 38 sets. Within each cell, positive genes
are ranked from highest to lowest across the full RNA feature universe, using
average ranks for ties and a rank cutoff of 1,500. Breadth is the fraction of
signature genes detected; intensity is the sum of linear rank weights divided
by the full signature size; and the combined score uses equal weights
(`alpha=0.5`). Each profile orders its UMAPs and matching 5%-window median/IQR
trajectories as native rank-mean, combined, breadth, and intensity. All 26,723
cells are scored, while trajectory plots use only the 12,642 cells with defined
existing `pt`; the 14,081 missing values remain unchanged. Results are under
`intermediate_objects/bdev_reference/reports/all_cytokine_breadth_intensity_v2_20260929/`.

Following review, the Combined flattened-rank family was promoted as the
accepted default Bdev gene-set score. Its scWorkbench artifact is
`cell_vector_family_20260929_120047_i11wgdku`. Breadth and Intensity remain
available as the separately promoted component families
`cell_vector_family_20260929_121111_0lx2kuzv` and
`cell_vector_family_20260929_122255_gny6iw1j`. Promotion used native identity
materialization (`cell_vector_family_calibrate`, `zscore` with centering and
scaling disabled), and every promoted value was checked exactly against the
reviewed score report. Recipient validation completed with only the two
expected external-source trust-boundary warnings. The promotion receipt and
validation are under
`intermediate_objects/bdev_reference/reports/all_cytokine_breadth_intensity_promotion_v1_20260929/`.

## Cytokine signatures along Bdev pseudotime

The latest review repeats the integration with all 33 passing competitor-aware
full signatures. Native UCell scores cover all 26,723 Bdev cells, and pooled
trajectory profiles use the 12,642 cells with defined existing `pt`. The former
five-gene/80% Bdev mapping gate is diagnostic in this branch: 24 signatures meet
it and nine remain plotted with 3–11 mapped genes. The profiles use 10% windows
with 1% steps and IQR bands, plus 5% and 15% sensitivity. Exact source
pseudotime and all 14,081 missing values are preserved. The review report is
`intermediate_objects/bdev_integrated/reports/competitor_aware_pt_all33_20260925/trajectory_review/REPORT.html`.
The session and scores remain unpromoted pending review.

The earlier integration scored 22 of the 23 original accepted Parse signatures on
all 26,723 Bdev cells and summarized them along the 12,642 existing `pt` values.
IFN-gamma up is excluded from scoring because only 35/50 source genes map to
Bdev, below the prespecified 80% coverage threshold. Exact source `pt` values
and all 14,081 missing entries are preserved. Directional signatures remain
separate. That earlier analysis pools all samples and cell-cycle states, with
continuous sliding-window median curves and interquartile-range bands. Windows
span 10% of the observed pt range and advance by 1%, with 5% and 15% width
sensitivity. Window membership, SD and cell counts are exported. Low-coverage
centers are marked; connecting segments across them are for display. The earlier
fixed-bin pooled report is preserved; sample and phase splits are historical only.

The integration is **awaiting figure-report review**, remains unaccepted, and
has not completed full recipient validation. The original 22-signature branch
promoted analysis artifacts before report review. The new 33-signature branch
is staged and unpromoted. Reproducible sources are in
`processing/integration/` and `figure_scripts/bdev_integrated/`.

## Workflow documentation

Execution phases and validation boundaries are documented in:

- `processing/parse_cytokines/README.md`
- `processing/scmulti_flt3ko_il7ko/README.md`
- `processing/integration/README.md`

Source selections and operation parameters are recorded after inspecting each
input and confirming the applicable scWorkbench interfaces inside the supplied
image. Future KO mapping parameters remain unset.
