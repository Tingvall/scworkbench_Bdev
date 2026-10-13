# scWorkbench B-cell development use case

This repository contains the B-cell development workflows for the scWorkbench
paper. It has four independently processed datasets and one derived integration
bundle:

1. `bdev_reference`: the reference B-cell development Seurat/Signac state and
   its accepted `pt` pseudotime.
2. `parse_cytokines`: Parse RNA cytokine-stimulation data used to derive
   cytokine response gene sets.
3. `scmulti_flt3ko_il7ko`: IL7KO/FLT3LKO RNA and ATAC multiome data used as a
   query and validation dataset.
4. `gse259355_tcf3_prob`: bulk pro-B ATAC-seq after acute E2A/TCF3 depletion,
   used to identify TCF3-dependent accessible regions.
5. `bdev_integrated`: a derived integration bundle with `bdev_reference` as
   the reference. Its current cytokine phase uses the accepted Bdev and Parse
   parents; a separate staged KO-mapping derivative uses accepted Bdev and KO
   parents without modifying either source bundle.

The source bundles remain independent and read-only during integration.
All joint state is written to `bdev_integrated`.

## Repository organization

```text
config/                       Bundle registry and project-specific settings
processing/
  bdev_reference/             Reference import and validation
  parse_cytokines/            Parse inspection, processing, and validation
  scmulti_flt3ko_il7ko/       KO multiome inspection, processing, and validation
  gse259355_tcf3_prob/         TCF3-depletion bulk ATAC processing and DESeq2 DA
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

## GSE259355 TCF3-depletion bulk ATAC bundle

The six selected bone-marrow pro-B libraries comprise three untreated controls
and three acute E2A/TCF3-depleted samples (two 6-hour and one 8-hour). After
paired-end mm10 alignment, duplicate/blacklist filtering, MACS3 peak calling,
and featureCounts quantification, the merged universe contains 84,088 peaks.
The replicate-aware DESeq2 comparison tested 84,083 peaks and identified 7,432
at BH-adjusted P < 0.05: 6,096 had lower accessibility after depletion and are
classified as TCF3-dependent, while 1,336 had higher accessibility.

The raw 84,088-by-6 consensus-peak count matrix and sample metadata are in the
accepted `gse259355_tcf3_prob` scWorkbench bundle. Its columns are explicitly
recorded as bulk ATAC libraries rather than single cells. Recipient validation
passed with only the expected external-source trust-boundary warnings. The
DESeq2 report is under
`intermediate_objects/gse259355_tcf3_prob/reports/tcf3_depleted_vs_control_da_v1/`,
and bundle validation evidence is under
`intermediate_objects/gse259355_tcf3_prob/reports/source_bundle_v1/`.

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
image. KO fragments have now been quantified over 97,788 exact-reciprocal
Bdev-reference peaks for the accepted 23,252-cell query set. PeakVI projection
was contract-tested on deterministic 256-cell reference and 256-cell query
subsets: the frozen-model projection, model-bound latent union, and display
UMAP all passed without promotion. The full GPU reference fit and frozen-model
projection of all 23,252 KO cells are now complete in an unpromoted durable
session; the 49,975-cell joint latent space and display UMAP passed dimension,
ordering, finite-value, and no-refit checks. Mapping QC and any downstream
label or pseudotime transfer were kept as a separate review stage.

The subsequent Bdev-restricted mapping QC and transfer is now complete and
unpromoted. Exact k=15/30/50 RANN searches were used because the registered
scWorkbench neighbor operation cannot perform the required reference-restricted
search at this scale. Five-fold Bdev calibration gave k=30 cluster accuracy
0.800 (balanced 0.787) and pseudotime Spearman correlation 0.990. Of 23,252 KO
cells, 20,595 are supported, 1,876 borderline, and 781 out-of-reference by the
held-out Bdev distance calibration. The report is
`intermediate_objects/bdev_integrated/reports/ko_mapping_qc_v1_20260930/REPORT.md`;
no results have been promoted into a bundle.

A subsequent report-only visualization places supported and borderline KO
query cells onto the unchanged native Bdev `umap.atac` coordinates using the
distance-weighted mean of their 30 exact Bdev-restricted PeakVI neighbors.
k=15/50 coordinates and neighbor-UMAP dispersion are retained as sensitivity
diagnostics; the 781 out-of-reference cells are withheld from the primary
display. Results are under
`intermediate_objects/bdev_integrated/reports/ko_mapping_reference_umap_v1_20261001/`.

The two-image positive-signal follow-up combines the accepted Bdev score and
driver trajectories with that reference-anchored projection and query WT/KO
score trajectories for Flt3L and IL-7. Query RNA counts use the same accepted
combined flattened-rank formula, with complete 45/45 and 20/20 feature coverage
respectively. Every trajectory panel includes dashed BB1/BB2 boundaries, and
low-support query windows remain explicit gaps. The report is under
`intermediate_objects/bdev_integrated/reports/ko_cytokine_positive_signal_v1_20261001/`.

The v2 display revision keeps the accepted positive-score UMAP at a larger
display size with its scale on the right and retains the Bdev score trajectory,
removes the driver-gene panel, and places
the WT/KO mapped-pseudotime cell distribution directly below the Bdev
trajectory. That distribution is the fraction of supported/borderline cells
within each genotype, displayed as a centered three-bin rolling mean over 100
underlying bins, not a cytokine score. The two rows have equal height; panels B
and D share the same horizontal span, with deliberate whitespace around A, and
both retain dashed BB1/BB2 boundaries. C restores the v1 point styling and
omits its redundant genotype legend; A uses prominent colored points so its
score remains visible. The two revised
figures are under
`intermediate_objects/bdev_integrated/reports/ko_cytokine_positive_distribution_v2_20261001/`;
the v1 report remains preserved for provenance.

The v3 report preserves v2 and enlarges panel C's projected points with
v1-like styling. Panel D corrects for the nonuniform Bdev pseudotime density:
within each centered three-bin window, the query fraction is divided by the
fraction of the 12,642 Bdev cells with defined `pt` in the same bins. A ratio of
1 therefore means that query and Bdev reference densities match. Results are
under
`intermediate_objects/bdev_integrated/reports/ko_cytokine_positive_distribution_bdev_normalized_v3_20261001/`.

The v4 six-panel revision adds `Flt3` or `Il7r` expression on native Bdev
`umap.atac` and its mean `SCT_data` trajectory in the exact saved 5%-wide,
1%-step windows. The v3 panels move to C--F. All trajectory panels B, D, and F
share identical horizontal placement and pseudotime limits and retain dashed
BB1/BB2 boundaries. Results are under
`intermediate_objects/bdev_integrated/reports/ko_cytokine_gene_expression_bdev_normalized_v4_20261001/`;
v3 remains preserved.

The follow-up region-wise RNA comparison extracts the maximal WT intervals in
v4 panel F whose Bdev-normalized rolling fraction is strictly above 1. It marks
WT and KO query cells in one Flt3L interval and two separate IL-7 intervals,
then runs native scWorkbench Wilcoxon DE on staged log1p-CP10K RNA. The result
is unpromoted and available under
`intermediate_objects/bdev_integrated/reports/ko_wt_enriched_peak_rna_de_v1_20261001/`.
Because there is one sample per genotype for each cytokine, these are
descriptive cell-level contrasts and not biological-replicate-level inference.
The report highlights the exact Flt3L and two IL-7 pseudotime intervals used
for DE as shaded regions over the WT/KO normalized trajectory distributions.
The report also contains region-specific WT/KO violins for significant genes
with absolute effect above 0.2. WT-higher and KO-higher genes are clearly
separated and capped at the 10 most significant genes per direction and region.
It also includes one two-panel violin figure comparing all four samples using
B- and T-lineage gene sets from supplementary `Table S1B`. The aggregate score
is mean log1p-CP10K expression across exactly mapped genes (246 B-lineage and
233 T-lineage genes); source membership and mapping coverage are exported.
The same report shows all four samples for the accepted Parse-derived Flt3L and
IL-7 UP/DOWN signatures in a single four-panel violin figure. All 105 signature
genes map exactly, and the per-cell score is mean log1p-CP10K expression across
each profile's members.
All-four-sample trajectory profiles are included for B/T lineage scores and
Flt3L-UP/IL-7-UP signals. They retain unsupported intervals as gaps and mark
BB1/BB2 with dashed lines.

The Ebf1 enhancer follow-up uses the exact two IL-7 enriched-trajectory cell
assignments from that report and separately compares IL-7 WT with KO across
all 259 Bdev peaks overlapping ±3 Mb around the GRCm38/mm10 Ebf1 TSS. Native
scWorkbench cell-level Wilcoxon tests on the Bdev peak-count matrix identify
23 significant peaks in IL-7 peak 1 (21 WT-higher, 2 KO-higher) and 128 in
IL-7 peak 2 (127 WT-higher, 1 KO-higher) at within-region BH-adjusted P < 0.05.
Four and one supplied peaks, respectively, had no finite test statistic. The
unpromoted report and complete tables are under
`intermediate_objects/bdev_integrated/reports/ebf1_enhancer_il7_enriched_pt_da_v1_20261001/`.
These remain descriptive cell-level results because each genotype has one
sample.

A follow-up combined ATAC track figure displays those ten most significant
IL-7 peak-2 loci at ±5 kb resolution. Flt3L WT/KO and IL-7 WT/KO are all
restricted to the same IL-7 peak-2 pseudotime interval. Fragment-multiplicity-
weighted Tn5 insertions are shown in 25-bp bins with a centered 125-bp rolling
mean and sample-specific total-fragment normalization. The report is under
`intermediate_objects/bdev_integrated/reports/ebf1_enhancer_il7_peak2_atac_tracks_v1_20261001/`.

The differential-accessibility calculation was subsequently repeated after
native scWorkbench per-cell log1p-CP10K normalization across the complete
97,788-peak Bdev matrix. IL-7 peak 1 retains all 23 raw-count significant peaks;
IL-7 peak 2 retains 121 of 128, with no normalized-only additions. Region 01
remains the most significant normalized peak, although its effect estimate is
reduced from 0.145 to 0.070. Full normalized and raw-versus-normalized results
are under
`intermediate_objects/bdev_integrated/reports/ebf1_enhancer_il7_enriched_pt_lognorm_da_v1_20261001/`.

The ATAC track figure was then revised to show the ten normalized,
BH-significant IL-7 peak-2 regions with higher accessibility in IL-7 WT than
IL7KO. Selection uses only the IL-7 comparison; all four samples are displayed
for context. The normalized-ranked v2 figure is under
`intermediate_objects/bdev_integrated/reports/ebf1_enhancer_il7_peak2_lognorm_atac_tracks_v2_20261002/`,
and the raw-ranked v1 is preserved.

The current v3 ATAC tracks correct the earlier display-normalization mismatch:
each deduplicated fragment record is counted once, the multiplicity/read-support
column is ignored, and coverage is normalized within each cell by its complete
Bdev peak-matrix depth before log1p transformation and averaging across all
cells. The matrix-aligned report is under
`intermediate_objects/bdev_integrated/reports/ebf1_enhancer_il7_peak2_matrix_aligned_atac_tracks_v3_20261002/`.
The v1/v2 figures are retained for provenance but superseded for biological
interpretation.

The subsequent three-way Ebf1 enhancer screen evaluates all 259 Bdev peaks in
the mm10 Ebf1 TSS ±3 Mb window. It combines GSE259355 TCF3-depletion peaks
(`abs(log2FC) > 1`, BH-adjusted P < 0.05) mapped to Bdev peaks, Ebf1
expression/accessibility partial Spearman tests in the exact 971 Bdev cells in
`Il7_shared_peak2`, and GSE162858 230-238 H3K4me3 PLAC-seq promoter contacts
with distal anchors expanded by 20 kb. Eight Bdev peaks have TCF3 support, 63
have promoter-contact support, and none pass the correlation BH threshold. The
distal union contains 64 putative enhancers; four have both TCF3 and PLAC-seq
support. The complete 259-peak audit and one-to-many mappings are under
`intermediate_objects/bdev_integrated/reports/ebf1_putative_enhancers_three_way_v2_20261005/`.

Normalized GSE259355 coverage tracks were generated across the same Ebf1
±3 Mb interval by normalizing each filtered paired-end library to fragment
coverage CPM and then averaging the three replicates per condition with equal
weight. The full-region figure includes all 28 230238 H3K4me3 PLAC-seq
Ebf1-promoter links; a second figure shows ±10-kb views of all eight passing
TCF3 peaks. Figures, 50-bp condition-average bedGraphs, replicate tracks, and
normalization audits are under
`intermediate_objects/bdev_integrated/reports/tcf3_ebf1_coverage_placseq_tracks_v4_20261005/`.

The saved RNA log1p-CP10K matrix and complete 33,696-gene scWorkbench
Wilcoxon tables for the same exact IL-7 peak-2 cells were reused without
repeating normalization or differential testing. Four separate violin panels
show the top ten WT-higher and KO-higher genes for FLT3L WT versus FLT3LKO and
IL-7 WT versus IL7KO, ranked within direction by genome-wide BH-adjusted P.
The unpromoted report and selection table are under
`intermediate_objects/bdev_integrated/reports/il7_peak2_top_de_genes_v1_20261002/`.
These remain descriptive cell-level comparisons because each group contains
one biological sample.

A subsequent effect-size-filtered report replaces the native scWorkbench
effect statistic with standard
`log2((mean linear CP10K WT + 1)/(mean linear CP10K KO + 1))` for filtering,
while retaining the saved native Wilcoxon P values and genome-wide BH
adjustment. Requiring adjusted P < 0.05 and at least a 1.5-fold change retains
14 IL-7 WT-higher, 13 IL7KO-higher, 83 FLT3L WT-higher, and 96 FLT3LKO-higher
genes in the same exact IL-7 peak-2 cells. Paginated violin plots and the full
206-gene table are under
`intermediate_objects/bdev_integrated/reports/il7_peak2_fc1p5_bh05_de_genes_v1_20261002/`.

The current shared-enrichment analysis defines one region per cytokine from
maximal contiguous pseudotime bins whose Bdev-normalized rolling cell fraction
is strictly above 1 in both WT and KO. FLT3L uses bins 19--26
(`pt` -25.397 to -22.919); IL-7 is restricted to its original peak 2 and uses
bins 21--31 (`pt` -24.777 to -21.370). The saved log1p-CP10K RNA matrix was
reused and new native scWorkbench Wilcoxon tests evaluated all 33,696 genes in
each changed cell set. With genome-wide BH-adjusted P < 0.05 and
`|log2FC| >= log2(1.5)`, 187 FLT3L genes pass (92 WT-higher, 95 KO-higher) and
24 IL-7 genes pass (13 WT-higher, 11 KO-higher). The report includes the
selected-region plot, new volcano plot, full normalized-count/DEG tables, and
up to 20 most significant cutoff-passing genes in each direction as separate
WT-higher and KO-higher violin figures:
`intermediate_objects/bdev_integrated/reports/ko_shared_enriched_peak_rna_de_v1_20261002/`.

The full detected-gene rankings from those shared-coverage comparisons were
also analyzed by user-approved preranked GSEA outside scWorkbench, whose
available enrichment operation implements over-representation analysis rather
than GSEA. KO-over-WT standard log2FC ranks were tested against the 50
mouse-native MSigDB 2026.1.Mm Hallmark sets using `fgseaMultilevel`; positive
NES denotes KO-higher enrichment. The frozen gene-set memberships, complete
rankings, leading edges, tables, and figures are under
`intermediate_objects/bdev_integrated/reports/flt3l_il7_ko_vs_wt_hallmark_gsea_v1_20261010_run2/`.

The project default for GO trajectory gene-set analysis is now the completed
RNA-feature-residualized UCell run
`go_ucell_rna_feature_residual_dynamic_range_perm1000_v1_20261012`. It uses the
unchanged Bdev metadata `pt`, 5%-width/1%-step windows, the median residualized
UCell score, 1,000 peak-to-trough dynamic-range permutations, groupwise BH
q < 0.01, and trajectory-wide positive-region p < 0.01. This is the analysis
referred to informally as residual-UCell-based GSEA; statistically it is a
trajectory gene-set permutation analysis rather than preranked GSEA. The
default pattern summary uses only supported regions and the reviewed k-means
resolution BP=25, MF=20, and CC=20. The bundle was not mutated. Older raw-UCell,
z-score, and lower-k reports remain comparison records rather than defaults.

The same residual-UCell protocol has now been completed for the 38 accepted
cytokine signatures: residualize each per-cell native UCell score for
`log1p(nFeature_RNA)` using a four-degree-of-freedom natural spline without
pseudotime, then rerun the 1,000-permutation trajectory dynamic-range and
positive-region tests on those residual scores. The completed default run is
`cytokine_ucell_rna_feature_residual_dynamic_range_perm1000_v1_20261012`:
32/38 signatures are dynamically significant, with 56 supported positive
peaks merged into 34 regions. Median absolute score correlation with
`nFeature_RNA` decreased from 0.187 to 0.0483, and the bundle fingerprint was
unchanged. The earlier cytokine run
`cytokine_ucell_dynamic_range_perm1000_v1_20261012` tested raw UCell scores.
Its residual-UCell follow-up is a curve/complexity diagnostic only and states
that permutation calls were not recomputed; it and the raw-UCell calls are now
historical comparisons rather than the default inference.

For WT-versus-KO preranked GSEA, the selected test intervals are the mapped
shared-enrichment regions, not the later Bdev residual-UCell activity regions.
They retain substantial coverage in both genotypes: FLT3L has 4,394 WT and
3,342 KO cells from `pt` -25.397 to -22.919; IL-7 has 1,912 WT and 1,355 KO
cells from -24.777 to -21.370. By contrast, the early Bdev FLT3L residual
region contains 672 WT but only 35 KO cells, and the Bdev IL-7 region contains
274 WT but only 20 KO cells; the later FLT3L region has only 11 WT and 6 KO
cells. These depleted later regions remain biological evidence consistent with
blocked KO progression, but they are unsuitable for a balanced within-stage
WT-versus-KO expression ranking. The completed analysis
`flt3l_il7_ko_vs_wt_pseudotime_matched_gsea_v1_20261012` therefore uses
pseudotime-bin-matched KO-minus-WT ranks within the mapped shared-enrichment
intervals. Each bin contributes its KO-versus-WT log2 mean-CP10K effect with a
Bdev-density weight. It tests mouse Hallmark, GO BP/MF/CC, and the 38 accepted
cytokine signatures with `fgseaMultilevel` and collection-specific BH control.

FLT3L yields 21 significant Hallmarks (19 WT-higher, 2 KO-higher); WT retains
MYC/E2F, oxidative-phosphorylation, mTORC1, DNA-repair, and related biosynthetic
programs. IL-7 yields seven significant Hallmarks (five WT-higher, two
KO-higher); WT retains MYC/E2F/G2M and DNA-repair programs, while KO is higher
for coagulation and KRAS-signaling-up annotations. Focused cytokine GSEA places
IL-7-UP on the WT side in both comparisons and IL-7-DOWN on the IL7KO side.
GO terms are more redundant and include broad shared-gene annotations, so they
are secondary to Hallmark and leading-edge interpretation. Excluding bins with
fewer than 20 cells in either genotype gives Hallmark NES correlations of 0.988
for FLT3L and 0.991 for IL-7. Results remain descriptive because each genotype
has one biological sample.

The dynamic follow-up
`flt3l_il7_ko_vs_wt_dynamic_change_gsea_v1_20261012` asks a different question:
which programs change more from the early to late portion of the mapped region
in WT than in KO. Genes are ranked by `(WT late - WT early) - (KO late - KO
early)` using the first and last two bins with at least 20 cells in each
genotype. The sparse terminal bin is reserved for sensitivity analysis.
FLT3L has strict failed/attenuated-induction examples in E2F targets and the
unfolded-protein response, with additional DNA-replication, repair, electron
transport, and respiratory-chain GO terms. IL-7 has strict examples in TNFA
signaling via NFKB and the early estrogen-response Hallmark. Many other
positive differential-change enrichments reflect a stronger decline in KO
rather than an absolute rise in WT; the report and curve labels keep those
categories separate. Results are descriptive because there is one biological
sample per genotype.

All GO BP, MF, and CC terms were retained in the dynamic test (10,540
comparison-specific tests). The focused top-GO report shows at most five
strict WT-induction/KO-attenuation terms per GO group and comparison rather
than presenting the full redundant inventory. It is stored under
`flt3l_il7_ko_vs_wt_dynamic_change_top_go_v1_20261012`.
Curves for the eight strict GO terms that also pass the sparse-terminal-bin
sensitivity analysis are under
`flt3l_il7_ko_vs_wt_fdr_supported_go_curves_v1_20261012`.
The reciprocal failure-to-silence screen requires a WT decrease and a more
upward KO change. Its compact top-five-per-GO-group report and sensitivity-
supported curves are under
`flt3l_il7_ko_vs_wt_go_failure_to_silence_v1_20261012`.

The matched individual-gene screen is under
`flt3l_il7_ko_vs_wt_individual_gene_divergence_v1_20261012`. It reports both
WT-induction/KO-attenuation and failure-to-silence directions, requires the
direction to survive the sparse-terminal sensitivity policy, and ranks genes
by the smaller divergence across endpoint definitions. Curves show the
mean-expression trajectory with cell-level IQR ribbons; the ribbons describe
heterogeneity rather than inferential uncertainty. The gene ranking is
descriptive because each genotype has one biological sample.

A preliminary Bdev modality-timing follow-up uses the accepted 30-cell windows
with 5-cell stride and the robust-LOESS-smoothed log10 RNA/ATAC progression
slope. Centered derivative sign changes produced 24 candidate extrema; requiring
at least 1.10-fold local prominence retained four exact 30-cell turning windows
centered at `pt` -24.586, -15.130, -12.458, and -6.524. These regions are
highlighted consistently on the slope, smoothed Ebf1/Pax5/Ezh2 expression, and
the existing Bdev `umap.atac` embedding. The read-only, unregistered report is
under
`intermediate_objects/bdev_reference/reports/preliminary_slope_turning_regions_w30_s5_span008_20261008/`.

The compact follow-up display omits Pax5; vertically aligns the slope, Ebf1,
and Ezh2 panels; and adds the RNA-versus-ATAC trajectory beside a UMAP of the
exact turning-window cells. It is under
`intermediate_objects/bdev_reference/reports/preliminary_turning_summary_ebf1_ezh2_w30_s5_v2_20261008/`.

The next preliminary display adds Nipbl, Mau2, and Wapl and colors trajectory
cells on `umap.atac` by their membership-weighted mean
`log2(smoothed delta-RNA/delta-ATAC)`. Red denotes RNA-faster progression, blue
ATAC-faster progression, and white matched progression; cells outside the
trajectory are gray. The report is under
`intermediate_objects/bdev_reference/reports/preliminary_modality_ratio_umap_five_genes_20261008/`.

The expanded nine-gene version additionally shows Ctcf, Eed, Tcf7, and Lef1:
`intermediate_objects/bdev_reference/reports/preliminary_modality_ratio_umap_nine_genes_v3_20261008/`.
Tcf7 alone required a recorded Gaussian-LOESS fallback after robust LOESS
failed numerically on its low, sparse expression trajectory.

A separate preliminary figure replaces the individual-gene curves with
B-lineage, T-lineage, Flt3L-UP, and IL-7-UP arithmetic-mean `SCT_data` scores
under the same slope panel and T1–T4 bands. Bdev mapping retains 236, 214, 45,
and 20 genes, respectively. The figure and complete membership/coverage tables
are under
`intermediate_objects/bdev_reference/reports/preliminary_slope_gene_set_signals_v2_20261008/`.

An LSI1 sensitivity check compares the accepted full LSI1–50 ATAC path with
LSI2–50 while holding RNA coordinates and all windows fixed. Excluding LSI1
changes relative ATAC distance only slightly (median absolute shift 0.00074;
maximum 0.00224), and the smoothed log-slope curves correlate at 0.99995. The
comparison is under
`intermediate_objects/bdev_reference/reports/preliminary_lsi1_sensitivity_20261008/`.

A separate Tsix figure shows `SCT_data` expression on the existing Bdev
`umap.atac` beside its robust-LOESS-smoothed 30-cell/5-cell pseudotime profile.
The exact T1–T4 cells are outlined and labeled on the UMAP, and the same fixed
turning windows are banded and labeled on the trajectory. The read-only,
unregistered report also contains a standalone white-to-crimson expression
UMAP without T1–T4 overlays and a standalone zero-baseline pseudotime plot
using the source BB1/BB2 ranges in place of T1–T4. It is under
`intermediate_objects/bdev_reference/reports/preliminary_tsix_expression_v2_20261009/`.

The local modality-balance curve is also available explicitly as
`log2(delta RNA / delta ATAC)`, centered at zero, with the exact inverse
`log2(delta ATAC / delta RNA)` on the secondary axis. T1–T4 remain fixed and
unchanged. The report and a proposed autocorrelation-aware gene-pattern screen
are under
`intermediate_objects/bdev_reference/reports/preliminary_log2_modality_ratio_v1_20261009/`.

The two proposed gene-pattern screens were implemented over 11,745 detected,
variable Bdev `SCT_data` genes using the accepted 30-cell/5-cell windows and a
fixed centered 8%-window smoother. The state-following screen retains 4,884
genes (486 RNA-faster and 4,398 ATAC-faster/inverse), consistent with broad
developmental-state structure. The change-following screen is much more
selective and retains four RNA-faster genes—Angpt1, Adgrl4, Gcnt2, and Mecom—and
no inverse genes at pooled-null BH < 0.05 and |correlation| >= 0.3. Circular
shifts preserve trajectory autocorrelation; results remain descriptive because
feature-fold cross-fitting of the RNA trajectory was not performed. The report
is under
`intermediate_objects/bdev_reference/reports/preliminary_modality_gene_pattern_screen_v2_20261009/`.

Separate figures show the top 20 RNA-faster state genes and top 20 ATAC-faster
inverse state genes as individual smoothed trajectories with within-window
cell-level IQR ribbons and their matching modality target. A companion figure
shows the four selected change-associated genes (Angpt1, Adgrl4, Gcnt2, and
Mecom), with expression/IQR panels above the derivative profiles used by the
screen. These read-only outputs are under
`intermediate_objects/bdev_reference/reports/preliminary_state_top20_separate_figures_v1_20261009/`.

Because whole-trajectory correlation can favor genes that barely change over
long quiet portions, an event-focused replacement now ranks genes by up/down
agreement over the strongest RNA/ATAC changes and explicitly across both
flanks of T1-T4. A p95-p05 `SCT_data` range of at least 0.50 removes nearly
flat profiles; selection additionally requires six of eight matching flanks
and at least two complete turns. Four direct genes (Clint1, Arhgap15, Camk2d,
and Baz2b) and six inverse genes (Cox4i1, Rplp2, Rps6, Impdh2, Rpl37, and
Rpl12) pass. Separate top-20 direct and inverse figures include within-window
IQR, while complete scores and per-turn calls remain available under
`intermediate_objects/bdev_reference/reports/preliminary_modality_turn_pattern_screen_v2_20261009/`.

Pass/fail trajectory plots for the chromatin-regulator candidates are under
`intermediate_objects/bdev_reference/reports/preliminary_epigenetic_pattern_candidates_v1_20261009/`.
They contain separate direct and inverse panels plus a focused cohesin/PRC2
audit. Each panel overlays standardized expression and within-window IQR on
the relevant modality target and reports the expression range, matched flanks,
and complete turns. Baz2b is the sole strict passing epigenetic candidate.

The cohesin set is also plotted specifically against the inverse ATAC/RNA
target under
`intermediate_objects/bdev_reference/reports/preliminary_inverse_cohesin_turn_audit_v2_20261009/`.
Green/red trajectory markers and a companion pass/fail heatmap identify every
T1-T4 pre-flank, post-flank, and complete-turn call. Wapl passes both flanks at
T2 and T3 and the T4 pre-flank, but remains a strict near-hit at five of eight
flanks and p95-p05 range 0.461.

An unbiased continuous replacement excludes T4 and ranks all 11,745 detected
genes without regulator preselection or binary turn cutoffs. The raw T1-T3
score combines continuous derivative similarity (60%) and six-flank vector
similarity (40%); amplitude and within-window IQR/SD contribute only smooth
evidence weights. Nipbl has raw direct shape score 0.788 and ranks 36th after
evidence weighting, while Smarca4 ranks 5th. Separate direct and inverse
top-20 IQR figures and full genome-wide score tables are under
`intermediate_objects/bdev_reference/reports/preliminary_t1_t3_continuous_gene_screen_v1_20261009/`.

For the more specific question of concordance around every T1-T3 event, a
turn-local balanced screen removes full-pseudotime correlation entirely and
scores all six pre/post flanks separately. It balances the two flanks within
each turn and then T1, T2, and T3 equally, with consistency penalties, so a
single strong peak cannot conceal an opposite response at another turn. This
demotes Igll1 from a visually T2-driven hit to direct rank 11,723 because its
T3 score is strongly opposite (-0.254), while Smarca4 ranks first and Nipbl
ranks 21st. The figures retain the entire pseudotime range and shade only the
scored +/-100-window T1-T3 regions; T4 is dashed unscored context. Outputs are
under
`intermediate_objects/bdev_reference/reports/preliminary_turn_local_balanced_gene_screen_v2_20261009/`.

A matched depth diagnostic compares the inverse top 20 with total RNA
(`nCount_RNA`) and total ATAC (`nCount_peaks`) using the same windows,
smoothing, complete pseudotime display, and shaded +/-100-window T1-T3
regions. Eleven of the 20 are Rpl/Rps genes. Their median absolute
scored-region correlation is 0.837 with total RNA but only 0.222 with total
ATAC; non-ribosomal hits are almost as correlated with total RNA (0.815).
Moreover, the ribosomal median cell-level RNA-depth correlation falls from
0.455 to 0.074 after removing an 8-df pseudotime trend. This supports some
RNA-depth or compositional contribution, but not a ribosomal-specific artifact
explained solely by raw library size. Figures and the full correlation table
are under
`intermediate_objects/bdev_reference/reports/preliminary_inverse_hits_depth_diagnostic_v1_20261009/`.

Total RNA and ATAC jointly explain 87.7% of the broad RNA/ATAC ratio level
in-sample and 78.7% under 10 contiguous-fold prediction. This should not be
interpreted as the same fraction of the local turning pattern: total-signal
changes explain only 23.8% of the full-path ratio derivative and 7.4% within
the scored T1-T3 regions. A model trained on two turns fails to predict the
centered derivative shape of the held-out turn (pooled R2 -0.043). Thus total
signal is a major correlate of the broad baseline but explains little of the
local signal used in the turn-focused gene ranking. Fitted, residual, and
derivative-residual curves are under
`intermediate_objects/bdev_reference/reports/preliminary_ratio_total_signal_explained_v1_20261009/`.

The preferred depth-removal sensitivity analysis is upstream and conditional
on preserved Bdev pseudotime: regress each RNA-PC on spline(pt) plus log-total
RNA and each ATAC-LSI dimension on spline(pt) plus log-total ATAC, retain the
spline(pt) component plus residual, and then recompute window centroids,
cumulative distances, the ratio, and gene ranks. This targets within-stage
depth dependence while avoiding automatic removal of biological changes in
cell size or transcriptional output. Direct subtraction of the fitted ratio
curve is retained only as a more aggressive post-hoc sensitivity analysis.

The post-hoc residual ratio and residual ratio-change curves are also shown in
separate vertically aligned full-pseudotime figures above Ebf1, Ezh2, Nipbl,
and Wapl expression with within-window IQR. These displays replace the T1-T3
bands with the exact source `bb_up` BB1/BB2 pseudotime ranges and median lines.
They are under
`intermediate_objects/bdev_reference/reports/preliminary_residual_ratio_bb_gene_stack_v4_20261009/`.

The same full-path residual ratio is mapped onto the existing Bdev
`umap.atac` by linear interpolation from window-center rank to every cell's
exact trajectory rank. This increases distinct mapped values from 2,530 to
12,614 among 12,642 cells and reduces adjacent exact ties from 80.0% to 0.22%;
the original membership-mean projection is retained in the per-cell table. A
zero-centered blue-white-red scale shows residual ATAC-faster, expected, and
residual RNA-faster balance, respectively; cells outside the trajectory are
gray. Both the absolute-p99-clipped symmetric scale and an unclipped symmetric
scale using the observed maximum absolute value are exported. The figure and
per-cell values are under
`intermediate_objects/bdev_reference/reports/preliminary_residual_ratio_umap_v3_20261009/`.

A genome-wide descriptive screen compares depth-adjusted gene-expression
up/down patterns with this corrected ratio without using BB or T-region
definitions. The ratio and genes receive identical smoothing, and their slopes
are normalized by their own scales before comparison, allowing magnitudes to
differ. Every window contributes, while ratio magnitude, slope, and curvature
continuously upweight informative parts. Eight equal-interest trajectory
segments receive balanced emphasis. A separate full-pattern audit retains only
genes with the requested direction in all eight segments: 15 direct and 71
inverse candidates, although the similarities remain modest rather than exact.
Because the RNA representation used
to construct the ratio contains the genes being screened, these ranks remain
descriptive; feature-fold cross-fitting is required for unbiased inferential
ranks. The direction screen and full-pattern audit are under
`intermediate_objects/bdev_reference/reports/preliminary_residual_ratio_direction_gene_screen_v4_20261009/`
and
`intermediate_objects/bdev_reference/reports/preliminary_residual_ratio_full_direction_matches_v1_20261009/`.

The same conditional RNA-depth adjustment and matched smoothing are shown for
Ebf1 together with the focused cohesin-associated set, and separately for Ebf1
with Ezh2, Eed, and Suz12. The full pseudotime is displayed, ribbons are
cell-level adjusted-expression IQR, and each label reports direct versus
inverse agreement across the eight balanced trajectory segments. Results are
under
`intermediate_objects/bdev_reference/reports/preliminary_adjusted_cohesin_prc2_profiles_v1_20261009/`.

A compact landscape version restricted to Ebf1, Tcf3, Stat5a, Stat5b, Ezh2,
Nipbl, and Wapl places the corrected ratio alongside standardized,
RNA-depth-adjusted expression. It shades the exact source BB1/BB2 ranges and
median positions in every panel and is under
`intermediate_objects/bdev_reference/reports/preliminary_adjusted_selected_genes_bb_v1_20261009/`.

Nipbl, Wapl, and Ezh2 `SCT_data` expression is also shown on the existing Bdev
`umap.atac`. Gene-specific white-to-crimson scales run from zero to each
positive-cell p99, with higher values clipped only for display and high cells
drawn last. The figure and per-cell values are under
`intermediate_objects/bdev_reference/reports/preliminary_nipbl_wapl_ezh2_umap_v2_20261009/`.

Cox6a1, Cox6a2, and Ube2c `SCT_data` expression is shown similarly on the
unchanged Bdev `umap.atac`, using independent white-to-crimson scales capped at
each gene's positive-cell p99 for display. Uncapped cell-level values and scale
audits are under
`intermediate_objects/bdev_reference/reports/preliminary_cox6a1_cox6a2_ube2c_umap_v1_20261010/`.

The FLT3L/IL-7 follow-up includes separate motif-accessibility heatmaps for the
exact `Flt3l_peak1` and `Il7_peak2` intervals. Motifs are selected solely by
their temporal variability in Bdev: Bdev also fixes accessible peaks, chromVAR
background/expectation, row order, centering, and scaling. WT and knockout
cells enter only afterward using mapped Bdev pseudotime, with sparse bins
(fewer than 20 cells) left blank and a descriptive WT-minus-KO panel. The
unpromoted report is under
`intermediate_objects/bdev_integrated/reports/bdev_selected_motif_pt_heatmaps_v1_20261010_run5/`.
This run also exports trajectories for all 746 Bdev-eligible motifs.

A full-universe follow-up retains every motif with a strong Bdev up/down
pattern, defined by at least 1.5 Bdev temporal SD early-to-late change,
absolute pseudotime Spearman correlation of at least 0.6, and concordant
directions. It compares WT and KO within matched mapped-pseudotime bins and
weights those contrasts to the Bdev cell distribution. The complete score
table and paginated figures are under
`intermediate_objects/bdev_integrated/reports/bdev_strong_directional_motif_trajectory_adjusted_wt_ko_v1_20261010/`.

Companion heatmaps show the 20 strongest adjusted WT-higher and 20 strongest
adjusted KO-higher motifs separately for Flt3L peak 1 and IL-7 peak 2. Bdev,
WT, KO, and matched-bin WT-minus-KO panels are aligned across mapped
pseudotime; stars identify motifs passing the strong descriptive difference
flag. Results are under
`intermediate_objects/bdev_integrated/reports/top_adjusted_wt_ko_motif_heatmaps_v1_20261010/`.

A single-cell-axis version retains the same selected motif rows but uses every
Bdev cell in each interval as an ordered focal column. Each panel is smoothed
independently from its 30 nearest pseudotime cells, and WT-minus-KO is evaluated
at identical Bdev focal positions. The displayed interval is restricted to the
largest contiguous span where both query groups have a 30-neighbor radius no
greater than 0.25 pseudotime units; full-axis values remain exported. Results
are under
`intermediate_objects/bdev_integrated/reports/single_cell_30nn_top_adjusted_motif_heatmaps_v1_20261010_run4/`.

SNAI2, TCF3, RUNX1, SREBF2, and STAT5A::STAT5B motif activity is also shown
on the unchanged Bdev `umap.atac`. A companion violin figure compares WT and
KO cells within Flt3L peak 1 and IL-7 peak 2 using per-cell chromVAR deviations
from region-specific Bdev-defined calibration. The raw violin distributions
are annotated with matched-bin WT-minus-KO effects weighted to Bdev. Outputs
are under
`intermediate_objects/bdev_integrated/reports/selected_motif_bdev_umap_query_violins_v2_20261010/`.

A row-aligned composite repeats this display for TCF3, EBF1, ZEB-family,
STAT5A::STAT5B, and RUNX1, with all Bdev UMAPs on the left and the corresponding
Flt3L/IL-7 WT/KO violins on the right. ZEB1 is explicitly labeled as a proxy
because the motif collection contains no ZEB2 entry. Results are under
`intermediate_objects/bdev_integrated/reports/tcf3_ebf1_zeb_stat5_runx1_composite_v1_20261010/`.

A focused trajectory view uses the same seven motifs in both cytokine regions:
EBF1, PAX5, TCF3, CEBPG as the most variable C/EBP-family motif, HLF, BACH2,
and LEF1. Each panel overlays Bdev in black, WT in red, and the corresponding
knockout in blue; the y-axis is chromVAR deviation centered and scaled using
Bdev, with motif-specific ranges. Results are under
`intermediate_objects/bdev_integrated/reports/selected_tf_motif_pt_curves_v1_20261009/`.

The corresponding unbiased RNA screen selects the 30 most temporally variable
eligible genes independently within Bdev `Flt3l_peak1` and `Il7_peak2`, then
displays Bdev, WT, KO, and WT-minus-KO heatmaps. Both datasets use native
scWorkbench log1p-CP10K normalization, and all selection and scaling is fixed
from Bdev before query display. Every nonempty bin is calculated; black
borders/crosses flag bins with fewer than 20 cells and gray marks zero-cell
bins. Results are under
`intermediate_objects/bdev_integrated/reports/bdev_variable_gene_expression_pt_heatmaps_v1_20261009_run4/`.

Companion difference heatmaps compare WT and KO separately with Bdev in the
same pseudotime bins, retaining the source top-30 genes and Bdev scaling.
Positive values mean higher query expression and negative values mean higher
Bdev expression. Results are under
`intermediate_objects/bdev_integrated/reports/query_vs_bdev_gene_expression_diff_heatmaps_v1_20261009_run2/`.

The KO-specific screen first limits its candidate pool, independently per
region, to the top quartile by Bdev temporal SD after basic Bdev-expression
eligibility. It then ranks those genes for large KO-minus-Bdev and small
WT-minus-Bdev divergence across matched pseudotime bins. Sparse bins remain
included with reduced reliability weight. The
unfiltered result contains many ribosomal/translation genes, particularly in
the IL-7 comparison, alongside lineage-related Flt3L hits; no technical gene
class was removed. Each figure also shows the corresponding Bdev expression
beside the WT-minus-Bdev and KO-minus-Bdev panels. Results are under
`intermediate_objects/bdev_integrated/reports/ko_specific_bdev_variable_gene_expression_diff_heatmaps_v1_20261009/`.
