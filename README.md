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
4. `bdev_integrated`: a derived integration bundle that consumes the three
   accepted parent bundles with `bdev_reference` as the reference.

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

## Parse cytokine bundle

The accepted `parse_cytokines` bundle is built from the checksummed
pre-filtered Seurat RDS at
`intermediate_objects/parse_cytokines/prepared/parse_cytokines_prefiltered_min500counts_seurat.rds`.
It contains the raw `RNA_counts` matrix for 251,599 cells and 78,298 features.
Ingestion retained only the required `sample_id`, `replicate`, and `condition`
metadata for all 96 samples; the reviewed treatment step later added the
derived `treatment_category` and `cytokine_count` fields. It contains no
normalized matrix or dimensional reduction.
The imported-cell universe remains the unchanged `set:parse_cytokines` set.

`processing/parse_cytokines/03_create_parse_cytokines.R` performs atomic
scWorkbench ingestion, compares source and imported identifiers, metadata, and
distributed count values, finalizes the bundle, and runs recipient validation.
The accepted bundle and run evidence are under
`intermediate_objects/parse_cytokines/bundle/` and
`intermediate_objects/parse_cytokines/reports/bundle_import_20260919T233232Z/`,
respectively.

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

The score commit is complete, but scWorkbench's enclosing JSON for its 1.398-GB
provenance row-order value exceeds R's 2^31-1-byte character-string limit
during a writable open or full recipient validation. Bounded checks passed and
cell memberships remain unchanged. Final feature-set registration is therefore deferred, and
the central accepted fingerprint remains the last fully validated Parse parent
state. No direct database rewrite or artifact deletion was used to work around
the limitation.

## Workflow documentation

Execution phases and validation boundaries are documented in:

- `processing/parse_cytokines/README.md`
- `processing/scmulti_flt3ko_il7ko/README.md`
- `processing/integration/README.md`

Source selections and operation parameters remain unset until each input has
been inspected and the applicable scWorkbench interfaces have been confirmed
inside the supplied image.
