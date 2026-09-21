# Parse cytokine workflow

The first step prepares the user-selected unfiltered Parse matrix outside the
bundle. It uses `Seurat::ReadParseBio()` in the supplied image, creates a Seurat
object without cell or feature filtering, preserves the delivered metadata,
and adds explicit `sample_id`, `condition`, and `replicate` fields. The output
directory contains the checksummed RDS, input inventory, sample summary,
environment manifest, and validation result.

The delivered metadata files are gzip-compressed, whereas Seurat 5.5.1 expects
the literal filenames `cell_metadata.csv` and `all_genes.csv`. The script
temporarily decompresses those two files and links the compressed count matrix
in a staging directory. It removes that staging directory immediately after
`ReadParseBio()` returns; the original delivery remains read-only.

Run it from the repository root:

```bash
processing/parse_cytokines/run_01_readparsebio_unfiltered.sh
```

The default output is outside the bundle at
`intermediate_objects/parse_cytokines/prepared/`. The
wrapper mounts the Parse source read-only and runs R from the supplied image,
which currently provides Seurat 5.5.1, SeuratObject 5.4.0, Matrix 1.7-4, and
data.table 1.18.4. A different input or empty output directory can be supplied
as the first and second positional arguments.

The unfiltered matrix declares 13,413,901 cells, 78,298 features, and
1,422,994,774 stored entries. The wrapper defaults to four computational
threads and a 160 GiB virtual-memory ceiling. On this 503 GiB shared server that
caps the process below one third of RAM. The wrapper refuses any requested
ceiling above 40% of installed server RAM. Override the lower limits with
`PARSE_THREADS` and `PARSE_MAX_MEMORY_GIB` if needed; the 40% guard remains.

The completed outside-bundle preparation is recorded in
`config/parse_cytokines.yml`. It retained all 13,413,901 cells and 78,298
features across 96 samples and passed the source, dimension, identifier-order,
metadata, and no-filtering checks. The checksummed RDS and validation evidence
remain in the external work directory named above.

The second step applies the user-selected prefilter in Seurat:

```bash
processing/parse_cytokines/run_02_prefilter_min500counts.sh
```

It calculates each cell's total directly from the raw RNA `counts` layer and
retains cells whose total is at least 500. It applies no sample, feature, or
additional cell filter. The completed output is:

`intermediate_objects/parse_cytokines/prepared/parse_cytokines_prefiltered_min500counts_seurat.rds`

The filter retained 251,599 of 13,413,901 cells, all 78,298 features, and all
96 samples. The minimum retained raw RNA count is exactly 500. The script
verified that `nCount_RNA` exactly matched the raw counts-layer column sums,
retained cell order is the input subset order, feature order is unchanged,
cell and feature metadata are unchanged, every retained cell meets the
threshold, and every removed cell is below it. The output SHA-256 is
`1e4da7149190e6c93f71bc28fc47d1f58a20e22ee9cbec6894f0824394e9bc40`.
The unfiltered source remained unchanged.

Validation, per-sample counts, environment versions, and session information
are written beside the filtered RDS with the same filename stem. The previously
reviewed single/combination/control classification was not used as a filter in
this count-only prefilter.

The third step creates the separate `parse_cytokines` scWorkbench bundle from
that accepted pre-filtered object:

```bash
processing/parse_cytokines/run_03_create_parse_cytokines.sh
```

It imports only the raw RNA `counts` layer as `RNA_counts` and the required
`sample_id`, `replicate`, and `condition` cell metadata. It does not import a
normalized layer, reductions, or unrelated Seurat metadata. The runner mounts
the input RDS read-only, refuses to overwrite an existing bundle or report
directory, and writes run-specific source comparison and recipient-validation
evidence under `reports/bundle_import_<UTC timestamp>/`.

The import verifies the accepted RDS checksum, exact cell and feature order,
exact required metadata values, all 96 sample mappings, and distributed raw
count values across every row and column of the source matrix. It then
finalizes the bundle and runs the scWorkbench recipient validation gate.

The fourth step calculates factual RNA QC metrics through the registered
scWorkbench `rna_qc_metrics` operation:

```bash
processing/parse_cytokines/run_04_calculate_rna_qc.sh
```

It uses every cell in the imported `set:parse_cytokines` universe and promotes
total RNA counts, detected features, and explicitly matched mitochondrial and
ribosomal-prefix counts and percentages. Feature matching is case-sensitive:
`^mt-` and `^Rp[sl]` in the mouse-symbol namespace. The exact membership is an
audit artifact. The ribosomal-prefix metric is not asserted to be a curated
structural-ribosome set.

The runner uses 400-cell matrix chunks because the current BPCells getter's
dense-fallback guard limits any one materialized slice to 256 MiB. This changes
only bounded-memory batching, not metric values or the analyzed cell universe.

This step does not select a threshold, create an exclusion set, or derive an
accepted analysis population. The imported cell set remains unchanged. Render
the distributions and neutral empirical-tail sensitivity previews with:

```bash
figure_scripts/parse_cytokines/run_02_rna_qc_report.sh
```

The user selected the following inclusive acceptance criteria after reviewing
the distributions:

- total RNA counts >= 1,000;
- detected RNA features >= 500;
- mitochondrial percentage <= 10%.

Apply them with:

```bash
processing/parse_cytokines/run_05_apply_rna_qc_filter.sh
```

The fifth step writes the combined proposed-exclusion list, criterion overlap,
and per-sample/per-condition effects before it promotes any QC-derived cell
set. It then records `set:parse_qc_excluded` as the accepted exclusion set and
`set:parse_qc_accepted` as the accepted analysis universe. The two sets are
validated as a disjoint, exhaustive partition of the imported cells. It never
replaces or narrows the imported `set:parse_cytokines` membership or the raw
`RNA_counts` matrix.

The completed filter retained 194,215 analysis cells (77.1923%) and marked
57,384 cells (22.8077%) as excluded. All 96 samples retain accepted cells. The
pre-promotion preview, complete compressed exclusion list, overlap tables,
cell-set recipes, promotion records, and recipient validation are under:

`intermediate_objects/parse_cytokines/reports/rna_qc_filter_20260921T035031Z/`

The sixth step applies the reviewed treatment classification:

```bash
processing/parse_cytokines/run_06_classify_treatments_and_select.sh
```

The explicit mapping is `config/parse_cytokines_treatment_categories.tsv`.
It defines `control` as zero cytokine stimuli, `single` as exactly one, and
`combined` as more than one. The map covers all 32 conditions and 96 samples;
`Csf_1` and `Csf_2` remain reviewed single-stimulus labels despite their
underscores.

Before bundle mutation, the step writes the complete condition and sample
classification, per-category counts, and all proposed treatment-scope
exclusions. It then imports `treatment_category` and `cytokine_count` for all
251,599 imported cells, promotes the categorical `parse_treatment_category`
cell vector, and derives `set:parse_qc_control_single` as the exact intersection
of `set:parse_qc_accepted` with the `control` and `single` categories.

The accepted downstream set contains 154,209 cells across 81 samples and 27
conditions: 4,816 control cells and 149,393 single-cytokine cells. The 40,006
QC-accepted combined-stimulation cells are outside that downstream set. They
remain in the bundle, as do all QC-excluded cells; `set:parse_cytokines`,
`set:parse_qc_accepted`, and `set:parse_qc_excluded` retain their exact prior
memberships. Classification, provenance, cell-set recipes, and recipient
validation are under:

`intermediate_objects/parse_cytokines/reports/treatment_selection_20260921T145825Z/`

Workflow phases are:

1. `01_readparsebio_unfiltered.R`
2. `02_prefilter_min500counts.R`
3. `03_create_parse_cytokines.R`
4. `04_calculate_rna_qc.R`
5. QC distribution review and explicit user threshold selection
6. `05_apply_rna_qc_filter.R`, which writes the exclusion preview before
   deriving the durable exclusion and accepted-analysis sets
7. `06_classify_treatments_and_select.R`, which classifies every sample/cell,
   writes the treatment-selection preview, and derives the QC-accepted
   control/single downstream analysis set

The accepted independently validated output is the `parse_cytokines` bundle at
`intermediate_objects/parse_cytokines/bundle/`. Reviewed cytokine gene sets may
later be consumed by `processing/integration/`.
