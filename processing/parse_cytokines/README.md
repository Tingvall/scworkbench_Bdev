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
8. `07_build_signature_states.R`, which learns 2,000 HVGs from control cells,
   builds a 30-PC representation for the control/single set, and promotes a
   fixed-resolution Parse-only Leiden state vector
9. `08_discover_specific_signatures.R`, the explicitly approved outside-
   scWorkbench replicate-aware step: exact raw-count pseudobulk aggregation,
   limma-voom response/specificity contrasts, and three-fold leave-one-
   replicate-out signature-size validation
10. `09_score_and_report_signatures.R`, which uses scWorkbench UCell on raw
    counts and checks score transfer at the sample-median level
11. `09_finalize_committed_ucell_report.R`, the bounded report-finalization
    path for the committed score artifact described below
12. `10_recover_parse_bundle.R`, which preserves the original bundle, repairs
    the engine runtime, registers the accepted sets in a copy, and completes
    recipient validation; see `recovery/README.md`

Run the three signature stages in order:

```bash
processing/parse_cytokines/run_07_build_signature_states.sh
processing/parse_cytokines/run_08_discover_specific_signatures.sh
processing/parse_cytokines/run_09_score_and_report_signatures.sh
```

The separate competitor-aware review is specified in
`COMPETITOR_AWARE_SIGNATURE_WORKFLOW.md`. It evaluates sizes 5, 10, 15, 20,
30, 40, and 50 without a gene-level FDR cutoff, selects size from full-data
target-to-strongest-competitor gaps, and reports selected-size cross-fitted
stability and core membership without promoting results. Run it with:

```bash
processing/parse_cytokines/run_14_competitor_aware_signature_review.sh
```

The installed scWorkbench build committed all 31 UCell score keys, their
summary, and both QC tables in commit
`commit_20260921_155948_s32blwqt`. Its generic artifact serializer then failed
while checkpointing the durable-session receipt because the promoted family
contains a 1.398-GB provenance row-order value; serializing the enclosing
manifest/receipt exceeds R's 2^31-1-byte character-string limit. The score
payload itself is complete in `artifact_pseudotime`: every key has exactly
154,209 distinct downstream cells.

The completed report was therefore finalized read-only, one signature at a
time, with:

```bash
processing/parse_cytokines/run_09_finalize_committed_ucell_report.sh
```

It accepted 23 signatures (9 induced and 14 repressed) spanning 14 of the 26
single cytokines. Eight discovery candidates failed UCell transfer and 21
cytokine/direction pairs had not advanced from discovery. The report and final
TSV/GMT exports are under
`intermediate_objects/parse_cytokines/reports/signature_report_20260921T221910Z/`.
All imported and downstream cell-set memberships were rechecked unchanged;
SQLite `quick_check` and recipient-readiness checks passed.

That original build could not complete a writable open or full recipient
validation with the oversized provenance record. The stage-10 recovery has now
resolved this in a separate scWorkbench `0.1.0.9004` runtime built from the
image's source. It retains the original identity scheme and score payload,
streams oversized registry JSON for hashing, accelerates equivalent character
row-key serialization, preserves imported read-only flags, and avoids sorting
the large provenance field in SQLite. The unchanged TSV was round-tripped to
CSV for the documented feature-set table parser.

All 23 accepted sets are registered read-only in `parse_cytokine_signatures`.
Full recipient validation passed, including the existing 31-signature score
family's content identity. Five missing materialization files from the original
interrupted commit were restored from its saved session, each with exactly the
raw-file SHA-256 recorded in the registry. The original bundle's 41 files are byte-for-byte
unchanged. The accepted copy and its runtime are recorded in `config/bundles.yml`;
processing and figure runners resolve them automatically. Evidence is under
`intermediate_objects/parse_cytokines/reports/signature_recovery_20260923T200125Z/finalization/`.
The earlier bounded report remains historical evidence. No direct database
rewrite or artifact deletion was used for recovery.

The checked-in stage-09 implementation now keeps per-cell score families in
the durable session and promotes only compact summaries/QC tables before the
feature-set import. This prevents a clean rerun from recreating the oversized
row-order provenance artifact. Stage 09 now supplies CSV to the table parser
and checks every staged set before promotion.

The exclusion audit and all-candidate visualization were added with the
read-only figure command:

```bash
figure_scripts/parse_cytokines/run_03_signature_sample_treatment_violins.sh
```

Its 31-page PDF contains every scored discovery gene set, including the eight
that failed UCell transfer. For each treatment it draws three separate
cell-level violins, one per biological sample, and marks each sample median.
The companion exclusion table accounts for all 52 cytokine/direction pairs:
21 lacked the minimum five discovery-eligible genes, eight failed the
sample-median UCell thresholds, and 23 were accepted. Directions without a
gene set are listed rather than assigned an undefined score or plot.

Signature discovery is restricted to `set:parse_qc_control_single`. Combined
stimulation cells remain retained but are not used. Bdev is not read by these
stages. The primary statistical unit is the biological sample: each cytokine
is tested both against control and against the equal-weight mean of the other
25 single-cytokine conditions. Induced and repressed programs are separate.
Unsupported directions are reported rather than forced into the final library.

The planned outside-scWorkbench boundary was limited to replicate-aware
statistical work that the installed build does not provide: BPCells exact
sample/state count aggregation and limma-voom modelling. State construction,
per-cell UCell scoring, and score promotion ran in scWorkbench. UCell uses raw
`RNA_counts`, does not graph-smooth scores, and treats sample medians as the
primary confirmation unit. Raw scores are compared within a signature, not
directly across signatures. The read-only finalizer queries scWorkbench's
durable score payload table in bounded per-signature slices, as required when
the original build encountered the provenance limitation.

The accepted independently validated output is the `parse_cytokines` bundle at
`intermediate_objects/parse_cytokines/bundle_recovered_20260923T200125Z/`.
Its required runtime is registered alongside it. Reviewed cytokine gene sets
may later be consumed by `processing/integration/`.

## Separate rest-only cytokine collection

The requested omission of treatment-versus-control selection is configured in
`config/parse_cytokines_rest_only.yml`. It uses a separate working review bundle,
`parse_cytokines_rest_only`, with `parse_cytokines` as its accepted parent.
The original parent bundle and its 23 registered signatures remain unchanged.

```bash
bash processing/parse_cytokines/run_rest_only.sh discover
bash processing/parse_cytokines/run_rest_only.sh prepare
bash processing/parse_cytokines/run_rest_only.sh score
bash processing/parse_cytokines/run_rest_only.sh report
```

Stage 11 reselects membership and candidate size from the saved specificity
fits, using only target-versus-other-single-stimuli effects, FDR, fold signs,
specificity-only ranking/state concordance, and rest AUC. The control AUC,
combined AUC and control-based dominance term are absent. The saved models
still estimate variance using the original 81 samples, including control;
they are not refitted after deleting the control samples. Existing Parse
states also retain their original control-derived HVG construction.

Stage 12 imports the complete 154,209-cell analysis population and all 78,298
RNA features into a separate bundle through scWorkbench. Its ingest adapter
requires a sparse `Matrix` payload rather than the backed BPCells object;
counts are unchanged. Ordered identifiers, metadata, all cell totals and every
gene-by-sample pseudobulk count are verified against the accepted source.

Stage 13 scores the new memberships with native scWorkbench UCell and tests
sample-median discrimination among the 78 single-stimulus samples only.
Passing sets are imported into `parse_cytokine_rest_only` in the private durable
session overlay. No operation outputs are promoted: the user's report-before-
promotion requirement takes precedence over the generic engine recommendation
to promote before reporting. The working input bundle contains no published
rest-only namespace. The durable workspace and staged bundle are recorded in
`scoring.json`; terminal score primitives are also retained separately as RDS.

Reports and exports are under
`intermediate_objects/parse_cytokines_rest_only/reports/rest_only_20260924/`.
The runner mounts the accepted parent and repository read-only and permits
writes only under the separate derived storage root. Scoring resumes completed
durable checkpoints without repeating the expensive UCell operation.

The completed review run has 51 directions with sufficient eligible genes,
47 passing discovery, and 37 passing UCell confirmation (13 up, 24 down;
24 stimuli). The final collection has 1,140 memberships covering 756 distinct
genes. Compared with the original accepted directions, 16 additional directions
pass and two fail after reselection: c-kit up and TNF-alpha up each select five
genes and have sample-median UCell rest AUC 0.50. All 47 AUCs were independently
recomputed from the exported sample medians; staged memberships match the GMT
and TSV exports exactly. No new result has been promoted.
