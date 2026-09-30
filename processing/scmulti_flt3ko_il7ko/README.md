# FLT3LKO/IL7KO multiome workflow

This workflow processes the independently delivered knockout multiome data.
`FLT3LKO` denotes knockout of the ligand gene *Flt3l*; it is not an *Flt3*
receptor knockout. The stable technical bundle key remains
`scmulti_flt3ko_il7ko` for compatibility.

Implemented unfiltered QC stage:

1. `01_calculate_unfiltered_qc.R`: validate identities and calculate RNA/ATAC
   metrics for every filtered-matrix cell;
2. `02_render_unfiltered_qc_report.R`: render the Markdown report;
3. `run_unfiltered_qc.sh`: run both stages in the supplied scWorkbench image.

No cell filtering or threshold selection occurs in this stage. Later phases
will create and validate the source bundle, apply only user-approved cell
selection, preprocess RNA/ATAC, and prepare the separate Bdev integration
bundle.

The approved non-destructive QC selection is implemented by:

4. `03_mark_qc_analysis_set.R`: retain every cell and mark the six individual
   criteria plus the combined `qc_use_downstream` flag;
5. `04_render_qc_selection_report.R`: report selection counts and criterion
   failures by library;
6. `run_qc_selection.sh`: reproduce the selection and report.

The source-bundle stage is implemented by:

7. `05_create_source_bundle.R`: import all 28,541 cells, raw `RNA_counts` and
   `ATAC_aggregate_peaks`, all QC metadata, and the accepted
   `set:qc_use_downstream` analysis universe through public scWorkbench APIs;
8. `run_create_source_bundle.sh`: run ingestion, exact source comparisons,
   finalization, and recipient validation in the supplied image.

The independently validated output is the accepted `scmulti_flt3ko_il7ko`
bundle. No normalization, dimensionality reduction, or mapping has been run.
It is the query bundle for later mapping to the Bdev reference.
