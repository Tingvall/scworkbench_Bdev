# Parse cytokine reports

`01_report_prepared_seurat.py` summarizes the validated unfiltered Seurat
preparation and supplied the initial `single`, `combination`, and `control`
classification review without loading the large RDS:

```bash
python3 figure_scripts/parse_cytokines/01_report_prepared_seurat.py \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines/prepared \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines/reports/object_inventory_20260915
```

The generated sample and condition TSV files are historical review inputs. The
accepted canonical map now uses the category name `combined` and is recorded at
`config/parse_cytokines_treatment_categories.tsv`; it was applied by
`processing/parse_cytokines/06_classify_treatments_and_select.R`.

After `processing/parse_cytokines/run_04_calculate_rna_qc.sh` has promoted the
factual scWorkbench RNA-QC artifacts, render their distributions with:

```bash
figure_scripts/parse_cytokines/run_02_rna_qc_report.sh
```

The report covers all imported cells, includes overall and per-sample
distributions, and supplies empirical-tail sensitivity previews. Those
previews are not threshold recommendations: the renderer creates no exclusion
or accepted analysis cell set.

After reviewing that report, the user selected total counts >=1,000, detected
features >=500, and mitochondrial percentage <=10%. The combined exclusion
preview was written before the durable sets were derived and is available at:

`intermediate_objects/parse_cytokines/reports/rna_qc_filter_20260921T035031Z/exclusion_preview.md`

The distribution report still represents all imported cells. Filtering did
not replace those figures or remove imported cells. `set:parse_qc_accepted`
retains all QC-accepted treatment categories; the narrower downstream analysis
universe is `set:parse_qc_control_single`, containing only QC-accepted control
and single-cytokine cells. Treatment-selection evidence is under
`intermediate_objects/parse_cytokines/reports/treatment_selection_20260921T145825Z/`.

The individual-cytokine signature workflow generates its report directly from
the promoted scWorkbench UCell score family:

```bash
processing/parse_cytokines/run_07_build_signature_states.sh
processing/parse_cytokines/run_08_discover_specific_signatures.sh
processing/parse_cytokines/run_09_score_and_report_signatures.sh
```

For the current already-committed 31-key UCell family, the report was completed
with the bounded read-only finalizer:

```bash
processing/parse_cytokines/run_09_finalize_committed_ucell_report.sh
```

The final run directory contains `REPORT.md`, an oriented condition-by-
signature specificity heatmap, multi-page sample-median UCell distributions,
exact sample medians, descriptive cell-level quantiles, all candidate/final
statuses, the final feature-set table and GMT, and the committed UCell
promotion audit. The completed report is:

`intermediate_objects/parse_cytokines/reports/signature_report_20260921T221910Z/REPORT.md`

It contains 23 accepted signatures across 14 cytokines. The report
distinguishes the unchanged imported set from the control/single analysis
population and records that neither combined-stimulation cells nor Bdev were
used for signature discovery. `validation.json` records bounded payload,
membership, SQLite, and recipient-readiness checks. Full recipient validation
and bundle feature-set registration remain deferred because the committed
31-key family has a 1.398-GB provenance value whose enclosing JSON exceeds R's
2^31-1-byte character-string limit; no artifact or registry row was rewritten
directly.

To explain every exclusion and render all scored gene sets—not only the 23
accepted sets—as cell-level, per-biological-sample violins grouped by
treatment, run:

```bash
figure_scripts/parse_cytokines/run_03_signature_sample_treatment_violins.sh
```

This read-only figure step queries each committed UCell key in a bounded slice.
It writes a 31-page PDF and 31 PNGs (23 accepted and 8 UCell-transfer failures),
plus a manifest and a complete table covering all 52 cytokine/direction
possibilities, to the same signature report directory. The 21 directions with
fewer than five discovery-eligible genes are documented but cannot be plotted
because no gene set or UCell score exists. The explanatory supplement is:

`intermediate_objects/parse_cytokines/reports/signature_report_20260921T221910Z/SIGNATURE_EXCLUSIONS_AND_VIOLINS.md`
