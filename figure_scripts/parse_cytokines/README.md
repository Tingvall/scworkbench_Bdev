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
