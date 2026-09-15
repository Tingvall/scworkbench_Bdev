# Parse cytokine reports

`01_report_prepared_seurat.py` summarizes the validated unfiltered Seurat
preparation and proposes `single`, `combination`, and `control` stimulus
classes without loading the large RDS:

```bash
python3 figure_scripts/parse_cytokines/01_report_prepared_seurat.py \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines/prepared \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines/reports/object_inventory_20260915
```

The generated sample and condition TSV files are proposed annotations. Review
the noted source-label ambiguities before using them as canonical metadata or
importing them into a scWorkbench bundle.
