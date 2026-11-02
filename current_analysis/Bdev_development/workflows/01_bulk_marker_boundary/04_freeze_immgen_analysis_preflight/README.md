# Freeze the ImmGen analysis preflight

This step freezes and validates the exact count matrix, complete GENCODE vM25
annotation, gene-symbol mapping policy, state-gene exclusions, and immutable
outside-analysis software environment. It creates complete annotation and
state-flag tables in paired storage.

It does not normalize counts, fit an ImmGen model, select marker genes, invoke
scWorkbench, or create a scWorkbench bundle.

```bash
python current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/04_freeze_immgen_analysis_preflight/01_freeze_immgen_analysis_preflight.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/03_immgen_analysis_preflight.yaml \
  --software-lock current_analysis/Bdev_development/config/01_bulk_marker_boundary/04_immgen_bulk_analysis_software_lock.yaml
```

Generated output is written under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/04_freeze_immgen_analysis_preflight/
```

The validator uses Python 3, PyYAML, and only standard-library modules for GTF,
CSV, tar, and SQLite processing. The edgeR environment is separately smoke
tested as recorded in the software lock.
