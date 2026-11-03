# Qualify 64-population ImmGen labeling on a bounded Bdev fixture

This development acceptance workflow uses the complete 64-population working
default and the native scWorkbench signal-family labeling chain. It never
modifies the archived Bdev bundle. A copied execution clone is migrated in
active storage, and only a deterministic 1,746-cell fixture is scored.

From the repository root:

```bash
python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/15_qualify_immgen_bdev_labeling/01_freeze_fixture.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/16_immgen_bdev_labeling_acceptance.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/15_qualify_immgen_bdev_labeling/02_prepare_execution_clone.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/16_immgen_bdev_labeling_acceptance.yaml

/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --bind /mnt/sata2:/mnt/sata2 \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif \
  Rscript current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/15_qualify_immgen_bdev_labeling/03_score_and_label_fixture.R \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/16_immgen_bdev_labeling_acceptance.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/15_qualify_immgen_bdev_labeling/04_finalize_manifest.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/16_immgen_bdev_labeling_acceptance.yaml
```

Generated tables, figures, reports, logs, and manifests are written to paired
development storage. The execution clone and durable sessions are stored under
the mirrored development intermediate-object path.
