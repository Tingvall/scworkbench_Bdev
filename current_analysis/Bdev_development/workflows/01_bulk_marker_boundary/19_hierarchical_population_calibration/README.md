# Test hierarchical population-specific calibration

This diagnostic workflow tests the global-local population programs using a
broad-state gate, local population competition, and population-specific
thresholds derived from held-out ImmGen predictions. It then compares no
margin, fixed 0.10 margin, and calibrated policies on the bounded Bdev fixture.

The population thresholds are deliberately marked optimistic: they are chosen
and summarized on the same small held-out inventory. No label is accepted or
promoted.

```bash
python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/19_hierarchical_population_calibration/01_calibrate_reference.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/20_hierarchical_population_calibration.yaml

/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --bind /mnt/sata2:/mnt/sata2 \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif \
  Rscript current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/19_hierarchical_population_calibration/02_score_bdev_fixture.R \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/20_hierarchical_population_calibration.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/19_hierarchical_population_calibration/03_finalize_manifest.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/20_hierarchical_population_calibration.yaml
```
