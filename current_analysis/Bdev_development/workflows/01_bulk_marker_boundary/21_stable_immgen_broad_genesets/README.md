# Discover stable ImmGen broad-lineage programs

This exploratory step gives every ImmGen population equal weight, selects
positive markers that are consistently rank-visible within a lineage and low
in every competing lineage, and evaluates selection with nested
leave-one-population-out predictions. The final programs are then scored on the
bounded Bdev fixture without using Bdev in gene selection.

The Bdev reducer preserves borderline multi-state evidence and uses the fixed
default-compatible score and margin gates. All outputs are diagnostic.

```bash
python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/21_stable_immgen_broad_genesets/01_discover_and_validate.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/22_stable_immgen_broad_genesets.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/21_stable_immgen_broad_genesets/02_prepare_session_clone.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/22_stable_immgen_broad_genesets.yaml

/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --bind /mnt/sata2:/mnt/sata2 \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif \
  Rscript current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/21_stable_immgen_broad_genesets/03_score_bdev_fixture.R \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/22_stable_immgen_broad_genesets.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/21_stable_immgen_broad_genesets/04_finalize_manifest.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/22_stable_immgen_broad_genesets.yaml
```
