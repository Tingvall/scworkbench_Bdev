# Test lineage-anchored narrow ImmGen programs

This exploration gives all 64 ImmGen populations an anchor module specific
against incompatible lineages and a local module specific against neighbouring
populations in the same lineage. No broad label is assigned or used as a gate.
The narrow composite penalizes the strongest incompatible population anchor.

Discovery, module weighting, and target-membership holdout evaluation use only
ImmGen. The frozen modules and rule are then scored on the bounded Bdev fixture.

```bash
python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/22_lineage_anchored_narrow_genesets/01_discover_and_validate.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/23_lineage_anchored_narrow_genesets.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/22_lineage_anchored_narrow_genesets/02_prepare_session_clone.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/23_lineage_anchored_narrow_genesets.yaml

/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --bind /mnt/sata2:/mnt/sata2 \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif \
  Rscript current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/22_lineage_anchored_narrow_genesets/03_score_bdev_fixture.R \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/23_lineage_anchored_narrow_genesets.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/22_lineage_anchored_narrow_genesets/04_finalize_manifest.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/23_lineage_anchored_narrow_genesets.yaml
```
