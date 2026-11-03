# Qualify hierarchical ImmGen labeling on the bounded Bdev fixture

This development prototype uses nine globally contrasted broad-lineage panels
for a multi-state first gate, then restricts population competition to the
eligible broad states. It compares the pinned 64-set collection with a variant
that removes genes shared across at least three broad lineages. The parent
step-15 bundle is copied to a session-only clone because scWorkbench requires a
writable bundle handle for durable sessions. Native operations then chain
through session-local outputs; no imports, registrations, or artifacts are
promoted, and the clone registry fingerprint must remain unchanged.

From the repository root:

```bash
python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/16_qualify_immgen_bdev_hierarchy/01_prepare_collections.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/17_immgen_bdev_hierarchical_labeling.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/16_qualify_immgen_bdev_hierarchy/02_prepare_session_clone.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/17_immgen_bdev_hierarchical_labeling.yaml

/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --bind /mnt/sata2:/mnt/sata2 \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif \
  Rscript current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/16_qualify_immgen_bdev_hierarchy/03_score_hierarchical_fixture.R \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/17_immgen_bdev_hierarchical_labeling.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/16_qualify_immgen_bdev_hierarchy/04_finalize_manifest.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/17_immgen_bdev_hierarchical_labeling.yaml
```

Outputs are diagnostic development evidence only. The workflow does not
register feature sets, promote scWorkbench artifacts, register winner cell
sets, emit accepted labels, or unlock implementation.
