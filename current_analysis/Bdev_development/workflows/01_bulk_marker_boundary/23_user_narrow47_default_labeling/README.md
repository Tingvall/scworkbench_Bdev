# User-supplied narrow47 labeling with the default scWorkbench chain

This diagnostic development workflow applies the 47 user-supplied, 20-gene
population sets to the frozen 1,746-cell Bdev fixture. It uses native
scWorkbench feature-set import and overlap QC, UCell scoring, score-family QC,
percentile calibration, and argmax reduction with the default 0.80 score and
0.10 margin thresholds. No broad or lineage gate is used.

Only the feature-set import is promoted, and only into the isolated execution
clone because the native overlap-QC operation requires a registered namespace.
Score and label artifacts are not promoted, and the source bundle is unchanged.

Run from the repository root with the pinned SIF:

```bash
SIF=/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif
CONFIG=current_analysis/Bdev_development/config/01_bulk_marker_boundary/24_user_narrow47_default_labeling.yaml
/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec --bind /mnt/sata2:/mnt/sata2 "$SIF" python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/23_user_narrow47_default_labeling/01_prepare_inputs_and_clone.py --config "$CONFIG"
/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec --bind /mnt/sata2:/mnt/sata2 "$SIF" Rscript current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/23_user_narrow47_default_labeling/02_score_with_scworkbench.R --config "$CONFIG"
/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec --bind /mnt/sata2:/mnt/sata2 "$SIF" python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/23_user_narrow47_default_labeling/03_finalize_manifest.py --config "$CONFIG"
```
