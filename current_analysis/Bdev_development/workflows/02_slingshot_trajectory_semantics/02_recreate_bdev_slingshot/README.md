# 02 Recreate Bdev Slingshot trajectory

This bounded development workflow recreates the published Bdev trajectory on
the authenticated `umap.atac` space with cluster 5 as the root, clusters 9, 15,
18, 19, 22, and 23 excluded, `stretch = 0`, `omega = TRUE`, and
`omega_scale = 5`. It uses `scwbTrajectory` 0.1.3 with the development core.

The fitted lineage family is the primary output. The preserved Bdev `pt` field
is read without modification and used only to compare the published `-L4`
axis. The fit remains in its durable scWorkbench session and is not promoted.

Run from the repository root:

```bash
python3 current_analysis/Bdev_development/workflows/02_slingshot_trajectory_semantics/02_recreate_bdev_slingshot/01_prepare_execution_clone.py \
  --config current_analysis/Bdev_development/config/02_slingshot_trajectory_semantics/02_bdev_trajectory_recreation.yaml

/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --bind /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif \
  env R_LIBS_USER=/home/jtg/.local/lib/R/scworkbench-development \
  Rscript current_analysis/Bdev_development/workflows/02_slingshot_trajectory_semantics/02_recreate_bdev_slingshot/02_fit_and_compare.R \
  --config current_analysis/Bdev_development/config/02_slingshot_trajectory_semantics/02_bdev_trajectory_recreation.yaml

python3 current_analysis/Bdev_development/workflows/02_slingshot_trajectory_semantics/02_recreate_bdev_slingshot/03_finalize_manifest.py \
  --config current_analysis/Bdev_development/config/02_slingshot_trajectory_semantics/02_bdev_trajectory_recreation.yaml
```

Generated state and results mirror this workflow beneath paired active
storage. The archived source bundle remains read-only.
