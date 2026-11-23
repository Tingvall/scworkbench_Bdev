# Frozen PeakVI mappings on the fixed ATAC UMAP

This read-only development visualization consumes the completed step-11
four-model benchmark. It does not reopen or modify a PeakVI model, recompute
neighbors, tune features, or promote a candidate.

`01_plot_fixed_umap_comparison.R` places each ImmGen library at the
distance-kernel-weighted barycenter of its 30 already-recorded Bdev neighbors.
All four candidate panels use the same approved `umap.atac` coordinates and
the same coordinate limits. Weighted neighbor ellipses show how broadly each
bulk library is supported on the fixed reference; the detailed views also draw
spokes to the five nearest recorded neighbors.

The authenticated fitted Slingshot `Lineage4` curve is read from the approved
trajectory model and drawn in its native orientation. The separate L4 panel
shows the existing provisional post-freeze values, including explicit missing
calls. These values remain `diagnostic_unaccepted` and are not trajectory
assignments.

Run inside the pinned development SIF:

```bash
R_LIBS='<colon-separated runtime.r_libraries>' \
  apptainer exec <sif> Rscript 01_plot_fixed_umap_comparison.R \
  --config ../../../config/04_reference_mapping/19_frozen_peakvi_fixed_umap_comparison.yaml
python 02_finalize_run.py \
  --config ../../../config/04_reference_mapping/19_frozen_peakvi_fixed_umap_comparison.yaml
```
