# Step 4 ImmGen PeakVI UMAP

This derived development workflow fits `scwbrepresentation::umap_fit` on the
authenticated model-bound union of the 1,382-cell qualified Bdev-only 10D
PeakVI reference and all 28 projected ImmGen ATAC library profiles.

Run the numbered actions with
`config/04_reference_mapping/10_immgen_peakvi_joint_umap.yaml`:

1. `01_prepare_derived_bundle.py` verifies the completed mapping inputs and
   creates a copy-on-write derived bundle.
2. `02_fit_transform_and_plot.R` fits and promotes the joint UMAP model,
   exports coordinates, and draws population and native-L4 views.
3. `03_finalize_run.py` verifies required outputs and writes the checksummed
   run inventory.

The completed mapping bundle is read-only input. This workflow is a
development visualization, does not build a SIF, and does not unlock
`Bdev_implementation`.

Run `immgen_peakvi_umap_20261111T025344Z`, declared by configuration 07, is
retained as a failed diagnostic: the reference fit completed, but the query
transform correctly stopped at the output cell-order grant because the query
cell set was implicit. Configuration 08 declares the exact 28-library query
cell set and passes it to the transform contract.

Run `immgen_peakvi_umap_20261111T030153Z`, declared by configuration 08, is
also retained as a failed preflight diagnostic: it demonstrated that direct
cell-set registry mutation is prohibited on a committed bundle. Configuration
09 instead derives and promotes the exact query cell set through the
registered `cell_set_derive` operation before fitting or transforming UMAP.

That provenance-bound set passed, but the disjoint query transform still
stopped before publication because the installed writer grant is bound to the
1,382-cell fitted source-space axis. Configuration 10 therefore uses the
already authenticated 1,410-object model-bound latent union for a supported
joint UMAP fit. No core or operation-pack code is changed.

The joint-UMAP visualization was rejected for the requested trajectory view on
2026-11-11 because it does not use the established `umap.atac` coordinates.
Retain its diagnostics for traceability, but do not use them as the ImmGen
trajectory overlay. Step `05_plot_immgen_atac_umap` is the selected replacement.
