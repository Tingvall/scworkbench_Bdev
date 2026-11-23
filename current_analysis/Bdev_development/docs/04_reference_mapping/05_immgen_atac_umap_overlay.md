# ImmGen overlay on the approved Bdev ATAC UMAP

Run `immgen_atac_umap_20261111T032028Z` completed a Step 4 diagnostic
visualization. It uses the established `umap.atac` coordinates from the Bdev
bundle used for fitted Slingshot L4; no UMAP was refit and no source coordinate
was changed.

Each of the 28 ImmGen ATAC libraries was placed at the kernel-weighted
barycenter of its primary 30 authenticated PeakVI neighbors. The weights use
`exp(-(distance / kth_distance)^2)`, matching the reference-mapping kernel.
This is a neighbor-based placement on an existing visualization, not a UMAP
transform. The figure labels every point by population and GSM, uses circles
for 25 supported mappings, and uses orange triangles for the 3 borderline
mappings.

The figure retains all 26,723 Bdev cells as context and colors the 13,451
authenticated L4-family cells by native fitted L4. The original mapping and
Bdev bundle manifest checksums remained unchanged. The run is a development
visualization only: it did not build a replacement SIF, promote paper results,
or unlock `Bdev_implementation`.

## Scientific status

This visualization was not accepted on 2026-11-11. The placements do not make
biological sense, especially for later-stage populations. A technically valid
neighbor-weighted barycenter is therefore not sufficient evidence that the
ImmGen libraries have been mapped correctly onto the trajectory UMAP. Treat
the run as a retained diagnostic, not a completed Step 4 result. Do not use its
coordinates or visual interpretation downstream; the mapping/visualization
problem remains unresolved and will be revisited.

The reproducible definition is
`config/04_reference_mapping/11_immgen_atac_umap_overlay.yaml`, and the numbered
actions are in
`workflows/04_reference_mapping/05_plot_immgen_atac_umap/`.
