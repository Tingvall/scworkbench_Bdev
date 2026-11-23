# Provisional Bdev development default: scwbATAC depth-adjusted HVP 20k

**Selected:** 2026-11-17  
**Status:** superseded as the model default on 2026-11-22 by the qualified 20D
reference and gate; retained as feature-selection and sensitivity evidence

The depth-adjusted HVP20k feature-selection decision remains active. The frozen
10-dimensional model described below is no longer the later-implementation
model default; that role is now assigned to the exact qualified 20-dimensional
model and limits recorded in
`19_qualified_20d_reference_mapping_default.md`.

The current working default for bulk-to-Bdev reference mapping is the frozen
`hvp_depth20k` model from run
`hvp20k_depth_comparison_20261117T044105Z`. Its 20,000 peaks were selected from
all 191,433 peaks and all 26,723 Bdev cells by the official ATAC-specific
`scwbatac::hvp_select` operation, using residual binary accessibility variance
after adjustment for `log1p(detected peaks per cell)`. PeakVI used the frozen
source feature order, ten latent dimensions, and no ImmGen information during
selection or fitting.

## Why this is the default

The user explicitly selected the official depth-adjusted HVP model after
reviewing its same-query comparison with raw HVP and the earlier `overall20k`
model. Depth-adjusted HVP was better than raw HVP on five of eight bulk-mapping
diagnostics and had lower peak-saturation association than `overall20k`. The
choice also establishes the reusable ATAC-specific operation as the intended
method path.

This was not an automatic best-score selection. In the same external evidence,
`overall20k` retained better population separation, replicate agreement,
support distance, neighbor-sample entropy, and depth association. It therefore
remains the required sensitivity comparator. Model discordance is reported as
feature-selection uncertainty; coordinates are not averaged.

## Working policy

Use `hvp_depth20k` for the first-stage whole-reference Bdev neighborhood and
fixed-UMAP mapping during development. Trajectory interpretation remains a
separate second stage: use the authenticated fitted Slingshot L4 trajectory in
native orientation only when the whole-reference neighborhood supports L4
membership. Current bulk L4 outputs remain diagnostic and unaccepted, and
non-L4 mappings are not forced onto the trajectory.

## Qualification boundary

The frozen model and its external projection are valid development evidence,
but the operation is not yet eligible for a replacement SIF. The real Bdev
matrix exposed a repeated 2.9 GiB sparse-to-dense coercion for each 2,048-cell
chunk. A sparse-safe source fix, a raw-mode bypass, and a large sparse-axis
regression contract are now present in the `scwbATAC` working tree. They still
require an immutable commit and repeat Bdev qualification. The post-fix
`scwbATAC` commit—not the unfixed `e9fd2da` base—must be installed in the final
implementation SIF.

This selection does not promote the development model into
`Bdev_implementation`, authorize paper analysis, accept trajectory assignments,
or waive the requirement for a newly built, qualified, checksummed, and
explicitly accepted immutable SIF. Exact artifacts, checksums, metrics,
sensitivity policy, and boundaries are pinned in
`config/04_reference_mapping/25_scwbatac_hvp_depth20k_default.yaml`.
