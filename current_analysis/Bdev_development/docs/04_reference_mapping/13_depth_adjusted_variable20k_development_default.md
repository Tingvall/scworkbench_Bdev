# Provisional Bdev development default: depth-adjusted variable 20k PeakVI

**Selected:** 2026-11-16  
**Status:** superseded on 2026-11-17 by the official scwbATAC depth-adjusted
HVP 20k development default; never promoted to implementation

This document preserves the historical `overall20k` decision. The active
declaration is now config 25 and development-default document 16.

The working default selected at that time for bulk-to-Bdev reference mapping was the frozen
`overall20k` candidate from run
`full_bdev_peakvi_immgen_4way_20261116T202157Z`. It uses all 26,723 Bdev cells
and the 20,000 depth-adjusted overall-variable peaks in their frozen source
order. The user selected it after review of the four-model external ImmGen
benchmark and the common-coordinate fixed-ATAC-UMAP comparison.

## Selection evidence

Among the four frozen candidates, overall20k had the strongest external bulk
geometry: a 0.176 within-population to between-population latent-distance
ratio, perfect replicate cluster agreement, 477 distinct reference neighbors,
a 0.0107 maximum neighbor-hub fraction, and a median normalized reference
distance of 1.565. The fixed-UMAP comparison showed ordered population mapping
and much tighter support than the all-peak model.

The principal limitation is its 0.736 maximum absolute association between a
latent coordinate and selected-peak saturation. That limitation remains
visible rather than being folded into a post hoc composite score.

## Working policy

Overall20k is the primary model for whole-reference Bdev neighborhood and
fixed-UMAP mapping. Cluster20k remains the required sensitivity model because
it produced nearly identical normalized support distance and perfect replicate
agreement with less saturation association. Concordance supports the broad
mapping; discordance is reported as uncertainty. Coordinates from the two
models are not averaged.

Trajectory interpretation remains a separate second stage. The authenticated
fitted Slingshot L4 trajectory is used in native orientation only when the
whole-reference neighborhood supports L4 membership. Current bulk L4 outputs
remain diagnostic and unaccepted, and non-L4 mappings are not forced onto the
trajectory.

This selection applies only to ongoing `Bdev_development`. It does not promote
the existing development model artifact into `Bdev_implementation`, authorize
paper analysis, or replace the requirement for a newly built, qualified,
checksummed, and explicitly accepted immutable SIF. Exact artifacts, checksums,
metrics, sensitivity policy, and boundaries are pinned in
`config/04_reference_mapping/20_depth_adjusted_variable20k_default.yaml`.
