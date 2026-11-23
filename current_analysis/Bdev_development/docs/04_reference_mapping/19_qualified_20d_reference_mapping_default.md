# Qualified 20D reference-mapping default

**Selected:** 2026-11-22  
**Status:** selected for later Bdev implementation; activation awaits the new
qualified immutable SIF

The qualified 20-dimensional depth-adjusted HVP20k PeakVI reference is now the
default reference-mapping model for the later Bdev implementation. This
supersedes the earlier 10-dimensional frozen HVP20k model as the working model
default while preserving the earlier run as historical and sensitivity
evidence.

The selected frozen artifact is
`model_ref_20261122_194532_hka11fvs` in qualification bundle manifest
`7f6e5fb2968ec86c69513be3ebd382884a23815b787215ad6326fa7fe6eea506`.
It uses 20,000 depth-adjusted HVP peaks, 20 PeakVI latent dimensions, and the
authenticated fitted Slingshot L4 trajectory in native orientation.

## Reference and pseudotime limits

The implementation must preserve the exact qualified gate:

- primary neighborhood size: 30;
- supported when the calibrated distance statistic is at most
  `2.02620492264996`;
- borderline when it is above the supported limit and at most
  `2.56846898217638`;
- out of reference above `2.56846898217638`;
- local trajectory annotation must be `in_domain`;
- annotation confidence must be at least `0.51`; and
- the accepted native-L4 estimate must be finite.

Raw native-L4 estimates and explicit rejection reasons remain required for
every query. An accepted value is published only when the complete
support-and-annotation gate passes.

## Bulk mapping

The same exact frozen 20D model and calibration are also the default for bulk
ATAC mapping. Bulk libraries must be projected without refitting on the exact
HVP-depth20k feature order and evaluated with the same support and native-L4
limits. A separate bulk-specific PeakVI reference fit or a relaxed bulk gate is
not the default. Fixed-UMAP display uses the primary-30 neighbor
distance-kernel barycenter on the unchanged Bdev `umap.atac` coordinates.

The previous 28-library frozen-model benchmark provides a timing precedent:
projection took about 86 seconds per model and fixed-UMAP plotting took about
37 seconds. For the already fitted 20D model, expected execution is roughly
10–20 minutes including cloning, input validation, projection, gated mapping,
plotting, and finalization. Building and validating a new traceable workflow is
expected to take about 30–45 minutes end to end.

Development run `immgen_bulk_gated_20d_20261123T001821Z` subsequently applied
this exact policy to all 28 ImmGen ATAC libraries. All 28 passed calibrated
distance support, while four passed the complete native-L4 domain gate and 24
were withheld by an `out_domain` local vote. The result and fixed-UMAP figure
are documented in `20_immgen_bulk_gated_20d_mapping.md`.

This selection is recorded in
`config/04_reference_mapping/32_qualified_20d_reference_mapping_default.yaml`.
It does not activate `Bdev_implementation`, authorize a paper analysis, or
replace the requirement to build, qualify, checksum, and explicitly accept the
new immutable SIF.
