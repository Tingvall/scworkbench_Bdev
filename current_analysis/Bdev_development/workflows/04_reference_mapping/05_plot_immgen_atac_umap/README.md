# Step 4 ImmGen overlay on the approved ATAC UMAP

This derived development workflow displays all 28 ImmGen ATAC libraries on
the established Bdev `umap.atac` coordinates used for fitted Slingshot L4.
It reproduces the authenticated 10D PeakVI neighbor ledger, takes the primary
30 neighbors, and places each library at their kernel-weighted ATAC-UMAP
barycenter. It does not refit or change the approved UMAP or L4 orientation.

Every library is labeled by population and GSM using repelled labels. Run the
numbered actions with
`config/04_reference_mapping/11_immgen_atac_umap_overlay.yaml`.

Run `immgen_atac_umap_20261111T032028Z` completed technically, but its
biological placement was not accepted on 2026-11-11. In particular, the
later-stage populations are implausibly positioned. Retain the checksummed run
as a diagnostic only; Step 4 remains unresolved and this overlay must not be
used downstream.
