# Test median-depth ImmGen pseudo-cells

This bounded development diagnostic samples 30 pseudo-cells from each of the
28 ImmGen bulk ATAC libraries at the median qualified Bdev reference matrix
depth of 18,952 peak counts. It fits the unchanged deterministic Bdev-only
PeakVI recipe, requires the qualified reference-latent checksum, projects all
840 pseudo-cells without query retraining, and compares their L4 and neighbor
reuse with the full-depth projection.

Run the numbered actions with
`config/04_reference_mapping/13_immgen_median_depth_downsample_test.yaml`:

1. `01_run_median_depth_test.R` creates and projects the pseudo-cells, derives
   the primary-30 neighbor ledger, and writes comparison tables and figures.
2. `02_finalize_run.py` validates and inventories the generated result tree.

The original mapping bundle is read-only. This is a diagnostic test, not a
paper analysis or accepted replacement mapping.

Run `immgen_median_depth_downsample_20261115T235830Z` completed all 840
pseudo-cell projections. Median-depth consensus retained 228 unique reference
neighbors versus 230 at full depth, median full-depth/consensus Jaccard was
0.968, and median absolute L4 shift was effectively zero. Thus matrix depth
alone does not explain the shared-neighbor hubs.
