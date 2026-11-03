# Provisional Bdev development default: simple ImmGen mean contrast

**Selected:** 2026-11-02  
**Status:** active provisional development default; not promoted to implementation

The current working default for Bdev development exploration is the 64-population
ImmGen-only collection from run `simple_mean_contrast_20261102T235008Z`.
The user selected this collection after reviewing the Bdev ATAC UMAP label
distribution.

## Gene-set discovery

The method converts ImmGen counts to `log2(CPM + 0.5)`, takes the median across
libraries in each population, and ranks each gene by the target-population
median minus the mean of the other 63 population medians. Genes must have target
log2 CPM of at least 1 and a positive contrast; the top 20 genes are retained
per population. Bdev and the supplied comparison sets were not used to select
genes.

## Labeling

The current diagnostic labeling rule uses the RNA count matrix and native UCell
scores with `maxRank = 1500`. It selects the raw-score argmax without percentile
calibration or a broad-lineage gate. Calls require a raw winner-versus-runner-up
margin of at least 0.01; ties and lower-margin cells are `unassigned`.

## Evidence and boundary

The supplied-set post hoc comparison placed 519 of 580 supplied genes within
the generated top 100 and 224 within the generated top 20. On the bounded Bdev
fixture, 1,416 of 1,746 cells were assigned. Samples 3 and 4 were B-heavy and
the previous DP/ISP failure was absent apart from three ISP calls in sample 1.

This selection is the default for ongoing development exploration only. It does
not accept the cell labels as truth, promote an artifact into scWorkbench,
modify core or an operation pack, or unlock `Bdev_implementation`. The exact
artifacts and checksums are pinned in
`config/01_bulk_marker_boundary/15_default_immgen_geneset_collection.yaml`.
