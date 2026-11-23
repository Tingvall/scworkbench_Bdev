# Diagnose ImmGen L4 and fixed-UMAP placement drivers

This development diagnostic explains the rejected fixed-`umap.atac` overlay
without changing any source bundle. For every ImmGen population it shows the
native fitted L4 distribution of the primary 30 Bdev PeakVI neighbors and the
fixed-UMAP coordinates that contribute to the recorded library placement.

The workflow independently reconstructs two distinct summaries from the same
neighbor ledger:

- mapped L4 is the distance-kernel-weighted median of neighbor native-L4;
- fixed-UMAP placement is the distance-kernel-weighted arithmetic mean of the
  neighbors' pre-existing `umap.atac` coordinates.

Run with
`config/04_reference_mapping/12_immgen_l4_umap_driver_diagnostic.yaml`:

1. `01_build_diagnostic.R` validates checksums, reconstructs both outputs,
   writes driver ledgers and summaries, and renders population-level figures.
2. `02_finalize_run.py` checks required outputs and inventories the diagnostic.

The source mapping has not been scientifically accepted. These figures are
diagnostic only and must not be used downstream.

Run `immgen_l4_umap_drivers_20261115T233506Z` completed for all 28 libraries
and 16 populations. It reconstructed mapped L4 with maximum absolute error
`4.62e-14` and the fixed-UMAP coordinates with maximum Euclidean error
`7.91e-15`; the Bdev source-bundle manifest remained unchanged.
