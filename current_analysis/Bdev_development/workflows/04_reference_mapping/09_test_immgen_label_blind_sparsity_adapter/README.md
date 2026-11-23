# Test the label-blind ImmGen sparsity adapter

This bounded development diagnostic converts each of the 28 dense ImmGen bulk
peak-count profiles into 100 sparse pseudo-cells before PeakVI projection. The
adapter uses only the global Bdev reference median peak-count depth and the
global sparsity distribution among reference cells within 10% of that depth.
Within each ImmGen library, peaks are sampled in proportion to their count
rank. No population, RNA label, cluster, pseudotime, sample identity, or
expected position enters pseudo-cell construction or parameter selection.

The current scWorkbench runtime does not provide an authenticated cross-bundle
model-import primitive. The workflow therefore reconstructs the already
qualified deterministic Bdev-only PeakVI recipe, requires the exact qualified
reference-latent checksum, and projects queries without retraining. The fixed
ATAC-UMAP display remains the existing primary-30 PeakVI-neighbor-weighted
barycentric overlay; it is not `umap_transform` and does not refit UMAP.

Run the numbered actions with
`config/04_reference_mapping/15_immgen_label_blind_sparsity_adapter.yaml`:

1. `01_run_sparsity_adapter.R` constructs, projects, and audits the pseudo-cell
   ensemble, then writes fixed-map, native-L4, convergence, rank-preservation,
   and neighbor-reuse diagnostics.
1b. `01b_audit_hub_concentration.R` performs a strictly post-projection audit
   of reference-cluster, reference-sample, and recurrent-cell hub usage.
2. `02_finalize_run.py` validates the label-blind construction declaration and
   inventories the generated results.

This is a diagnostic test, not an accepted mapping or paper analysis. It does
not modify the source bundles, operation packs, scWorkbench core, or PyBench.

Run `immgen_sparsity_adapter_20261116T012301Z` completed all 2,800 projections.
The adapter exactly matched the intended sparsity distribution and preserved
bulk rank, but it produced a stronger shared-hub collapse: only 75 reference
cells appeared anywhere in the 84,000 primary-30 neighbor assignments, and the
top hub appeared for 2,775 of 2,800 pseudo-cells. The median library neighbor
weight was 0.463 for cluster 6 and 0.725 for reference sample 1. Most libraries
collapsed near native L4 25.8. The rank-sparse adapter is therefore rejected as
a replacement mapping; matching marginal sparsity without cell-level
co-accessibility is insufficient.
