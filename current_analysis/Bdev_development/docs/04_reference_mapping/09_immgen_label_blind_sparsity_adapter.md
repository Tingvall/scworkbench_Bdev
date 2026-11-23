# ImmGen label-blind sparsity adapter

## Question

Can the existing qualified Bdev PeakVI recipe yield more informative ImmGen
bulk placements if each dense bulk profile is first converted into an ensemble
with single-cell-like fragment depth and sparsity, without using any biological
label or expected position?

## Construction boundary

The adapter may use only:

- the global median peak-count depth of the qualified Bdev reference;
- the global distribution of nonzero peak counts among reference cells within
  10% of that depth; and
- the within-library rank of nonzero ImmGen peak counts.

It may not use population labels, provisional RNA labels, reference clusters,
pseudotime, reference sample identity, or an expected UMAP location. Each
pseudo-cell receives exactly the global median fragment count and an empirical
sparsity draw. Peaks are selected by label-free rank-weighted successive
sampling, and fragments are allocated using the same within-library ranks.

## Interpretation boundary

The pseudo-cells preserve marginal bulk peak ordering but cannot reconstruct
cell-level co-accessibility that was lost during bulk preparation. The output
therefore represents uncertainty in population-level placement, not inferred
single cells or biological heterogeneity. Any post-projection comparison with
the earlier full-bulk location is audit-only and cannot select adapter
parameters. The resulting mapping remains diagnostic and unaccepted.

## Result

Run `immgen_sparsity_adapter_20261116T012301Z` completed 100 pseudo-cells for
each of 28 libraries. All pseudo-cells had exactly 18,952 fragments. Their
nonzero-peak distribution (median 8,144; range 7,053–9,401) matched the
unlabelled reference calibration pool (median 8,149.5; range 7,053–9,401),
whereas the earlier depth-only pseudo-cells retained a median of 15,914 nonzero
peaks. Median bulk-rank preservation was 0.806 for accessibility frequency and
0.868 for allocated fragments. The seeded ensemble was numerically stable:
the median 50-versus-100 pseudo-cell fixed-map difference was 0.0039.

Despite satisfying those technical targets, the biological geometry worsened.
The median fixed-map displacement from the full-bulk projection was 2.25 and
the median absolute native-L4 shift was 2.78. Most libraries converged near
native L4 25.8. Across all 84,000 primary-30 neighbor assignments, only 75
distinct reference cells were used; the most recurrent cell occurred for
2,775 of 2,800 pseudo-cells. Post-projection auditing found a median cluster-6
neighbor weight of 0.463 and median reference-sample-1 weight of 0.725.

This rejects independent rank-weighted peak sampling as the adapter. It shows
that matching marginal fragment count and sparsity is not enough: independently
selecting peaks destroys cell-level co-accessibility and makes different bulk
libraries resemble the same generic sparse profile. A further fixed-model
adapter should preserve multipeak structure, for example by fitting a
label-blind mixture of existing PeakVI decoder profiles to the quantitative
bulk ranks. Merely changing the rank exponent is not justified by this result.
