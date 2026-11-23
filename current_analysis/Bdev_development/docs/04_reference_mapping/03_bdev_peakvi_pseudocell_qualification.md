# Bdev PeakVI pseudobulk and pseudo-cell qualification

## Outcome

Run `bdev_peakvi_pseudocell_10d_20261110T180851Z` passed the predeclared
ten-dimensional scientific gate on the exact 191,433-feature Bdev ATAC peak
universe. PeakVI was fit only on 1,382 Bdev reference cells. A frozen retained
model then projected 428 disjoint held-out Bdev single cells, 32 raw Bdev
pseudobulks, and 96 depth-matched Bdev pseudo-cells without refitting.

| Query type | Accepted | Native-L4 MAE | Native-L4 Spearman |
|---|---:|---:|---:|
| held-out single cell | 1.000 | 2.087 | 0.940 |
| raw pseudobulk | 1.000 | 1.310 | 0.988 |
| depth-matched pseudo-cell | 1.000 | 1.327 | 0.979 |

All three query classes passed their frozen thresholds. Pseudo-cell depths were
exactly 9,781, 18,952, and 29,738 counts, corresponding to the reference-depth
quartiles. Native fitted Slingshot L4 was used without reorientation. The
active 03c source manifest remained
`cc8b92e8b55bb1b12b7f818c3bc272193c2139abec81b9777d48bfe2dd022a05`.

## Contract acceptance

The run used scWorkbench 0.1.2 commit
`2e7cbbdf50e615504412c5f9bbda76ceabdbd148` and scwbLatent 0.1.3 commit
`cc93bb713d9d66a32ab697a9fc2d949025d0c678`. `peakvi_project` retained all ten
model dimensions in exact query order, and `model_bound_latent_union` produced
the authenticated 1,938-cell joint space without direct attachment or refit.

All structural checks passed: exact feature axis, source-manifest preservation,
disjoint reference and query cells, projection without refit, requested latent
width, authenticated union, exact pseudo-cell depths, and native fitted L4.
The final inventory contains 18 checksummed files. The bounded Step 4 gate is
qualified; real ImmGen projection is now eligible as a separate controlled
development run. No SIF was built and implementation remains locked.
