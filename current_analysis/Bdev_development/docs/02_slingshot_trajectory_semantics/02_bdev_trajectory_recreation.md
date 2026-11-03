# Bdev Slingshot trajectory recreation

**Run:** `bdev_slingshot_recreation_20261103T044735Z`  
**Status:** complete development validation; staged outputs not promoted

The published Bdev recipe was recreated through
`scwbtrajectory::slingshot_fit` from `scwbTrajectory` 0.1.3 using the ATAC UMAP,
cluster 5 as the root, exclusion of clusters 9, 15, 18, 19, 22, and 23,
`stretch = 0`, `omega = TRUE`, and `omega_scale = 5`.

The fit retained 23,710 of 26,723 cells and produced eight lineages. Lineage
L4 followed clusters `5-1-2-16-13-12-7-3`. Its negated pseudotime overlapped
the preserved Bdev `pt` field for 12,510 cells and had Spearman correlation
0.999817 and Pearson correlation 0.999942 on that shared finite set. Membership
was not identical: 941 cells were defined only in the recreated L4 and 132 only
in preserved `pt`. Values were also not byte-identical (mean absolute
difference 0.5732; maximum 1.0209), consistent with a near-exact recreation of
ordering rather than exact recovery of the original numeric vector.

The preserved `pt` vector retained its exact binary hash
`2a9e9b9b8cd85daa00a85c3897cb2706844a306513efe70956bda7cd4745d6ee`,
including 12,642 defined and 14,081 missing values. The archived source
manifest hash was unchanged. The lineage family remains the primary trajectory
result. Based on this evidence, fitted L4 is now the selected default Bdev
pseudotime trajectory. Use native L4 orientation for new analyses and negate it
only for explicit alignment with the legacy `pt` orientation. Preserved `pt`
remains a validation reference and fallback.

Generated tables, ATAC-UMAP figures, the scWorkbench session receipt, and the
full result report are under the paired active-storage run directory declared
in
[`02_bdev_trajectory_recreation.yaml`](../../config/02_slingshot_trajectory_semantics/02_bdev_trajectory_recreation.yaml).
