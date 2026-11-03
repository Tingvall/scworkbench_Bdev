# L4 default trajectory decision

The updated scWorkbench runtime worked: core commit `da9477d` together with
`scwbTrajectory` 0.1.3 commit `080a192` fitted the published Bdev ATAC-UMAP
Slingshot recipe and generated an almost identical trajectory ordering.
Negated fitted L4 and preserved `Bdev$pt` had Spearman correlation 0.999817 and
Pearson correlation 0.999942 across 12,510 shared finite cells.

Fitted lineage `L4` is therefore the current default pseudotime trajectory for
Bdev analyses. It follows clusters `5-1-2-16-13-12-7-3`. Use native Slingshot
L4 orientation by default. Negate L4 only when explicitly aligning with or
reproducing the legacy scalar `pt` orientation.

The original `Bdev@meta.data$pt` field remains unchanged. It is retained as a
validation reference, a fallback if authenticated fitted L4 is unavailable,
or an explicit legacy source when requested. It must not be overwritten or
silently substituted for fitted L4.

This selects the development and proposed later-implementation strategy. It
does not activate the reserved paper implementation before a new immutable SIF
is qualified and accepted.

The machine-readable decision is
[`03_l4_default_trajectory.yaml`](../../config/02_slingshot_trajectory_semantics/03_l4_default_trajectory.yaml).
