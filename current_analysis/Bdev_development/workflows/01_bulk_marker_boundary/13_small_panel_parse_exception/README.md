# Apply the small-panel Parse exception

This step applies a user-authorized lower gene-count allowance to two-fold
stable candidates that failed only the latest Parse 20-gene floor. The lower
floor is five genes, matching the earlier archived Parse discovery minimum.

Small panels retain their complete two-fold top-100 intersection and are not
greedily reduced further. Their `moderate_twofold_small_panel` label must remain
attached in downstream evaluation. This diagnostic does not accept or import
gene sets and reads no single-cell expression data.
