# Simple mean-contrast winner-margin sensitivity

**Run:** `simple_mean_contrast_margin005_20261106T230556Z`  
**Status:** complete diagnostic; selected 0.01 default unchanged

This analysis reapplied the label-reduction rule to the same 64 raw UCell
scores and 1,746-cell Bdev fixture using a winner–runner-up margin of 0.05
instead of 0.01. It did not regenerate gene sets or UCell scores.

The stricter threshold retained 443 of the 1,416 baseline assignments
(31.3%), made 973 additional cells unassigned, and reduced the represented
labels from 31 to 17. Total unassigned cells increased from 330 to 1,303.
Sample 1 retained 137/555 baseline assignments (24.7%); samples 2–4 retained
34–37%.

The largest absolute losses were `B.FrE.BM` (326), `MPP4.135+.BM` (170), and
`STHSC.150-.BM` (139). `proB.FrBC.BM` retained 93/150 calls (62%), whereas
`B.T3.Sp` retained none of its 37 calls. The assignment-loss UMAP shows that
removed calls occur broadly across the embedding rather than at one isolated
boundary.

The result root is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/01_bulk_marker_boundary/29_simple_mean_contrast_margin_sensitivity/simple_mean_contrast_margin005_20261106T230556Z/
```

All ten manifest-listed outputs reproduced their recorded sizes and SHA-256
values. The final manifest SHA-256 is
`b9210aeabab2a04d95800d1a9bb3345e854dccd3a07bf900dbd3037093933a0e`.
The 0.05 rule is retained as a high-specificity sensitivity analysis and does
not replace the provisional 0.01 development default.

A subsequent comparison against pure argmax is documented in
[`27_simple_mean_contrast_three_rule_comparison.md`](27_simple_mean_contrast_three_rule_comparison.md).
