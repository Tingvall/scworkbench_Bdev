# Simple mean-contrast three-rule comparison

**Run:** `simple_mean_contrast_margins_20261106T234706Z`  
**Status:** complete diagnostic; selected 0.01 default unchanged

This analysis reduced the same authenticated 64-program raw UCell score family
with three nested rules: pure argmax, winner margin at least 0.01, and winner
margin at least 0.05. Exact top-score ties remained unassigned rather than
being broken arbitrarily. The independently recomputed 0.01 decisions exactly
matched the recorded baseline labels, scores, runner-up scores, and margins.

| Rule | Assigned cells | Unassigned cells | Represented labels |
|---|---:|---:|---:|
| Pure winner | 1,709 | 37 | 33 |
| Margin >= 0.01 | 1,416 | 330 | 31 |
| Margin >= 0.05 | 443 | 1,303 | 17 |

The nested cell strata are 37 exact ties, 293 cells assigned only by pure
argmax, 973 cells passing 0.01 but not 0.05, and 443 cells passing 0.05. Thus,
pure argmax maximizes coverage but forces 293 very small-margin calls, while
0.05 removes most assignments. The 0.01 rule remains the provisional balance
between coverage and winner separation.

The result root is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/01_bulk_marker_boundary/30_simple_mean_contrast_margin_comparison/simple_mean_contrast_margins_20261106T234706Z/
```

All 13 manifest-listed outputs reproduced their recorded sizes and SHA-256
values. The final manifest SHA-256 is
`a3fefdf93acbf1bb99dc1b97e9b613126b534ccf8da65f77509013655c8f2710`.

