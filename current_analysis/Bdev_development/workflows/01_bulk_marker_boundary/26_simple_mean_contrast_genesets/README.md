# Simple ImmGen population-balanced mean-contrast gene sets

This exploration generates 20 genes per ImmGen population using one simple
rule: median target log2-CPM minus the mean of the other population-median
log2-CPM values. Genes must have target log2-CPM at least 1 and positive
contrast. The supplied sets are used only for post hoc similarity assessment.

The 64 generated sets are scored on the bounded Bdev fixture using native
scWorkbench UCell and direct raw-score reduction at margin 0.01. No percentile
calibration, broad gate, label promotion, or Bdev-informed selection is used.
