# Flt3L/IL-7 competitor-refinement pilot

## Scope

This unpromoted pilot analyzes four Parse cytokine directions: Flt3L up,
Flt3L down, IL-7 up, and IL-7 down. All other single-cytokine conditions and
control remain competitors. Existing accepted and review analyses are inputs
only and are not modified.

## Candidate construction

The feature universe is restricted to genes present in the accepted Bdev
`SCT_counts` matrix. For each target and direction, genes are ranked
independently in each of the three leave-one-replicate-out Parse specificity
fits. A gene is eligible when its direction-oriented specificity log2 fold
change is strictly greater than 0.5. Eligible genes with FDR below 0.05 are
ranked first, followed by direction-oriented moderated t statistic,
direction-oriented log2 fold change, FDR, and feature identifier.

The reusable contribution pool is the union of the top 100 eligible genes
from every fold and all four directions. The candidate tested in this pilot is
the intersection of the three fold-specific top-100 lists. Construction passes
when this intersection contains at least 20 genes; it is not padded.

## Algebraic UCell contributions

Parse `RNA_counts` are ranked with UCell 2.12.0 using average ties and
`maxRank = 1500`. The exact UCell capped rank for each pool gene is aggregated
to its mean within each of the 81 samples and cached. For a set containing
`n` genes with capped ranks `r`, its sample mean score is reconstructed as:

```text
1 - (sum(mean(r)) - n(n + 1)/2) /
    (n * 1500 - n(n + 1)/2)
```

This works because the UCell statistic is affine in the sum of capped ranks
for a fixed set size. It permits exact recalculation of sample mean scores
after removing a gene without reranking cells. Native scWorkbench UCell scores
for the initial and final sets are used to verify the reconstructed sample
means and to calculate cell-level responder fractions.

## Full-data refinement

For each candidate, sample mean scores are calculated for all 81 samples.
Each condition is represented by the median of its three sample means. Down
scores are multiplied by -1 so that larger values always indicate a stronger
response in the requested direction. The target is compared with the strongest
of the other 25 cytokines and control. The optimization metric is:

```text
(target condition median - strongest competitor median) /
  SD(the 27 condition medians)
```

At each iteration, every one-gene deletion is tested. The gene producing the
largest strictly positive improvement is removed, with original candidate
rank breaking an exact tie. Contributions are recomputed after every removal.
Refinement stops when no deletion improves the metric or 20 genes remain.

## Final tests

The mean route passes when the target condition has a strictly larger median
sample mean response than every competitor. For the responder route, matched
control cells define replicate-specific thresholds: the 95th percentile for
up sets and the 5th percentile for down sets. The responder fraction is
calculated per sample, summarized as the median of three samples per condition,
and passes when the target is strictly largest. Overall pass is mean or
responder.

Candidate construction uses held-out fits, but refinement and final testing
reuse the complete three-replicate dataset. The final pass is therefore an
internal discrimination result, not independent validation.

## Bdev expression figures

Final genes are displayed on the existing Bdev `umap.atac` coordinates using
values read directly from `SCT_data`. Eligibility is based on `SCT_counts`, so
genes absent from `SCT_data` are retained in the gene-set tables and explicitly
reported as unavailable for plotting.
