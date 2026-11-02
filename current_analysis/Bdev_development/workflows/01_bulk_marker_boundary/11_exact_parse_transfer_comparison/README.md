# Compare the exact applicable Parse transfer on ImmGen

This diagnostic transfers the latest archived Parse construction to ImmGen:
three cross-fitted top-100 rankings, their intersection, UCell-style capped-rank
scoring, and greedy deletion against the strongest competitor. ImmGen lacks
matched three-replicate blocks and control-cell distributions, so target units
are deterministically partitioned into three folds and only the Parse mean
route is applicable. These structural adaptations are recorded explicitly.

The workflow uses only ImmGen counts and metadata for construction. It compares
the transferred memberships with the prior production, positive-pilot, and
nested-stable memberships under the same score. No single-cell expression data
are read, and no set is accepted or exported for scWorkbench import.
