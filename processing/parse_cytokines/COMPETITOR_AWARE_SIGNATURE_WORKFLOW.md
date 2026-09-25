# Competitor-aware cytokine signature workflow

## Purpose and scope

This review analysis produces one induced (`up`) and one repressed (`down`)
candidate for each of the 26 individual cytokine conditions. It evaluates all
52 cytokine-direction combinations without a discovery pass/fail gate. The
analysis is separate from the accepted `parse_cytokines` bundle and from the
existing `parse_cytokines_rest_only` staged collection. It does not promote or
replace either result.

The final pass decision describes discrimination within the complete Parse
dataset. A cross-fitted analysis of the selected size is reported separately
as a replicate-stability check. It does not control the final pass decision.

## 1. Gene eligibility and ranking

For each cytokine, the specificity contrast is the target individual stimulus
minus the equal-weight mean of the other 25 individual stimuli. Control samples
do not enter this gene-level contrast.

Separate up and down rankings are constructed. A gene is eligible when:

- it is in the recorded expressed-gene universe allowed for signature
  membership; and
- its specificity log2 fold-change is at least `0.25` in the requested
  direction.

No p-value or FDR cutoff is applied. Eligible genes are ordered by the
direction-oriented moderated t statistic, which prioritizes lower p-values and
FDR values within a contrast, followed by the larger oriented log2 fold-change
and feature identifier for deterministic ties. FDR is retained in the evidence
tables but is ranking evidence rather than an eligibility gate.

The same rule is applied to the full-data model and independently to each of
the three saved leave-one-replicate-out training models.

## 2. Full-data candidate sets

The full-data ranking is truncated at each feasible size in the prespecified
grid:

`5, 10, 15, 20, 30, 40, 50`

Every candidate set is scored per cell with the native scWorkbench `ucell`
operation on `RNA_counts`, using `max_rank = 1500`, no smoothing, and all
154,209 analysis cells.

For each biological sample, two summaries are calculated:

- `mean`: the mean cell-level UCell score;
- `responder_fraction`: for an up signature, the fraction of cells strictly
  above the 95th percentile of the matched control sample; for a down
  signature, the fraction strictly below the 5th percentile of the matched
  control sample.

The mean response is oriented so that larger always means a stronger requested
response: raw mean UCell for up signatures and negative mean UCell for down
signatures. Responder fraction is already oriented so that larger means a
stronger response.

## 3. Competitors, gaps, and route passes

Within each replicate, the target is compared with the other 25 cytokines and
control. The strongest competitor is the non-target condition with the largest
oriented response. All exactly tied strongest competitors are retained.

For each size, replicate, and route, the standardized gap is:

```text
(target response - strongest competitor response) /
  SD(response across all 27 conditions)
```

A replicate gap greater than zero means that the target ranks strictly first.
A zero gap is a tie and does not pass. If the across-condition standard
deviation is zero, the replicate is nondiscriminatory and does not pass.

The full-data discriminatory result for a size and route is the median of its
three replicate gaps. A route passes when this median is greater than zero,
which is equivalent to a strict target rank of one in at least two of three
replicates.

## 4. Full-data size selection

Mean separation is the primary selection route:

1. If any size passes the mean route, select the size with the largest median
   mean gap. This rule also applies to sizes that pass both routes.
2. Otherwise, if any size passes the responder route, select the size with the
   largest median responder gap.
3. If no size passes either route, retain the size with the largest median gap
   across the two routes for reporting and mark the final signature failed.
4. The smaller size is preferred only when the relevant median gaps are
   exactly equal.

The top genes at the selected size in the full-data ranking form the final gene
set. Its mean, responder, and overall pass statuses are inherited from the
full-data calculations. Overall pass requires either the mean or responder
route to pass. Cross-fitted stability never changes this status.

## 5. Holdout stability at the selected size

Only the selected size is tested by cross-fitting. For each held-out replicate:

1. rank genes using the saved model fitted to the other two replicates;
2. take the top genes at the already selected size;
3. score that fold-specific set with native scWorkbench UCell; and
4. calculate the mean and responder gaps in the held-out replicate.

The three held-out gaps are summarized exactly as in the full-data analysis.
The report includes cross-fitted mean and responder pass statuses, median gaps,
strongest competitors, and fold membership overlap. These are stability
annotations for the gene-selection procedure at the full-data-selected size.
Because the size was selected using the complete dataset, this is not an
independent validation of model size.

## 6. Core sets

A core set is reported for every final gene set that passes the full-data mean
or responder rule. The core contains only genes that are:

- members of the final full-data gene set; and
- selected at the same fixed size in each of the three leave-one-replicate-out
  training rankings.

Equivalently, the core is the intersection of the final membership and all
three fold-specific memberships. Core membership is a stability definition,
not an individual-gene significance test. Core sets may contain fewer than
five genes or be empty; this is reported and does not alter the full-set pass.
Failed final signatures do not receive a core set.

## 7. Required outputs

The review report contains:

- one row for every cytokine-direction combination, with the selected size,
  full-data mean and responder gaps and passes, overall pass, strongest
  competitors, holdout stability results, and core size;
- one row per candidate size, with both route gaps, passes, and replicate-level
  competitors;
- one row per selected-size holdout fold;
- complete final membership and evidence for all evaluable signatures;
- complete core membership for all passing signatures;
- TSV and GMT exports for the passing full sets and their core sets; and
- provenance, input hashes, runtime information, and confirmation that no
  result was promoted to a bundle.

## Interpretation and limitations

The final pass status is an internally optimized description of this dataset,
not an independent validation. Searching seven sizes makes the reported
full-data separation optimistic. The holdout analysis tests replicate
stability of gene selection at the chosen size but cannot remove the optimism
from full-data size selection. With only three replicates, each fold ranking is
estimated from two samples per condition and may be unstable. The responder
route can also be uninformative for down signatures when scores accumulate at
zero. Independent data are required to validate the exact exported membership.
