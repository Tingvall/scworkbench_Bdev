# Competitor-refined cytokine gene sets and Bdev reporting

## Status and next step

This document freezes the method tested on Flt3L and IL-7 and defines how it
was extended to all Parse cytokine signals. The unpromoted all-cytokine review
completed on 2026-09-25.

The next analysis will attempt an up and a down signature for each of the 26
single-cytokine conditions. Control is a competitor and threshold reference,
not a target signature. This gives 52 attempted cytokine directions. A
direction can fail during candidate construction or fail the final Parse
discrimination tests; these failures must remain in the final inventory with
an explicit reason.

The workflow has two parts:

1. Define and test competitor-refined gene sets in the Parse dataset.
2. For final passing sets, calculate native UCell scores in Bdev and report
   them on the existing ATAC UMAP and along the existing Bdev pseudotime.

Existing analyses and bundles remain unchanged. The new analysis must use new,
versioned output directories and must not be promoted until it has been
reviewed.

This is a fixed top-100 stability-intersection method followed by gene removal.
It does not select among several set sizes, does not run the earlier
`full`/`top25` dual strategy, and does not choose a set merely because it has
the smallest allowed size. Final size is determined by whether another
single-gene deletion improves the competitor-aware gap, with a floor of 20
genes.

## Part 1: define the cytokine gene sets in Parse

### 1. Fix the eligible gene universe

Only genes present in the accepted Bdev `SCT_counts` feature universe are
eligible. This ensures that every selected gene can be represented in the Bdev
reference. Do not substitute genes that are absent from this universe.

The Bdev `SCT_counts` matrix defines eligibility only. Parse differential
expression and Parse `RNA_counts` define and validate the signatures.

### 2. Rank genes separately in three cross-fitted analyses

For every cytokine and direction, use the three existing leave-one-replicate-out
specificity fits. Each fit is trained on two replicates and leaves the third
replicate out. Rank genes independently in each fit.

For an up set, orient the target-versus-other-cytokines specificity coefficient
and moderated t statistic in the positive direction. For a down set, multiply
both by -1 so that a larger oriented value always means a stronger response in
the requested direction.

A gene is eligible for ranking when:

- it is present in Bdev `SCT_counts`;
- its oriented specificity log2 fold change is strictly greater than 0.5; and
- its coefficient and moderated t statistic are finite.

FDR is used for ranking but is not a hard inclusion requirement. Rank genes in
this order:

1. FDR < 0.05 before FDR >= 0.05 or unavailable;
2. larger oriented moderated t statistic;
3. larger oriented specificity log2 fold change;
4. smaller FDR; and
5. feature identifier as the deterministic final tie-breaker.

The specificity contrast ranks the target cytokine against the other 25
single-cytokine conditions. Control is introduced later during competitor
testing and responder-threshold calculation.

### 3. Construct the stable candidate

Keep the top 100 eligible genes from each of the three fold-specific rankings.
The candidate set is their three-way intersection. Order candidate genes by
mean fold rank, then worst fold rank, then feature identifier.

Construction passes when the intersection contains at least 20 genes. Do not
pad a smaller intersection with less stable genes. A direction with fewer than
20 genes fails as `FAIL_FEWER_THAN_20_IN_THREE_FOLD_TOP100_INTERSECTION`.

Every candidate gene has therefore appeared in the top 100 in all three
cross-fitted analyses. Because refinement only removes candidate genes, every
gene in the final set retains this property. The final set is itself the
cross-fitted core; this version does not define a separate smaller core set.

### 4. Build one reusable UCell contribution cache

Take the union of all fold-specific top-100 lists across all 52 attempted
directions. Rank Parse `RNA_counts` with UCell 2.12.0 using average ties and
`maxRank = 1500`. Cache each union gene's mean capped-rank contribution in each
of the 81 Parse samples.

For a set of `n` genes with sample-level mean capped ranks `r`, reconstruct the
sample mean UCell score exactly as:

```text
1 - (sum(mean(r)) - n(n + 1)/2) /
  (n * 1500 - n(n + 1)/2)
```

This cache is reusable when genes are removed because the UCell statistic is
affine in the capped-rank sum at a fixed set size. The all-cytokine run must
verify reconstructed sample means against native scWorkbench UCell results.

### 5. Measure separation from the strongest competitor

For each candidate set, calculate one mean UCell score per Parse sample. For a
down set, multiply sample scores by -1 when calculating discrimination metrics,
so a larger oriented response always means a stronger response in the requested
direction.

Summarize each of the 27 conditions by the median of its three sample means.
For the target cytokine, identify the largest median among all other cytokines
and control. This is the strongest competitor. Define the optimization metric
as:

```text
(target condition median - strongest competitor median) /
  SD(the 27 condition medians)
```

Record the target value, strongest competitor and its value, standardized gap,
and target rank at every step.

### 6. Refine the candidate by greedy gene removal

Starting with the complete candidate, test deletion of every remaining gene.
Remove the gene whose deletion produces the largest strictly positive increase
in the standardized mean gap. Recalculate all deletion choices after each
removal.

There is no arbitrary minimum improvement. Any strictly positive improvement
qualifies. If two deletions give exactly the same improvement, remove the gene
with the earlier original candidate rank.

Stop when no single-gene deletion improves the gap or 20 genes remain.

Keep the full optimization trace, including every tested deletion, the gap
before and after deletion, the strongest competitor, the selected deletion,
and the stop reason.

### 7. Recalculate native UCell scores and apply the final tests

Run native scWorkbench UCell on Parse `RNA_counts` for the initial candidate
and final refined set, using all Parse cells and `max_rank = 1500`. Confirm that
native sample means agree with the algebraic reconstruction to numerical
precision.

Apply two final tests to the refined set.

**Mean test**

- Calculate the native mean UCell score in each sample.
- Orient down scores by multiplying by -1.
- Take the median of the three samples for every condition.
- Pass only when the target median is strictly larger than every other cytokine
  and control.

**Responder test**

- Within each replicate, use matched control cells to define the threshold.
- For an up set, the threshold is the control-cell 95th percentile and a cell
  responds when its score is strictly above the threshold.
- For a down set, the threshold is the control-cell 5th percentile and a cell
  responds when its score is strictly below the threshold.
- Calculate the responding-cell fraction in each sample.
- Take the median of the three responder fractions for every condition.
- Pass only when the target median responder fraction is strictly larger than
  every other cytokine and control.

The final set passes when the mean test and/or responder test passes. Report
the two routes separately as well as the combined result. A tie with a
competitor is a failure.

Refinement optimizes only the mean standardized gap. The responder result is an
additional final test and does not influence gene removal.

### 8. Required Parse outputs

The all-cytokine result must include:

- one row for every attempted cytokine direction, including construction and
  validation failures;
- all three fold-specific ranked top-100 lists;
- candidate and final gene memberships with fold ranks and final ranks;
- the union-gene capped-rank cache and its tabular export;
- the complete greedy-removal trace;
- native initial and final UCell summaries;
- sample-level mean scores, responder thresholds and responder fractions;
- condition-level mean and responder metrics, strongest competitors and ranks;
- explicit mean-pass, responder-pass and overall-pass columns; and
- input, parameter, runtime, checksum and bundle-fingerprint records.

## Part 2: report passing signatures in Bdev

Only final sets that pass the Parse mean and/or responder test are treated as
passing cytokine signals in the primary Bdev figures. Construction failures and
nonpassing refined sets remain in the inventory and may receive clearly labeled
diagnostic figures, but they must not be presented as passing signals.

### 9. Calculate Bdev UCell scores

For every final passing set, run native scWorkbench UCell on Bdev raw
`RNA_counts` for all 26,723 cells, using `max_rank = 1500`, no graph smoothing,
and an error on missing features. Preserve the exact memberships selected in
Parse; do not replace missing genes.

Save raw cell-level UCell scores and a signature summary containing the input,
matched and missing gene counts. The Bdev bundle must remain unchanged and the
scoring session must remain unpromoted during review.

### 10. Make UCell-score UMAP figures

Plot the raw UCell scores on the existing Bdev `umap.atac` coordinates. Create:

- an overview containing all passing signatures;
- grouped figures that are practical to inspect by cytokine and direction; and
- separate high-resolution figures for each signature.

Use the established bright blue-cyan-green-yellow-orange-magenta scale. Preserve
the raw UCell values in the exported table. For display, scale each panel to the
99th percentile among cells with a positive score and draw zero-score cells in
gray. State prominently that panel-specific color normalization makes colors
unsuitable for comparing absolute UCell magnitude between signatures.

Each panel must show the cytokine, direction, final gene count and percentage
of cells with UCell score greater than zero. Use the unchanged ATAC-derived UMAP
only as a display coordinate; do not refit the embedding.

### 11. Make 5% sliding-window pseudotime figures

Use the existing `Bdev@meta.data$pt` values exactly. Include only the 12,642
cells with finite `pt` when calculating trajectories, and preserve the 14,081
missing entries as missing. Do not refit, reverse, rescale or impute pseudotime.

For every passing signature:

- use a window width equal to 5% of the observed `pt` range;
- advance windows in steps of 1% of the observed range;
- include cells satisfying `window_left <= pt <= window_right`;
- clip endpoint windows to the observed range;
- place each plotted point at the midpoint of its actual window,
  `(window_left + window_right) / 2`;
- require at least 20 cells in a window;
- plot the median UCell score as the line; and
- plot the type-7 25th to 75th percentiles as the IQR band.

Windows with fewer than 20 cells are omitted from the curve and shown as gray
unsupported regions. Do not add circles, hollow points or other point markers
for either supported or omitted windows. With the current Bdev cell density,
this retains the windows near `pt` -18.5 containing 23 and 25 cells while
omitting the terminal windows containing only 2–3 cells.

Do not add curve smoothing. Do not connect across unsupported windows. The IQR
describes cell-to-cell spread and is not a confidence interval.

Add dotted vertical lines at the median original `pt` of each nonmissing group
in the Bdev `bb_up` metadata field. In the current reference these groups are
`up1` and `up2`. Export their cell counts, `pt` ranges and exact median line
locations rather than hard-coding the numerical positions.

Create an overview across all passing signatures, cytokine-level up/down
figures, and tabular exports of every window's bounds, midpoint, support,
median and IQR.

## Interpretation and limitations

- Cross-fitting stabilizes candidate membership, but it is not an independent
  validation test: refinement and final pass decisions use the complete three-
  replicate Parse dataset.
- Greedy backward removal finds a locally improved set and need not find the
  globally optimal subset.
- The same data define the strongest competitor, guide gene removal and test
  the final set, so final gaps can be optimistic.
- Three replicates provide limited information about between-replicate
  variability. The strict largest-condition rule is descriptive rather than a
  formal significance test.
- Overall pass uses mean **or** responder. The passing route must always be
  shown; the combined label must not hide failure of the other route.
- FDR prioritizes robust genes but does not filter all FDR >= 0.05 genes.
  Stability across fold rankings and multigene discrimination provide the
  additional constraints.
- Shared genes and developmental programs can make different cytokine
  signatures correlated. A Bdev pattern describes similarity to a Parse-derived
  response program, not actual cytokine exposure or causality.
- Panel-normalized UMAP colors show spatial pattern, not comparable absolute
  activity between signatures.
- Sliding windows overlap and are therefore statistically dependent. Their IQR
  bands describe the cells in each window and do not quantify uncertainty of
  the median trajectory.

## Executed workflow

1. Generalize the Flt3L/IL-7 refinement script from two cytokines to all 26
   single-cytokine conditions and both directions.
2. Build the three fold-specific top-100 lists and one reusable union-gene
   contribution cache.
3. Construct each three-way intersection, mark sets with fewer than 20 genes as
   construction failures, and greedily refine every constructible set to a
   minimum of 20 genes.
4. Run native Parse UCell for initial and final memberships, verify the
   algebraic scores, and calculate the mean and responder pass routes.
5. Produce the complete 52-direction inventory and review the gene memberships,
   competitors, optimization traces and pass decisions.
6. After that review, score final passing sets in Bdev with native UCell.
7. Generate the UCell-score UMAP figures and the 5% sliding-window median/IQR
   pseudotime figures with `bb_up` vertical lines.
8. Validate cell and feature alignment, exact source `pt`, missingness, bundle
   fingerprints, output checksums and the absence of promotions.

All steps in this action list were executed for the version-1 unpromoted
review. The Parse stage attempted 52 directions, constructed 43 candidates and
retained 38 final passing signatures. The Bdev stage scored those 38 signatures
and generated 23 trajectory and 61 UMAP figures. The registered report paths
are in `config/parse_cytokines_rest_only.yml` and
`config/bdev_reference.yml`.

## Reference implementation from the completed pilot

The method is grounded in the completed Flt3L/IL-7 top-100 pilot and its Bdev
reports:

- `processing/parse_cytokines/17_flt3_il7_competitor_refinement_pilot.R`:
  cross-fitted ranking, top-100 intersection, capped-rank cache, greedy
  refinement and final Parse tests;
- `processing/parse_cytokines/COMPETITOR_REFINEMENT_PILOT_WORKFLOW.md`:
  concise pilot method;
- `figure_scripts/bdev_reference/05_flt3_il7_refined_ucell_trajectory.R`:
  native Bdev UCell scoring and 5% sliding-window trajectories; and
- `figure_scripts/bdev_reference/06_flt3_il7_refined_ucell_umap.R`:
  UCell-score display on Bdev `umap.atac`.

The all-cytokine implementation should generalize these scripts rather than
altering the completed pilot outputs.
