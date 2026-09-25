# Bdev-compatible dual cytokine signature workflow, version 1

## Scope

This workflow is a new review branch. It does not replace, modify, or promote
over the existing competitor-aware analysis or its accepted feature-set
namespaces. It creates separate `full` and `top25` candidates for every
cytokine-direction combination and leaves all results unpromoted.

The Bdev reference contributes only its `RNA_counts` feature identifiers. Bdev
expression values, metadata, pseudotime, and trajectory results do not enter
gene ranking, construction, or validation.

## Gene eligibility and ranking

For each of 26 single-cytokine conditions, the specificity contrast is the
target condition minus the equal-weight mean of the other 25 cytokines.
Controls do not enter this gene-level contrast. Up and down programs are
constructed separately.

A gene is effect eligible when it:

- belongs to the recorded Parse signature-membership universe;
- is present in the accepted Bdev `RNA_counts` feature universe; and
- has direction-oriented specificity log2FC strictly greater than `0.5`.

An effect-eligible gene is called strictly significant when its specificity
FDR is strictly less than `0.05`. Ranking places strictly significant genes
first, then orders genes by direction-oriented moderated t statistic,
direction-oriented log2FC, FDR, and feature identifier. The same deterministic
rule is applied to the full-data model and independently to each of the three
leave-one-replicate-out training models.

## Full and top25 construction

Let `n_significant` be the number of full-data effect-eligible genes with
FDR below `0.05`.

- If `n_significant >= 50`, `full` contains the top 50 significant genes.
- If `25 <= n_significant < 50`, `full` contains every significant gene.
- If `n_significant < 25` but at least 25 effect-eligible genes exist, `full`
  contains all significant genes and the highest-ranked non-significant,
  effect-eligible genes needed to reach 25.
- If fewer than 25 effect-eligible genes exist, construction fails.

`top25` contains the first 25 genes in the same ordered list and also fails
when fewer than 25 effect-eligible genes exist. `full` and `top25` retain
separate identifiers and results even when their memberships are identical.
Fallback genes are effect-size eligible; they are not described as
FDR-significant.

## Full-data discrimination

Both variants are scored independently with native scWorkbench UCell on raw
`RNA_counts`, with `max_rank = 1500` and no graph smoothing. Scores are
summarized separately in each of the three biological replicates.

The mean route uses the sample mean UCell score, oriented so larger always
means a stronger requested response. The responder route uses the fraction of
cells strictly above the matched-control 95th percentile for up sets or below
the matched-control 5th percentile for down sets.

Within each replicate and route, the target is compared with the strongest of
the other 25 cytokines and control. The standardized gap is:

```text
(target response - strongest competitor response) /
  SD(response across all 27 conditions)
```

A route passes when the median of its three replicate gaps is greater than
zero, equivalent to strict target rank one in at least two replicates. Overall
full-data pass requires the mean or responder route to pass. Results for
`full` and `top25` remain separate; neither is selected as a winner.

## Cross-fitted stability

The full-data `full` size is held fixed during cross-fitting. For each held-out
replicate, genes are re-ranked using the model trained on the other two
replicates. The fold-specific `full` set takes the first `N_full` genes and the
fold-specific `top25` set takes the first 25 genes. Significant genes remain
prioritized and effect-eligible fallback genes fill the fixed size when
needed. A fold is not evaluable when too few effect-eligible genes exist.

Each fold-specific set is scored only in its omitted replicate. Mean and
responder gaps are summarized separately for `full` and `top25`. These are
stability annotations and do not alter the full-data pass status.

## Core sets

For each constructed variant, its core is the ordered intersection of the
full-data membership and all three fold-specific memberships:

```text
variant_core = variant_full_data ∩ variant_fold1 ∩ variant_fold2 ∩ variant_fold3
```

Core sets are calculated for every constructed variant, including variants
that fail full-data discrimination. Cores are not padded to a minimum size.

## Interpretation

The minimum size and Bdev feature filter improve technical portability but do
not guarantee expression in Bdev cells. The number of significant genes is
power dependent, and two-replicate training fits can require fallback genes.
The two variants are nested and should not be treated as independent evidence.
Independent data remain necessary to validate exact memberships.
