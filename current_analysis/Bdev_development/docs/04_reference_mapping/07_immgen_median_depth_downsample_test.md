# ImmGen median-depth downsample test

Run `immgen_median_depth_downsample_20261115T235830Z` tested whether the
full-depth ImmGen neighbor hubs were caused primarily by library depth and peak
saturation. Each of the 28 libraries produced 30 seeded multinomial
pseudo-cells at exactly 18,952 peak counts, for 840 frozen-recipe PeakVI
projections. The deterministic Bdev-only reference fit reproduced its qualified
latent checksum, and the source mapping manifest remained unchanged.

## Result

Downsampling did not remove the shared-neighbor structure.

| Diagnostic | Full depth | Median-depth consensus |
|---|---:|---:|
| Unique cells among 28 × 30 neighbors | 230 | 228 |
| Cells selected by at least 2 libraries | 174 | 173 |
| Cells selected by at least 5 libraries | 79 | 80 |
| Cells selected by at least 10 libraries | 3 | 2 |
| Maximum library occurrence | 10 | 11 |

The median Jaccard similarity between each full-depth neighbor list and its
30-replicate median-depth consensus was 0.968. The median absolute shift in
mapped native L4 was `4.62e-14`; the largest shift was only 0.273. Twenty-one
libraries were supported in all 30 replicates. The remaining support behavior
mostly reproduced the original borderline calls; no pseudo-cell was called
out-of-reference.

Median pseudo-cell sparsity improved from 178,848 nonzero peaks at full depth
to 15,914, but remained above the Bdev reference median of 8,145. Replicate
neighbor stability was nevertheless high (median pairwise Jaccard 0.771), and
the consensus was nearly identical to full depth.

## Interpretation

Raw depth and saturation are therefore not the primary cause of neighbor reuse.
Multinomial sampling preserves each bulk library's relative peak composition,
and the resulting profiles still behave as population mixtures rather than
real single cells. Remaining explanations include shared bulk composition,
PeakVI/reference-space hub geometry, the documented reference-sample bias, and
the absence of an explicit post-FrE negative/domain gate. A forced 30-neighbor
query will continue to return the same dense reference neighborhoods even for
biologically ineligible populations.

The test is diagnostic only. It does not rehabilitate the original mapping or
authorize downstream use.
