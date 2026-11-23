# ImmGen position-driver and technical-bias diagnostic

Run `immgen_position_drivers_20261116T003136Z` attributes the rejected ImmGen
mapping to provisional RNA-label neighborhoods, PeakVI latent dimensions,
empirical peak contrasts, and technical covariates. It is diagnostic only. The
provisional ImmGen RNA labels remain unaccepted, the authenticated native L4
was not modified, and neither source bundle was changed.

## RNA-label comparison

The 1,746-cell provisional RNA-label fixture was projected through the same
deterministic Bdev-only PeakVI recipe without query retraining. Labelled cells
were never added to the 1,382-cell mapping reference. Fifteen provisional
labels had at least five cells and produced sample-balanced, stability-weighted
latent centroids. Seventeen of the 28 ImmGen libraries had a sufficiently
populated centroid for their exact source label; sparse mature labels were not
treated as targets.

For those 17 libraries, the distance to the observed 30-neighbor centroid was
only 0.262 times the distance to the corresponding provisional RNA-label
centroid at the median. The median neighbor weight assigned to the expected
provisional label was zero. `proB.CLP.BM` was the exception: both libraries had
expected-label neighbor weight 1. This demonstrates that most evaluated
libraries are being placed close to a coherent but different latent
neighborhood, rather than merely fluctuating around their expected label.

## Latent-vector driver

`latent_3` was the dominant displacement component for all 17 evaluable
libraries and accounted for a median 31.0% of absolute latent pull. Within
replicated ImmGen populations, `latent_3` correlated with peak-overlap fraction
at Spearman 0.840. Within Bdev clusters it also correlated with TSS enrichment
(`0.443`), ATAC feature count (`-0.422`), fragment count (`-0.414`), ATAC count
(`-0.391`), and nucleosome signal (`-0.355`). The same relationships persisted
after centering the labelled Bdev diagnostic cells within provisional labels.
This dimension therefore mixes biology with a strong technical accessibility
and library-complexity axis.

## Technical bias

After centering within replicated ImmGen populations, peak-overlap fraction
was strongly associated with fixed UMAP 1 (`0.894`), median neighbor distance
(`-0.868`), native L4 (`-0.856`), and fixed UMAP 2 (`-0.787`). More mapped
fragments (`0.574`) and more nonzero peaks (`0.535`) were associated with
greater cluster-6 neighbor weight.

Reference-cell reuse was also technically structured after centering within
the original ATAC clusters: reuse correlated negatively with TSS enrichment
(`-0.283`) and positively with fragment count (`0.254`), ATAC feature count
(`0.250`), nucleosome signal (`0.243`), and ATAC count (`0.238`). Reference
sample 2 supplied 44.8% of mean neighbor weight despite containing 17.1% of the
reference, a 2.63-fold enrichment; sample 4 was depleted to 0.16-fold.

## Peak-axis mismatch

The 191,433-peak audit found 68,570 peaks detected in all 28 ImmGen bulk
libraries. Of these, 48,618 were accessible in at most 10% of Bdev reference
cells. Conversely, 2,709 peaks were absent from all bulk libraries, including
36 accessible in at least half of the reference cells. Reference-versus-bulk
peak prevalence correlation was only 0.574.

Sixty-five peak candidates recurred among the top 100 drivers for at least five
libraries; 62 were all-or-none across the full ImmGen panel. These peaks are
strong candidates for the shared bridge pull, but the ranking is empirical and
does not prove encoder causality.

## Interpretation and next test

The current evidence supports a combined assay-granularity, peak-saturation,
reference-sample, and QC-axis explanation for the cluster-6 hubs. Raw depth
alone was already excluded by the median-depth experiment. The next bounded
test should mask the recurrent all-or-none candidate peaks in the ImmGen
queries, project without retraining, and compare the change in `latent_3`,
expected-label distance, cluster-6 weight, fixed-UMAP position, and native L4
against prevalence-, width-, and count-matched random peak sets. Only that
perturbation can establish which candidate peaks causally drive placement.
