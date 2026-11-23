# Full-Bdev PeakVI feature and mapping contract

## Scientific decomposition

The reference representation and trajectory assignment are separate. PeakVI
learns a whole-Bdev accessibility space from all qualified Bdev cells. Mapping
first reports support and Bdev neighborhood. L4 eligibility is then evaluated
from authenticated fitted Slingshot L4 weights; native L4 pseudotime is
transferred only for eligible mappings and only from L4-compatible neighbors.
Cells omitted from the fitted lineage-family weight axis because their clusters
were excluded from Slingshot fitting are represented as weight zero for the
secondary L4-eligibility task. Supplied weights and pseudotimes are unchanged.

## Feature contract

The primary feature axis is selected without external queries and without a
sample-replication requirement. For binary accessibility, depth is represented
by Bdev depth quantiles. The cluster branch ranks depth-adjusted one-versus-rest
binomial deviance, stabilized by deterministic within-cluster split halves.
Cluster quotas are proportional to the square root of cluster cell count so
large clusters cannot dominate while very small clusters do not receive an
unsupported equal quota. The overall branch ranks residual binomial deviance
without cluster labels. It contributes only peaks not already selected by the
cluster branch, yielding exactly 10,000 peaks from each route.

The final matrix retains the source feature order. The selection table records
both ranks, the selected route, cluster target where applicable, prevalence,
and split-half stability. Cluster-only, overall-only, and complete-axis models
are Bdev-only qualification controls; they are not alternative external-query
tuned models.

## Freeze boundary

Candidate choice, support calibration, L4 eligibility, and pseudotime transfer
must be completed from Bdev validation data. Once the selected recipe is fit on
all Bdev cells, the ordered feature list, contained model, reference latent
coordinates, thresholds, configuration, runtime sources, and checksums are
immutable. ImmGen is absent from every design input and can be evaluated only
after that freeze boundary.

## Qualification result

The exact 50/50 axis was constructed and all four candidates completed on
21,389 reference cells, 5,334 disjoint held-out cells, and 104 Bdev-derived
pseudobulks. No external query was read.

| Candidate | Macro cluster accuracy | Pseudobulk accuracy | Maximum absolute depth-latent correlation | L4 AUC | L4 MAE |
|---|---:|---:|---:|---:|---:|
| hybrid20k | 0.640 | 0.798 | 0.543 | 0.929 | 0.635 |
| cluster20k | 0.698 | 0.769 | 0.390 | 0.935 | 0.548 |
| overall20k | 0.517 | 0.779 | 0.715 | 0.886 | 0.787 |
| full191433 | 0.660 | 0.837 | 0.375 | 0.939 | 0.585 |

The primary hybrid exceeded its allowed 0.03 macro-accuracy loss and its
allowed depth-correlation excess over the complete axis. The gate therefore
failed and the all-cell freeze step was not run. Candidate artifacts are
preserved. A cluster-majority hybrid is the evidence-supported next bounded
test, but its ratio must be declared in a new run rather than changed here.
