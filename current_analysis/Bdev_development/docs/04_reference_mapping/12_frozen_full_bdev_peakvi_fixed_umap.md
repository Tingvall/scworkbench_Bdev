# Four frozen full-Bdev PeakVI mappings on the fixed ATAC UMAP

Development step 12 is a read-only visualization of the completed step-11
external benchmark. It addresses the comparison question directly: all four
frozen PeakVI models are shown on the same approved Bdev `umap.atac` reference,
using the same 28 ImmGen libraries and each model's already-recorded 30-neighbor
ledger.

For each model and library, the plotted location is the fixed-UMAP barycenter
of the 30 reference neighbors with weight
`exp(-(distance / kth_distance)^2)`. The ellipse is the corresponding weighted
reference-neighbor covariance. These are projections for visualization; the
fixed UMAP is not transformed and the PeakVI models are not reopened or fit.

The overview provides a common-scale four-panel comparison. Four detailed
figures label every library and show the nearest five reference-neighbor
spokes. A separate trajectory diagnostic shows all provisional native-L4
values, including missing calls. The UMAP panels draw the authenticated fitted
Slingshot `Lineage4` curve in native orientation.

The trajectory values remain `diagnostic_unaccepted`. The figures neither
promote a feature model nor assign the final trajectory position of any bulk
library.

## Observed comparison

The common-scale view supports the quantitative external benchmark. All three
20k feature models place the ImmGen populations in a substantially more
ordered geometry than the all-191,433-peak model. The complete-peak mapping has
broad, strongly overlapping neighbor ellipses through the lower-right portion
of the reference and mixes several mature populations. Cluster20k and
overall20k give the tightest visual population groupings; hybrid20k retains the
same broad ordering but has wider support for several B populations.

The native-L4 panel also makes the current limitation visible. Finite
provisional values are population-structured where present, but the 20k models
leave many mature and stem/progenitor libraries without a finite L4 value.
Full191433 yields more finite values, but its poor fixed-UMAP bulk geometry
means that coverage alone is not evidence for a better trajectory mapping.
