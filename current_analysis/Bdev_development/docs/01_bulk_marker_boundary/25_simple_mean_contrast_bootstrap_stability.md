# Simple mean-contrast bootstrap stability

**Run:** `simple_mean_contrast_bootstrap_stability_20261103T010000Z`  
**Status:** complete diagnostic qualification; provisional default unchanged

Thirty seeded bootstrap rounds resampled ImmGen libraries with replacement
within each of the 64 populations. Every round regenerated the top-20
population-balanced mean-contrast sets. Bdev was then scored with native UCell
and the selected raw argmax/margin-0.01 labeling rule. Bdev was not used to
select genes.

The batched native UCell operation produced 3,352,320 scores (1,920 signatures
by 1,746 cells). The documented reducer contract was applied in memory for
efficiency and matched a completed native round-1 reducer exactly for labels,
abstention reasons, scores, runner-up scores, and margins.

## Result

- Overall cell-round agreement with the baseline was 71.2%.
- Median per-cell agreement was 76.7%.
- At an 80% agreement threshold, 859/1,746 cells (49.2%) were stable.
- Among the 1,416 baseline-assigned cells, 57.1% were stable.
- The strongest abundant labels were `B.FrE.BM` (93.8% stable cells),
  `proB.FrBC.BM` (76.0%), and `proB.CLP.BM` (59.6%).
- The weakest abundant labels were `B.T3.Sp` (8.1%), `STHSC.150-.BM`
  (11.7%), `MPP3.48+.BM` (18.4%), and `LTHSC.34+.BM` (22.0%).

The ATAC UMAP localizes most instability to the adjacent stem/progenitor states
on the left manifold. The B-cell-heavy right manifold is substantially more
stable. Most changes are local biological alternatives or abstentions, such as
STHSC to MPP4 and MPP3 to MPP4, rather than broad lineage jumps.

Eight singleton ImmGen populations cannot have target-replicate stability
assessed by within-population resampling. They remain explicitly structurally
unassessable even though resampling of competitor populations can change their
selected genes and scores.

This result does not replace the provisional default. It supports carrying
bootstrap agreement as an uncertainty field in subsequent exploratory use and
avoiding hard interpretation of unstable stem/progenitor sublabels. No artifact
or label was promoted and the source bundle fingerprint was unchanged.
