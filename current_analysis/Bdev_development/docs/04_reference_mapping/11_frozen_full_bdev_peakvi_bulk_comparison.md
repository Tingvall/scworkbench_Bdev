# Frozen full-Bdev PeakVI bulk comparison contract

Development step 10 answered a same-domain question and therefore cannot by
itself choose the peak axis for bulk transfer. Its cluster-directed features
were selected using Bdev cluster labels from the complete reference, so better
Bdev cluster recovery by those features is expected and is not an independent
external test.

Step 11 freezes four models before external data access:

1. hybrid20k: 10,000 cluster-discriminating plus 10,000 additional
   depth-adjusted overall-variable peaks;
2. cluster20k: 20,000 cluster-discriminating peaks;
3. overall20k: 20,000 depth-adjusted overall-variable peaks;
4. full191433: every source ATAC peak.

Every model uses all 26,723 Bdev cells, the same ten-dimensional PeakVI
architecture and seed, and an exact source-ordered feature axis. The freeze
receipt records the model artifact, payload checksum, feature-order identity,
reference latent checksum, runtime, and source bundle checksum. Config 17 has
no ImmGen data path or library identity. The bundle reserves 28 anonymous slots
only so query matrices can be attached after the freeze boundary.

Config 18 then exposes the external stage to the already-processed ImmGen bulk
matrix and its identities. The same 28 libraries are projected into all four
models without refitting. Comparison emphasizes population-replicate coherence,
neighbor diversity and hub concentration, and latent association with common
bulk depth and candidate-axis peak saturation. Cluster mapping against the
whole Bdev reference is primary. Native-L4 neighborhood weight and pseudotime
are retained only as diagnostic, unaccepted second-stage outputs.

This benchmark does not automatically promote a model. Any selection is a
separate scientific decision made after reviewing the frozen external results;
it cannot trigger feature changes or refitting within this run.

## Completed benchmark

The four frozen models and their external ImmGen comparison completed under run
`full_bdev_peakvi_immgen_4way_20261116T202157Z`. The result makes the limitation
of the earlier same-Bdev gate explicit: all three informative 20k axes mapped
the bulk queries closer to the Bdev reference, separated population replicates
more coherently, and used reference neighbors less concentratively than the
all-peak model.

Overall20k was best on the external mapping geometry, while hybrid20k was best
at limiting the maximum association of latent coordinates with bulk depth and
selected-peak saturation. Consequently the benchmark does not justify choosing
the all-peak model, but it also does not collapse mapping quality and saturation
robustness into one post hoc score. The raw comparison is retained for an
explicit scientific choice between those objectives.

The subsequent explicit decision selected overall20k as the provisional
development default for whole-reference bulk mapping. Cluster20k is retained
as a required sensitivity model, and the saturation association of overall20k
remains an explicit limitation. The decision is recorded separately in
`20_depth_adjusted_variable20k_default.yaml` and
`13_depth_adjusted_variable20k_development_default.md`; it does not alter this
benchmark or retroactively turn its comparison into automatic promotion.
That decision was later superseded on 2026-11-17 by the explicit selection of
the official scwbATAC depth-adjusted HVP 20k model. Config 25 is the active
declaration, and overall20k is now its required sensitivity comparator.
