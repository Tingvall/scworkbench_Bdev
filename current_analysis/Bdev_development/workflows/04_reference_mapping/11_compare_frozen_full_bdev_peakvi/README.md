# Frozen full-Bdev PeakVI bulk comparison

This workflow corrects the selection logic from development step 10. Same-Bdev
cluster recovery is retained as a diagnostic, but it cannot eliminate a
candidate before the intended external bulk test. In particular, the cluster
labels used for feature selection make held-out Bdev cluster recovery an
optimistic measure of the cluster-directed axes.

The freeze stage fits four independent PeakVI references on all 26,723 Bdev
cells: the 10k/10k hybrid, cluster-only 20k, depth-adjusted overall 20k, and all
191,433 peaks. It has no ImmGen path or library identity in its configuration;
only 28 anonymous query slots are reserved in the derived bundle. Each model,
feature order, reference latent space, and model payload checksum is recorded
before external data are opened.

The external stage then applies the identical already-processed 28-library
ImmGen matrix to every frozen model without refitting. It compares replicate
coherence, reference-neighbor use, hub concentration, and latent association
with bulk depth and peak saturation. Whole-reference cluster mapping is the
first output. Native-L4 neighborhood weight and pseudotime are diagnostic
second-stage outputs and are not accepted trajectory assignments here.

Run `01_install_runtime.py` and `02_freeze_full_bdev_candidates.R` with config
17, then run `03_project_immgen_all_candidates.R` and `04_finalize_run.py` with
config 18. No candidate is automatically promoted from this benchmark.

## Observed outcome

Run `full_bdev_peakvi_immgen_4way_20261116T202157Z` completed all four
full-Bdev fits and post-freeze ImmGen projections. Every model payload checksum
was identical before and after projection. The reduced axes all improved the
external bulk comparison over all 191,433 peaks: their median query/reference
distance ratios were 1.56--1.88 rather than 4.37, and their within-population to
between-population latent-distance ratios were 0.176--0.210 rather than 0.332.

Overall20k had the strongest bulk geometry and neighbor use (distance ratio
1.565, population-distance ratio 0.176, perfect replicate cluster agreement,
477 reference neighbors, and maximum hub fraction 0.0107), but it also had the
largest maximum latent association with peak saturation (0.736). Hybrid20k had
the lowest maximum saturation association (0.429) and lowest maximum depth
association (0.470), at the cost of weaker replicate cluster agreement (0.692)
and a larger distance ratio (1.881). Cluster20k lay between those choices. The
complete peak model was weakest for bulk mapping and replicate coherence.

Thus the external result supports informative-peak reduction, but leaves a
real multi-objective choice between overall20k for the strongest bulk mapping
geometry and hybrid20k for the most saturation-resistant representation. L4
outputs remain diagnostic and unaccepted.

## Subsequent development-default decision

After reviewing this benchmark and the separate fixed-ATAC-UMAP comparison,
the user selected overall20k (the depth-adjusted variable 20k model) as the
provisional Bdev development default. Cluster20k is retained as the required
sensitivity model, not as a coordinate-averaging ensemble. The explicit
selection, exact artifact checksums, limitations, and implementation boundary
are recorded in config 20 and the corresponding development-default document.
That decision was subsequently superseded on 2026-11-17 by the user's explicit
selection of the official scwbATAC depth-adjusted HVP 20k model. Config 25 is
the active development-default declaration; `overall20k` remains its required
sensitivity comparator.
