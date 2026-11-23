# Map all FLT3LKO/IL7KO cells and display them on the fixed ATAC UMAP

This development workflow maps all 23,252 QC-selected cells from the four
FLT3LKO/IL7KO libraries. It reuses the qualified reciprocal liftOver map,
quantifies the immutable GRCm39 fragments on the exact frozen HVP-depth20k
axis, projects through the frozen PeakVI model without refitting, assigns
calibrated support and native L4, and places every query on the unchanged Bdev
`umap.atac` by a 30-neighbor distance-kernel barycenter.

The run is descriptive development analysis. It does not refit PeakVI, UMAP,
or Slingshot; it does not perform inferential WT/KO testing; and it does not
build a SIF or unlock implementation. A new SIF remains deferred until all
remaining development analyses are complete.

Run in order:

1. `01_prepare_full_inputs.py --config <config>`
2. `02_quantify_full_axis.R --config <config>`
3. install the pinned development core and packs into the isolated run library
4. `03_map_and_plot_fixed_umap.R --config <config>`
