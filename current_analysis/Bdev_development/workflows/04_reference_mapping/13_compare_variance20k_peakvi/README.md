# Non-depth-adjusted variance20k versus depth-adjusted overall20k

**Historical status:** this workflow retained the then-current `overall20k`
default. The user's later explicit decision selected the official scwbATAC
depth-adjusted HVP 20k model; config 25 is now the active declaration.

This development workflow tests the ATAC-appropriate non-depth-adjusted mode
of the existing `scwbrna::hvg_select` operation. The blind stage clones the
source bundle, selects the top 20,000 peaks with `method = "variance"`, creates
a source-ordered PeakVI axis, and freezes a full-Bdev model before any ImmGen
path or identity is available.

The external stage projects the same 28 already-processed ImmGen libraries
without refitting and compares the result with the checksum-pinned overall20k
evidence from step 11. The comparison uses the same geometry, replicate,
neighbor-use, depth, saturation, fixed-UMAP, and diagnostic native-L4 outputs.

This workflow does not test the RNA-oriented VST/dispersion modes, modify
`scwbRNA`, or automatically replace the then-current overall20k development
default. A default change requires a separate explicit decision after review.

## Observed outcome

Run `variance20k_vs_depth_adjusted_20261117T031433Z` completed. The authenticated
raw-count selector scored all 191,433 ATAC peaks across all 26,723 Bdev cells,
selected exactly 20,000 peaks, and froze the PeakVI model before opening the 28
ImmGen libraries. Projection used the frozen model without refitting.

The two axes share 19,174 peaks (95.87% of either 20k axis). Raw variance had a
lower within/between-population distance ratio and lower support-distance ratio,
but it reduced replicate cluster agreement from 1.0 to 0.846, increased maximum
saturation association from 0.736 to 0.808, increased reference-hub use, and
reduced neighbor-sample entropy. The fixed-UMAP comparison showed the largest
changes in mature B-cell mappings and replicate-discordant calls for `B.FrE.BM`
and `B.T2.Sp`.

The result supported retaining depth-adjusted `overall20k` as the default at
that time. Raw-count `variance20k` is preserved as a focused
non-depth-adjusted sensitivity model. No scWorkbench or operation-pack source
change was required.
