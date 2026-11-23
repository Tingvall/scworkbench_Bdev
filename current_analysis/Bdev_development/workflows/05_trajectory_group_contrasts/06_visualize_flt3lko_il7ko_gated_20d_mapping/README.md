# Visualize the accepted full gated 20D mappings

This development-only visualization consumes the finalized Step 05 mapping
table and its 47 checksummed neighbor checkpoints. It colors only cells with
`accepted == TRUE` and finite `accepted_native_l4`; all other query cells are
shown only as a neutral background on the native query UMAPs.

The workflow produces:

1. the authenticated Bdev L4 reference and accepted query placements on the
   unchanged Bdev `umap.atac`;
2. accepted query cells on the source Cell Ranger ARC ATAC UMAP;
3. accepted query cells on the source Cell Ranger ARC GEX UMAP; and
4. descriptive accepted native-L4 distributions by library.

The fixed-reference placements use the same 30-neighbor distance-kernel
barycenter contract as the prior development visualization. No PeakVI, UMAP,
or Slingshot model is refit. The source UMAP CSV files are read from their
declared archived paths and are not modified.

Run `01_build_mapping_figures.R --config <config>` inside the configured SIF
with the configured R library.

Action `02_build_all_cells_fixed_umap.R` creates a separately identified
follow-up view in which every query cell is displayed on the fixed Bdev ATAC
UMAP. Complete-gate rejected cells are dark grey and accepted cells are drawn
on top with the native-L4 color scale. It consumes the checksummed placement
table from action 01 and does not overwrite the accepted-only run.

## Result

Run `gated_20d_mapping_visualization_20261122T235003Z` completed and passed.
It aligned all 23,252 mapped query cells to both source query UMAPs and
colored the 2,974 complete-gate accepted cells by native L4 pseudotime. Four
numbered PNG/PDF figure pairs, the combined coordinate/acceptance table, and
the library-level accepted-pseudotime summary were checksummed in the final
inventory. No source bundle, archived input, UMAP, or Slingshot fit changed.
