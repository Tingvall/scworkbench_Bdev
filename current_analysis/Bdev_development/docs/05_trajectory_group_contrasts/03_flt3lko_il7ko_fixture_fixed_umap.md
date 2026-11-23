# FLT3LKO/IL7KO bounded fixture on the fixed Bdev ATAC UMAP

**Run:** `flt3lko_il7ko_fixture_fixed_umap_20261117T190116Z`  
**Status:** bounded diagnostic visualization

The 28 balanced projection-fixture cells are displayed on the established Bdev
`umap.atac` without changing its 26,723 reference coordinates. Query positions
are the distance-kernel-weighted barycenters of their 30 exact PeakVI
neighbors. This is the same fixed-coordinate diagnostic convention used for
prior external mappings; it is not a learned or transformed query UMAP.

All cells remain visible whether or not native L4 is finite. The figure is not
a genotype comparison: seven deterministic fixture cells cannot estimate each
library's full distribution, and the source contains no independent biological
replicate field.

The completed overlay shows most finite-L4 fixture cells on the yellow/orange
later-L4 arm. The single low-L4 Flt3l-KO fixture cell lies near the upper
dark-purple L4 region. Most Il7 cells without finite L4 occupy the lower gray
spur, while the remaining missing-L4 cells appear off the fitted L4 curve in
the left compartment. All 28 placements passed the mapping-support gate; an
off-curve location therefore means that L4 is undefined, not that the mapping
was classified out of reference.

The checksummed PNG is
`02_figures/01_flt3lko_il7ko_fixture_fixed_atac_umap.png` under the run result
root. Its SHA-256 is
`57fbfa90c115a5405c2067ef213a3153e66a0844ba4f2ea697fb3998be52d5a4`.
