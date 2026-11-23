# Plot the bounded FLT3LKO/IL7KO fixture on the fixed ATAC UMAP

This development visualization places all 28 balanced projection-fixture cells
on the unchanged 26,723-cell Bdev `umap.atac`. Each query position is the
distance-kernel-weighted barycenter of its 30 exact PeakVI neighbors in the
authenticated frozen-model joint space. It is not a UMAP transform, and it
does not refit the UMAP, PeakVI model, or native Slingshot L4 trajectory.

The four facets contain seven cells each for Flt3l WT, Flt3l KO, Il7 WT, and
Il7 KO. Cyan-edged circles have finite transferred native L4. Magenta crosses
are supported mappings whose reference neighborhoods did not yield finite L4;
they remain visible and are not imputed.

This is a bounded diagnostic only. It does not estimate a population
distribution, compare genotypes, or unlock the full dataset or implementation.
