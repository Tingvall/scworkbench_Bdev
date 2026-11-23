# ImmGen native-L4 and fixed-UMAP driver diagnostic

Run `immgen_l4_umap_drivers_20261115T233506Z` exposes the reference cells that
drive the rejected ImmGen fixed-`umap.atac` overlay. It is a development
diagnostic only; it does not change the scientific disposition of Step 4.

## Exact reconstruction

Both reported outputs were independently reconstructed from the primary 30
PeakVI neighbors using the kernel
`exp(-(distance / kth_distance)^2)`:

- native L4 is the weighted median of the neighbors' authenticated native-L4
  values;
- fixed-UMAP location is the weighted arithmetic mean of those same neighbors'
  pre-existing `umap.atac` coordinates.

The maximum absolute L4 reconstruction error was `4.62e-14`, and the maximum
Euclidean coordinate reconstruction error was `7.91e-15`. The driver ledger
contains 230 unique Bdev reference cells. The Bdev source manifest remained
unchanged.

## Population findings

`B.FrE.BM` is the largest replicate disagreement. Its supported library maps to
L4 13.90 at `(4.59, -5.41)`, while its borderline library maps to L4 23.34 at
`(-3.26, -4.08)`. Its pooled neighbor distribution is consequently broad and
multimodal; a single population location is not representative.

Other notable within-population L4 ranges are `proB.FrBC.BM` 3.47, `B.Sp` 3.40,
`B.T3.Sp` 2.35, and `B.T2.Sp` 2.32. In contrast, the two `B.MZ.Sp` libraries
are nearly identical (range 0.04).

The five single-library stem/progenitor populations—`MPP3.48+.BM`,
`LTHSC.34+.BM`, `LTHSC.34-.BM`, `MPP4.135+.BM`, and `STHSC.150-.BM`—have
weighted neighbor medians of only 26.05–26.15 and occupy nearly the same
lower-left UMAP region. Their driver composition is also very similar and is
dominated by Bdev numeric clusters 6, 3, and 0. The mapping therefore provides
little discrimination among these biologically distinct populations.

Many mature/transitional populations are also dominated by Bdev cluster 6:
84% of `B.Fo.Sp`, 85% of `B1b.PC`, 76% of `B.T3.Sp`, 69% of `B.Sp`, 67% of
`B.T2.Sp`, and 58% of `B.MZ.Sp` population-normalized neighbor weight. This
shared neighbor pull explains why their fixed-UMAP placements concentrate on
the same lower trajectory arc. It does not establish that the placement is
biologically correct.

## Interpretation boundary

L4 does not determine either UMAP coordinate. PeakVI distance selects and
weights the Bdev neighbors. Their existing UMAP coordinates determine the
displayed barycenter, while their authenticated L4 values separately determine
the mapped weighted median. A bulk library can therefore have a plausible L4
summary while being pulled between disconnected or broad UMAP neighborhoods.

The reproducible definition is
`config/04_reference_mapping/12_immgen_l4_umap_driver_diagnostic.yaml`; the
numbered actions are under
`workflows/04_reference_mapping/06_diagnose_immgen_l4_umap_drivers/`.
