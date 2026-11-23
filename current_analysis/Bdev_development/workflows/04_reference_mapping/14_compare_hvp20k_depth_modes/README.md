# Official scwbATAC HVP 20k depth-mode comparison

**Decision update (2026-11-17):** after this benchmark, the user explicitly
selected `hvp_depth20k` as the provisional development default. Config 25 and
development-default document 16 record the active decision. The benchmark's
resource-qualification finding remains open.

This workflow runs `scwbatac::hvp_select` on all Bdev cells in two externally
blind modes: log-depth-adjusted residual binary variance and raw binary
variance. Each mode selects exactly 20,000 peaks after the same inclusive
detection filter, restores source peak order for PeakVI, and freezes a
full-Bdev model before ImmGen is opened.

The external stage projects the same 28 ImmGen libraries without refitting and
plots both frozen mappings on the same fixed Bdev ATAC UMAP with authenticated
native L4 shown only as a diagnostic second-stage trajectory reference. The
existing `overall20k` default is not replaced automatically.

## Observed outcome

Run `hvp20k_depth_comparison_20261117T044105Z` completed. Both 20k PeakVI
models were frozen before the 28 ImmGen libraries were opened, and projection
did not refit either model. The HVP axes shared 19,948 peaks and differed by 52
peaks each.

Depth-adjusted HVP was better than raw HVP on population separation, support
distance, unique-neighbor use, hub concentration, and saturation association;
replicate agreement was tied. Raw HVP had higher neighbor-sample entropy and a
lower maximum depth-latent association. The fixed-UMAP figure showed broader
mature-B support ellipses for raw HVP. Both modes had one discordant
`B.FrE.BM` replicate pair.

Compared with the then-current depth-bin `overall20k` default, adjusted HVP improved
saturation association and used more distinct neighbors, but had worse
population separation, replicate agreement, support distance, sample entropy,
and depth association. The benchmark therefore did not trigger an automatic
default replacement; the subsequent explicit selection did.

The real-data selector run exposed a pack implementation problem: every
2,048-cell chunk caused a 2.9 GiB sparse-to-dense coercion. The depth-adjusted
operation completed and its output was promoted into the frozen evaluation
bundle. Because its complete authenticated
ledger contains the raw binary variance for every peak, the exact raw-mode axis
was reranked from that ledger rather than repeating the same costly
cross-product. This is valid scientific comparison evidence but not a clean
qualification of the raw operation execution path. `scwbATAC::hvp_select`
requires a sparse-safe implementation update and repeat qualification before
container acceptance.
