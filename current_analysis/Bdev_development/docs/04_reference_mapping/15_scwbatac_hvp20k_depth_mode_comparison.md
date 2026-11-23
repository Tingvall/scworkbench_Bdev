# Official scwbATAC HVP 20k depth-mode comparison

## Outcome

Depth-adjusted HVP is preferable to raw-binary HVP for this bulk-mapping task,
although this benchmark alone did not justify replacing the then-current
depth-bin `overall20k` development default.
The official selector also does not yet pass real-data resource qualification
because its sparse input is coerced to a 2.9 GiB dense object per 2,048-cell
chunk during the depth cross-product.

Run identity: `hvp20k_depth_comparison_20261117T044105Z`.

## Subsequent explicit default selection

After reviewing this evidence, the user explicitly selected the official
`scwbatac::hvp_select` depth-adjusted HVP 20k candidate as the provisional
development default on 2026-11-17. This is a deliberate method/default choice,
not an automatic benchmark winner. `overall20k` is retained as the required
sensitivity comparator because it remains stronger on several external
mapping diagnostics. The selector's resource-qualification failure remains
open, so the selection does not qualify a new SIF or unlock implementation.

## Frozen comparison

Both models used all 26,723 Bdev cells, exactly 20,000 source-ordered peaks,
the same PeakVI architecture and seed, and 28 anonymous external slots. Models
were frozen and checksummed before ImmGen access. The same 28 libraries were
then projected without refitting and placed on the fixed Bdev ATAC UMAP from
their frozen 30-neighbor ledgers.

The depth-adjusted selector used residual binary accessibility variance after
linear adjustment for `log1p(detected peaks per cell)`. The raw axis used the
operation's authenticated `raw_variance` ledger and exact documented raw-mode
ordering. This reranking is numerically equivalent to
`depth_adjustment = FALSE`; it avoided repeating an unnecessary depth
cross-product after the real-data dense-coercion defect was observed.

The two HVP axes shared 19,948 peaks and differed by 52 peaks each (Jaccard
0.9948). Adjusted HVP shared 19,873 peaks with the existing `overall20k` axis;
raw HVP shared 19,833.

## Mapping diagnostics

| Diagnostic | HVP depth-adjusted | HVP raw binary | Existing overall20k |
|---|---:|---:|---:|
| Within/between population distance | 0.1873 | 0.1980 | 0.1761 |
| Replicate cluster agreement | 0.9231 | 0.9231 | 1.0000 |
| Unique reference neighbors | 508 | 474 | 477 |
| Maximum reference-hub fraction | 0.0107 | 0.0131 | 0.0107 |
| Median support-distance ratio | 1.6010 | 1.8295 | 1.5650 |
| Median neighbor-sample entropy | 0.2841 | 0.3463 | 0.3905 |
| Maximum absolute depth-latent Spearman | 0.5424 | 0.4576 | 0.4915 |
| Maximum absolute saturation-latent Spearman | 0.6691 | 0.7123 | 0.7361 |

Depth adjustment improves five of the eight paired HVP diagnostics: population
geometry, reference-neighbor breadth, hub use, support distance, and saturation
association. Raw HVP improves neighbor-sample entropy and depth association;
replicate agreement is tied. Both HVP modes split the two `B.FrE.BM`
replicates across clusters 2 and 15. Their individual cluster calls differ for
only two stem/progenitor libraries.

Against the then-current default, adjusted HVP has lower saturation association and
more unique reference neighbors, but worse population separation, replicate
agreement, support distance, sample entropy, and depth association. That did
not justify an automatic default change from the benchmark alone.

## Operation-pack qualification finding

The source contract test passed, but it did not expose the real-matrix resource
behavior. On the 191,433-by-26,723 Bdev matrix, each 2,048-cell chunk warned of
a 2.9 GiB sparse-to-dense coercion. The implementation computes the depth
cross-product even when `depth_adjustment = FALSE`, so running raw mode would
repeat the same avoidable cost.

Before the operation is accepted in a new immutable SIF, `scwbATAC` should:

1. implement the binary peak-by-log-depth cross-product with a sparse-safe
   Matrix method;
2. skip all depth cross-products when depth adjustment is disabled;
3. add a sparse large-axis contract or integration fixture that rejects dense
   coercion; and
4. repeat both modes on the Bdev matrix and reproduce the selected axes and
   downstream model results.

This is an operation-pack change, not a scWorkbench or PyBench core change.

## Evidence

Results root:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/04_reference_mapping/14_compare_hvp20k_depth_modes/hvp20k_depth_comparison_20261117T044105Z/
```

- `01_tables/08_hvp_external_metrics.tsv`: SHA-256
  `5be8f57eb31b8e1f4f1d738ba89ce92e12eb72f92471a03e703908db1c8f3205`;
- `02_figures/01_hvp_depth_vs_raw_fixed_atac_umap.png`: SHA-256
  `e6dd1e6bab759a8a3e0fc50facda9a439d3f20c717f092ebea0513812df439e1`;
- `02_figures/04_hvp_external_metric_comparison.png`: SHA-256
  `91967d8774bcdbcd8a05ad3733c25010d38fbb6696acde92af0efab4009eada1`;
- `06_manifests/05_run_inventory.json`: SHA-256
  `342acb12dc299f341aef945c18f8b6d9009e40e0da62b575d7139ac316298238`.

The benchmark did not change the default automatically. The later explicit
selection is declared separately in config 25. Implementation remains locked,
and no new SIF was built.
