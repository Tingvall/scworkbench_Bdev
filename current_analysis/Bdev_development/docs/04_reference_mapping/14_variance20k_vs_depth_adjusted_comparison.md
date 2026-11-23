# ATAC raw-variance 20k versus depth-adjusted variable 20k

**Historical status:** this comparison supported the then-current `overall20k`
default. An explicit user decision on 2026-11-17 subsequently selected the
official `scwbatac::hvp_select` depth-adjusted HVP 20k model as the development
default. The results below remain valid comparison evidence.

## Decision

Retain depth-adjusted `overall20k` as the provisional Bdev development default.
Preserve raw-count `variance20k` as a non-depth-adjusted sensitivity model; do
not promote it automatically.

This is a development-stage decision only. It does not unlock implementation,
accept a bulk trajectory assignment, or replace the requirement for a newly
qualified immutable SIF.

## Comparison contract

Run identity: `variance20k_vs_depth_adjusted_20261117T031433Z`.

- The non-depth-adjusted axis used the existing `scwbrna::hvg_select`
  operation with `method = variance`, `feature_type = peak`, and
  `n_features = 20000` on `ATAC_counts`.
- All 191,433 peaks and all 26,723 Bdev cells were eligible. No batch/sample
  consistency constraint was applied.
- The full-Bdev PeakVI model used the same architecture, training cap, seed,
  and 10-dimensional latent size as `overall20k`.
- Feature selection and model fitting were completed and checksummed before
  the external ImmGen bundle was opened.
- The same 28 ImmGen libraries were projected without refitting. The
  `overall20k` comparison evidence was checksum-pinned to step 11.
- Whole-reference mapping was evaluated first. Native-L4 results remain a
  diagnostic, unaccepted second-stage trajectory output.

The raw-variance and depth-adjusted axes share 19,174 peaks. Each axis contains
826 peaks absent from the other, giving a Jaccard overlap of 0.921.

## External bulk-mapping results

| Diagnostic | Raw variance 20k | Depth-adjusted 20k | Preferred direction |
|---|---:|---:|---|
| Within/between population latent distance | 0.1666 | 0.1761 | Lower |
| Replicate predicted-cluster agreement | 0.846 | 1.000 | Higher |
| Unique reference neighbors | 494 | 477 | Higher |
| Maximum reference-hub fraction | 0.0131 | 0.0107 | Lower |
| Median support-distance ratio | 1.4379 | 1.5650 | Lower |
| Median neighbor-sample entropy | 0.2005 | 0.3905 | Higher |
| Maximum absolute depth-latent Spearman | 0.4641 | 0.4915 | Lower |
| Maximum absolute saturation-latent Spearman | 0.8083 | 0.7361 | Lower |

Raw variance improves the two latent-geometry diagnostics, uses more unique
reference neighbors, and modestly lowers the maximum depth association. Those
gains do not translate into a uniformly more robust mapping: raw variance loses
replicate coherence, concentrates the neighbor samples more strongly, has a
higher maximum hub fraction, and is more associated with peak saturation.

The raw-variance replicate failures occur for `B.FrE.BM` and `B.T2.Sp`; the
depth-adjusted model has no population with discordant replicate cluster calls.
Seven of 28 individual library cluster calls differ between the models. On the
fixed ATAC UMAP, the largest paired changes include the two `B.Sp` libraries,
`B1b.PC`, and one `B.T2.Sp` library. The overall biological layout remains
broadly concordant.

## Interpretation

For the intended bulk-mapping use, depth adjustment remains the better default.
The modest raw-variance geometry improvement is outweighed by its loss of
replicate stability and its stronger saturation association. This is also the
less biased default with respect to library depth/saturation because it avoids
selecting peaks solely because their unadjusted counts vary strongly.

The raw-variance model is still useful as a sensitivity check. Agreement with
`overall20k` supports a stable whole-reference placement; disagreement should
be reported as feature-selection uncertainty rather than resolved by averaging
latent coordinates or forcing a trajectory assignment.

## Traceable evidence

Generated results root:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/04_reference_mapping/13_compare_variance20k_peakvi/variance20k_vs_depth_adjusted_20261117T031433Z/
```

Key outputs:

- `01_tables/11_variance_vs_depth_adjusted_metrics.tsv`, SHA-256
  `67e2aef52cfb6a67f46361ca470a794cbbb9de28d61dde27ce6147931afbbb4b`;
- `02_figures/01_variance_vs_depth_adjusted_fixed_atac_umap.png`, SHA-256
  `94e271092f9da9cc53b4777e74be2bf225504bfde40ca5115f0163d9985473d8`;
- `02_figures/04_external_metric_comparison.png`, SHA-256
  `ae1e6fbe59824bd39ca52d416c512317786de46d70db77eafa6627a9b987095c`;
- `06_manifests/05_run_file_inventory.json`, SHA-256
  `947093e9e0c28a518f0d8e03ec00ba6a2ddc917ced9132fab6553192ec57ddab`.

No `scworkbench`, PyBench, `scwbRNA`, or `scwbLatent` source was modified.
The existing ATAC-compatible variance mode was sufficient.
