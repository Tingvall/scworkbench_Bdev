# 03c Depth-adjusted Parse cytokine trajectory qualification

**Status:** bounded Bdev development qualification passed; immutable SIF build
and acceptance still pending

## Qualified runtime

Run `bdev_cytokine_depth_adjusted_20261106T145817Z` used these development
sources in an isolated run library:

| Component | Version | Commit or loaded manifest identity |
|---|---:|---|
| scworkbench | 0.1.1 | `e6289d5353fbdb1a0e14a7d3a806ce4b958f56ff` |
| scwbRNA | 0.1.0 | `151dd29211b219343dce6247c5e48fc932b4c83d` |
| scwbTrajectory | 0.1.8 | `f5b6863afc99ca7118b63fe915b759f9e38f6d2b` |
| scwbTrajectory loaded manifest | — | `cc66ffd389761a115d24788be77ad0b0ba8655a42577d413fdb5466642e18d9c` |

The core and trajectory updates authenticate a multi-feature-set score output
as an exact ordered feature-set-by-effective-cell Cartesian product. This is
the infrastructure contract that blocked the preceding failed 03c attempts.
All signaling, depth-adjustment, interval, GSEA, Slingshot, and publication
contracts passed in the source qualification sweep. The broader
scwbTrajectory suite passed 111/112 contracts; its sole failure was the
pre-existing Palantir numeric-oracle comparison, which was unchanged by the
03c fix commits and is outside this qualification. The pack therefore has
complete 03c evidence but does not yet have a completely green pack-wide test
claim.

## Bdev qualification

The run reused the checksummed fitted Slingshot `L4` pseudotime, contained
model, and exact deterministic fixture from 03b. It did not refit Slingshot.
Analysis retained native L4 orientation; figures alone display
`LSK (0) -> CD19 (1)` as one minus the native-L4 fraction.

The fixture contained 1,810 cells in 20 bins across four `orig.ident` samples,
selected from 13,451 cells with finite L4 pseudotime. RNA-count library sizes
matched `nCount_RNA` exactly. `scwbrna::matrix_log_normalize` created
`RNA_log1p_cp10k` as `log1p(counts / library_size * 10000)`.

All 38 frozen Parse38 signatures passed overlap screening on both
`RNA_log1p_cp10k` and `SCT_data`. The reusable score operation adjusted each
gene for centered `log1p(nCount_RNA)` while preserving an
`ns(L4, df=8)` component, standardized the adjusted gene values, and averaged
genes with equal weight. All 753 fitted genes were adjusted and standardized.
The authenticated output contains 68,780 values: 38 signatures by 1,810 cells.

The matched downstream methods produced:

| Method | Configuration | Result size |
|---|---|---:|
| Dynamics | 20 median bins | 38 summaries; 722 curve rows |
| Adjusted-score intervals | width 0.05, step 0.025, 99 unblocked permutations | 72 supported intervals; 32/38 global tests supported |
| SCT preranked GSEA | width 0.05, step 0.025, 99 unblocked permutations, no added depth fit | 38 global results; 1,406 local results; 46 supported intervals |

Every interval adjustment record was `already_adjusted`, so depth was not
removed twice. Every GSEA fit record retained `adjust_depth = FALSE` and
`permutation_mode = unblocked`. fgsea emitted non-fatal tied-rank warnings and
one warning that a true P value was below the configured `eps = 1e-8` floor;
the warnings are preserved in the run log.

## Preservation and output audit

The archived source bundle was not modified. Preserved `pt` values and
missingness had identical before/after digests. No result was promoted to the
source bundle or to `Bdev_implementation`.

The run produced five PNG/PDF figure pairs: dynamics, interval, and local-GSEA
overview figures plus focused FLT3L-up and IL-7-up four-view summaries. All
five rendered successfully. The final inventory contains 53 files, and an
independent audit reproduced every recorded size and SHA-256.

The inventory is stored at:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/03c_parse_cytokine_depth_adjusted_trajectory/01_qualify_depth_adjusted_cytokine_trajectory/bdev_cytokine_depth_adjusted_20261106T145817Z/06_manifests/03_run_file_inventory.json
```

Its SHA-256 is
`f5f0c4804c9d144ee29584fd09b9359cd9c7542d21535c24c1bb50147d329452`.

These outputs are bounded development evidence only. They do not unlock the
implementation stage. The accepted paper workflow must rerun from the frozen
inputs after the replacement immutable SIF has been built, qualified,
checksummed, and explicitly accepted.

## Post-qualification UMAP views

A separate checksummed visualization run created one ATAC-UMAP for each of the
21 `up` signatures. It colors only fitted-L4 fixture cells within supported
local-GSEA intervals by local NES and retains all Bdev cells as gray context.
Ten signatures have supported regions; eleven are explicit no-region maps.
See `02_up_geneset_gsea_umaps.md` for the exact result root, interpretation,
and inventory identity. The parent qualification run remained unchanged.
