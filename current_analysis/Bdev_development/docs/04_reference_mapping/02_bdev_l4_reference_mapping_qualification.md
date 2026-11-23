# Bdev L4 reference-mapping qualification

## Outcome

Run `bdev_l4_reference_mapping_20261110T042541Z` completed the bounded
reference-mapping semantics qualification without building a SIF or creating
implementation-stage outputs. The run exercised the pinned
`scwbrepresentation::knn_neighbors`,
`scwbreferencemap::reference_map_calibrate`, and
`scwbreferencemap::reference_map_transfer` operations.

All 1,810 cells in the authenticated 03b native-L4 fixture were evaluated once
as query cells. Each of the four `orig.ident` samples was held out in turn, and
only the other three samples supplied reference neighbors. Mapping used ATAC
LSI dimensions 2--30; dimension 1 was excluded. Calibration used k values 15,
30, and 50, primary k 30, and q95/q99 support thresholds. Exact neighbor search
used the RANN backend at commit
`c49da6aea8fcb8ace7a2ee36d51dc754d1fdac78`.

## Contract acceptance

The following contract checks passed:

- all four samples were held out exactly once;
- all 1,810 fixture cells had exactly one query result;
- categorical and native-L4 values were withheld for every out-of-reference
  query;
- calibration models were retained and consumed by transfer;
- native fitted Slingshot L4 was used without reorientation;
- the active 03c source manifest remained exactly
  `cc8b92e8b55bb1b12b7f818c3bc272193c2139abec81b9777d48bfe2dd022a05`;
- no outputs were promoted into the 03c source; and
- no replacement SIF was built and `Bdev_implementation` remained locked.

The final fresh fixture-bundle manifest is
`42c8ba5189ff89d9c3ad294136b6baae0931a8819f63725404a49447405a9e1a`.
The checksummed run inventory contains 18 files.

## Scientific result and decision

The LSI qualification is a successful operation-contract test but is not an
accepted biological projection recipe.

| Metric | Aggregate result |
|---|---:|
| accepted query fraction | 0.9862 |
| cluster accuracy among accepted queries | 0.4762 |
| native-L4 MAE among accepted queries | 5.5252 |
| native-L4 RMSE among accepted queries | 9.3291 |
| native-L4 Spearman among accepted queries | 0.2712 |

The strongest warning is sample 1: its accepted-query native-L4 Spearman was
`-0.5809` and its cluster accuracy was `0.0564`. Samples 2 and 3 performed much
better, showing that the operation can preserve useful structure, but the
cross-sample instability rules out adopting direct LSI nearest-neighbor mapping
for the bulk-ATAC analysis.

The next Step 4 gate is therefore the planned Bdev-only PeakVI validation with
depth-matched Bdev pseudobulks/pseudo-cells. Real ImmGen bulk libraries must not
be projected until that validation establishes acceptable recovery, support,
and out-of-reference behavior.

## ImmGen ATAC preflight

The preflight froze 16 selected B/Stem&Prog populations covering 28 GSE100738
ATAC libraries. The exact tables are:

- `01_tables/09_selected_immgen_atac_populations.tsv`; and
- `01_tables/10_selected_immgen_atac_libraries.tsv`.

GSE100738 is mm10. Its population-level normalized bigWigs are retained only
for provenance or visualization; they are not count matrices and are not valid
PeakVI inputs. Future model input requires the raw SRA libraries to be
reprocessed and counted against the exact Bdev peak universe.

## Infrastructure observation

Calibration promotion on the migrated legacy-derived execution clone failed
because its `artifact_model_ref` table lacks the current `linkage_json` column.
No core edit or direct database repair was attempted. The completed
qualification instead used a fresh fixture-only scWorkbench bundle containing
the real 10-feature chromVAR axis slice, real LSI coordinates, authenticated
metadata, and native L4. That fresh bundle supported ordinary retained-model
promotion and transfer.

This observation should be resolved or explicitly waived before a runtime is
accepted if migrated bundles must support current model promotion. It does not
authorize a scWorkbench core change.
