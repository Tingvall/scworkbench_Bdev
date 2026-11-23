# ImmGen bulk mapping with the selected gated 20D reference

Run `immgen_bulk_gated_20d_20261123T001821Z` completed as a development-stage
scWorkbench validation of the selected later-implementation default. All 28
checksummed ImmGen ATAC libraries were projected through exact frozen model
`model_ref_20261122_194532_hka11fvs` without model refitting. The run reused
calibration `model_ref_20261122_195229_ewshintf`, the exact HVP-depth20k
feature order, primary `k = 30`, and the qualified support/native-L4 gate.

All 28 libraries were supported by the calibrated reference-distance limit;
none was borderline or out of reference. Four libraries passed the complete
L4-domain gate and received accepted native-L4 values:

- `GSM2692318`, `MPP4.135+.BM`: `26.1768611032085`;
- `GSM2692236`, `proB.FrA.BM`: `22.5991908103007`;
- `GSM2692238`, `proB.FrBC.BM`: `15.7253960116066`; and
- `GSM2692239`, `proB.FrBC.BM`: `15.8529915113917`.

The other 24 libraries were withheld because their local reference-domain vote
was `out_domain`, not because they failed the calibrated distance-support
limit. Their raw estimates, domain votes, confidences, and explicit rejection
reasons remain in the mapping table.

The figure uses the unchanged Bdev `umap.atac` coordinates and the native
fitted L4 curve. Accepted libraries are colored by accepted native L4 and all
withheld libraries are dark grey. Neither UMAP nor Slingshot was refitted.

Outputs are under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/04_reference_mapping/17_map_immgen_atac_gated_20d/immgen_bulk_gated_20d_20261123T001821Z/
```

Key files are:

- `01_tables/02_immgen_bulk_gated_20d_mapping.tsv`;
- `01_tables/03_immgen_bulk_gated_20d_neighbors.tsv.gz`;
- `01_tables/04_immgen_bulk_fixed_umap_placements.tsv`;
- `02_figures/01_immgen_bulk_gated_20d_fixed_atac_umap.{png,pdf}`;
- `03_reports/01_immgen_bulk_gated_20d_mapping.json`; and
- `06_manifests/03_run_inventory.json`.

The first plotting attempt established the projection, authenticated latent
union, and gated mapping but correctly stopped when the qualification clone did
not contain the visualization-only `umap.atac` space. The workflow then pinned
the checksummed upstream Bdev bundle as a separate read-only fixed-UMAP input
and resumed the existing projection and union. The selected model payload,
qualified source manifest, bulk source manifest, and fixed-UMAP source manifest
all remained unchanged. This run is development evidence only; implementation
and paper analysis remain locked pending the qualified immutable SIF.
