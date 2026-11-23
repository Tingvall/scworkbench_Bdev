# Reference-mapping working-state cleanup

Cleanup authorized on 2026-11-17 for large reproducible working state that is
not consumed by the current full FLT3LKO/IL7KO mapping. All repository
definitions and results-side tables, figures, reports, receipts, and
inventories are retained.

## Removed paths

Paths are relative to
`intermediate_objects/04_reference_mapping/`.

| Relative path | Allocated bytes | Reason |
|---|---:|---|
| `14_compare_hvp20k_depth_modes/hvp20k_depth_comparison_20261117T044105Z/05_objects/02_selection_bundle/` | 13,164,937,216 | Reproducible selection clone; selected frozen model bundle retained. |
| `01_qualify_bdev_l4_reference_mapping/bdev_l4_reference_mapping_20261110T042541Z/05_objects/04_sessions/` | 13,804,765,184 | Recomputable qualification sessions; execution and fixture bundles plus results retained. |
| `15_requalify_scwbatac_hvp_sparse/scwbatac_hvp_sparse_source_requalification_20261118T004129Z/05_objects/02_sessions/` | 6,901,837,824 | Recomputable requalification sessions; ledgers, reports, receipts, and adjudication retained. |
| `03_map_immgen_atac_trajectory/immgen_atac_l4_20261110T202521Z/05_objects/03_raw_fastq/` | 16,765,227,008 | Re-downloadable preprocessing working input; ENA resolution and preprocessing receipts retained. |

Total allocated size before removal: 50,636,767,232 bytes.

## Critical retained artifacts

- selected HVP-depth20k frozen bundle manifest SHA-256:
  `1f9e399600a9ab74b6037301af8a6b521bf4669b8cfd4b76f80f2d07ba437b98`
- HVP external execution receipt SHA-256:
  `2c69e2c9395a3bd611f0121f523d98a0077b67b8bc9ab16549d9568170dab0c2`
- sparse requalification summary SHA-256:
  `36cf4b19ad382a705705d96b60ec437e980831b7055b90f9ca4ad0591ddca25a`
- ImmGen ENA resolution, peak-axis export, preprocessing receipt, execution
  receipt, and run-file inventory remain in the results manifest directory;
  the run-file inventory SHA-256 is
  `caacd244497cc3aed180cd9129683b725084dc1ba367acab69c1ef8a05c1a13f`

The selected frozen bundle, exact feature-axis tables, corrected sparse-core
qualification evidence, canonical cytokine bundle, and all inputs required by
the later implementation remain unchanged.

## Execution

All four approved working-state paths were removed on 2026-11-17. Retained
checksums were reverified, and every input named by the current full-mapping
workflow was confirmed present. Combined with the cytokine cleanup, available
space on `/mnt/sata2` increased from 7.0 GiB to 121 GiB.
