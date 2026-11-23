# Earlier-step working-state cleanup

Cleanup authorized on 2026-11-17 for abandoned and recomputable working state
from steps 1 and 3 of `05_trajectory_group_contrasts`. Qualified bundles,
results, audits, and execution receipts are retained.

## Removed paths

| Step | Relative path below `intermediate_objects/05_trajectory_group_contrasts/` | Allocated bytes | Reason |
|---|---|---:|---|
| 1 | `01_qualify_flt3lko_il7ko_mapping_fixture/.flt3lko_il7ko_mapping_fixture_20261117T174408Z.tmp-4114292/` | 544,768 | Abandoned temporary run directory. |
| 1 | `01_qualify_flt3lko_il7ko_mapping_fixture/flt3lko_il7ko_mapping_fixture_20261117T175427Z/05_objects/03_sessions/` | 204,853,248 | Recomputable projection session workspace; promoted artifacts retained. |
| 1 | `01_qualify_flt3lko_il7ko_mapping_fixture/flt3lko_il7ko_mapping_fixture_20261117T175427Z/04_logs/sessions/` | 478,949,376 | Recomputable calibration and transfer session workspaces; qualified outputs retained. |
| 3 | `03_plot_flt3lko_il7ko_fixture_fixed_umap/flt3lko_il7ko_fixture_fixed_umap_20261117T190116Z/04_logs/01_sessions/` | 4,096 | Empty session directory. |

Total allocated size before removal: 684,351,488 bytes.

## Explicitly retained

- step 1 qualified projection bundle, manifest SHA-256
  `1b1935f7af3721a29fcf748a01c077a5b4628bbeae74cdb0ba1aa0e89a37a9d0`
- step 3 completed visualization execution bundle, manifest SHA-256
  `97ef2bbdfab6e344267b318083a135601d55feb014a365ed71b9625014533c89`
- finalized quantified matrices, tables, figures, reports, manifests, and
  repository definitions

## Execution

All four approved paths were removed on 2026-11-17. The retained bundle
manifest checksums were reverified after removal. Post-cleanup free space on
`/mnt/sata2` was 7.0 GiB.
