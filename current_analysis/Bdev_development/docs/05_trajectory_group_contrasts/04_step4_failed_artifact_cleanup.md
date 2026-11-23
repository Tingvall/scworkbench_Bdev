# Step 4 failed-artifact cleanup

Cleanup authorized on 2026-11-17 for failed and superseded working state from
`04_map_flt3lko_il7ko_full`. The qualified quantification and the validated
corrective baseline bundle are retained.

## Removed paths

| Run | Relative path below the step intermediate root | Allocated bytes | Reason |
|---|---|---:|---|
| `20261118T022951Z` | `04_logs/01_sessions/` | 41,967,616 | Failed durable-session workspace; results-side failure logs retained. |
| `20261118T022951Z` | `05_objects/02_ATAC_hvp_depth20k_counts.tmp-736175/` | 119,943,168 | Abandoned temporary quantification; published matrix retained. |
| `20261118T022951Z` | `05_objects/03_mapping_bundle/` | 947,073,024 | Superseded cross-bundle-model mapping attempt. |
| `20261118T031500Z` | `04_logs/01_sessions/` | 3,127,476,224 | Recomputable operation workspaces; promoted baseline artifacts retained in the validated bundle. |
| `20261118T031500Z` | `05_objects/04_mapping_bundle_primary_k_only/` | 2,350,469,120 | Incomplete primary-k clone stopped during ingestion; must be recreated cleanly. |

Total allocated size before removal: 6,586,929,152 bytes.

## Explicitly retained

- `20261118T022951Z/05_objects/02_ATAC_hvp_depth20k_counts/`
- `20261118T031500Z/05_objects/03_mapping_bundle/`
- all results-side audits, logs, checksums, and failure diagnostics
- repository configurations, inputs, and workflow scripts

The retained baseline bundle manifest SHA-256 is
`3e68752163dfa31fe28adef5fc9354ac62357e3eeafb63ec8ed0d58d99ae7ce6`,
matching the corrective configuration.

## Execution

All five approved paths were removed on 2026-11-17. Post-cleanup free space on
`/mnt/sata2` was 6.4 GiB. Removal is permanent; the incomplete primary-k bundle
must be recreated from its retained, checksummed sources.

## Follow-up

On 2026-11-18 the primary-`k` bundle was recreated cleanly from the retained,
checksummed sources. Corrective run
`flt3lko_il7ko_full_mapping_fixed_umap_20261118T031500Z` then completed all 47
bounded transfer batches and the fixed-UMAP outputs. Completion and validation
are recorded in `06_step4_full_mapping_completion.md`; this cleanup record is
retained as provenance for the removed failed state.
