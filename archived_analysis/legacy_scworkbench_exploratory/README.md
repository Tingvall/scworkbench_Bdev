# Archived analysis — exploratory legacy scworkbench

- **Status:** frozen on 2026-10-26
- **Software generation:** older, pre-module-framework scworkbench
- **Role:** provenance and reference only

This directory preserves the repository side of the exploratory B-cell
development analysis. Its scripts, configuration, shared source, and local
figure outputs were moved here without rewriting their historical path strings.

## Contents

| Path | Description |
|---|---|
| `config/` | Legacy bundle registry and analysis configuration |
| `processing/` | Legacy dataset and integration workflows |
| `figure_scripts/` | Legacy report and figure generators |
| `figure_outputs/` | Local generated outputs present when archived; ignored by Git |
| `src/` | Shared legacy helpers |
| `LEGACY_ANALYSIS.md` | Detailed README from the exploratory analysis |
| `LEGACY_REPOSITORY_GUIDANCE.md` | Repository guidance that governed the exploratory analysis; retained as a non-active historical document |
| `LEGACY_GITIGNORE` | Root ignore rules in effect when the exploratory analysis was archived |
| `ARCHIVE_MANIFEST.tsv` | Top-level repository snapshot recorded during the move |

The associated data, bundles, sessions, runtime materials, and other large
artifacts are preserved at:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/archived_analysis/legacy_scworkbench_exploratory/
```

New work belongs in [`../../current_analysis/`](../../current_analysis/).

## Archive boundary

This archive is read-only by project policy. Historical commands and absolute
paths are preserved as provenance and may no longer be runnable from the new
repository root. Reuse an archived artifact by declaring its exact storage path
under `../../current_analysis/inputs/`; do not modify it in place.
