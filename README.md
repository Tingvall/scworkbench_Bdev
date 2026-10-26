# B-cell development analysis project

This project was reorganized on 2026-10-26 to separate the completed
exploratory analysis from the next analysis based on the updated,
module-oriented scworkbench.

## Project map

| Location | Status | Purpose |
|---|---|---|
| [`current_analysis/`](current_analysis/) | **ACTIVE** | Configuration, documentation, input declarations, custom modules, workflows, and output contracts for the updated scworkbench analysis |
| [`archived_analysis/`](archived_analysis/) | **ARCHIVED / READ-ONLY** | Legacy exploratory code, configuration, documentation, and local generated figures retained for provenance |

These are the only two analysis roots. Legacy root-level `config/`,
`processing/`, `figure_scripts/`, `figure_outputs/`, and `src/` directories
have been moved into the archive.

## Mirrored storage

Large inputs, containers, bundles, working state, and results live under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/
```

That workspace has the same two analysis roots:

- `current_analysis/` is paired with this repository's active scaffold.
- `archived_analysis/legacy_scworkbench_exploratory/` contains the data and
  runtime artifacts associated with the archived repository workflow.

No new analysis should use or recreate the former storage-root paths
`datasets/`, `intermediate_objects/`, or `results/`.

## Current transition state

The legacy analysis is frozen. The current analysis directory is intentionally
a clean scaffold: the exact analyses to reproduce and the updated scworkbench
module configuration have not yet been specified. The next step is to pin the
new scworkbench version, select the legacy steps to reproduce, and map each
selected step to the new module framework.
