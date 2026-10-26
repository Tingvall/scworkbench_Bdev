# Current analysis — module-based scworkbench

**Status:** active scaffold; analysis selection pending

This is the only active project area for the redo using the updated scworkbench
module framework. It is deliberately separate from the archived exploratory
run so legacy configuration and new-framework outputs cannot be confused.

## Layout

| Directory | Purpose |
|---|---|
| `config/` | Pinned scworkbench version and project/run configuration |
| `docs/` | Reproduction plan and old-to-new module mapping |
| `inputs/` | Input declarations and provenance; large inputs remain in mirrored storage |
| `modules/` | Only project-specific modules not supplied by scworkbench |
| `workflows/` | Reproducible definitions for selected analysis steps |
| `intermediate_objects/` | Contract for generated working state; actual large artifacts live in mirrored storage |
| `results/` | Contract for final outputs; actual generated outputs live in mirrored storage |

The paired storage location is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/
```

## Before running analysis

1. Record the updated scworkbench release, immutable source revision, or
   container digest in `config/`.
2. Document the updated default modules and extension interface in `docs/`.
3. Select the legacy analysis steps to reproduce.
4. For each selected step, identify its archived inputs and expected outputs,
   then map it to default modules or a clearly scoped custom module.
5. Keep generated paths under the paired storage `current_analysis/` root.

The archived repository workflow is available at
[`../archived_analysis/legacy_scworkbench_exploratory/`](../archived_analysis/legacy_scworkbench_exploratory/).
