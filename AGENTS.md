# Repository operating mode

These instructions apply to the entire repository.

## Project

This repository develops the B-cell development use case for scWorkbench. It
contains one frozen exploratory generation and one active scaffold for a redo
with the updated, module-oriented scworkbench.

For the Bdev reference, use the existing `Bdev@meta.data$pt` field as
pseudotime, as selected by the user. Preserve its values and missing entries
when importing and using the reference.

## Archived analysis

Treat `archived_analysis/` and everything below it as read-only unless the user
explicitly requests a legacy repair. Do not use archived paths for new outputs.

Preserve legacy filenames, manifests, timestamps, and internal references. If
a legacy artifact is needed by the new analysis, declare its exact archived
storage path in `current_analysis/inputs/` and document its provenance instead
of modifying the archive.

The archived repository code is paired with storage at:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/archived_analysis/legacy_scworkbench_exploratory/
```

## Active analysis

Put all new work under `current_analysis/`:

- framework and run configuration in `current_analysis/config/`
- migration plans and analysis documentation in `current_analysis/docs/`
- input manifests and provenance in `current_analysis/inputs/`
- project-specific scworkbench modules in `current_analysis/modules/`
- reproducible analysis definitions in `current_analysis/workflows/`
- generated working state in the paired storage `current_analysis/intermediate_objects/`
- final reports and figures in the paired storage `current_analysis/results/`

The paired active storage root is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/
```

Do not put large inputs, containers, bundles, models, or generated results in
Git.

## Analysis framework

Use scWorkbench as the default analysis framework and create traceable
scWorkbench bundles for each selected dataset. Before implementing workflows,
pin the exact updated scworkbench release, commit, or container digest and
consult the documentation supplied with that version. Do not infer the new
module interface from archived runtime source or duplicate framework
documentation in this repository.

Use analysis outside scWorkbench only when the user explicitly requests it, or
when scWorkbench cannot perform the task and the user approves the proposed
alternative. Existing explicit authorization for the same scope is sufficient.

For authorized outside analysis, use reproducible scripts with recorded
dependencies, inputs, parameters, and random seeds where applicable. Preserve
raw data and scWorkbench bundles, save outside results separately, and validate
cell/sample identifiers and metadata alignment.

## Reproducibility

- Keep dataset identities, bundle locations, parent relationships, and source
  checksums explicit.
- Preserve source bundles during integration; write joint state only to a
  separate derived bundle.
- Record the configuration and workflow responsible for every generated result.
- Keep changes focused and report what changed and how it was validated.
