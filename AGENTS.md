# Repository operating mode

These instructions apply to the entire repository.

## Project

This repository develops the B-cell development use case for scWorkbench. It
contains one frozen exploratory generation, a development stage for extending
and qualifying scWorkbench, and an implementation stage for the paper analysis.

For the Bdev reference, use the existing `Bdev@meta.data$pt` field as
pseudotime, as selected by the user. Preserve its values and missing entries
when importing and using the reference.

## Archived analysis

Treat `archived_analysis/` and everything below it as read-only unless the user
explicitly requests a legacy repair. Do not use archived paths for new outputs.

Preserve legacy filenames, manifests, timestamps, and internal references. If
a legacy artifact is needed by the new analysis, declare its exact archived
storage path in the applicable development or implementation input manifest and
document its provenance instead of modifying the archive.

The archived repository code is paired with storage at:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/archived_analysis/legacy_scworkbench_exploratory/
```

## Active analysis

The current project phase is `Bdev_development`. Treat
`current_analysis/Bdev_implementation/` as a reserved, non-executing scaffold
until all approved framework and operation-pack updates are complete and a new
immutable SIF has been built, qualified, checksummed, and explicitly accepted.
Do not run paper analyses, create implementation results, or treat provisional
implementation configuration as active during development.

Put all new work under one of the two active stages:

- `current_analysis/Bdev_development/` develops missing reusable capabilities
  in scWorkbench core or official operation packs, qualifies them, and produces
  a new immutable SIF;
- `current_analysis/Bdev_implementation/` consumes only the qualified SIF for
  the frozen Bdev paper analysis.

Within each stage, use `config/`, `docs/`, `inputs/`, and `workflows/` for
repository-tracked definitions. Development may use `modules/` only for
temporary Bdev-owned prototypes or integration fixtures. Put generated working
state and results in the corresponding stage path under paired storage.

The paired active storage root is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/
```

Do not put large inputs, containers, bundles, models, or generated results in
Git.

## Analysis framework

Use scWorkbench as the default analysis framework and create traceable
scWorkbench bundles for each selected dataset. Development may use the current
candidate SIF as a toolchain baseline, but must implement reusable scientific
methods in the owning scWorkbench core or official operation-pack repository.

Treat `scworkbench` and PyBench as core infrastructure and leave them unchanged
during operation development by default. Add scientific operations to the
existing official pack whose domain and contracts fit, or create a new official
operation pack when none does. Modify either core repository only for a pressing,
documented infrastructure need, such as a missing generic primitive, authenticated
reader or writer, registry, provenance, or publication capability that cannot be
implemented correctly within an operation pack.

Never modify `scworkbench` or PyBench without the user's explicit approval for
the specific proposed core change. If a core change appears necessary, stop
before editing, clearly identify and motivate the missing capability, show why
an existing or new operation pack cannot implement it correctly, describe the
smallest proposed core change and its qualification impact, and ask the user.
Do not infer approval from approval of the analysis, an operation, a pack, or a
previous unrelated core change.

Before paper workflows are implemented, build and qualify a new immutable SIF,
pin its exact digest and source commits, and consult its supplied documentation.
Do not infer the interface from archived runtime source or duplicate framework
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

## Organization and numbering

Maintain a stable hierarchical organization for all active work. Assign every
dataset or derived bundle a permanent two-digit `NN_slug` identifier in the
implementation dataset registry before adding workflows for it. Within each
dataset, use numbered `NN_step_slug` workflow directories and numbered
`NN_action_slug.R`, `.py`, or `.sh` scripts in execution order.

Generated storage must mirror the repository workflow hierarchy. Within each
step, use numbered output directories such as `01_tables/`, `02_figures/`,
`03_reports/`, `04_logs/`, `05_objects/`, and `06_manifests/`, omitting only
categories that do not apply. Number figure and report filenames within their
category. Do not write generated outputs into Git.

Once an identifier is referenced by a manifest, bundle, result, or report, do
not renumber or reuse it. Mark removed work as retired and allocate a new number
for replacements. Do not use ambiguous state names such as `final`, `latest`,
`new`, or ad hoc `v2`; distinguish reruns with recorded run identities,
configuration hashes, and timestamps while preserving the stable dataset and
step path.
