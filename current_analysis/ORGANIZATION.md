# Active-analysis organization and numbering

The current phase is development. Its active hierarchy is organized by numbered
framework capabilities and qualification steps:

```text
Bdev_development/
├── workflows/
│   └── NN_capability_slug/
│       └── NN_step_slug/
│           └── NN_action_slug.{R,py,sh}
└── results/                         # paired storage, not Git
    └── NN_capability_slug/
        └── NN_step_slug/
            └── NN_output_category/
```

The dataset hierarchy below is reserved for the later implementation phase. It
must not be populated with paper-analysis results until the new SIF is accepted.

The active repository and paired storage use the same stable hierarchy:

```text
stage/
├── workflows/
│   └── NN_dataset_slug/
│       └── NN_step_slug/
│           └── NN_action_slug.{R,py,sh}
└── results/                         # paired storage, not Git
    └── NN_dataset_slug/
        └── NN_step_slug/
            ├── 01_tables/
            ├── 02_figures/
            ├── 03_reports/
            ├── 04_logs/
            ├── 05_objects/
            └── 06_manifests/
```

## Dataset identifiers

Dataset and derived-bundle IDs are permanent two-digit identifiers followed by
a lowercase snake-case slug, for example `02_immgen_bulk_rna`. Reserve them in
`Bdev_implementation/config/dataset_registry.yaml` before creating a workflow.
The identifier is a stable key, not a mutable display title or an assertion that
all proposed datasets are approved.

Do not renumber or reuse an identifier after it appears in a manifest, bundle,
result, report, or external reference. Retired datasets remain in the registry
with their disposition and replacement ID when applicable.

## Steps and scripts

Within a dataset, assign steps in dependency order as `NN_step_slug`. Within a
step, number executable scripts in their intended order. One script should have
one primary responsibility. A workflow manifest records dependencies and remains
authoritative when execution is conditional or parallel; numeric prefixes make
the ordinary path legible but do not replace dependency declarations.

If a completed script or step is replaced, retain its original identity and add
the replacement under the next unused number. Do not create `final`, `latest`,
`new`, `redo`, or ad hoc version-suffixed paths.

## Generated outputs

Generated state belongs only in paired storage. Mirror the dataset and step path
from `workflows/`, then use the applicable numbered categories:

1. `01_tables/` — machine-readable result and QC tables;
2. `02_figures/` — numbered figures and their plotting-data manifests;
3. `03_reports/` — rendered reports and capsules;
4. `04_logs/` — execution, scheduler, and validation logs;
5. `05_objects/` — bundles, models, serialized objects, and other large state;
6. `06_manifests/` — run, input, configuration, checksum, and provenance records.

Number files within categories when more than one is produced, such as
`01_sample_qc.pdf` and `02_contrast_summary.pdf`. Every output manifest must name
the producing dataset ID, step ID, script, configuration, runtime, inputs, and
checksums. Use run identities and timestamps inside manifests or run-specific
subdirectories when reruns are retained; do not change the stable hierarchy.

## Development stage

Framework development follows the same principle, using numbered capability and
qualification-step identifiers rather than biological dataset identifiers. It
must remain distinct from implementation workflows and results. Reserve active
capability IDs in
`Bdev_development/config/capability_registry.yaml`. Development may use small,
explicit fixtures and Bdev acceptance fixtures, but those are qualification
evidence rather than paper-analysis results.
