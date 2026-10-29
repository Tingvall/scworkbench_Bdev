# Implementation workflows

This directory is currently a non-executing scaffold. Its numbered dataset
layout is reserved for later implementation and must remain inactive until the
replacement SIF is qualified and accepted.

Store the reproducible Bdev paper workflows here. Workflows must:

- verify the accepted SIF checksum before execution;
- activate the frozen official operation-pack profile explicitly;
- use qualified operation names and recorded parameters;
- preserve source bundles and write joint state only to derived bundles;
- preserve the authoritative Bdev `pt` values and missing entries; and
- write generated state only to the implementation storage paths.

Method development or runtime patching is out of scope for this directory.

Organize workflows first by the permanent dataset ID in
`../config/dataset_registry.yaml`, then by the numbered analysis step. Number
scripts within each step in execution order. For example:

```text
workflows/
└── 02_immgen_bulk_rna/
    ├── 01_input_audit/
    │   ├── 01_validate_manifest.R
    │   └── 02_render_input_qc.R
    └── 02_marker_discovery/
        ├── 01_fit_contrasts.R
        └── 02_export_marker_sets.R
```

Do not renumber a dataset, step, or script after it has produced referenced
state. Retire it in place and allocate the next unused number.
