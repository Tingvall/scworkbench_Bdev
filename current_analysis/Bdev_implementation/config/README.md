# Implementation configuration

The implementation runtime is frozen. Record:

- the accepted SIF path, filename, architecture, and SHA-256;
- scWorkbench, PyBench, and operation-pack source commits and package versions;
- activated official packs and their installed identities;
- qualified operation names and definition identities;
- workflow parameters, random seeds, and compute profiles; and
- repository and paired-storage roots.

Do not point implementation workflows at a development branch, mutable package
library, or unqualified candidate SIF.

[`reference_mapping.yaml`](reference_mapping.yaml) records the selected
`gated_20d_shared_reference_mapping` proposed method, its exact 20-dimensional
HVP20k model, support/native-L4 limits, fixed-UMAP reporting convention, and
successful development bulk validation. The same exact model is required for
single-cell and bulk ATAC mapping. The method remains inactive until the
qualified implementation SIF is pinned and explicitly accepted.

[`runtime.yaml`](runtime.yaml) is deliberately unaccepted and incomplete until
development produces and qualifies the replacement SIF.

[`dataset_registry.yaml`](dataset_registry.yaml) reserves the permanent numbered
dataset and derived-bundle identifiers. Its entries remain proposed until their
scientific inputs and bundle roles are approved, but an identifier must not be
renumbered or reused after it is referenced.
