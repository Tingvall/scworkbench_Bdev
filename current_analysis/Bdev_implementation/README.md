# Bdev paper implementation

**Purpose:** run the Bdev paper analysis with the immutable, qualified SIF
produced by `Bdev_development`.

**Current status: locked planning scaffold.** Do not execute workflows or create
implementation state or results until development has completed all accepted
framework updates, built and qualified the replacement SIF, and the user has
explicitly accepted that runtime for implementation.

This stage consumes released capabilities. It does not modify scWorkbench core,
PyBench, operation-pack source, or the accepted SIF. If a required operation is
missing or scientifically inadequate, record the gap and return it to the
development stage; do not patch the implementation environment in place.

## Layout

| Directory | Purpose |
|---|---|
| `config/` | Accepted SIF and operation identities, analysis parameters, seeds, and compute profiles. |
| `docs/` | Approved scientific plan, provenance mapping, decisions, and paper-analysis documentation. |
| `inputs/` | Checksummed raw inputs, source bundles, external preprocessing outputs, and the qualified runtime manifest. |
| `workflows/` | Reproducible paper-analysis definitions using qualified operation names. |
| `intermediate_objects/` | Repository contract for generated bundles and working state stored outside Git. |
| `results/` | Repository contract for final tables, reports, figures, and capsules stored outside Git. |

Every result must resolve to a workflow, configuration, input manifest,
qualified operation catalog, and SIF checksum. Preserve the exact values and
missing entries of `Bdev@meta.data$pt` whenever the reference is imported or
used.
