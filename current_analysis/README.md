# Current analysis — Bdev development and implementation

The active work is divided into two gated stages. Development extends and
qualifies scWorkbench for the Bdev paper use case. Implementation then runs the
paper analysis with the immutable runtime produced by development.

**Current phase: `Bdev_development`.** The implementation tree is a reserved
planning scaffold only. No paper workflow is executed until framework updates
are complete and the replacement SIF is qualified and explicitly accepted.

```text
Bdev_development
  current candidate SIF + source branches
  -> approved core changes, if any, and official operation-pack changes
  -> contract and scientific validation
  -> new immutable qualified SIF

Bdev_implementation
  qualified SIF manifest
  -> frozen analysis configuration and workflows
  -> scWorkbench bundles
  -> paper results and figures
```

## Stage boundaries

| Directory | Responsibility |
|---|---|
| `Bdev_development/` | Develop missing reusable operations in scWorkbench core or the appropriate official operation packs, qualify them, and build the replacement SIF. |
| `Bdev_implementation/` | Run the approved Bdev paper analysis using only the qualified replacement SIF and frozen analysis inputs/configuration. |

Preprocessing that is not appropriately represented as a reusable scWorkbench
operation may remain outside the framework, but its dependencies, inputs,
parameters, and outputs must be declared before import.

Development and implementation must not share writable bundles or result
directories. Implementation may start only after its runtime manifest points to
a development-produced SIF with a verified checksum and accepted qualification
record.

The paired active storage root remains:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/
```

Use stage-specific paths below that root:

```text
Bdev_development/{intermediate_objects,results}/
Bdev_implementation/{intermediate_objects,results}/
```

All datasets, workflow steps, scripts, and generated output categories use the
stable numbered hierarchy defined in [`ORGANIZATION.md`](ORGANIZATION.md).
Repository workflow paths and paired-storage result paths must mirror one
another so that every output has an unambiguous producing dataset and step.

The archived exploratory generation remains read-only at
[`../archived_analysis/legacy_scworkbench_exploratory/`](../archived_analysis/legacy_scworkbench_exploratory/).
