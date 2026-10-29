# Implementation results

This result area is locked during `Bdev_development`. Do not create paper tables,
figures, reports, bundles, or other implementation outputs before the qualified
replacement SIF is explicitly accepted.

Final paper tables, reports, figures, exports, and capsules belong under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_implementation/results/
```

Every result must link to the producing workflow, configuration, input and
bundle identities, qualified operation definitions, and accepted SIF checksum.

Mirror the numbered workflow hierarchy in storage. Within each step, use
numbered output categories and numbered files, for example:

```text
results/
└── 02_immgen_bulk_rna/
    └── 02_marker_discovery/
        ├── 01_tables/
        ├── 02_figures/
        ├── 03_reports/
        ├── 04_logs/
        ├── 05_objects/
        └── 06_manifests/
```

Intermediate objects use the same dataset and step identifiers beneath the
implementation `intermediate_objects/` root.
