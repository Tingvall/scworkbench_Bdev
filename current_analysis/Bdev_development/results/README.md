# Development results

Retain qualification reports, operation-catalog snapshots, source/manifest
identity records, and immutable SIF builds under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/
```

A SIF becomes eligible for implementation only after its checksum, exact source
commits, installed pack identities, test coverage, exits, skips, and known
limitations are recorded and approved.

Mirror the numbered development workflow and step paths under this result root.
The cytokine work uses the stable split labels
`03a_parse_cytokine_signal_generation/` and
`03b_parse_cytokine_trajectory_analysis/`; do not recreate the provisional
`03_parse_cytokine_signatures/`, `06_parse_cytokine_signatures/`, or
`07_trajectory_signal_regions/` roots.
Within each step, use the applicable numbered `01_tables/`, `02_figures/`,
`03_reports/`, `04_logs/`, `05_objects/`, and `06_manifests/` categories.

The checksummed run
`03b_parse_cytokine_trajectory_analysis/01_qualify_cytokine_trajectory_methods/bdev_cytokine_trajectory_20261105T153001Z`
is retained as the preliminary cytokine trajectory analysis for comparison with
the later implementation rerun. Its development outputs are not promoted paper
results.
