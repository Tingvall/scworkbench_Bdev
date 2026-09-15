# Shared source code

This directory contains reusable project code shared by more than one pipeline
or report. Dataset-specific orchestration remains in its `processing/` or
`figure_scripts/` directory.

- `bundle_config.R` reads and validates the project bundle registry.
- `bundle_validation.R` checks an accepted registry entry against a read-only
  scWorkbench bundle and can run recipient validation.
- `cell_identity.R` provides reusable ordered and cross-bundle cell-ID checks.
- `report_helpers.R` provides small output and table helpers.
- `reference_validation.R` contains the existing Bdev import comparisons.
