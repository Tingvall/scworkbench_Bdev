# Shared source code

This directory contains reusable project code shared by more than one pipeline
or report. Dataset-specific orchestration remains in its `processing/` or
`figure_scripts/` directory.

- `bundle_config.R` reads and validates the project bundle registry.
- `bundle_validation.R` checks an accepted registry entry against a read-only
  scWorkbench bundle and can run recipient validation, using a registered
  runtime library when the bundle requires one.
- `parse_runtime.sh` resolves the active Parse bundle and optional recovery
  runtime for processing and figure runners from the central registry.
- `scworkbench_registry_json.R` is installed into the separate recovery engine
  build; it preserves registry identity while streaming large JSON for hashing.
  Its build and validation scripts are in `processing/parse_cytokines/recovery/`.
- `scworkbench_object_row_keys.R` supplies the recovery engine's equivalent
  vectorized character row-key serializer for full score-family validation.
- `cell_identity.R` provides reusable ordered and cross-bundle cell-ID checks.
- `report_helpers.R` provides small output and table helpers.
- `reference_validation.R` contains the existing Bdev import comparisons.
