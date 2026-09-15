# Bdev reference workflow

1. `01_inspect_bdev_reference.R` inventories the source without modifying it.
2. `02_create_bdev_reference.R` creates and validates a new bundle.
3. `03_validate_bdev_reference.R` validates an existing import without
   recreating its data.

The accepted bundle and its exact source, pseudotime, and validation state are
recorded in `config/bdev_reference.yml`.
