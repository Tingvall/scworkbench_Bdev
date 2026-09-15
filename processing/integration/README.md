# Bdev-centered integration workflow

Integration creates the separate `bdev_integrated` bundle. It must never add
derived integration state to any of the three accepted parent bundles.

Expected phases are:

1. `01_preflight_inputs.R`: validate bundle identities, accepted fingerprints,
   source/sample metadata, cell identifiers, and feature/genome compatibility.
2. `02_create_bdev_integrated.R`: initialize the derived integration bundle and
   record all three parents, with `bdev_reference` as the reference.
3. `03_map_parse_to_bdev.R`: attach reviewed Parse-derived cytokine signals to
   the explicit Bdev reference trajectory.
4. `04_map_multiome_to_bdev.R`: map the KO query state into the selected Bdev
   reference representation.
5. `05_joint_analysis.R`: run analyses that require the linked datasets.
6. `06_validate_bundle.R`: perform source, cell, coordinate, provenance, and
   recipient validation.

Do not create an integration script until the applicable scWorkbench operation
interfaces have been inspected in the supplied image and all required parent
states have been accepted.
