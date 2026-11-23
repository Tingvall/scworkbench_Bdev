# Step 4 Bdev L4 reference-mapping qualification

This bounded development workflow freezes the selected ImmGen ATAC preflight
and qualifies held-out-sample reference mapping on the authenticated 1,810-cell
native-L4 Bdev fixture. It is not a bulk-ATAC projection, paper analysis, or SIF
build.

The workflow:

1. clones the completed active 03c bundle for traceability and freezes the 16
   selected ImmGen populations and 28 GSE100738 ATAC libraries;
2. installs the pinned core and operation-pack commits into an isolated run
   library;
3. assembles a fresh fixture-only scWorkbench bundle from real 03c chromVAR,
   LSI, metadata, and authenticated native-L4 inputs;
4. runs four leave-one-`orig.ident`-out calibration and transfer folds using
   LSI dimensions 2--30 and exact RANN neighbors; and
5. inventories every generated result with SHA-256.

Run from the repository root with the commands encoded by scripts `01` through
`04` and configuration
`config/04_reference_mapping/02_bdev_l4_reference_mapping_qualification.yaml`.
Generated state is written only to paired `Bdev_development` storage.

The contract qualification passed, but the aggregate LSI mapping performance
did not justify direct bulk projection. The next gate is depth-matched Bdev
pseudobulk/pseudo-cell validation through a Bdev-only frozen PeakVI model.
