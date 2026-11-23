# Step 4 Bdev PeakVI pseudocell qualification

This bounded development workflow fits PeakVI only on Bdev reference cells and
projects disjoint held-out Bdev single cells, raw Bdev pseudobulks, and
depth-matched Bdev pseudo-cells without refitting. It then calibrates support
and transfers native fitted Slingshot L4 in the shared PeakVI latent space.

The workflow does not download or project ImmGen libraries, build a SIF, write
implementation-stage results, or modify scWorkbench or its operation packs.
Run scripts `01` through `03` with configuration
`config/04_reference_mapping/05_bdev_peakvi_pseudocell_10d_qualification.yaml`.
Generated state is written only to paired `Bdev_development` storage.

Run `bdev_peakvi_pseudocell_10d_20261110T180851Z` passed the bounded scientific
and operation-contract gate using the fixed scWorkbench and scwbLatent sources.
Real ImmGen projection remains a separate controlled development action.
