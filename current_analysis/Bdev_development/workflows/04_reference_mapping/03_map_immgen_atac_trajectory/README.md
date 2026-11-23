# Step 4 ImmGen ATAC trajectory mapping

This development workflow downloads the 28 frozen GSE100738/SRP110978 ImmGen
B-cell ATAC libraries, reprocesses them against mm10, counts deduplicated and
blacklist-filtered Tn5-shifted fragments on the exact 191,433-feature Bdev peak
axis, and projects the resulting library profiles through the qualified 10D
Bdev PeakVI recipe.

The scientific output is reference support plus native fitted Slingshot L4
pseudotime. Categorical Bdev cluster transfer is intentionally omitted. The
workflow preserves raw FASTQs, writes generated state only to paired
`Bdev_development` storage, and treats the declared archived mm10 index and
blacklist as read-only legacy inputs.

Run numbered scripts in order with
`config/04_reference_mapping/06_immgen_atac_trajectory_mapping.yaml`.

1. `01_install_runtime.py` installs the pinned core and operation-pack sources
   into the isolated run library.
2. `02_prepare_immgen_inputs.py` resolves the frozen GSM panel against the ENA
   study report.
3. `03_export_bdev_peak_axis.R` exports the exact ordered Bdev peak axis.
4. `04_preprocess_immgen_atac.py` downloads and preprocesses the 28 libraries;
   it is resumable and supports bounded parallel workers.
5. `05_map_immgen_trajectory.R` performs the scWorkbench projection, support
   calibration, and L4-only transfer.
6. `06_finalize_run.py` inventories and checksums the result tree.

## Completed run

Run `immgen_atac_l4_20261110T202521Z` completed all six actions. All 28
libraries passed preprocessing on the exact 191,433-peak axis. The qualified
10D reference checksum was reproduced, frozen-model projection and the
authenticated latent union completed, and support calls were 25 supported, 3
borderline, and 0 out-of-reference. The final inventory contains 48 files.
No categorical cluster labels were transferred, no SIF was built, and
implementation remains locked.
