# ImmGen bulk ATAC mapping through the selected gated 20D reference

This development workflow projects the checksummed 28-library ImmGen ATAC
panel through the exact selected HVP-depth20k PeakVI model without refitting.
It reuses the qualified calibration and applies the same reference-support,
L4-domain, confidence, and finite-value gate used for single-cell queries.

The fixed-UMAP figure keeps the original Bdev `umap.atac` coordinates. Bulk
libraries that pass the complete gate are colored by accepted native L4;
libraries that fail are shown in dark grey. Neither UMAP nor Slingshot is
refitted. Outputs are development validation evidence and do not activate the
implementation stage.

Run in order:

1. `01_prepare_clone.py`
2. `02_project_map_and_plot.R`
3. `03_finalize.py`
