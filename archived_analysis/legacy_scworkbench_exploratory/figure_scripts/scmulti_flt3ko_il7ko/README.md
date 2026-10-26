# FLT3LKO/IL7KO multiome reports

`01_unfiltered_qc_report.Rmd` reports per-library RNA and ATAC QC for all cells
in the Cell Ranger ARC filtered matrix. It applies no thresholds or filtering.

`02_qc_selection_report.Rmd` reports the user-approved downstream-use marker.
All cells remain retained; cells passing every strict threshold are marked for
downstream analysis.

The stable technical key is retained for compatibility; the biological
knockout is the ligand *Flt3l* (`FLT3LKO`).
