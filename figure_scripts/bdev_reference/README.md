# Bdev reference reports

`01_bdev_reference_report.Rmd` renders the accepted Bdev reference dataset
report from the promoted bundle. Its UMAP panels use `umap.atac`, the peaks
assay embedding, and the accepted pseudotime artifact recorded in the bundle.

`02_il7_flt3_top25_core_atac_umap.R` is a read-only sanity check for the
Bdev-compatible Parse signature review. It plots log1p(CP10k) Bdev RNA
expression for every IL-7 and Flt3L top25 core gene on the accepted
`umap.atac` coordinates, with separate combined figures for up and down sets.
Run it with `run_02_il7_flt3_top25_core_atac_umap.sh`; outputs and provenance
are recorded under the report path in `config/bdev_reference.yml`.

`03_flt3_il7_refined_sct_umap.R` plots the final top-100 refinement-pilot
memberships as four separate figures: Flt3L up/down and IL-7 up/down. It reads
expression directly from Bdev `SCT_data` and uses the existing `umap.atac`
coordinates. Non-positive cells are gray; positive cells are overlaid as
larger, fully opaque points on a saturated blue-to-magenta scale. Run it with
`run_03_flt3_il7_refined_sct_umap.sh`. The final v2 output and provenance are
registered in `config/bdev_reference.yml`.

`04_flt3_il7_refined_aggregate_sct_umap.R` creates one four-panel UMAP for
the same final sets. A cell's aggregate expression score is the arithmetic
mean of its `SCT_data` values across all genes in the set. The raw per-cell
scores and panel summaries accompany the figure. Run it with
`run_04_flt3_il7_refined_aggregate_sct_umap.sh`; its report path is registered
in `config/bdev_reference.yml`.

`05_flt3_il7_refined_ucell_trajectory.R` reproduces the earlier native-UCell
trajectory method for the four refined sets. It scores raw Bdev `RNA_counts`
with `max_rank = 1500`, preserves the existing `pt` values and missingness, and
uses 5% sliding windows with 1% steps, median lines, and IQR bands. Each curve
point is placed at the midpoint of its actual window, including clipped endpoint
windows. Dotted vertical lines mark the median source `pt` for the `up1` and
`up2` groups in Bdev metadata field `bb_up`. Run it with
`run_05_flt3_il7_refined_ucell_trajectory.sh`; the unpromoted report and score
session are registered in `config/bdev_reference.yml`.

`06_flt3_il7_refined_ucell_umap.R` displays the saved native UCell scores from
the refined trajectory analysis on Bdev `umap.atac`. It creates a four-panel
overview and four separate figures. Raw UCell values are retained in the
output table; the bright color scale is normalized within each panel to the
positive-score 99th percentile. Run it with
`run_06_flt3_il7_refined_ucell_umap.sh`; its read-only report is registered in
`config/bdev_reference.yml`.

`07_flt3_il7_ucell_trajectory_all_windows_diagnostic.R` is a sensitivity
display of the refined Flt3L/IL-7 trajectories. It retains the six windows
below the primary 30-cell support threshold and marks their unchanged native
medians with hollow points. It does not recalculate UCell scores or window
summaries. Run it with
`run_07_flt3_il7_ucell_trajectory_all_windows_diagnostic.sh`; the diagnostic
report is registered in `config/bdev_reference.yml`.

# All-cytokine refined UCell report

After reviewing a completed all-cytokine Parse refinement report, run:

```bash
figure_scripts/bdev_reference/run_08_all_cytokine_refined_ucell_report.sh \
  /path/to/all_cytokine_competitor_refinement_report
```

The unpromoted report scores only final Parse sets passing the mean and/or
responder route. It uses native scWorkbench UCell on Bdev `RNA_counts`, the
unchanged `umap.atac` coordinates, and the unchanged `pt` metadata. Pseudotime
curves use 5% windows, 1% steps, a 20-cell minimum, median lines and type-7 IQR
bands. Missing `pt` values remain missing.
