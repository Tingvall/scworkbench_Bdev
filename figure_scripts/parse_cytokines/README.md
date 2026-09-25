# Parse cytokine reports

`01_report_prepared_seurat.py` summarizes the validated unfiltered Seurat
preparation and supplied the initial `single`, `combination`, and `control`
classification review without loading the large RDS:

```bash
python3 figure_scripts/parse_cytokines/01_report_prepared_seurat.py \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines/prepared \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines/reports/object_inventory_20260915
```

The generated sample and condition TSV files are historical review inputs. The
accepted canonical map now uses the category name `combined` and is recorded at
`config/parse_cytokines_treatment_categories.tsv`; it was applied by
`processing/parse_cytokines/06_classify_treatments_and_select.R`.

After `processing/parse_cytokines/run_04_calculate_rna_qc.sh` has promoted the
factual scWorkbench RNA-QC artifacts, render their distributions with:

```bash
figure_scripts/parse_cytokines/run_02_rna_qc_report.sh
```

The report covers all imported cells, includes overall and per-sample
distributions, and supplies empirical-tail sensitivity previews. Those
previews are not threshold recommendations: the renderer creates no exclusion
or accepted analysis cell set.

After reviewing that report, the user selected total counts >=1,000, detected
features >=500, and mitochondrial percentage <=10%. The combined exclusion
preview was written before the durable sets were derived and is available at:

`intermediate_objects/parse_cytokines/reports/rna_qc_filter_20260921T035031Z/exclusion_preview.md`

The distribution report still represents all imported cells. Filtering did
not replace those figures or remove imported cells. `set:parse_qc_accepted`
retains all QC-accepted treatment categories; the narrower downstream analysis
universe is `set:parse_qc_control_single`, containing only QC-accepted control
and single-cytokine cells. Treatment-selection evidence is under
`intermediate_objects/parse_cytokines/reports/treatment_selection_20260921T145825Z/`.

The individual-cytokine signature workflow generates its report directly from
the promoted scWorkbench UCell score family:

```bash
processing/parse_cytokines/run_07_build_signature_states.sh
processing/parse_cytokines/run_08_discover_specific_signatures.sh
processing/parse_cytokines/run_09_score_and_report_signatures.sh
```

For the current already-committed 31-key UCell family, the report was completed
with the bounded read-only finalizer:

```bash
processing/parse_cytokines/run_09_finalize_committed_ucell_report.sh
```

The final run directory contains `REPORT.md`, an oriented condition-by-
signature specificity heatmap, multi-page sample-median UCell distributions,
exact sample medians, descriptive cell-level quantiles, all candidate/final
statuses, the final feature-set table and GMT, and the committed UCell
promotion audit. The completed report is:

`intermediate_objects/parse_cytokines/reports/signature_report_20260921T221910Z/REPORT.md`

It contains 23 accepted signatures across 14 cytokines. The report
distinguishes the unchanged imported set from the control/single analysis
population and records that neither combined-stimulation cells nor Bdev were
used for signature discovery. Its `validation.json` records the bounded
payload, membership, SQLite, and readiness checks available at that time.
The provenance limitation was subsequently resolved in a separate recovery
runtime: all 23 sets are registered, and full recipient validation has passed.
See
`intermediate_objects/parse_cytokines/reports/signature_recovery_20260923T200125Z/finalization/REPORT.md`.
The original report is preserved; figure runners now use the accepted recovered
bundle and runtime from `config/bundles.yml`.

To explain every exclusion and render all scored gene sets—not only the 23
accepted sets—as cell-level, per-biological-sample violins grouped by
treatment, run:

```bash
figure_scripts/parse_cytokines/run_03_signature_sample_treatment_violins.sh
```

This read-only figure step queries each committed UCell key in a bounded slice.
It writes a 31-page PDF and 31 PNGs (23 accepted and 8 UCell-transfer failures),
plus a manifest and a complete table covering all 52 cytokine/direction
possibilities, to the same signature report directory. The 21 directions with
fewer than five discovery-eligible genes are documented but cannot be plotted
because no gene set or UCell score exists. The explanatory supplement is:

`intermediate_objects/parse_cytokines/reports/signature_report_20260921T221910Z/SIGNATURE_EXCLUSIONS_AND_VIOLINS.md`

For a complete inventory of all 52 individual-stimulus directions, including
excluded candidates, run:

```bash
bash figure_scripts/parse_cytokines/run_04_cytokine_signal_inventory.sh
```

The optional first argument is a new, empty output directory. By default,
reports go to the ignored `figure_outputs/parse_cytokines/signal_inventory_<UTC>/`
directory. Both the accepted Parse bundle and source reports are mounted
read-only; only this output directory is writable in the container.

The offline HTML report provides a searchable, sortable inventory, individual
signature records, three figures, and CSV/TSV downloads. It summarizes saved
gene-level effects/FDR, gene counts at each selection step, discovery and UCell
AUCs against control and other single stimuli, sample-median IQRs, off-target
condition profiles, native score correlations, and gene-membership overlap.
The separate scope table also includes the five combined treatments and control.
No DE, scoring, or promotion is rerun. Rest means the other 25 individual
stimuli, excluding control and combinations; missing metrics for unscored
directions remain NA. The report's metric dictionary explains raw versus
response-oriented scores and distinguishes gene-level FDR from signature AUC.

The inventory verifies the accepted bundle fingerprint, all registered gene
memberships, sample metadata, gene effects, and score coverage. Input/output
checksums, the runtime session, and `validation.json` accompany each report.

`CHECKS.html` is the concise view of the same inventory: pass/fail/not-evaluated
for each acceptance check, eligible and selected gene counts, and the raw
p-value and BH-FDR ranges across the selected genes. `check_summary.csv` and
`individual_check_results.tsv` retain the decisions; `selected_gene_pvalues.tsv`
gives exact values and both contrast filter results for every selected gene.
P-values belong to genes. The AUC and stability checks do not have saved
signature-level p-values. Directions without selected gene sets retain NA in
the selected-set summaries; their few eligible genes, where present, have a
separate p-value download.

To count the consequences of omitting checks from the saved results, run:

```bash
bash figure_scripts/parse_cytokines/run_05_check_omission_report.sh
```

This report separates acceptance of existing gene sets with unchanged
memberships from candidate eligibility after gene filters are removed.
Single-check and grouped scoring omissions use the saved AUC/stability results;
missing downstream checks remain unassessed. Gene-filter counts reuse saved
full-data and training-fold coefficients, verify the original 52 eligible-gene
counts, and report newly eligible directions. They are not final acceptance
counts for revised signatures: no reselection, cross-validation, or UCell
rescoring is performed. The original bundle and discovery outputs stay
read-only. Results, per-direction details, checksums and validation go to a
new ignored output directory; an optional first argument chooses that directory.

The new collection that actually omits treatment-versus-control selection is
rendered by `06_rest_only_signal_report.R`, after its native scoring stage:

```bash
bash processing/parse_cytokines/run_rest_only.sh report
```

Its `REPORT.html` includes all 52 direction outcomes, pass/fail checks, gene
counts, old/new membership comparisons, a three-page overview figure PDF, and
sample-distribution plots for every scored candidate. Passing sets are exported
as GMT/TSV with exact specificity p-values and FDR in a separate gene-evidence
table. Up/down here means expression relative to other individual stimuli;
it does not establish a change from unstimulated control. The report records
all thresholds, the reuse of original models/states, input hashes, and the
unpromoted durable-session location. Both the accepted Parse parent and Bdev
integration remain unchanged.

To render cell-level UCell violin maps for all 51 rest-only directions with at
least five eligible genes, including control and every individual cytokine
condition, run:

```bash
bash figure_scripts/parse_cytokines/run_07_rest_only_all_candidate_violin_maps.sh
```

The combined 51-page PDF and four PDFs split by UCell pass/fail and up/down are
written under the rest-only report. The 47 discovery-passing memberships reuse
their completed scores. The four discovery failures are assigned their
prespecified best-sized exploratory memberships and scored in a separate
unpromoted scWorkbench session solely to complete the requested figure. Every
page shows three cell-level violins per condition, one per biological sample,
with sample medians overlaid. Control is displayed but remains excluded from
the rest-only UCell AUC. No result is promoted.

To render the same violin layout for all 52 gene sets selected by the newer
competitor-aware analysis, run:

```bash
bash figure_scripts/parse_cytokines/run_08_competitor_aware_selected_violin_maps.sh
```

This writes a combined 52-page PDF and four PDFs split by full-data pass/fail
and direction inside the competitor-aware review directory. It reuses the
completed native UCell scores, shows the three biological samples separately,
and records the selected size, full-data gaps, holdout stability, and core size
on each page. Control is displayed and is part of the competitor-aware test.
No result is promoted.
