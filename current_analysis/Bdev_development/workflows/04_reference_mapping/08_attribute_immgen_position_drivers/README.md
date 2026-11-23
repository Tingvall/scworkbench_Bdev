# Attribute ImmGen position drivers

This diagnostic asks why an ImmGen bulk ATAC library selects its observed Bdev
neighbors instead of cells carrying its corresponding provisional ImmGen RNA
label. It preserves the authenticated native L4 and does not accept the
provisional labels as truth.

The workflow projects the labelled Bdev diagnostic fixture through the same
deterministic Bdev-only PeakVI model, without adding those cells to the mapping
reference. It then reports expected-label neighbor weight, latent-dimension
pull, empirical peak contrasts, reference-cell hub reuse, donor/sample
enrichment, and associations with ATAC/RNA QC variables. Peak rankings are
candidate drivers and require masking/reprojection for causal confirmation.

Run the numbered actions with
`config/04_reference_mapping/14_immgen_position_driver_attribution.yaml`:

1. `01_prepare_attribution.R` creates the isolated diagnostic bundle, performs
   the frozen-recipe projections, and exports peak candidates and coordinates.
2. `02_analyze_attribution.py` builds the label, latent-vector, and technical
   bias diagnostics.
3. `03_finalize_run.py` validates and inventories the result tree.

Run `immgen_position_drivers_20261116T003136Z` completed for all 28 libraries.
Seventeen had an exact provisional RNA-label centroid supported by at least
five diagnostic cells. The dominant unexpected-position component was
`latent_3`, median cluster-6 neighbor weight was 0.539, reference sample 2 was
enriched 2.63-fold among neighbor weights, and 62/65 recurrent peak candidates
were all-or-none across the ImmGen panel. Candidate peaks remain non-causal
until matched masking and frozen-model reprojection are completed.
