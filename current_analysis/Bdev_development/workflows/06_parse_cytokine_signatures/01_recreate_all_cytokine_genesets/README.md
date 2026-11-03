# 01 Recreate the 38 Parse cytokine gene sets

This development workflow reconstructs the selected 38-set all-cytokine Parse
collection without loading or calling scWorkbench. It preserves the archived
scientific method:

1. rank genes independently in each of three saved limma specificity fits;
2. intersect each direction's three top-100 lists after restricting to the
   Bdev `SCT_counts` feature universe;
3. greedily remove genes when doing so improves the target-to-strongest-
   competitor standardized gap, stopping at 20 genes; and
4. retain the 38 sets whose target has the strongest median sample response.

The frozen full-precision sample-level capped-rank cache replaces the expensive
cell-level UCell calculation. It is an input statistic, not a copied final gene
set. The script independently rebuilds candidate and final membership and then
requires exact agreement with the archived 38-set result.

Run from the repository root:

```bash
Rscript current_analysis/Bdev_development/workflows/06_parse_cytokine_signatures/01_recreate_all_cytokine_genesets/01_recreate_all_cytokine_genesets.R \
  --fits /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/archived_analysis/legacy_scworkbench_exploratory/intermediate_objects/parse_cytokines/reports/signature_discovery_20260921T221102Z/limma_voom_fits.rds \
  --bdev-sct-count-features /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/archived_analysis/legacy_scworkbench_exploratory/intermediate_objects/bdev_reference/bundle/mat/mat_frag_SCT_counts_20260914_162226/row_names \
  --bdev-sct-data-features /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/archived_analysis/legacy_scworkbench_exploratory/intermediate_objects/bdev_reference/bundle/mat/mat_frag_SCT_data_20260914_162228/row_names \
  --rank-cache /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/archived_analysis/legacy_scworkbench_exploratory/intermediate_objects/parse_cytokines_rest_only/reports/all_cytokine_competitor_refinement_top100_v1_20260925/top100_union_sample_rank_cache.rds \
  --sample-metrics /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/archived_analysis/legacy_scworkbench_exploratory/intermediate_objects/parse_cytokines_rest_only/reports/all_cytokine_competitor_refinement_top100_v1_20260925/sample_metrics.tsv \
  --reference-report /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/archived_analysis/legacy_scworkbench_exploratory/intermediate_objects/parse_cytokines_rest_only/reports/all_cytokine_competitor_refinement_top100_v1_20260925 \
  --output /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/06_parse_cytokine_signatures/01_recreate_all_cytokine_genesets/parse38_standalone_reconstruction_20261103T052911Z
```

Generated outputs use numbered table, report, log, and manifest directories.
The archived analysis is read-only and no paper-analysis implementation is run.
