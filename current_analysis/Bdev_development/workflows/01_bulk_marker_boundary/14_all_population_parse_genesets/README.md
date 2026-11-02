# All-population ImmGen Parse genesets

This diagnostic workflow extends the replication-adaptive Parse construction from the previously approved 27 hierarchy targets to every population in the frozen ImmGen vM25 matrix.

Population-specific contrasts are performed within `table_s1_lineage`. Discovery targets use discovery libraries only. For the three formerly held-out broad populations, the target is added only to its lineage-specific discovery pool; these three derived sets are labeled and do not retain an independent validation claim.

Evidence tiers follow available RNA libraries: three-fold cross-fitting for populations with at least three libraries, two-fold leave-one-library-out intersections for two-library populations, and provisional full-data sets for singleton populations. Empty cross-fold intersections use an explicitly provisional ranked-union fallback. Nonempty intersections below 20 genes are retained as small diagnostic panels.

The workflow runs outside scworkbench, reads no single-cell data, and emits no accepted or directly importable genesets. Its GMT is a diagnostic candidate artifact.

Run from the workspace root:

```bash
/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --env R_LIBS_USER=/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/07_execute_immgen_vm25_markers/05_objects/01_r_library \
  current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif \
  Rscript current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/14_all_population_parse_genesets/01_fit_population_rankings.R \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/14_all_population_parse_genesets.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/14_all_population_parse_genesets/02_refine_and_export.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/14_all_population_parse_genesets.yaml
```
