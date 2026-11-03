# Discover globally and locally specific ImmGen population programs

This ImmGen-only exploratory workflow reconstructs the 64 population programs
for flat scWorkbench competition. It reuses the pinned cross-fitted local edgeR
rankings from step 14, then requires a population-balanced cross-lineage log2-CPM
effect before global refinement. Target-membership holdouts emulate the default
percentile, score, margin, and abstention rules.

Run from the repository root:

```bash
python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/17_global_local_population_genesets/01_prepare_joint_rankings.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/18_global_local_population_genesets.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/17_global_local_population_genesets/02_refine_and_qualify.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/18_global_local_population_genesets.yaml
```

No Bdev single-cell data are read. Both GMT files remain diagnostic candidates;
nothing is imported or promoted in scWorkbench.

An earlier full-atlas edgeR-per-target formulation was stopped before its first
fit completed because repeating the 61-population dispersion model for every
target fold was computationally disproportionate. Its isolated fold-plan output
is retained under the unversioned step-17 intermediate root as retired evidence.
