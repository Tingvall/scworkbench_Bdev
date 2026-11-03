# Percentile-margin ImmGen population refinement

**Step:** `18_percentile_margin_population_genesets`  
**Run:** `immgen_vm25_20261102_percentile_margin_population_refinement`  
**Status:** exploratory negative result; no flat-label programs accepted

This follow-up tested whether the remaining incompatibility was caused by
optimizing raw UCell separation while the default flat reducer consumes
empirical percentiles. Starting from the global-local fold rankings, it used
repeated greedy coordinate refinement to maximize each target program's mean
percentile margin over the strongest competing program on its target
libraries. The same optimization was repeated inside each target-membership
holdout. The default 0.80 score and 0.10 margin gates were unchanged.

The optimization converged in two complete-data coordinate passes and removed
36 genes across the 64 complete programs. It did not produce a
default-compatible label collection:

- 103/126 complete-atlas libraries had the correct unthresholded winner, but
  0/126 were assigned; 100 failed the margin gate and 26 tied.
- 77/118 target-membership holdouts had the correct winner, but 0/118 were
  assigned correctly.
- The median complete-atlas margin was 0.015873 and the maximum was 0.079365,
  so no library reached the required 0.10 margin.
- Zero of 64 programs qualified as flat-label candidates.

Direct percentile-margin pruning therefore does not rescue flat 64-way
labeling and is worse for winner identity than the preceding global-local
programs. Do not adopt the pruned programs. Retain the original 64-set
collection as the score/lookup default and treat either a coarser ontology or
a genuinely joint multiclass discovery method as a separate future
exploration.

No Bdev data or scWorkbench bundle was read, no feature set was registered, no
artifact was promoted, and no accepted set was emitted. The checksum manifest
is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/01_bulk_marker_boundary/18_percentile_margin_population_genesets/immgen_vm25_20261102_percentile_margin_population_refinement/06_manifests/01_run_manifest.json
SHA-256: b3f8dc5f29e634a5ec4f6779aa5d8522a4ece0b3596c0e9aef4afd77a57343ba
```
