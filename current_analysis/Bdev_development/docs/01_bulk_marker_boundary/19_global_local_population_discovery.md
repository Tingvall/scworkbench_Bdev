# Global-local ImmGen population discovery

**Step:** `17_global_local_population_genesets`  
**Run:** `immgen_vm25_20261102_global_local_population_parse_adaptive_20261102T1910Z`  
**Status:** exploratory qualification complete; no flat-label programs accepted

This ImmGen-only exploration rebuilt the 64 population programs so that every
retained gene had both cross-lineage and within-lineage evidence. It reused the
pinned target-library folds from the earlier Parse analysis, added an
equal-population/equal-lineage cross-lineage contrast, and refined each program
against its strongest population competitor. No Bdev expression data entered
gene discovery or qualification.

All 64 programs remained scoreable. On the complete 126-library atlas, the
correct population was the unthresholded winner for 126/126 libraries. In
target-membership holdout, the correct population won for 85/118 libraries.
This is a substantial improvement in winner identity over the original
population programs.

It did not make the programs compatible with the frozen default flat reducer.
Only 1/126 complete-atlas libraries and 1/118 held-out libraries cleared both
the 0.80 percentile score and 0.10 percentile margin gates. Consequently, zero
of 64 programs qualified as discrete flat-label candidates. The median
complete-atlas winner margin was 0.015873, far below 0.10.

The result isolates the remaining problem: gene membership can recover the
right winner, but 64 empirical percentile-calibrated programs yield a very
close runner-up. The collection therefore remains diagnostic score evidence,
not a replacement default and not a label collection.

No scWorkbench bundle was read, no feature set was registered, no artifact was
promoted, and no accepted set was emitted. The checksum manifest is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/01_bulk_marker_boundary/17_global_local_population_genesets/immgen_vm25_20261102_global_local_population_parse_adaptive_20261102T1910Z/06_manifests/01_run_manifest.json
SHA-256: 94edc9dce30b1be977498003fbcfcc2acbdf336a29f1855c7aba7d8848b7ea09
```
