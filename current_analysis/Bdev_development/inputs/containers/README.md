# Development containers

Container binaries are stored outside Git. The corrected `20261029` candidate
SIF is used as the development toolchain baseline for the scWorkbench update.
It supersedes the `20261022` candidate, which did not provide the required CPU
PeakVI runtime. Both images remain immutable at their declared storage paths;
do not modify or overwrite either one.

Develop and test source changes against the corrected SIF, then build a new,
uniquely named SIF by extending that immutable base with the qualified updated
packs. Reconstructing the corrected baseline from its source commits is not a
development requirement.

New development builds should use unique names under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/containers/
```

For every build, retain the definition, resolved Git commits, dependency locks,
build log, test report, filename, architecture, and SHA-256. Only a qualified
immutable build may be declared as an input to `Bdev_implementation`.
