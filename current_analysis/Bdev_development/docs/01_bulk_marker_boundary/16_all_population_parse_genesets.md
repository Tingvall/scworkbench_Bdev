# Working default ImmGen genesets

As of 2026-11-02, the complete 64-population Parse collection from run `immgen_vm25_20261102_all_population_parse` is the working default ImmGen geneset version.

The default is intended for Bdev integration planning, integration evaluation, and population-level geneset lookup. It replaces the earlier 27-target collection for working use, but it remains a diagnostic candidate collection. It has not been accepted for production and has not been imported into scworkbench.

The default selection is declared in `current_analysis/Bdev_development/config/01_bulk_marker_boundary/15_default_immgen_geneset_collection.yaml`. The declaration pins the run configuration, GMT, membership table, target summary, and run manifest by SHA-256.

Evaluation is deferred to Bdev integration. That evaluation must separately review singleton-provisional and held-out-derived populations, assess specificity and abstention, and inspect tissue, laboratory, state, and progenitor-program confounding. Any approved revision must be frozen as a new version rather than modifying this collection in place.
