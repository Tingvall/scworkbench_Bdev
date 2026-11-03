# Former working default: all-population Parse genesets

The complete 64-population Parse collection from run `immgen_vm25_20261102_all_population_parse` was the initial working default on 2026-11-02. It was subsequently superseded for Bdev development exploration by the simple ImmGen mean-contrast collection documented in [`24_simple_mean_contrast_development_default.md`](24_simple_mean_contrast_development_default.md).

This former default was intended for Bdev integration planning, integration evaluation, and population-level geneset lookup. It remains preserved as a diagnostic candidate collection and historical provenance. It has not been accepted for production or imported into scWorkbench.

The current selection is declared in `current_analysis/Bdev_development/config/01_bulk_marker_boundary/15_default_immgen_geneset_collection.yaml`.

Evaluation is deferred to Bdev integration. That evaluation must separately review singleton-provisional and held-out-derived populations, assess specificity and abstention, and inspect tissue, laboratory, state, and progenitor-program confounding. Any approved revision must be frozen as a new version rather than modifying this collection in place.
