# Map all FLT3LKO/IL7KO cells with the qualified 20D gate

This development workflow clones the passed step-16 qualification bundle so
the qualified 20D PeakVI model and five-fold reference calibration remain
bundle-local. It ingests the already quantified 23,252-cell HVP20k query
matrix, projects without model refitting, and applies the qualified generic
`support_and_annotation` native-L4 gate in bounded, checksummed batches.

The workflow produces cell-level support, raw native-L4 estimates, accepted
values, domain votes, confidence, and explicit acceptance reasons. It performs
no genotype contrast, does not refit UMAP or Slingshot, does not write to
`Bdev_implementation`, and does not build a SIF. It is a newly numbered run
and does not overwrite the earlier 10D mapping.

Run in order:

1. `01_prepare_runtime_and_clone.py --config <config>`
2. `02_project_and_map_full.R --config <config>` inside the configured SIF and run library
3. `03_finalize.py --config <config>`

## Result

Run `flt3lko_il7ko_full_gated_20d_20261122T202618Z` passed. It projected
all 23,252 query cells through the qualified 20-dimensional model without
refitting and completed 47 checksummed transfer batches. General reference
support classified 21,885 cells as supported, 1,186 as borderline, and 181
as out of reference.

The complete gate published accepted native-L4 values for 2,974 cells and
withheld values for 20,278 cells. The withheld reasons were 19,955 local
`out_domain` votes, 181 out-of-reference calls, and 142 annotation-confidence
failures. Every cell has a raw-estimate field; 8,161 are explicit missing
values, all on rejected rows, principally where local L4 neighbors were not
available after an `out_domain` vote. No accepted row has a missing estimate,
and no rejected row publishes an accepted value. The qualified source bundle
manifest remained unchanged at
`7f6e5fb2968ec86c69513be3ebd382884a23815b787215ad6326fa7fe6eea506`.
The run performed no genotype contrast, paper analysis, model refit, or SIF
build.
