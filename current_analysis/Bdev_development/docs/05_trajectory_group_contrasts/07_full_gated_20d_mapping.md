# Full gated 20D KO mapping

The bounded 20D continuous-gate qualification passed on 2026-11-22. The next
development action maps all 23,252 already quantified FLT3LKO/IL7KO cells
through the exact qualified model and calibration. This replaces neither the
earlier 10D development evidence nor any implementation-stage artifact.

The full run retains a raw native-L4 estimate field for every query and
publishes an accepted value only when general reference support, local
`in_domain` membership, 0.51 confidence, and finite-value checks all pass.
Only descriptive mapping counts are allowed; genotype inference remains
prohibited.

## Result

Run `flt3lko_il7ko_full_gated_20d_20261122T202618Z` passed on 2026-11-22.
All 23,252 query cells were projected through the frozen qualified 20D model
without refitting and were accounted for across 47 bounded, checksummed
transfer batches.

General reference support called 21,885 cells supported, 1,186 borderline,
and 181 out of reference. The complete acceptance gate retained native-L4
values for 2,974 cells and withheld 20,278 values: 19,955 because the local
annotation vote was `out_domain`, 181 because the query was out of reference,
and 142 because annotation confidence was below 0.51. Every cell retained a
raw-estimate field and acceptance-ledger row. Of those fields, 8,161 contain
an explicit missing value, all on rejected rows and principally following an
`out_domain` vote. No accepted row has a missing estimate, and no rejected row
publishes an accepted value.

Descriptive accepted/withheld counts were 1,129/4,950 for FLT3LKO BM WT,
449/5,199 for FLT3LKO BM KO, 1,033/6,586 for IL7KO BM WT, and 363/3,543 for
IL7KO BM KO. These are mapping-QC counts, not genotype contrasts or biological
effect estimates.

The qualified source bundle manifest remained exactly
`7f6e5fb2968ec86c69513be3ebd382884a23815b787215ad6326fa7fe6eea506`.
The finalized run inventory is complete, and the mapping performed no paper
analysis, genotype contrast, model refit, or SIF build. The next development
decision is whether this coverage is sufficient for the planned downstream
trajectory occupancy contract; implementation remains locked until all
framework work is accepted and a new immutable SIF is built and qualified.
