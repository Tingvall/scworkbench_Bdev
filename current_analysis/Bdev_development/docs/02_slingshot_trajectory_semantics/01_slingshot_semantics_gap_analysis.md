# Slingshot trajectory-semantics gap analysis

**Capability:** `02_slingshot_trajectory_semantics`  
**Date:** 2026-11-02  
**Status:** gap analysis complete; no source or bundle changed

## Result

The published Bdev Slingshot fit cannot currently be reproduced exactly through
`scwbtrajectory::slingshot_fit`. The operation already supports the ATAC UMAP,
cluster labels, cluster exclusion, start cluster, `stretch = 0`, a contained
fit, and lineage-specific pseudotime. It does not expose `omega` or
`omega_scale`, although the installed Slingshot 2.16.0 method accepts both.

The published recipe filters clusters 9, 15, 18, 19, 22, and 23, starts at
cluster 5, uses `stretch = 0`, `omega = TRUE`, and `omega_scale = 5`, selects
lineage 4, aligns it back to the full cell axis, and stores its negative as
`Bdev$pt`. The source is the upstream
[`scMulti_analysis.Rmd`](https://github.com/Tingvall/Bdev_epigenetics/blob/main/script/scMulti/scMulti_analysis.Rmd).
An exact source commit must be pinned before qualification rather than relying
on the moving `main` URL.

There is a second, independent interface gap. The authoritative `pt` field is a
scalar `cell_vector`, while the explicit trajectory-dynamics operations consume
a lineage-keyed `cell_vector_family`. No current `scwbTrajectory` operation
adapts a preserved scalar vector into a one-lineage family. The existing
downstream `pseudotime_direction = "invert"` option computes
`max(pseudotime) - pseudotime`; this preserves the reverse ordering but is not
numerically identical to the published `-pseudotime` transformation.

## Selected trajectory-source policy

Fitted lineage-family pseudotime is the default. Downstream analyses should use
the authenticated fitted model or its lineage-specific pseudotime family and
select the required lineage key explicitly.

The preserved scalar `pt` vector is not the default trajectory source. It is a
fallback when fitted lineage output is unavailable, a source when specifically
requested, and a validation reference for the published Bdev trajectory. Its
values and missing entries remain unchanged in every case.

## Smallest reusable change

Keep both changes in the existing official `scwbTrajectory` pack:

1. Extend `slingshot_fit` with typed `omega` and `omega_scale` parameters and
   pass them unchanged to `slingshot::slingshot`. Omitting them must preserve
   current behavior.
2. Add fallback operation `trajectory_pseudotime_from_cell_vector`, which wraps one authenticated
   continuous `cell_vector` as a one-key `cell_vector_family`. It must not
   refit, orient, negate, rescale, shift, or impute the values. It should retain
   source-artifact provenance and report total, finite, and missing counts.

These are separate operations because refitting a trajectory and adapting an
already selected fallback scalar axis have different scientific meanings. No
scWorkbench core or PyBench change is indicated.

## Why existing adoption is insufficient

Core operation `slingshot_adopt` can authenticate and import a complete
external Slingshot fit without refitting. The archived Bdev inspection found
the scalar `pt` field but no retained complete Slingshot fit. Adoption therefore
cannot recover the original fit, and it does not adapt the preserved scalar
field for lineage-family consumers. The existing core operation remains
unchanged.

## Qualification boundary

The fit change needs interface, forwarding, model-parameter, deterministic
outgroup, negative-validation, replay, and atomic-promotion tests. The scalar
adapter needs exact-value, missingness, lineage-key, source-immutability,
downstream-consumer, and failure-without-mutation tests.

Only after those pack tests pass should a bounded Bdev validation fit use the
published ATAC UMAP and cluster recipe. That validation compares `-L4` with the
preserved `pt` on the exact shared finite cell set and reports cell membership,
missingness, numerical differences, and rank correlation. It must never
overwrite `pt`.

The machine-readable audit is
[`01_slingshot_semantics_gap_analysis.yaml`](../../config/02_slingshot_trajectory_semantics/01_slingshot_semantics_gap_analysis.yaml).
