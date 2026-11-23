# Trajectory group-occupancy operation contract

**Decision date:** 2026-11-17  
**Status:** contract defined; operation implementation pending

The next reusable capability is a support-aware group-occupancy operation for
an authenticated trajectory. It will consume explicit query pseudotime, group,
mapping-support, and reference-pseudotime primitives. It will preserve missing
trajectory values, exclude out-of-reference cells with complete accounting,
and publish fixed-bin occupancy curves, a declared reference-relative envelope,
candidate intervals, and an accounting summary.

The reference envelope is descriptive. It asks where a same-sized draw from
the Bdev reference trajectory would lie; it is not a WT-versus-KO null model.
Every random draw count, seed, bin boundary, smoothing window, denominator, and
threshold must be retained.

An inferential mode is reserved for data with independent biological samples.
It must reject fewer than two samples per group before computing statistics.
The current FLT3LKO/IL7KO metadata contain one library per target/genotype and
no mouse or donor field, so cells cannot be substituted as replicates. The
current dataset may qualify structural accounting, but it cannot qualify or
run an inferential genotype contrast without additional provenance.

The operation belongs in `scwbTrajectory` unless implementation reveals a
coherent broader contrast domain that warrants a dedicated official pack. No
scWorkbench or PyBench core change is indicated.
