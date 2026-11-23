# Define and qualify trajectory group occupancy

This step follows the passed FLT3LKO/IL7KO bounded mapping fixture. Its
configuration is:

```text
current_analysis/Bdev_development/config/05_trajectory_group_contrasts/02_trajectory_occupancy_contract.yaml
```

Implement the generic operation in the owning official operation pack and
qualify it on synthetic fixtures before using any biological dataset. The
qualification must demonstrate exact cell accounting, deterministic
reference-relative envelopes, missing-pseudotime preservation,
out-of-reference exclusion, and hard rejection of an inferential request with
only one biological sample per group.

The current FLT3LKO/IL7KO data may be used only as a bounded structural
fixture. Full mapping, interval inference, genotype testing, motif analysis,
and paper figures remain locked until the replacement SIF is qualified,
checksummed, and explicitly accepted.
