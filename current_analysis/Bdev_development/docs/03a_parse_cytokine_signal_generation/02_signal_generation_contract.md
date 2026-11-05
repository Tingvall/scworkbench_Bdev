# 03a Parse cytokine signal-generation contract

**Status:** planned development qualification; no paper result is active

## Purpose

Stage 03a generates the authenticated cell-level cytokine signals consumed by
the trajectory methods in 03b. It has two parts:

1. reconstruct and checksum the selected 38 Parse cytokine gene sets; and
2. import those sets against `RNA_counts` and generate one UCell
   `cell_vector_family` using `scwbgeneprograms::ucell`.

The completed standalone reconstruction is documented in
`01_standalone_parse38_reconstruction.md`. Its selected GMT has SHA-256
`22f041528a1e882c0bec7bd920b7638646b20589020669bf0f42b59880e0bdb7`.

## Trajectory and bundle boundary

The authenticated fitted Slingshot lineage `L4` from the approved ATAC-UMAP
recipe is the trajectory source. Native L4 orientation is retained. Preserved
`Bdev@meta.data$pt` is only a validation reference or documented fallback.

The existing L4 recreation result remains unchanged as the parent evidence.
Create at most one writable derived analysis bundle for the complete 03a/03b
chain. Generate the UCell family in that bundle and run all three 03b methods
against the same artifact identities. Do not create a bundle per method.

A writable child is needed because feature-set import, UCell output, sessions,
and operation provenance change bundle state. Writing those records into the
already reported L4 qualification bundle would make its recorded fingerprint
and result no longer reproducible. This child is an execution state derived
from the same Bdev data, not a new biological dataset or another normalization.

## Signal-generation parameters

| Item | Primary contract |
|---|---|
| Matrix | `RNA_counts` |
| Feature-set namespace | versioned Parse38 namespace bound to `RNA_counts` |
| Import source | checksummed 38-set GMT from step 01 |
| Overlap gate | report all requested, overlapping, and absent genes per set before scoring |
| Scoring operation | `scwbgeneprograms::ucell` |
| `maxRank` | `1500` |
| Output | all 38 raw scores in one authenticated `cell_vector_family` |
| Reduction | none; do not select a winning cytokine program |
| Cell axis | cells eligible for native fitted L4, with exclusions recorded |

The same GMT must also be imported as a separate namespace bound to
`SCT_data` for 03b preranked GSEA. That registry import does not create a
second UCell signal family.

## Required records

- parent and derived bundle fingerprints;
- fitted L4 pseudotime and membership artifact identities;
- GMT path and SHA-256;
- matrix identities and feature-overlap audit;
- feature-set import manifests for `RNA_counts` and `SCT_data`;
- UCell operation identity, parameters, seed/backend information, and output
  artifact IDs; and
- exact cell inclusion and exclusion counts.

