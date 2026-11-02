# Small-panel exception for the adaptive Parse transfer

**Capability:** `01_bulk_marker_boundary`  
**Step:** `13_small_panel_parse_exception`  
**Run:** `immgen_vm25_20261031_small_panel_parse_exception`  
**Status:** complete 27-target diagnostic collection; no accepted gene sets  
**Executed:** 2026-10-31

## Scope and policy

Following explicit user authorization to allow fewer genes while preserving
the limitation, this step adds a small-panel exception to the completed
replication-adaptive collection. It does not change or overwrite prior runs.

The machine-readable definition is
[`13_small_panel_parse_exception.yaml`](../../config/01_bulk_marker_boundary/13_small_panel_parse_exception.yaml),
and the workflow is
[`13_small_panel_parse_exception`](../../workflows/01_bulk_marker_boundary/13_small_panel_parse_exception/).

The exception applies only when a moderate two-fold candidate:

- is the complete intersection of two leave-one-target-unit-out top-100 lists;
- contains at least five but fewer than 20 genes; and
- failed the adaptive run only because of the 20-gene floor.

Five genes is the prespecified minimum used by the earlier archived Parse
discovery workflow. The complete stable intersection is retained without
further greedy deletion. Every affected set is labeled
`moderate_twofold_small_panel`.

No single-cell expression value, label, embedding, or pseudotime was read.

## Results

Both previously missing B-branch panels meet the exception:

| Target | Units | Genes | Discovery gap | Target units above strongest-competitor median | Strongest competitor |
|---|---:|---:|---:|---:|---|
| Committed marrow B | 2 | 13 | 0.229 | 1/2 | Early lymphoid progenitor |
| Mature conventional B | 2 | 9 | 1.999 | 2/2 | Splenic transitional B |

Both aggregate target medians rank first on the complete discovery data. The
committed-marrow-B result is weak: only one of its two populations exceeds the
strongest-competitor median and its standardized gap is small. Its 13-gene
panel must therefore remain conspicuously lower confidence. The nine-gene
mature-conventional-B panel separates the two available target populations
more strongly, but still has only two-fold evidence and a small membership.

The resulting diagnostic collection contains all 27 targets:

- 14 high three-fold sets;
- eight standard moderate two-fold sets;
- three provisional singleton sets; and
- two moderate two-fold small panels.

All 27 aggregate discovery medians rank their target first. The broad-lineage
memberships are unchanged from step 12 and continue to call all three
validation-only broad populations correctly.

## Decision

The lower-count allowance makes the collection complete, but it does not make
the evidence equal. In particular:

- the two small panels must keep their gene counts and evidence label;
- the three singleton panels still have no replication evidence;
- committed marrow B is fragile even within the two available populations;
  and
- all full-data score results remain descriptive and optimistic.

No set is accepted, no import artifact is created, and scWorkbench is not
invoked. Biological/confounding review and a frozen evaluation policy remain
required before any single-cell evaluation.

## Generated outputs

Generated artifacts are outside Git under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/13_small_panel_parse_exception/
```

| Artifact | SHA-256 |
|---|---|
| `01_tables/01_small_panel_membership.tsv` | `ed66cff882c4d741a84a13ced69c0d5d505fec09880578d3cade734aa7bc4c8a` |
| `01_tables/02_complete_membership.tsv.gz` | `c06bc05410c9fc30800b77bc23f21080ceb5aaeee753dae81a16735d0610b633` |
| `01_tables/03_target_summary.tsv` | `3687600a7777ae342c0491f0336b484f7265ec6dfb5018ce30db845aeacfdcd2` |
| `01_tables/04_method_target_comparison.tsv` | `4f060f7f2d629fb8afa272f24c3edbe64fb0ccf1e1782ea47ef759c954fdfa4d` |
| `01_tables/05_heldout_broad_comparison.tsv` | `800eb32fd1e3a193f9c597c5354fb14c42556e3cd1526a1117de1071400037ca` |
| `01_tables/06_method_summary.tsv` | `42b4228a40e00dcac084722bf14538960cb1a6332ee99d1acff930a782cdef2e` |
| `03_reports/01_summary.json` | `1b628c24fbef75066d208b18312bb820915119fa2224a9d88ef72643c1093b94` |
| `06_manifests/01_small_panel_manifest.json` | `f495b728f653a8c37480667f099d0756c9a9f9a2d5ad0893a64011167bda16cf` |

The closure audit verified all inputs, the workflow script, and all outputs.
The two small memberships exactly equal their previously generated two-fold
stable intersections, and the other 25 memberships exactly match step 12.
