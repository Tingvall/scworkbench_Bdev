# Replication-adaptive Parse transfer

**Capability:** `01_bulk_marker_boundary`  
**Step:** `12_replication_adaptive_parse_transfer`  
**Run:** `immgen_vm25_20261031_replication_adaptive_parse_transfer`  
**Status:** diagnostic adaptive collection complete; no accepted gene sets  
**Executed:** 2026-10-31

## Scope

This step extends the Parse-style ImmGen construction to the available
replication without changing or overwriting the completed exact-transfer run.
All membership selection uses ImmGen only.

The machine-readable definition is
[`12_replication_adaptive_parse_transfer.yaml`](../../config/01_bulk_marker_boundary/12_replication_adaptive_parse_transfer.yaml),
and the workflow is
[`12_replication_adaptive_parse_transfer`](../../workflows/01_bulk_marker_boundary/12_replication_adaptive_parse_transfer/).

## Frozen evidence tiers

| Independent target units | Evidence tier | Construction |
|---:|---|---|
| 3 or more | `high_threefold` | Reuse the completed three-fold top-100-intersection and refinement membership from step 11. |
| 2 | `moderate_twofold` | Intersect two leave-one-target-unit-out top-100 lists, require at least 20 genes, then refine on complete discovery profiles. |
| 1 | `provisional_singleton` | Use a full-data top-100 list and refine, with no reproducibility claim. |

All tiers use a positive target-versus-equal-weight-all-other effect greater
than 0.5, FDR as ranking evidence, the UCell capped-rank score with
`max_rank = 1500`, strongest-competitor refinement, and a 20-gene floor.

No single-cell expression value, label, embedding, or pseudotime was read.

## Results

The adaptive method produced panels for 25 of 27 targets:

- 14 high-tier three-fold sets, unchanged from step 11;
- eight moderate two-fold sets; and
- three singleton provisional sets.

Every constructed set ranked its target first on the complete discovery
profiles. This is a descriptive full-data result, not independent validation.

| Family | Target | Units | Tier | Stable candidate | Refined set | Result |
|---|---|---:|---|---:|---:|---|
| Broad lineage | monocyte | 2 | two-fold | 61 | 34 | constructed |
| B branch | early lymphoid progenitor | 2 | two-fold | 57 | 55 | constructed |
| B branch | committed marrow B | 2 | two-fold | 13 | 0 | below 20-gene minimum |
| B branch | mature conventional B | 2 | two-fold | 9 | 0 | below 20-gene minimum |
| B branch | B1b | 1 | singleton | 100 | 93 | provisional |
| Local B stage | `proB.CLP.BM` | 2 | two-fold | 72 | 33 | constructed |
| Local B stage | `proB.FrA.BM` | 2 | two-fold | 64 | 23 | constructed |
| Local B stage | `proB.FrBC.BM` | 2 | two-fold | 79 | 55 | constructed |
| Local B stage | `B.FrE.BM` | 2 | two-fold | 81 | 35 | constructed |
| Local B stage | `B.T1.Sp` | 1 | singleton | 100 | 32 | provisional |
| Local B stage | `B.T2.Sp` | 2 | two-fold | 84 | 24 | constructed |
| Local B stage | `B.T3.Sp` | 1 | singleton | 100 | 21 | provisional |
| Local B stage | `B.MZ.Sp` | 2 | two-fold | 88 | 59 | constructed |

The two construction failures were not caused by absence of positive genes in
the full data. They failed the frozen cross-fold stability-and-size rule:
committed marrow B retained 13 genes in both fold-specific top-100 lists, and
mature conventional B retained nine. Neither was padded or passed by lowering
the minimum after observing the result.

`B.T2.Sp` passed but has the weakest discovery separation among constructed
sets: standardized gap 0.351 against `B.T3.Sp`. It should receive particular
attention in later independent evaluation.

## Method comparison

| Method | Sets available | Full-data mean passes | Broad panel sets | Validation-only broad calls |
|---|---:|---:|---:|---:|
| Strict production | 13/27 | 12 | 4 | 1/3 |
| Positive full-data pilot | 27/27 | 22 | 9 | 3/3 |
| Positive nested stable | 8/27 | 8 | 3 | 2/3 |
| Exact three-fold Parse transfer | 14/27 | 14 | 8 | 3/3 |
| Replication-adaptive Parse transfer | 25/27 | 25 | 9 | 3/3 |

The adaptive broad panel adds the two-fold monocyte set to the exact-transfer
panel and therefore contains all nine broad lineages. It still correctly calls
the three validation-only populations: `B.Sp` as B, `MF.PC` as myeloid, and
`Tgd.Sp` as gdT. Those held-out calls do not independently validate the new
monocyte set because none of the three validation populations is monocyte.

## Interpretation and decision

The replication-adaptive construction is feasible for the current ImmGen data
and substantially improves panel coverage while preserving evidence strength
as an explicit property. It does not make the tiers equivalent:

- three-fold sets have the strongest available membership-stability evidence;
- two-fold sets are supported by only one target unit per training fold; and
- singleton sets have no biological-replication evidence and remain
  provisional even when they separate the discovery profiles.

No set is accepted, and no scWorkbench import artifact is created. The two
sub-20 B-branch candidates require an explicit future decision: retain the
20-gene floor and leave the stages unlabeled, or prospectively qualify a
smaller-set policy without using single-cell performance to choose its cutoff.

## Generated outputs

Generated artifacts are outside Git under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/12_replication_adaptive_parse_transfer/
```

| Artifact | SHA-256 |
|---|---|
| `01_tables/01_adaptive_plan.tsv` | `b0f8147b865bc59482954b4e55d210d4e822c687c00df4d0128c9019b4289258` |
| `01_tables/02_low_replication_rankings.tsv.gz` | `540a420eba6bbf845e7e939a8442117e2f42c6f9f5df68707e08b5bb9345febb` |
| `01_tables/03_low_replication_candidates.tsv.gz` | `d021c224d813dee1d542022193993f753472e6f7a7092e45cdd0af022c264423` |
| `01_tables/04_model_audit.tsv` | `db015277a7ecdc77cc72b534d36fdcb2fb5ee32b8d77f2d6529d0acdc58dab18` |
| `01_tables/05_adaptive_final_membership.tsv.gz` | `4d7d54f55222bed29aed05956d4c4f5744c676bd0c204fd79260bb437d95f8dd` |
| `01_tables/06_refinement_trace.tsv.gz` | `69cc11ee65e96a09967f07e24d095478c9d8751ddd9f767d127d192d16fe63df` |
| `01_tables/07_adaptive_target_summary.tsv` | `19e1e6caaea852e816c36b1f551dc4b0000b6e0516df6bf75afd1b6a1aa44c01` |
| `01_tables/08_method_target_comparison.tsv` | `b85e27b1741816ed22869be920b142bd8b88a77eef353a0797ad125c8dee35f6` |
| `01_tables/09_heldout_broad_comparison.tsv` | `20cc533708da97f88e58e66149faf4a83163d3da4e1c2a2532f2b8aa0796067b` |
| `01_tables/10_method_summary.tsv` | `f839da2ad8ab38e68b571ef05f3ee038027ca7e0619ee9d3453aa606562c71d0` |
| `03_reports/01_model_summary.json` | `a982d9c25a6cedd8a948e07c0d3ba7f6c3fddfe5b58e224e8e766ae1f698b72d` |
| `03_reports/02_adaptive_summary.json` | `21d3dc67f5b808b09e952ee666ef8da61ec0b267c57d2cd99e5e1c44d2e9f107` |
| `06_manifests/01_adaptive_manifest.json` | `a966a7ae26d2be2b2294f0396e6df4dc9853ee435e01742792a2b84dc3b215d0` |

The closure audit verified all configured inputs, both workflow scripts, and
all generated outputs. Every two-fold candidate exactly equals the intersection
of its two fold-specific top-100 lists; singleton candidates equal their
full-data top-100 lists; every refined low-replication set is a subset of its
candidate; and all high-tier memberships exactly match step 11.
