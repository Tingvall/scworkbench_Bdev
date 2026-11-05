# 03a signal generation: standalone reconstruction of the selected 38 Parse cytokine gene sets

## Decision

The 38 passing signatures from the archived all-cytokine top-100 competitor
refinement are the selected Parse collection for later Bdev implementation.
This selection replaces the unresolved Parse-design decision in the proposed
implementation plan, but it does not activate the implementation stage.

## Development reconstruction

Run `parse38_standalone_reconstruction_20261103T052911Z` reconstructed the
collection outside scWorkbench. The workflow used:

- the frozen three-fold limma specificity fits;
- plain-text Bdev `SCT_counts` and `SCT_data` feature identifiers;
- the full-precision sample-level capped-rank cache produced by the archived
  discovery run; and
- the frozen sample-to-condition map.

The code independently rebuilt all three-fold top-100 intersections and reran
the greedy strongest-competitor deletion procedure. The archived final gene
membership and signature summary were used only as post-reconstruction
equivalence references.

## Result

| Check | Result |
|---|---:|
| Cytokine-direction attempts | 52 |
| Constructible candidates | 43 |
| Construction failures | 9 |
| Selected sets | 38 |
| Selected gene memberships | 934 |
| Distinct genes across selected sets | 753 |
| Candidate membership vs archive | exact |
| Final membership vs archive | exact |
| Selected-set identity vs archive | exact |
| Maximum absolute final-gap difference | 1.45e-13 |
| scWorkbench calls | 0 |

The full-precision rank-cache RDS is required for exact reconstruction. Using
the human-readable contribution TSV introduces enough decimal rounding in one
borderline IL-6-up deletion comparison to remove `Rtn4rl2`; this produces 31
rather than 32 genes for that signature. The validated workflow therefore
pins the RDS checksum and treats the TSV only as a display export.

## Outputs

The paired-storage result was relabeled from the provisional capability `06`
to the stable signal-generation capability `03a` without changing file
contents or checksums:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/03a_parse_cytokine_signal_generation/01_recreate_all_cytokine_genesets/parse38_standalone_reconstruction_20261103T052911Z/
```

Primary files are:

- `01_tables/05_passing_gene_sets.gmt` — the 38-set import collection;
- `01_tables/04_passing_gene_membership.tsv` — long-form membership;
- `01_tables/06_signature_summary.tsv` — all 52 attempted directions;
- `01_tables/08_archived_equivalence.tsv` — direction-level equivalence;
- `03_reports/01_reconstruction_summary.tsv` — concise result summary; and
- `06_manifests/01_input_manifest.tsv` and `02_output_manifest.tsv` — exact
  provenance and checksums.

This reconstruction supplies the gene-set definition but does not by itself
complete 03a. The remaining 03a operation is matrix-bound import followed by
UCell scoring on `RNA_counts` with `maxRank=1500`, producing the authenticated
cell-level signal family consumed by 03b. It does not justify copying or
promoting the archived bundle.
