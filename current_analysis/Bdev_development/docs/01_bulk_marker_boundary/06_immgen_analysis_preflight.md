# ImmGen analysis preflight freeze

**Capability:** `01_bulk_marker_boundary`  
**Bdev step:** freeze annotation, state-gene flags, and software  
**Status:** passed  
**Frozen:** 2026-10-31

## Outcome

The approved statistical contract now has executable, checksummed inputs. The
GSE109125 GENCODE vM25 count matrix and annotation remain in paired Bdev
development storage, outside scWorkbench. The preflight performed no
normalization, differential expression, resampling, or marker selection.

The machine-readable definitions are:

- [`03_immgen_analysis_preflight.yaml`](../../config/01_bulk_marker_boundary/03_immgen_analysis_preflight.yaml);
- [`04_immgen_bulk_analysis_software_lock.yaml`](../../config/01_bulk_marker_boundary/04_immgen_bulk_analysis_software_lock.yaml); and
- the reproducible
  [`04_freeze_immgen_analysis_preflight`](../../workflows/01_bulk_marker_boundary/04_freeze_immgen_analysis_preflight/README.md)
  workflow.

## Annotation compatibility

The complete GENCODE mouse M25 GTF for GRCm38.p6 was selected because it
contains every one of the 55,011 unique matrix symbols. Mapping first prefers a
unique record on `chr1`–`chr19`, `chrX`, `chrY`, or `chrM`, then permits a
single non-reference record only when no reference-chromosome record exists.

The frozen result contains:

- 54,828 unique reference-chromosome mappings;
- 76 unique non-reference mappings; and
- 107 ambiguous reference-chromosome symbols, which remain in complete
  statistics but cannot enter exported gene sets.

This resolves the apparent mismatch with the smaller reference-chromosome-only
GENCODE file, which omits 76 matrix symbols. The complete annotation is the
source file; the reference-only annotation is retained as a checksummed audit
input.

## State and identity exclusions

State flags are derived deterministically from `org.Mm.eg.db` 3.23.0, whose
Gene Ontology source is dated 2026-01-23, plus the exact 17-gene
warm-dissociation stress score reported by Denisenko et al. The union flags
2,245 matrix genes:

| Category | Matrix genes |
|---|---:|
| Cell cycle | 1,693 |
| Interferon response | 239 |
| Heat shock | 60 |
| Unfolded-protein or oxidative stress | 388 |
| Published dissociation-stress set | 17 |
| Immediate-early subset | 12 |

Categories overlap. These genes remain visible in complete statistics but are
not eligible for accepted core identity sets. This deliberately avoids using
the entire broad Gene Ontology `response to stress` branch, which would remove
many genuine lineage functions.

The separate mandatory identity filters flag Y-linked genes, `Xist`,
mitochondrial genes, GO-defined ribosomal and hemoglobin genes, and GENCODE
immunoglobulin/TCR V, D, and J segments. After mapping, protein-coding,
state-gene, and identity filters, 19,327 matrix symbols are eligible to compete
for core identity sets. Expression, effect, FDR, support, stability, and
validation requirements will reduce this number during marker analysis.

## Software freeze

The outside analysis uses the immutable development SIF only as an R toolchain;
it does not invoke scWorkbench. The locked runtime is R 4.5.2, Bioconductor
3.21, edgeR 4.6.3, limma 3.64.3, and locfit 1.5.9.12. The edgeR and locfit
library is stored as a deterministic archive with SHA-256:

```text
3aed57f0a07c4da25de1de689f1f42e3c22c17f90be1e930b298baca55b00087
```

Each analysis run must verify and extract this archive into its own run-owned
library directory. The mutable package-build directory is not an executable
input.

A synthetic-only smoke test confirmed package loading, robust `glmQLFit`, and
`glmTreat` with the contract's log2-fold-change threshold. No ImmGen model was
used for the smoke test.

## Generated artifacts

Generated output is under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/04_freeze_immgen_analysis_preflight/
```

| Artifact | SHA-256 |
|---|---|
| `01_tables/01_gene_annotation.tsv.gz` | `b7c909ecb1e2e83c749a26501415a9ea4d3c163f8861491d662ba9fe065c5f6f` |
| `01_tables/02_state_gene_flags.tsv` | `993197999bcb33aa2937c6991229963302140d249b8d69f0bc72e26acf2e5396` |
| `03_reports/01_preflight_summary.json` | `8a7d0966bc208a031242d4d0331025498e01c1d9b18494512b125e4594d150c7` |
| `06_manifests/01_preflight_run_manifest.json` | `d82ca2df2a1f25a17390fac3f4b9c36cc1d5c5822f53415d9602b7e0d4d312da` |

A repeat run reproduced all four hashes byte-for-byte.

## Next step

The model and deterministic-selection core passed the synthetic qualification
documented in
[`07_immgen_marker_workflow_core_qualification.md`](07_immgen_marker_workflow_core_qualification.md).
The 200 balanced stability resamples, leave-one-population/library-out checks,
tissue and laboratory confounding audits, and held-out scoring gates then
passed the qualification documented in
[`08_immgen_deferred_gate_qualification.md`](08_immgen_deferred_gate_qualification.md).
The real execution subsequently completed and is documented in
[`09_immgen_vm25_marker_execution.md`](09_immgen_vm25_marker_execution.md).
No gene-set import was produced; diagnostic and biological review is required.

## Sources

- GENCODE mouse M25 release:
  <https://www.gencodegenes.org/mouse/release_M25.html>.
- Bioconductor release compatibility:
  <https://www.bioconductor.org/about/release-announcements/>.
- Denisenko et al., *Systematic assessment of tissue dissociation and storage
  biases in single-cell and single-nucleus RNA-seq workflows*:
  <https://doi.org/10.1186/s13059-020-02048-6>.
