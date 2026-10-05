# GSE259355 TCF3/E2A pro-B ATAC-seq

This source workflow analyzes the six user-selected bone-marrow pro-B ATAC-seq
libraries from GSE259355: three untreated controls and three samples collected
after acute auxin-inducible E2A/TCF3 degradation. The depleted group contains
two 6-hour samples and one 8-hour sample, exactly as recorded in the manifest.

The raw SRA archives, mm10 reference, alignments, sample peaks, consensus peak
matrix, and differential results remain outside Git in mirrored project
storage. The repository contains only configuration, identifiers, parameters,
and reproducible scripts.

`run_tcf3_da.sh` performs:

1. SRA retrieval and validation.
2. Paired-end trimming and mm10 alignment.
3. MAPQ, proper-pair, duplicate, primary-chromosome, mitochondrial, and
   blacklist filtering.
4. Per-sample MACS3 peak calling and construction of a merged consensus peak
   universe.
5. Fragment counting with featureCounts.
6. User-approved outside-scWorkbench DESeq2 analysis with design
   `~ condition`, where log2 fold change is depleted/control.

The bulk replicate-aware DESeq2 step is outside scWorkbench because the image
does not provide a replicate-aware bulk differential-accessibility model.
Afterward, the raw peak-count matrix and sample metadata are ingested into a
separate scWorkbench source bundle; differential results remain a separately
recorded report artifact.

The workflow completed on 2026-10-05. The merged universe contains 84,088
peaks; DESeq2 tested 84,083 and identified 7,432 at BH-adjusted P < 0.05.
Of these, 6,096 have lower accessibility after depletion and are classified as
TCF3-dependent, while 1,336 have higher accessibility after depletion. Results
are under `reports/tcf3_depleted_vs_control_da_v1/` in mirrored storage.

`run_create_source_bundle.sh` completed the scWorkbench handoff. The accepted
bundle contains the raw 84,088-by-6 consensus-peak count matrix and metadata
that explicitly labels its columns as bulk ATAC libraries; scWorkbench's
generic column/cell APIs do not change their biological interpretation.
Recipient validation passed with only the expected external-source trust-
boundary warnings. Evidence is under `reports/source_bundle_v1/`.
