# Parse cytokine workflow

The first step prepares the user-selected unfiltered Parse matrix outside the
bundle. It uses `Seurat::ReadParseBio()` in the supplied image, creates a Seurat
object without cell or feature filtering, preserves the delivered metadata,
and adds explicit `sample_id`, `condition`, and `replicate` fields. The output
directory contains the checksummed RDS, input inventory, sample summary,
environment manifest, and validation result.

The delivered metadata files are gzip-compressed, whereas Seurat 5.5.1 expects
the literal filenames `cell_metadata.csv` and `all_genes.csv`. The script
temporarily decompresses those two files and links the compressed count matrix
in a staging directory. It removes that staging directory immediately after
`ReadParseBio()` returns; the original delivery remains read-only.

Run it from the repository root:

```bash
processing/parse_cytokines/run_01_readparsebio_unfiltered.sh
```

The default output is outside the bundle at
`intermediate_objects/parse_cytokines/prepared/`. The
wrapper mounts the Parse source read-only and runs R from the supplied image,
which currently provides Seurat 5.5.1, SeuratObject 5.4.0, Matrix 1.7-4, and
data.table 1.18.4. A different input or empty output directory can be supplied
as the first and second positional arguments.

The unfiltered matrix declares 13,413,901 cells, 78,298 features, and
1,422,994,774 stored entries. The wrapper defaults to four computational
threads and a 160 GiB virtual-memory ceiling. On this 503 GiB shared server that
caps the process below one third of RAM. The wrapper refuses any requested
ceiling above 40% of installed server RAM. Override the lower limits with
`PARSE_THREADS` and `PARSE_MAX_MEMORY_GIB` if needed; the 40% guard remains.

The completed outside-bundle preparation is recorded in
`config/parse_cytokines.yml`. It retained all 13,413,901 cells and 78,298
features across 96 samples and passed the source, dimension, identifier-order,
metadata, and no-filtering checks. The checksummed RDS and validation evidence
remain in the external work directory named above.

Bundle creation is paused until prefiltering has been implemented at the Seurat
object stage. The prefiltering workflow will use the unfiltered RDS above,
restrict its sample scope to the 78 single-stimulus samples and three controls,
calculate factual RNA QC metrics from raw counts, and show distributions and
proposed exclusions. It must not derive accepted cells until the user selects
the QC thresholds.

The intended accepted output is
`intermediate_objects/parse_cytokines/prepared/parse_cytokines_prefiltered_seurat.rds`.
That file does not exist yet. The unfiltered source RDS remains unchanged.

No Parse scWorkbench bundle or bundle-creation script is currently present.
After the prefiltered RDS is accepted and validated, implement a new creation
script from the current scWorkbench documentation and build the bundle cleanly
from that object.

Expected phases are:

1. `01_readparsebio_unfiltered.R`
2. implement prefiltering and calculate factual RNA QC metrics
3. review distributions and proposed exclusions
4. derive the accepted cells only after user-selected thresholds
5. write and validate `parse_cytokines_prefiltered_seurat.rds`
6. implement and run a clean scWorkbench bundle creation
7. validate the bundle before downstream analysis

The future independently validated output will be the `parse_cytokines` bundle.
Reviewed cytokine gene sets may later be consumed by `processing/integration/`.
