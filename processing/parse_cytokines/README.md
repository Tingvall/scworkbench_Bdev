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

The second step applies the user-selected prefilter in Seurat:

```bash
processing/parse_cytokines/run_02_prefilter_min500counts.sh
```

It calculates each cell's total directly from the raw RNA `counts` layer and
retains cells whose total is at least 500. It applies no sample, feature, or
additional cell filter. The completed output is:

`intermediate_objects/parse_cytokines/prepared/parse_cytokines_prefiltered_min500counts_seurat.rds`

The filter retained 251,599 of 13,413,901 cells, all 78,298 features, and all
96 samples. The minimum retained raw RNA count is exactly 500. The script
verified that `nCount_RNA` exactly matched the raw counts-layer column sums,
retained cell order is the input subset order, feature order is unchanged,
cell and feature metadata are unchanged, every retained cell meets the
threshold, and every removed cell is below it. The output SHA-256 is
`1e4da7149190e6c93f71bc28fc47d1f58a20e22ee9cbec6894f0824394e9bc40`.
The unfiltered source remained unchanged.

Validation, per-sample counts, environment versions, and session information
are written beside the filtered RDS with the same filename stem. The previously
reviewed single/combination/control classification was not used as a filter in
this count-only prefilter.

No Parse scWorkbench bundle or bundle-creation script is currently present.
Before bundle creation, consult the current documentation inside the supplied
image, implement a new creation script, and build the bundle cleanly from the
prefiltered object.

Workflow phases are:

1. `01_readparsebio_unfiltered.R`
2. `02_prefilter_min500counts.R`
3. implement and run a clean scWorkbench bundle creation
4. validate the bundle before downstream analysis

The future independently validated output will be the `parse_cytokines` bundle.
Reviewed cytokine gene sets may later be consumed by `processing/integration/`.
