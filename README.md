Code related to the Bdevelopment use case for the scWorkbench paper. 

More specifically this project is focused on:
1. Cytokine stimulation parse data preprocessing and analysis: Experiemnt run by Vathsa (LiU) on a broad population of progenitors. The parse data is used for extracting specific cytokine signals that should be tracked along the Bdev developmental trajectory.
2. IL7KO & FLT3KO scMulti data preprocessing and analysis: scMultiome (RNA & ATAC) on IL7KO and FLT3KO mice. Popualtion sorted on cells expressing respective receptor i.e. Il7R+ cells in Il7KO mice). Data will be autoencoder based mapped ontp the reference Bdev data set, as validation for the parse cytokine signal analaysis.

## Bdev reference bundle

The reference source is the `Bdev` Seurat/Signac object in
`/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/datasets/scmulti_bdev/Bdev.Rdata`.
Generated data remain in mirrored storage under
`/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/`:

- `bdev_reference/`: the scWorkbench bundle.
- `bdev_reference_inspection/complete/`: the read-only source inventory.
- `bdev_reference_import/`: source identity, import mappings, and validation reports.

`processing/01_inspect_bdev_reference.R` inspects the original object.
`processing/02_create_bdev_reference.R` creates a new bundle using scWorkbench
inside the supplied image; its arguments are the source file, bundle directory,
and report directory. Mount the source read-only and the bundle's parent
directory writable. The script checks the source checksum and refuses to
overwrite an existing bundle. General scWorkbench instructions remain inside
the image; see `AGENTS.md` for its location.

The import retains the nine available assay/layer matrices, all 68 metadata
columns, and `pca`, `lsi`, `umap.atac`, `umap.rna`, and `umap.wnn` coordinates.
Original assay/layer names, reduction names, feature and cell identifiers, and
the source checksum are recorded. Metadata factor labels are imported as text.
The existing `pt` metadata defines pseudotime: its 12,642 defined values form
the trajectory vector, while all 14,081 missing values remain missing. The
metadata retains the exact original numeric values. This image converts
cell-vector values to text internally, causing rounding of at most
`5e-14` in the derived pseudotime artifact; its ordering and missingness are
checked against the original. Use the metadata `pt` for exact source values.

`src/reference_validation.R` checks all metadata and reduction coordinates,
matrix dimensions and identifier ordering, and 128 distributed features per
matrix across every cell. The creation script also verifies the complete `pt`
vector and runs scWorkbench recipient validation. `completion.json` in the
import report directory records successful completion.
`processing/03_validate_bdev_reference.R` accepts the same three arguments to
validate and finalize an existing import without recreating its data or vector.
