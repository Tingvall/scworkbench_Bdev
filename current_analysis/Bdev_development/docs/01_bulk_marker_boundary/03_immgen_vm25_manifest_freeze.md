# ImmGen GENCODE vM25 manifest freeze

**Capability:** `01_bulk_marker_boundary`  
**Bdev step:** freeze the outside-analysis RNA input cohort  
**Status:** passed  
**Frozen:** 2026-10-31

## Outcome

The GSE109125 GENCODE vM25 integer-count matrix is accepted as the primary
expression source for the ImmGen steady-state gene-set analysis. Discovery is
restricted to the approved 64 Table S1 populations that have
population-matched ATAC data.

The frozen column manifest contains all 205 vM25 libraries:

- 126 libraries from the 64 approved populations are marked `discovery`;
- 39 libraries from excluded Table S1 populations are rejected because of
  perturbation, engineered background, or explicit sex-specific sampling; and
- 40 libraries from later populations outside paper Table S1 are rejected.

No excluded library is eligible for gene-set discovery. The manifest records
the exclusion reason for every rejected column.

## Validation

The freeze workflow verified that:

- all five source checksums match the approved policy;
- the matrix contains 55,011 unique gene-symbol rows and 205 unique library
  columns;
- all count values are non-negative integers, with no malformed rows;
- all 205 matrix columns map one-to-one to RNA GSM accessions;
- the selected set contains exactly 126 libraries and 64 populations;
- no unapproved population is selected; and
- every selected population has a Table S1 ATAC match with its ATAC GSM
  accessions recorded.

RNA and ATAC are matched at the population level only. Their replicate numbers
must not be interpreted as paired libraries.

## Frozen artifacts

- policy and source checksums:
  `../../inputs/01_bulk_marker_boundary/01_immgen_vm25_manifest_policy.json`;
- complete column manifest:
  `../../inputs/01_bulk_marker_boundary/02_immgen_vm25_library_manifest.tsv`;
- validation summary:
  `../../inputs/01_bulk_marker_boundary/03_immgen_vm25_manifest_summary.json`;
- reproducing workflow:
  `../../workflows/01_bulk_marker_boundary/01_freeze_immgen_vm25_manifest/`.

The manifest SHA-256 is
`34bebcdea2cc2e94c6e9da6febb351ecf0d73be9345e338902077c0ebb53dc71`.

## Remaining limitation

The vM25 matrix changes replicate coverage relative to paper Table S1. Five
selected populations gain one library and six lose one, producing 126 selected
vM25 libraries rather than 127 paper-listed libraries. These differences are
explicit in the validation summary and must be respected by the statistical
design.

## Next step

The replicate-aware statistical contract and analysis preflight are now frozen
and validated. The next step is to implement and qualify the marker-analysis
workflow without changing this manifest. No ImmGen differential-expression
model has yet been run.
