# FLT3LKO/IL7KO mapping and group-contrast gap analysis

**Decision date:** 2026-11-17  
**Status:** bounded mapping qualified; full paper analysis locked

## Confirmed identity

The query dataset contains bone-marrow `Flt3l` ligand knockout (`FLT3LKO`)
and `Il7` knockout (`IL7KO`) cells, each with its corresponding WT library.
`FLT3LKO` must never be relabeled as an *Flt3* receptor knockout. The permanent
dataset identifier is `05_scmulti_flt3lko_il7ko`.

The immutable archived source retains 28,541 cells. The previously approved QC
rule selects 23,252 cells: 6,079 FLT3L-WT, 5,648 FLT3L-KO, 7,619 IL7-WT, and
3,906 IL7-KO. These artifacts are inputs by exact archived path and checksum;
they are not copied back into or modified within the archive.

## Reference-axis incompatibility

The archived derived KO matrix has 97,788 rows on the legacy Bdev
`peaks_counts` axis. The selected current reference is the frozen
`hvp_depth20k` PeakVI model derived from the 191,433-peak active Bdev axis.
An exact feature-ID audit found zero overlap between the 20,000 selected model
features and the 97,788 legacy matrix rows. The legacy matrix is therefore not
a valid projection input and will not be padded, renamed, or silently adapted.

The first active step requantifies only a bounded four-library fixture from the
immutable GRCm39 fragments. The exact selected mm10 peaks are lifted to GRCm39
with unique, standard-chromosome, exact reciprocal mappings. The strict first
audit accepted 19,967/20,000 peaks and stopped before fragment access: 22 were
forward-multimapped, eight were unmapped non-primary contig peaks, two lacked a
unique reverse mapping, and one did not recover the exact mm10 interval. The
33 unavailable peaks (0.165%) are restored only as explicitly masked structural
zero rows after quantifying the 19,967 exact intervals. This retains the exact
20,000-feature frozen-model order without inventing ambiguous coordinates. The
contract rejects an unavailable fraction above 0.25%.

Run `flt3lko_il7ko_mapping_fixture_20261117T175427Z` then quantified the
bounded fixture successfully. The persisted matrix is exactly 20,000 by 400,
contains 1,633,008 nonzero entries, has no zero-count cells, and retains the
exact frozen-model feature and fixture-cell orders. No archived input was
modified and no genotype comparison was performed.

The frozen bundle exposes 28 preallocated external-query slots, so projection
qualification used a balanced 28-cell subfixture (seven cells per library).
All 28 cells were projected in ten dimensions with
`transfer_without_refit`; the frozen model payload checksum remained
`024d11e3...134a8`. Bdev-only leave-one-library-out calibration established
supported and out-of-reference distance thresholds of 2.48994 and 2.93105.
All 28 queries were supported. Eighteen had finite accepted native-L4 values;
the other ten were supported mappings whose reference neighborhoods did not
provide finite L4 lineage values, so they remain missing rather than imputed.

The calibration operation uses its matrix input to authenticate the reference
cell axis and metadata. A one-feature sparse authorizer with the exact 26,723
reference-cell axis was therefore attached to avoid materializing the unrelated
20,000-peak count matrix. The full reference space and native-L4 values, not
the authorizer values, drive calibration. This eliminated the rejected 4.0 GiB
sparse-to-dense preprocessing path without a core change.

## Existing reusable capabilities

- `scwblatent::peakvi_project` projects query counts through the frozen model
  without fitting or updating it.
- `scwblatent::model_bound_latent_union` authenticates the reference/query
  linkage.
- `scwbreferencemap::reference_map_calibrate` retains Bdev-only calibration.
- `scwbreferencemap::reference_map_transfer` produces supported, borderline,
  and out-of-reference calls and withholds accepted native-L4 values for
  out-of-reference queries.

These operations are sufficient for the bounded mapping fixture. No
scWorkbench or PyBench core change is indicated.

## Replication boundary

The available metadata identify one library for each target/genotype
combination and contain no mouse, donor, or other independent biological-sample
field. Cells are not biological replicates. Consequently, an inferential
FLT3LKO-versus-WT or IL7KO-versus-WT genotype test is not currently
identifiable. Cell-level permutation, Wilcoxon, or pseudobulk tests would
pseudoreplicate the single library per group and are prohibited for the final
analysis unless independent biological-replicate provenance is supplied.

A reference-relative occupancy curve and cell-count-matched descriptive
envelope can still be qualified, but it must not be labeled a genotype-effect
test and must not publish biological p-values. The explicit contract is
recorded in `02_trajectory_occupancy_contract.yaml`.

## Remaining reusable gap

There is no operation that summarizes supported group occupancy along an
authenticated trajectory with fixed bins, a declared reference-relative
envelope, and an explicit inference mode. That contract belongs in
`scwbTrajectory` or a dedicated official contrast pack. It must remain separate
from mapping and later motif-signal contrasts. An inferential genotype mode
must require independent biological samples and reject this dataset's current
one-library-per-group design.

The bounded mapping fixture has now established query support and the available
sample structure. The next development task is to implement and qualify the
descriptive/reference-relative occupancy contract on a synthetic replicated
fixture, with a hard rejection test for unreplicated inferential contrasts.
Full KO mapping, enrichment-interval inference, motif testing, and biological
figures remain implementation work and cannot run before the replacement SIF
is qualified and explicitly accepted.

## Runtime decision

The final implementation SIF must contain the immutable, qualified post-fix
`scwbATAC` commit providing sparse-safe `hvp_select`. The unfixed `e9fd2da`
base is not acceptable for final implementation.
