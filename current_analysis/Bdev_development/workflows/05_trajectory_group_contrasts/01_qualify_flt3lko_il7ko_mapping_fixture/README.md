# Qualify the bounded FLT3LKO/IL7KO mapping fixture

This development workflow prepares and maps 100 lexicographically selected
QC-passing cells from each of the four confirmed libraries:

1. audit immutable input checksums, extract the exact `hvp_depth20k` feature
   order, and freeze the 400-cell fixture;
2. reciprocally lift the exact selected mm10 peak axis to GRCm39, quantify only
   uniquely and exactly represented intervals, and retain a mask for the 33
   unavailable model features found by the strict first-pass audit;
3. quantify the immutable query fragments over those intervals into a new
   active-development matrix, restoring unavailable features as explicit zero
   rows in exact frozen-model order (maximum allowed unavailable fraction
   0.25%);
4. import the bounded matrix into a fresh derived bundle, project it through
   the frozen PeakVI model without refitting, and authenticate the model-bound
   latent union;
5. apply retained Bdev-only support calibration and transfer native L4 only
   for supported or borderline queries; and
6. inventory all outputs and verify that source manifests and archived inputs
   are unchanged.

The workflow is an acceptance fixture only. It does not map all knockout
cells, infer enrichment intervals, compare genotypes, generate paper figures,
or write to `Bdev_implementation`.

Status: passed on run
`flt3lko_il7ko_mapping_fixture_20261117T175427Z`. The 400-cell quantified
fixture passed exact-axis gates; the balanced 28-cell projection subfixture
passed frozen-model projection and calibrated support transfer. All 28 cells
were supported and 18 retained finite native-L4 values. No missing L4 value was
imputed.

Configuration:

```text
current_analysis/Bdev_development/config/05_trajectory_group_contrasts/01_flt3lko_il7ko_mapping_fixture.yaml
```

Input declaration:

```text
current_analysis/Bdev_development/inputs/05_trajectory_group_contrasts/01_flt3lko_il7ko_inputs.yaml
```
