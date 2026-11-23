# Development results

Retain qualification reports, operation-catalog snapshots, source/manifest
identity records, and immutable SIF builds under:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/
```

A SIF becomes eligible for implementation only after its checksum, exact source
commits, installed pack identities, test coverage, exits, skips, and known
limitations are recorded and approved.

Mirror the numbered development workflow and step paths under this result root.
The cytokine work uses the stable split labels
`03a_parse_cytokine_signal_generation/` and
`03b_parse_cytokine_trajectory_analysis/`, followed by the bounded
`03c_parse_cytokine_depth_adjusted_trajectory/` qualification; do not recreate the provisional
`03_parse_cytokine_signatures/`, `06_parse_cytokine_signatures/`, or
`07_trajectory_signal_regions/` roots.
Within each step, use the applicable numbered `01_tables/`, `02_figures/`,
`03_reports/`, `04_logs/`, `05_objects/`, and `06_manifests/` categories.

The checksummed diagnostic run
`01_bulk_marker_boundary/29_simple_mean_contrast_margin_sensitivity/simple_mean_contrast_margin005_20261106T230556Z`
compares the selected raw-UCell winner margin of 0.01 with 0.05 on the same
1,746-cell fixture. The 0.05 rule retained 443 assignments and made 1,303 cells
unassigned. It remains a sensitivity analysis and does not replace the
provisional 0.01 development default.

The checksummed diagnostic run
`01_bulk_marker_boundary/30_simple_mean_contrast_margin_comparison/simple_mean_contrast_margins_20261106T234706Z`
compares pure argmax, margin 0.01, and margin 0.05 using the authenticated raw
UCell scores. They assigned 1,709, 1,416, and 443 cells, respectively; exact
top-score ties remain unassigned. This comparison also leaves the provisional
0.01 default unchanged.

The checksummed run
`03b_parse_cytokine_trajectory_analysis/01_qualify_cytokine_trajectory_methods/bdev_cytokine_trajectory_20261105T153001Z`
is retained as the preliminary cytokine trajectory analysis for comparison with
the later implementation rerun. Its development outputs are not promoted paper
results.

The checksummed run
`03c_parse_cytokine_depth_adjusted_trajectory/01_qualify_depth_adjusted_cytokine_trajectory/bdev_cytokine_depth_adjusted_20261106T145817Z`
is the completed bounded qualification of log1p-CP10K, gene-level
depth-adjusted Parse38 scores, matched dynamics and intervals, and SCT-based
GSEA. It remains a development result and is not promoted to implementation.

The checksummed post-qualification visualization run
`03c_parse_cytokine_depth_adjusted_trajectory/02_plot_up_geneset_gsea_umaps/bdev_cytokine_up_gsea_umaps_20261106T210417Z`
contains one PNG/PDF ATAC-UMAP pair for each of the 21 Parse38 `up` signatures.
Only fitted-L4 fixture cells in supported local-GSEA intervals are colored, by
local NES; all Bdev cells remain as gray context. This run also remains a
development result and did not mutate its parent qualification.

The checksummed derived run
`03c_parse_cytokine_depth_adjusted_trajectory/03_combine_supported_up_geneset_gsea_umaps/bdev_cytokine_supported_up_gsea_grid_20261106T211733Z`
combines the ten `up` signatures having at least one supported local-GSEA
interval into one 2-by-5 PNG/PDF figure with a shared NES scale. It did not
modify either parent run and remains an unpromoted development visualization.

The checksummed run
`04_reference_mapping/02_qualify_bdev_peakvi_pseudocells/bdev_peakvi_pseudocell_10d_20261110T180851Z`
is the completed bounded Step 4 qualification of ten-dimensional PeakVI
projection, authenticated model-bound latent union, support calibration, and
native-L4 transfer. All structural and performance gates passed. It remains a
development result; no real ImmGen library was projected and implementation is
still locked.

The checksummed controlled-development run
`04_reference_mapping/03_map_immgen_atac_trajectory/immgen_atac_l4_20261110T202521Z`
completed real projection of the 28 frozen ImmGen ATAC libraries through the
qualified ten-dimensional Bdev-only PeakVI recipe. Preprocessing passed for all
28 exact-peak-axis profiles. Support calls were 25 supported, 3 borderline,
and 0 out-of-reference, so all libraries retained reportable native-L4
estimates. The 48-file inventory is complete and checksummed. This remains an
unpromoted development result. Scientific acceptance was subsequently reopened
because the approved-ATAC-UMAP diagnostic showed implausible later-stage
placement; the transferred L4 estimates must not be used downstream. No
replacement SIF was built and implementation is still locked.

The checksummed diagnostic visualization run
`04_reference_mapping/05_plot_immgen_atac_umap/immgen_atac_umap_20261111T032028Z`
places all 28 mapped ImmGen libraries on the established `umap.atac` used for
the fitted L4 trajectory. Placement is the kernel-weighted barycenter of each
library's primary 30 authenticated PeakVI neighbors; it is not a refitted UMAP
or a UMAP transform. All population and GSM labels are shown, with 25 supported
circles and 3 borderline orange triangles. The result was not scientifically
accepted because the placements, especially for later-stage populations, are
implausible. Retain it only as a diagnostic: Step 4 remains unresolved, it must
not be used downstream, and implementation remains locked. Source manifests
remained unchanged.

The checksummed diagnostic run
`04_reference_mapping/06_diagnose_immgen_l4_umap_drivers/immgen_l4_umap_drivers_20261115T233506Z`
shows, for all 16 ImmGen populations, the native-L4 distributions and fixed
`umap.atac` coordinates of the primary 30 Bdev PeakVI neighbors that drive each
library. It exactly reconstructed both the weighted-median L4 transfer and the
weighted-barycenter coordinates. The result exposes a major `B.FrE.BM`
replicate split and near-collapse of the five stem/progenitor populations onto
the same high-L4 reference neighborhood. It is diagnostic only; the mapping
remains scientifically unaccepted and must not be used downstream.

The checksummed median-depth diagnostic
`04_reference_mapping/07_test_immgen_median_depth_downsampling/immgen_median_depth_downsample_20261115T235830Z`
projected 30 seeded 18,952-count pseudo-cells for each of 28 ImmGen libraries.
Downsampling did not remove the neighbor hubs: the median full-depth versus
downsample-consensus neighbor Jaccard was 0.968, unique reference neighbors
were 230 versus 228, and median absolute L4 shift was effectively zero. This
rules out raw matrix depth alone as the primary explanation while leaving bulk
mixture composition, reference-space hubness/sample bias, and the missing
post-FrE domain gate unresolved. The run is diagnostic only.

The checksummed position-driver diagnostic
`04_reference_mapping/08_attribute_immgen_position_drivers/immgen_position_drivers_20261116T003136Z`
projected the provisional ImmGen RNA-labelled Bdev fixture through the same
frozen-recipe PeakVI model without adding it to the mapping reference. For the
17 libraries with sufficiently populated exact-label centroids, the observed
neighbor centroid was much closer than the expected-label centroid, and median
expected-label neighbor weight was zero. `latent_3` supplied a median 31% of
the displacement and tracked peak-overlap fraction at Spearman 0.84. Cluster-6
weight, fixed-UMAP coordinates, native L4, neighbor distance, reference sample,
and reference-cell reuse all showed technical associations. The whole-axis
audit also found 48,618 peaks present in every bulk library but accessible in
at most 10% of reference cells. These findings motivate matched peak-masking
and reprojection; they do not accept the provisional labels, candidate peaks,
or the ImmGen mapping for downstream use.

The checksummed label-blind sparsity-adapter diagnostic
`04_reference_mapping/09_test_immgen_label_blind_sparsity_adapter/immgen_sparsity_adapter_20261116T012301Z`
projected 100 rank-weighted sparse pseudo-cells for each of 28 ImmGen libraries.
Construction used only the global reference median of 18,952 fragments, the
global sparsity distribution among reference cells within 10% of that depth,
and within-library bulk peak ranks; it used no labels, reference metadata,
clusters, pseudotime, or expected positions. Sparsity and rank-preservation
targets passed, and the median 50-versus-100 pseudo-cell position difference
was 0.0039. Nevertheless, the adapter caused a stronger common-hub collapse:
only 75 reference cells occurred in all primary-30 neighborhoods, the top hub
occurred for 2,775 of 2,800 pseudo-cells, median cluster-6 neighbor weight was
0.463, and median reference-sample-1 weight was 0.725. Most libraries moved
toward native L4 25.8. Independent rank-sparse sampling is therefore rejected;
the mapping remains unaccepted and implementation remains locked.

The checksummed official-HVP comparison
`04_reference_mapping/14_compare_hvp20k_depth_modes/hvp20k_depth_comparison_20261117T044105Z`
froze depth-adjusted and raw-binary `scwbatac::hvp_select` 20k PeakVI models
before projecting the same 28 ImmGen libraries without refitting. The user
subsequently selected the depth-adjusted HVP model as the provisional
development default; `overall20k` remains its required sensitivity comparator.
The pack-side sparse fix at scwbATAC `3dce524` and authenticated typed-slice
fix at scWorkbench `4aca71d` subsequently passed their focused regressions and
all seven scwbATAC source contracts. Full-Bdev run
`04_reference_mapping/15_requalify_scwbatac_hvp_sparse/scwbatac_hvp_sparse_source_requalification_20261118T004129Z`
then completed both official HVP modes on all 26,723 cells and 191,433 peaks
without the previous dense-allocation warning. Both modes selected exactly
20,000 peaks and reproduced the prior selected membership. Tie-aware
adjudication `scwbatac_hvp_sparse_tie_adjudication_20261118T011000Z` passed;
the source scientific state was unchanged and no operation output was
promoted. No artifact was promoted to implementation, and no replacement SIF
has yet been accepted.

Per the development sequencing decision, no replacement SIF is to be built
until the remaining development analyses are complete. The premature
2026-11-18 build was interrupted during SIF creation and produced no `.sif`
artifact; its preflight record and build log are retained only as an audit of
the aborted attempt. Full FLT3LKO/IL7KO mapping and its fixed UMAP subsequently
completed in the development stage as recorded below.

The checksummed bounded knockout-mapping run
`05_trajectory_group_contrasts/01_qualify_flt3lko_il7ko_mapping_fixture/flt3lko_il7ko_mapping_fixture_20261117T175427Z`
requantified a balanced 400-cell FLT3LKO/IL7KO fixture on the exact frozen
HVP-depth20k axis and projected a balanced 28-cell subfixture through the
frozen model without refitting. All 28 projected cells passed Bdev-calibrated
support; 18 retained finite native-L4 values and ten remained missing rather
than being imputed. This is mapping qualification only. It performed no
genotype contrast, drew no biological conclusion, did not modify archived
inputs, and did not unlock implementation.

The checksummed derived visualization
`05_trajectory_group_contrasts/03_plot_flt3lko_il7ko_fixture_fixed_umap/flt3lko_il7ko_fixture_fixed_umap_20261117T190116Z`
places all 28 balanced projection-fixture cells on the unchanged 26,723-cell
Bdev `umap.atac`. Each point is the distance-kernel-weighted barycenter of its
30 exact frozen-PeakVI neighbors; it is not a UMAP transform. All 28 cells are
shown, including the ten supported mappings without finite native L4. This is
a bounded diagnostic, not a genotype comparison or population estimate.

The completed corrective full-cell run
`05_trajectory_group_contrasts/04_map_flt3lko_il7ko_full/flt3lko_il7ko_full_mapping_fixed_umap_20261118T031500Z`
maps all 23,252 QC-selected FLT3LKO/IL7KO cells through the frozen
HVP-depth20k PeakVI model without refitting and places every cell on the
unchanged 26,723-cell Bdev `umap.atac` using the primary 30-neighbor
distance-kernel barycenter. Support calls were 23,251 supported, one
borderline, and zero out-of-reference; 15,171 cells retained finite accepted
native-L4 values, while missing values remained missing. All 47 bounded
transfer batches were checksummed, both frozen source manifests remained
unchanged, and the PNG/PDF fixed-UMAP panel plus execution receipt were
written. This is descriptive development mapping only: it performs no
genotype contrast, draws no biological conclusion, does not build a SIF, and
does not unlock implementation.

The checksummed gated-transfer qualification
`04_reference_mapping/16_qualify_gated_continuous_transfer/gated_continuous_transfer_20261122T192735Z`
qualified `scwbReferenceMap` 0.2.0 on the fixed depth-adjusted HVP20k
reference. It used 20 latent dimensions, five deterministic folds balanced by
sample, domain, and native-L4 stratum, 30 neighbors, and a 0.51 local-domain
confidence threshold. Calibration accounted for all 26,723 reference cells;
7 of 28 balanced query-fixture cells passed the complete gate and the other 21
were withheld by an `out_domain` vote. The qualified bundle manifest is
`7f6e5fb2968ec86c69513be3ebd382884a23815b787215ad6326fa7fe6eea506`.

The checksummed full gated mapping
`05_trajectory_group_contrasts/05_map_flt3lko_il7ko_gated_20d/flt3lko_il7ko_full_gated_20d_20261122T202618Z`
then projected all 23,252 query cells through that frozen 20D model without
refitting and completed 47 bounded transfer batches. Support calls were
21,885 supported, 1,186 borderline, and 181 out of reference. The complete
gate published accepted native-L4 values for 2,974 cells and withheld values
for 20,278 cells. Every cell has a raw-estimate field and acceptance ledger;
8,161 raw fields are explicitly missing, all on rejected rows. The qualified
source manifest remained unchanged. This remains mapping QC only: no genotype
contrast or paper analysis was performed, no SIF was built, and implementation
remains locked.
