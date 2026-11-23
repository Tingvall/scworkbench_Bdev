# Full-Bdev PeakVI feature-axis qualification

This development workflow selects PeakVI features from Bdev alone. Its primary
20,000-peak axis contains exactly 10,000 peaks selected by depth-adjusted,
cluster-discriminating deviance and 10,000 additional peaks selected by
label-free depth-adjusted residual deviance. Sample agreement is not required:
sample-restricted Bdev states are retained as biology.

Candidate models are fit on a cluster-stratified Bdev reference split and
project disjoint held-out Bdev cells and Bdev-derived pseudobulks. The primary
hybrid model is compared with cluster-only, overall-only, and complete-peak
controls. L4 eligibility and native-L4 transfer are calibrated only after
whole-Bdev mapping, using the authenticated fitted Slingshot outputs.

If the predeclared hybrid gate passes, a fresh PeakVI model is fit on all
26,723 Bdev cells and frozen with its ordered feature list, model, latent space,
mapping thresholds, source identities, and checksums. ImmGen and every other
external query are prohibited throughout selection, fitting, and calibration.
They may be applied only in a later workflow after this reference is frozen.

Run the numbered scripts with
`config/04_reference_mapping/16_full_bdev_peakvi_feature_qualification.yaml`.
Generated matrices, models, logs, and reports are written only to paired
`Bdev_development` storage.

## Observed outcome

Run `full_bdev_peakvi_hybrid20k_20261116T174018Z` completed all four Bdev-only
candidate fits. The 50/50 hybrid failed its predeclared macro-accuracy and
depth-correlation gates, so no full-reference model was frozen. Cluster-only
20k was strongest for held-out single-cell cluster recovery, depth robustness,
L4 membership and L4 transfer; the complete peak axis was strongest for the
104 Bdev pseudobulks. The retained evidence supports testing a cluster-majority
hybrid in a new run rather than weakening this run's gate or retroactively
selecting a control.
