# Gated continuous-transfer qualification

The previous reference-map implementation accepted continuous trajectory
values whenever general latent-space support passed. That could retain a
native-L4 estimate for a query whose local neighbors did not predominantly
belong to fitted lineage L4.

`scwbReferenceMap` 0.2.0 adds a reusable `support_and_annotation` policy. This
development qualification tests that policy with a 20-dimensional HVP20k
PeakVI reference, five sample-balanced folds, complete `in_domain` versus
`out_domain` reference labels, 30 neighbors, and a 0.51 minimum vote
confidence. Whole-sample holdout is deliberately not used because Bdev sample
and trajectory position are strongly confounded; balanced folds retain each
sample and represented trajectory stratum in calibration training and test
sets.

The gate does not delete uncertain estimates. It retains the raw estimate and
records why an accepted value was withheld. The bounded 28-cell fixture must
pass before a newly numbered full KO mapping is allowed.

## Result

Run `gated_continuous_transfer_20261122T192735Z` passed. It used 20,000
depth-adjusted HVP features, 20 PeakVI dimensions, five balanced folds, and
30 neighbors. The reference contained 10,859 `in_domain` and 15,864
`out_domain` cells. Calibration evaluated all 26,723 cells, with domain-label
accuracy 0.821 and native-L4 median absolute error 0.568. The calibrated
supported and out-of-reference thresholds were 2.026 and 2.568.

Seven of 28 balanced KO fixture cells passed the complete gate. The other 21
were withheld because the local annotation vote was `out_domain`; none failed
general support or the 0.51 confidence threshold. Every query retained an
acceptance-ledger row and raw estimate field, including explicit missing raw
values where no finite local L4 estimate existed. The parent bundle manifest
remained exactly `1b1935f7...`, and no genotype contrast or full-query mapping
was performed.

This result qualifies the generic gate on the bounded fixture. The next step
is a newly numbered full KO mapping run using the same frozen 20D model,
calibration artifact, and acceptance policy; it must not overwrite the prior
10D mapping.
