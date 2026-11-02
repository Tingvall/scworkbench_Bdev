# Qualify deferred ImmGen marker gates

This step uses only the checksummed synthetic outputs from step 05. It
qualifies balanced stability resampling, target-population and target-library
holdouts, confounding audits, broad held-out scoring, and local-neighbor AUROC
validation.

It never reads ImmGen expression values. It cannot create `accepted_core` sets
or an importable GMT. Generated output belongs under the paired storage root
declared in `06_immgen_deferred_gate_qualification.yaml`.
