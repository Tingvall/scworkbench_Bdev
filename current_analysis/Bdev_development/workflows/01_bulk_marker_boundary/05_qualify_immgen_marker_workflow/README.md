# Qualify the ImmGen marker workflow core

This step uses synthetic integer counts with the exact frozen 126-library
cohort structure. It exercises population aggregation, four design matrices,
every declared contrast, edgeR filtering and quasi-likelihood testing, complete
statistics, annotation gates, deterministic ranking, family size rules, and
no-padding behavior.

It never reads ImmGen expression values. Stability resampling, confounding
audits, and held-out validation remain explicit gates for the next
qualification step. Therefore this workflow cannot emit `accepted_core` sets
or an importable GMT.

Generated output belongs under the paired storage path declared in
`05_immgen_marker_workflow_qualification.yaml`.
