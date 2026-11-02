# Nested qualification of positive competitor refinement

This workflow tests the positive-only set-level procedure with target-unit
holdouts. Candidate ranking and greedy refinement are rebuilt using training
units only. Targets with fewer than three independent units are inventoried but
not evaluated.

Execution uses the frozen development SIF and run-owned edgeR library for the
fold models, followed by an outside-scWorkbench Python refinement/evaluation
step. Outputs remain diagnostic and are never import artifacts.

