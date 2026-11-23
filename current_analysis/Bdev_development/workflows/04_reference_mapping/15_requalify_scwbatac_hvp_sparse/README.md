# Sparse-safe scwbATAC HVP requalification

This development qualification installs pinned scWorkbench commit `4aca71d`
and scwbATAC commit `3dce524` over the immutable baseline SIF, executes both
official `hvp_select` modes on all 26,723 Bdev reference cells and 191,433
peaks, and compares their complete ledgers and selected 20k axes with the prior
benchmark.

The source bundle is opened with the writable handle required for coherent
durable-session snapshots, but no outputs are promoted to it. Operation state
is retained only in run-specific session directories. Acceptance requires an
unchanged scientific-state fingerprint and permits only the framework's
`publication.published_at` refresh on writable close.

Attempts `20261117T193110Z` and `20261117T200929Z` stopped before selector
execution because a physical 13 GB copy exhausted paired storage and a
read-only source cannot create durable sessions, respectively. Attempt
`20261117T201433Z` reached the selector and was stopped after two 2.9 GiB
sparse-to-dense allocation warnings. A direct diagnostic of one exact
191,433-feature by 2,048-cell chunk showed that scwbATAC slicing,
binarization, row/column summaries, and sparse multiplication stayed sparse.
The remaining conversion was localized to scWorkbench's authenticated
typed-slice adapter.

scWorkbench commit `4aca71d` removed that upstream conversion and added a
regression that traps sparse typed-slice coercion. The regression and all seven
scwbATAC contracts passed against the fixed core.

Attempt `20261118T003258Z` correctly failed before operation execution because
a prior durable-session snapshot cannot be nested as a new durable-session
source. Fresh source-snapshot run `20261118T004129Z` then completed
depth-adjusted and raw HVP in 360 and 320 seconds without dense-allocation
warnings. Both selected exactly 20,000 peaks. The depth-adjusted result
reproduced exact ranks and all numerical columns within `1e-12`. Raw-mode
score-driving columns and selection reproduced; its 94 order/rank changes were
only tied-score permutations with maximum displacement three and no selection
change. Its unused depth diagnostics correctly changed to the tested
not-requested sentinels. Adjudication `20261118T011000Z` passed and the result
inventory is checksummed. Writable close again changed only
`publication.published_at`; scientific-state fingerprint `4d53e5f7...`
remained unchanged.

This workflow does not read or map the FLT3LKO/IL7KO query cells, create paper
results, or unlock implementation. A replacement SIF still must be built,
container-qualified, checksummed, and explicitly accepted.

For a new execution run, run in order:

1. `01_prepare_runtime_and_clone.py --config <execution-config>`
2. `02_requalify_hvp.R --config <execution-config>` inside the configured SIF and R library
3. `05_adjudicate_core_fixed_requalification.py --config <adjudication-config>`
4. `06_finalize_adjudicated_requalification.py --config <adjudication-config>`
