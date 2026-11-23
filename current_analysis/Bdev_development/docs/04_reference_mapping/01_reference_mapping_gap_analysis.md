# 04 Bdev reference-mapping gap analysis

**Status:** proposal only; no mapping executed and no operation-pack or core
source modified

## Required scientific contract

The Bdev reference must be fitted without query observations. Query cells are
then projected through the frozen reference model and mapped only to Bdev
reference cells. Mapping must publish calibrated support, uncertainty, and an
explicit out-of-reference result rather than forcing every query onto the
reference trajectory.

The trajectory source is the authenticated fitted Slingshot `L4` lineage in its
native orientation. The `LSK (0) -> CD19 (1)` convention remains a display-only
transform. Preserved `Bdev@meta.data$pt` values and missing entries remain
unchanged and are not the primary transfer target.

## Existing capabilities that should be reused

The pinned operation packs already provide most of the workflow:

- `scwblatent::peakvi_fit` fits Bdev only and publishes a contained model with
  its reference latent space;
- `scwblatent::peakvi_project` projects disjoint queries through the retained
  model without joint refitting;
- `scwblatent::model_bound_latent_union` joins reference and query coordinates
  without changing them;
- `scwbrepresentation::knn_neighbors` already has authenticated, separate query
  and candidate axes and publishes a query/rank neighbor ledger; and
- `scwbTrajectory` already publishes fitted L4 pseudotime and lineage weights.

No new latent model, trajectory method, scWorkbench core primitive, or PyBench
protocol is indicated.

## Comparison with existing operations

| Existing operation | What it already provides | What it does not provide for Bdev mapping | Decision |
|---|---|---|---|
| `scwblatent::peakvi_fit` | Fits Bdev-only PeakVI and publishes the retained model and reference latent coordinates. | It does not evaluate a query, calibrate support, transfer labels or L4, or identify out-of-reference observations. | Reuse unchanged. |
| `scwblatent::peakvi_project` | Applies the frozen PeakVI model to disjoint query observations without joint refitting. | Projection always yields coordinates when its technical contract succeeds; coordinates alone do not establish that a query is biologically represented by Bdev. | Reuse unchanged. |
| `scwblatent::model_bound_latent_union` | Row-binds authenticated reference and query coordinates under one model identity without changing them. | It intentionally performs no neighbors, calibration, transfer, or acceptance decision. | Reuse unchanged. |
| `scwbrepresentation::knn_neighbors` | Restricts query and candidate axes, returns deterministic neighbor identities/ranks/distances, and therefore has the right generic output contract. | Its current runner computes all query-candidate distances in an R loop. It also intentionally stops at neighbors: it has no held-out calibration, label confidence, L4 transfer, uncertainty, or out-of-reference semantics. | Keep the operation, add only a scalable exact backend, and leave scientific mapping to a separate pack. |
| `scwbtrajectory::slingshot_fit` and `slingshot_lineage_weights` | Publish authenticated fitted native-L4 pseudotime and lineage membership for reference cells. | They fit the reference trajectory; they do not transfer it to projected query cells. | Reuse unchanged as the authoritative reference target. |
| `scwbtrajectory::slingshot_project` | Projects fitted curves into a two-dimensional display space. | It projects curve geometry for visualization, not query cells onto pseudotime, and supplies no mapping support or uncertainty. | Reuse only for display if needed. |
| `scwbtrajectory::trajectory_pseudotime_from_cell_vector` | Wraps an already existing scalar vector without changing values. | It does not estimate query pseudotime or calibrate whether the estimate is supported. | Not a mapping method. |
| `space_path_define`, `space_path_score`, and `trajectory_space_path_segment` | Define a new checkpoint-anchored analytical path, project cells onto its arclength, and optionally filter by distance to that path. | The result is explicitly a distance-derived space-path trajectory, not transfer of the authenticated fitted Slingshot L4. Its distance cutoff is computed within the scored cells rather than from held-out reference samples, and it has no reference-label voting or out-of-reference calibration. | Do not substitute it for L4 transfer; it could be a separate sensitivity analysis only if later approved. |

For the planned single-cell mapping scale, the current brute-force neighbor
runner would evaluate approximately `n_query * n_reference` distances. Using
the archived dimensions as a scale example, 23,252 queries and 26,723 reference
cells require about 621 million 20-dimensional distance evaluations. The issue
is the runner, not the authenticated neighbor contract. An exact indexed
backend therefore avoids inventing a second neighbor artifact while making the
existing generic operation usable.

The proposed `scwbReferenceMap` pack is needed because none of the existing
operations turns authenticated neighbors into a calibrated scientific mapping.
Putting that logic in `scwbLatent` would incorrectly conflate model projection
with biological acceptance; putting it in `scwbRepresentation` would make a
generic neighbor operation responsible for labels and trajectories; and putting
it in `scwbTrajectory` would not cover non-trajectory reference mapping. A
small dedicated pack preserves those boundaries and lets both single-cell and,
if later validated, bulk queries use the same calibration contract.

## Why modifying an existing operation is not sufficient

The two proposed operations are new scientific stages, not additional output
columns from an existing stage. scWorkbench operation manifests bind inputs,
outputs, promotion closure, replay equivalence, and provenance to one operation
identity. Adding calibration and transfer to an existing operation would change
both its meaning and its artifact lifecycle.

### Extending `peakvi_project`

To perform reference mapping, `peakvi_project` would have to gain the reference
space, reference sample groups, annotations, fitted L4 pseudotime, neighbor
parameters, calibration thresholds, and mapping policies as new inputs. It
would also have to emit support calls, labels, uncertainty, and pseudotime in
the same promotion closure as the projected latent space.

That is undesirable for four reasons:

1. PeakVI projection is a model operation; successful projection is not a
   biological acceptance decision.
2. Calibration must be learned from the reference without seeing the query and
   should be reusable across queries. Coupling it to projection would repeat
   calibration and make leakage easier.
3. Mapping should also work with another authenticated shared space, not only
   PeakVI. Embedding it in `peakvi_project` would make the method unnecessarily
   model-specific.
4. A support-policy change should not require rerunning the expensive latent
   projection or changing the identity of the projected coordinates.

`peakvi_project` should therefore remain unchanged and terminate at an
authenticated query space plus projection summary.

### Extending `model_bound_latent_union`

The union operation is intentionally an exact row bind with strict producer
lineage. Adding mapping would make the same operation both preserve coordinates
unchanged and interpret their biological relationship. It would also prevent
creation of a reusable reference calibration before a query union exists.
The union should remain a coordinate and lineage operation only.

### Extending `knn_neighbors`

The neighbor operation could technically accept optional labels, pseudotime,
sample groups, and thresholds. Doing so would make a generic geometric
primitive conditional on reference-mapping semantics and would still not solve
calibration cleanly:

- one query-to-reference ledger cannot estimate held-out reference behavior;
- sample-held-out calibration requires multiple candidate-exclusion folds;
- label voting and pseudotime transfer have different missing-value and
  uncertainty policies;
- support thresholds must be retained as a reference-bound artifact, not
  inferred independently for every query; and
- many valid uses of `knn_neighbors` require no labels, trajectory, or support
  decision at all.

The correct modification to `knn_neighbors` is therefore limited to execution:
make the same neighbor contract scalable. Calibration and transfer consume its
output as separate stages.

### Extending `knn_pseudobulk`

`knn_pseudobulk` aggregates matrix features over neighbor membership using a
mean or sum. Reference mapping requires categorical weighted voting, confidence
and entropy, weighted pseudotime quantiles, cross-`k` sensitivity, calibrated
distance support, and missing-value withholding. Adding those behaviors would
turn a feature-aggregation operation into a multi-purpose mapping operation and
would not provide reference-only cross-validation.

### Extending a Slingshot operation

`slingshot_fit` owns reference trajectory fitting, and `slingshot_project` owns
curve geometry in a display space. Neither operation is given query-to-reference
neighbors. Adding query label and pseudotime transfer would couple Slingshot to
PeakVI/reference-mapping assumptions and would not support mappings that have no
trajectory. The fitted L4 output should instead be a typed input to a mapping
transfer operation.

### Extending the scalar pseudotime adapter

`trajectory_pseudotime_from_cell_vector` promises value preservation: it wraps
an existing scalar without estimating, transforming, or imputing values. Making
it compute neighbor-weighted query values would contradict that explicit
semantic contract and blur observed/reference pseudotime with transferred
pseudotime.

### Extending the analytical space-path operations

The space-path operations already provide path progress, distance-to-path, and
a distance-derived membership rule. They are the closest apparent alternative,
but they answer a different question. They define a new analytical path in the
shared latent space and return its arclength progress. The proposed transfer
instead estimates the authenticated reference L4 value from neighboring Bdev
cells. The two axes need not be numerically or biologically equivalent.

Their current threshold is also derived from the distribution of cells being
scored, so a shifted query can influence its own acceptance threshold. Bdev
mapping requires thresholds frozen from held-out reference samples before the
query is evaluated. Modifying the path operation to accept reference folds,
labels, fitted-model pseudotime, voting, uncertainty, and out-of-reference rules
would effectively recreate the proposed reference-mapping operation under a
trajectory-specific name.

## Why calibration and transfer are separate operations

`reference_map_calibrate` is reference-bound. It performs held-out reference
predictions, measures error, freezes thresholds, and publishes an artifact tied
to the exact reference space and caller-selected reference vector identities.

`reference_map_transfer` is query-bound. It consumes that frozen calibration
and one query-to-reference neighbor ledger, then publishes query-specific
support, labels, uncertainty, and mapped L4. This separation ensures:

- the query cannot influence its acceptance thresholds;
- one accepted calibration can be reused for multiple query datasets;
- alternate transfer policies can be compared without refitting PeakVI;
- calibration failure stops before any query interpretation; and
- provenance distinguishes reference validation from query inference.

## Pack ownership alternatives

The requirement is for two new operations; a new package is not technically
mandatory. There are two reasonable ownership choices:

| Owner | Advantage | Disadvantage |
|---|---|---|
| Add the two operations to `scwbRepresentation` | Fewer packages and direct proximity to the neighbor operation. | Expands a generic representation pack into sample-held-out biological calibration, label transfer, trajectory transfer, uncertainty, and out-of-reference policy. It would also make the pack depend conceptually on artifacts produced by multiple scientific domains. |
| Create `scwbReferenceMap` | Gives calibration artifacts, transfer policies, support categories, and cross-domain inputs one explicit owner; remains reusable with PeakVI or another shared space. | Adds one small official pack and its own manifest/test lifecycle. |

The recommended owner remains a small `scwbReferenceMap` pack because mapping
acceptance is a reusable scientific domain rather than a coordinate-construction
method. If minimizing package count is preferred, the same two operations could
instead be added to `scwbRepresentation`; they should still be new operations,
not modes added to `knn_neighbors` or `peakvi_project`.

## Generic operation-pack contract

The proposed scWorkbench operations must contain no Bdev-specific names,
trajectory keys, state labels, dataset identities, or hard-coded metadata
fields. Bdev is only one qualification fixture.

`reference_map_calibrate` should expose a categorical
`holdout_group_field` parameter. Its default may be `orig.ident`, because that
is a common single-cell sample field, but callers may select any authenticated
categorical metadata field, such as donor, sample, batch, animal, library, or
study. The operation should reject a missing or non-categorical selected field
and record the exact field, category values, and per-category coverage in its
calibration artifact.

The holdout policy is leave-one-category-out: each category is mapped using all
other categories as candidates. A state absent from a particular training fold
is a category/state coverage limitation, not an automatic failure of the full
calibration. The operation reports it explicitly and avoids claiming validated
transfer for that category/state combination.

`reference_map_transfer` should likewise accept optional typed reference
annotations and optional named continuous reference vectors or vector families.
It must not require a field called `cluster`, a trajectory called `L4`, or any
particular biological vocabulary. Output keys are caller-supplied and retain
the input artifact identities. The Bdev workflow will provide native L4 as one
continuous vector family, but the operation itself remains general.

## Gap 1: scalable exact reference-restricted neighbors

The current `scwbrepresentation::knn_neighbors` runner computes every
query-candidate distance in an R loop. Its contract is correct, but it is not
appropriate for tens of thousands of query and reference cells. The archived
exploration therefore used an outside exact `RANN::nn2` calculation.

The proposed update is to retain the existing operation and outputs while
adding a bounded exact `RANN` backend with `eps = 0` and query chunking. The
current brute-force backend should remain available for small-fixture numerical
equivalence and explicit fallback. Candidate restriction, distance definition,
tie/order policy, requested `k`, backend, chunk size, and dependency identity
must be retained in provenance.

This belongs in `scwbRepresentation`, not core.

## Gap 2: calibrated mapping and transfer

Neighbor finding alone does not define a scientifically acceptable mapping.
A new official `scwbReferenceMap` operation pack is proposed with two
operations:

### `scwbreferencemap::reference_map_calibrate`

This operation would perform category-held-out calibration using the
authenticated reference space. The categorical holdout field is configurable;
`orig.ident` is the default, not a hard-coded requirement. Each observed
category is held out in turn and mapped only to the other categories. It would
publish:

- held-out reference predictions and neighbor-distance diagnostics;
- categorical accuracy, balanced accuracy, confidence, margin, and entropy;
- continuous-vector transfer error and rank correlation where a selected
  continuous reference input is defined;
- the calibrated distance distribution and support thresholds; and
- a retained calibration artifact tied to the reference space, grouping field,
  neighbor counts, distance statistic, and selected input artifact identities.

For the Bdev qualification fixture only, the selected holdout field is
`orig.ident`. The preliminary method defaults are `k = 15, 30, 50`, primary
`k = 30`, and median first-30 distance. Support is `supported` at or below the
held-out-reference 95th percentile, `borderline` between the 95th and 99th
percentiles, and `out_of_reference` above the 99th percentile. These remain
proposal values until approved.

### `scwbreferencemap::reference_map_transfer`

This operation would consume the authenticated query-to-reference neighbor
ledger, retained calibration, optional reference annotations, and optional
named continuous reference vectors or vector families. It would publish, for
every query:

- support category and calibrated distance statistics;
- predicted label plus confidence, margin, entropy, and cross-`k` agreement;
- distance-weighted median values for each requested continuous input;
- weighted quartiles, IQR, and cross-`k` sensitivity;
- authenticated query vector families with caller-supplied keys; and
- a machine-readable mapping summary by sample and query group.

The provisional weight is `exp(-(distance / kth_distance)^2)`. Raw estimates
are retained for audit, but primary transferred values are missing for
out-of-reference queries. Missing reference values are excluded only for the
corresponding transfer calculation; they are never imputed or changed. In the
Bdev qualification configuration, native L4 is supplied as one such named
continuous input.

## Bulk ATAC boundary

No bulk-specific PeakVI operation is proposed yet. Projecting a bulk library
through a model trained on single cells is technically possible only after the
feature and count contracts align, but that does not establish biological
validity. Bdev pseudobulks and depth-matched pseudo-cells must first demonstrate
acceptable behavior. A bulk-specific operation should be considered only if
that proof succeeds and exposes a reusable contract not already handled by
`peakvi_project` plus `scwbReferenceMap`.

## Qualification required before real mapping

Before any full Bdev query mapping, the proposed generic operations should pass:

1. exact-neighbor equivalence against brute force on a small fixture;
2. deterministic query/candidate ordering and disjointness tests;
3. configurable categorical holdout-field tests, including a non-`orig.ident`
   field;
4. synthetic supported, borderline, and displaced out-of-reference queries;
5. missing-continuous-value preservation and out-of-reference withholding
   tests;
6. a bounded, `orig.ident`-held-out Bdev calibration and transfer fixture;
7. identifier-alignment and source-bundle immutability checks; and
8. operation-pack manifest, dependency, replay, promotion-closure, and negative
   tests.

The full single-cell knockout dataset, bulk libraries, and paper analyses must
not be run during this qualification step.

## Proposed decision

Approve only these two development changes:

1. a scalable exact backend for the existing
   `scwbrepresentation::knn_neighbors`; and
2. a new official `scwbReferenceMap` pack with calibration and transfer
   operations.

No scWorkbench core, PyBench, `scwbLatent`, or `scwbTrajectory` edit is proposed
for this capability at present. Source edits should begin only after this
boundary and the provisional calibration defaults are approved.
