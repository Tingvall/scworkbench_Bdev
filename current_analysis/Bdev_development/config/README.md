# Development configuration

Record both the baseline environment and the proposed replacement runtime.
Configuration should include:

- baseline SIF path and SHA-256;
- baseline scWorkbench, PyBench, and official pack commits;
- source repository URL, feature branch, base commit, and target commit for
  every changed component;
- operation manifest identities and compatibility ranges;
- dependency locks, architecture, seeds, and compute profiles;
- SIF build definition and build inputs; and
- required qualification suites and acceptance thresholds.

Do not treat a mutable branch name as a runtime pin. Record the resolved commit
for every candidate build.

The current source record is [`baseline_runtime.yaml`](baseline_runtime.yaml).
The external Git checkouts and their exact pull-time commits are recorded in
[`source_repositories.yaml`](source_repositories.yaml). Runtime commits and
checkout commits are separate identities and must not be substituted for one
another when they differ. Availability of every embedded baseline commit in the
public checkouts is not required: the corrected SIF itself is the immutable
toolchain and build base. Exact Git commits are required for components changed
and installed into the derived candidate SIF.

Populate [`target_runtime.yaml`](target_runtime.yaml) as branches are selected,
resolved, built, and qualified.

[`capability_registry.yaml`](capability_registry.yaml) assigns permanent numbers
to proposed Bdev framework gaps and records their likely owning pack or approved
outside boundary. A proposed entry is not authorization to implement it, and no
core-owned entry may proceed without the separate explicit approval gate.
