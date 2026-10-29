# Repository separation

The framework source and Bdev analysis are separate Git repositories with
different responsibilities:

```text
/home/jtg/msi_projects/
├── scworkbench/       # scWorkbench core source and its Git history
└── scworkbench_Bdev/  # Bdev analysis coordination, configuration, and workflows
```

The Bdev repository does not vendor the core, use it as a submodule, or store
core branches. It records the exact source URL, local checkout, branch, and
commit in `../config/source_repositories.yaml`. Generated libraries, test
bundles, image build contexts, and qualification results belong in paired
storage rather than either Git repository.

## Current core checkout

The authoritative standalone checkout is:

```text
path: /home/jtg/msi_projects/scworkbench
remote: https://gitlab.liu.se/hende58/scworkbench.git
branch: master
commit: d85353c3f8acfd181cc9d179a98f2c6339cc020d
```

The corrected development SIF contains core commit
`9f10e7c0a849a87249573a2ef78d8ff7643dccba`, which is not available from the
public remote. The SIF remains the immutable toolchain and image base; this
source mismatch is provenance, not a requirement to reconstruct the SIF.

No development branch has been created yet because its correct base must be
chosen explicitly. Creating a branch does not authorize core changes. Before
any scWorkbench edit, document the pressing core gap, explain why an operation
pack cannot solve it, motivate the smallest change and its qualification impact,
and obtain explicit user approval for that specific edit.

## Development flow

1. Define the Bdev capability and first route it to an existing or new operation
   pack.
2. Use the corrected SIF to run SDK checks and qualification against changed
   source.
3. If a missing core primitive is demonstrated, stop and request explicit
   approval before branching or editing scWorkbench.
4. Record every changed source commit in the Bdev target runtime manifest.
5. Extend the corrected SIF with the accepted updates, qualify the derived SIF,
   and pin it for implementation.
