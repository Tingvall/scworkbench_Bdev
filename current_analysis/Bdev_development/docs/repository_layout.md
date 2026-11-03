# Repository separation

The framework source and Bdev analysis are separate Git repositories with
different responsibilities:

```text
/home/jtg/msi_projects/
├── scworkbench_code/  # scWorkbench core, pybench, and operation-pack repositories
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
path: /home/jtg/msi_projects/scworkbench_code/scworkbench
remote: https://github.com/Tingvall/scworkbench.git
upstream: https://gitlab.liu.se/hende58/scworkbench.git
branch: development
commit: da9477d7ebe3010cedcaff6867e1e6eebbd617a3
```

The corrected development SIF contains core commit
`9f10e7c0a849a87249573a2ef78d8ff7643dccba`, which is not available from the
public remote. The SIF remains the immutable toolchain and image base; this
source mismatch is provenance, not a requirement to reconstruct the SIF.

The development branch is maintained in the user-owned GitHub repository.
The original GitLab repository is retained only as a read-only upstream. Before
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
