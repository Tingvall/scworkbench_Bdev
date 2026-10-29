# Development inputs

Declare the immutable baseline SIF, source repositories, SDK documentation,
test fixtures, and checksummed biological inputs required for scientific
acceptance tests. Large files remain in paired storage.

The currently inspected baseline is:

```text
path: /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif
sha256: 2c9aa71f69d5775166f7dd0625176a8795d00caa498500eae7c39c42c2a7695a
status: candidate development baseline; not an implementation runtime pin
scworkbench_version: 0.1.0
scworkbench_commit: 9f10e7c0a849a87249573a2ef78d8ff7643dccba
pybench_commit: 608066e77583a37d9187a7bb5c8a5c49c082b17d
architecture: amd64
```

This corrected image supersedes the immutable `20261022` image for active
development because the earlier image did not provide the required CPU PeakVI
runtime. It is the baseline for the scWorkbench update, not the qualified image
for the paper implementation. The superseded image remains at its declared
storage path and is not modified.

The corrected image's source commits and installed pack identities are recorded
in `../config/baseline_runtime.yaml`. A full CPU PeakVI lifecycle and the
operation-pack tutorial still require qualification for this exact image before
their status can be recorded as passed.

The authoritative scWorkbench core checkout is a sibling of this analysis
repository:

```text
/home/jtg/msi_projects/scworkbench
```

The remaining core and operation-pack checkouts currently remain under paired
development storage at
`/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/source_repositories/`.
Their explicit local paths, public remotes, branches, pull-time commits, and
relationship to the corrected SIF are pinned in
`../config/source_repositories.yaml`.

The corrected SIF's scWorkbench and `scwbLatent` commits are not present in the
currently published GitLab checkouts. This is recorded as baseline provenance
and does not block development because the corrected SIF is the toolchain and
immutable base, not an image that must be reconstructed. Only changed components
require pinned, available source commits for installation into the new derived
SIF.
