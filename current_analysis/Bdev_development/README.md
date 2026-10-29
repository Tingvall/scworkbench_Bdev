# Bdev scWorkbench development

**Purpose:** extend and qualify scWorkbench for the Bdev paper use case, then
produce a new immutable SIF for the implementation stage.

This is the current active stage. Work here identifies and closes reusable
capability gaps, qualifies every changed component, and produces the replacement
SIF. It does not execute the Bdev paper analysis. Biological data may be used
only as declared, bounded acceptance fixtures until implementation is unlocked.

The current candidate SIF is a development baseline. It supplies the installed
toolchain, dependencies, operation-pack SDK, and immutable base for the next
image; it must never be edited in place. The baseline SIF does not need to be
reconstructed from its embedded source commits. Source changes are made on
short-lived branches in the repository that owns each changed component, tested
against the baseline SIF, and then installed into a newly built derived SIF:

- scientific operations in the existing official pack that best fits their
  domain, or in a new official pack when no existing pack fits;
- generic engine capabilities in the scWorkbench core repository only when a
  pressing missing primitive or other infrastructure gap requires them;
- PyBench protocol/runtime changes only when a pressing core protocol gap
  requires them; and
- Bdev development coordination, fixtures, and acceptance workflows here.

`scworkbench` and PyBench are core and are frozen by default during operation-pack
development. A proposed core edit must document why an operation pack cannot
provide the capability. Before creating a core development branch or editing
either checkout, present the evidence, motivate the smallest proposed change,
state its qualification consequences, and obtain the user's explicit approval
for that specific change. Approval is never inferred from the broader analysis
or operation-pack work.

The scWorkbench core source is a separate sibling Git repository at
`/home/jtg/msi_projects/scworkbench`; this Bdev repository contains analysis
coordination and exact source pins only. See
[`docs/repository_layout.md`](docs/repository_layout.md).

Project-local operation packs are prototypes only. A method intended for the
paper analysis should be promoted into an official thematic pack, or into a new
official pack when no existing domain is suitable, before the implementation
runtime is frozen.

## Layout

| Directory | Purpose |
|---|---|
| `config/` | Baseline and target source pins, branch manifest, build configuration, and qualification policy. |
| `docs/` | Capability gap analysis, operation contracts, and development decisions. |
| `inputs/` | Declarations for the baseline SIF, externally stored source repositories, fixtures, and immutable reference inputs. |
| `modules/` | Temporary Bdev-owned pack prototypes or integration fixtures; do not vendor scWorkbench core or official pack repositories here. |
| `workflows/` | Reproducible SDK checks, package checks, integration qualification, and SIF build definitions. |
| `intermediate_objects/` | Repository contract for temporary libraries, build state, test bundles, and logs stored outside Git. |
| `results/` | Qualification reports and immutable SIF artifacts stored outside Git. |

Number current work by the permanent capability IDs in
[`config/capability_registry.yaml`](config/capability_registry.yaml), followed
by numbered qualification steps and scripts. Do not create active implementation
dataset workflows during this stage.

## Promotion gate

Development is complete only when every changed component commit and every
inherited baseline identity is immutable; all required checks are accounted
for; the new SIF checksum and operation catalog are recorded; and the
implementation stage has an accepted runtime manifest. The new SIF is an
artifact of the pinned baseline SIF plus accepted source changes, not a
reconstruction of the baseline image and not a substitute for changed
components' Git history.
