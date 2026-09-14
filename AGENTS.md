# Repository guidance

These instructions apply to the entire repository.

## Project

This repository develops the B cell development use case for scWorkbench.
The analyses described in `README.md` cover cytokine stimulation Parse data
and IL7KO/FLT3KO single-cell multiome data, including mapping to a reference
B cell development dataset.

## Default analysis: scWorkbench

- Use scWorkbench as the default for analysis and create a scWorkbench bundle
  for each dataset, following the documentation provided inside the `.sif`.
- The image is in the mirrored storage location:
  `/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif`.
- All scWorkbench documentation is inside that image. Consult it for commands,
  supported workflows, and bundle creation. Do not guess tool behavior or
  duplicate its documentation in this repository.
- Keep dataset identities and bundle locations explicit so each bundle can be
  traced to its inputs and analysis settings.
- Editing repository instructions alone does not require running the image.

## Analysis outside scWorkbench

Use analysis outside scWorkbench only when the user explicitly requests it,
or when scWorkbench cannot perform the task and the user approves an outside
approach. Do not switch automatically because another tool is more familiar
or convenient.

If scWorkbench cannot perform a requested task, explain the specific limitation
and the proposed alternative, then obtain approval before implementing or
running that outside analysis. Existing explicit authorization for the same
scope is sufficient; do not ask again.

For authorized outside analysis:

- Keep work within the requested or approved scope.
- Use reproducible scripts with recorded dependencies, inputs, parameters, and
  random seeds where applicable.
- Preserve raw data and existing scWorkbench bundles. Save outside results
  separately and record their source dataset or bundle.
- Verify cell/sample identifiers, metadata alignment, and relevant output
  checks. State assumptions and any validation that could not be performed.

## Repository structure

- `processing/`: numbered pipeline scripts, with numeric prefixes indicating
  execution order (for example, `01_inspect_reference.R`).
- `figure_scripts/`: final figure generation scripts.
- `src/`: reusable code shared by pipeline and figure scripts.

## Repository practices

- Keep container images, large datasets, bundles, models, and generated outputs
  out of Git; check `.gitignore` when adding output locations.
- Preserve sample metadata and record input provenance, analysis settings,
  and output locations needed to reproduce results.
- Keep changes focused and report what changed and how it was validated.
