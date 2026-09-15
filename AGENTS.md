# Repository guidance

These instructions apply to the entire repository.

## Project

This repository develops the B cell development use case for scWorkbench.
The analyses described in `README.md` cover cytokine stimulation Parse data
and IL7KO/FLT3KO single-cell multiome data, including mapping to a reference
B cell development dataset.

For the Bdev reference, use the existing `Bdev@meta.data$pt` field as
pseudotime, as selected by the user. Preserve its values and missing entries
when importing and using the reference.

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

- `config/`: the central bundle registry and one project configuration per
  source or derived bundle.
- `processing/<bundle_id>/`: numbered dataset pipeline scripts, with numeric
  prefixes indicating execution order within that workflow.
- `processing/integration/`: the Bdev-centered integration pipeline. It reads
  accepted parent bundles and writes only to the separate `bdev_integrated`
  bundle.
- `figure_scripts/<bundle_id>/`: final figure and report generation scripts.
- `src/`: reusable code shared by pipeline and figure scripts.

Keep `bdev_reference`, `parse_cytokines`, and `scmulti_flt3ko_il7ko`
independently reproducible. Do not write integration results into these parent
bundles. Record joint state in `bdev_integrated`, with `bdev_reference` as the
explicit reference parent.

## Repository practices

- Keep container images, large datasets, bundles, models, and generated outputs
  out of Git; check `.gitignore` when adding output locations.
- Preserve sample metadata and record input provenance, analysis settings,
  and output locations needed to reproduce results.
- Keep changes focused and report what changed and how it was validated.
