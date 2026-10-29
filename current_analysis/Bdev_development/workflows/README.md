# Development workflows

Store reproducible definitions for:

- SDK scaffold and static/installed/contract checks;
- package-level checks and scientific oracle tests;
- Bdev acceptance fixtures that do not modify authoritative source bundles;
- coordinated core, PyBench, and pack integration checks;
- replacement-SIF construction; and
- clean-process activation, catalog inspection, relocation, replay, and
  pack-free recipient validation.

Development workflows may use the baseline SIF as their toolchain, but final
qualification must execute inside the newly built candidate SIF.

These workflows develop and qualify framework capabilities; they do not run the
paper analysis. Use only synthetic fixtures, SDK fixtures, and explicitly
declared bounded Bdev acceptance fixtures. Paper-scale inputs, source-bundle
construction, terminal figures, and biological conclusions belong to the locked
implementation stage.

Use permanent numbered capability directories, numbered qualification steps,
and numbered scripts. Paired-storage build and test output must mirror those
identifiers and use numbered output categories as defined in
`../../ORGANIZATION.md`. Do not renumber referenced development evidence.
