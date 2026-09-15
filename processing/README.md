# Processing workflows

Each dataset has its own numbered pipeline. Numeric prefixes define execution
order within that dataset directory. The `integration` pipeline begins only
after all three parent bundles have an accepted registry state.

Run scientific steps in the supplied scWorkbench image and use the project
configuration under `config/` for stable dataset identities and locations.
Write run-specific inspection and validation evidence outside Git under
`intermediate_objects/reports/<bundle_id>/<run_id>/`.

The integration pipeline must open all three parent bundles read-only and write
only to `bdev_integrated`.
