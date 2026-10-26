# Parse signature bundle recovery

This directory records the project-specific repair for the large registry JSON
failure encountered after the 2026-09-21 UCell score promotion. The original
image and original Parse bundle are preserved. Recovery uses a separate bundle
copy and an isolated scWorkbench R package built from the source inside the
image. General scWorkbench interfaces remain documented inside that image.

The local package version is `0.1.0.9004`. `registry_json.patch` changes only the
registry-table hashing call site; `src/scworkbench_registry_json.R` supplies its
bounded JSON serializer inside the engine package. The serializer retains the
existing JSON bytes, 5,000-row hash boundaries, and identity scheme. It does not
remove, shorten, replace, or migrate the original provenance or score payload.
Large escaped JSON is written to a temporary file before hashing. The database
value itself is still read into memory, so this is a serialization-size fix,
not a claim of constant-memory bundle inspection.

`object_row_keys.patch` adds an equivalent vectorized path for plain character
row keys, retaining the original serializer for other types and missing values.
This makes verification of the existing multi-million-row score families
practical. `feature_set_readonly.patch` preserves the logical read-only flag
when a staged feature-set import is promoted. `artifact_listing.patch` orders
artifact row indices in R after fetching them, avoiding a SQLite temporary sort
of the large provenance record while retaining the public getter's ordering.
All changes are internal engine fixes; operation interfaces and scientific
algorithms are unchanged.

Run the phases inside the repository, using a new UTC run identifier:

```bash
bash processing/parse_cytokines/run_10_recover_parse_bundle.sh signature_recovery_YYYYMMDDTHHMMSSZ prepare
bash processing/parse_cytokines/run_10_recover_parse_bundle.sh signature_recovery_YYYYMMDDTHHMMSSZ recover
```

The runner uses the supplied SIF, a 96-GiB virtual-memory ceiling, and at most
four BLAS/OpenMP threads for recovery. It creates these separate locations
under the configured Parse storage root:

- `work/<run_id>/runtime/`: extracted source, installed R library, and build identity.
- `reports/<run_id>/runtime_tests/`: runtime tests and semantic audit.
- `reports/<run_id>/finalization/`: file hashes, operation evidence, and validation.
- `bundle_recovered_<timestamp>/`: the candidate recovered bundle.

Preparation verifies the exact inspected source-file checksum before patching.
Tests compare streamed and original jsonlite hashes, including escaping,
Unicode, missing values, large strings, hash boundaries, and tamper detection.
The image's manifest publication, object round-trip, deep integrity, and
feature-set import contract tests run against disposable fixtures, followed by its semantic
catalog audit. CSV registration also has a disposable-bundle test checking
exact gene order, read-only metadata, and full recipient validation. Artifact
listing is compared against the original SQL query on disposable fixtures with
ties, missing dates, filters, and lifecycle states.

An already tested runtime can be reused by passing its directory as the third
runner argument and its test-report directory as the fourth argument to
`recover`. If registration completed but validation stopped, `validate` resumes
recipient validation on that same copy using a tested runtime. It rechecks
signature memberships, read-only flags, metadata, cell sets, unaffected registry
identities, and public artifact listing before validation. It preserves the
initial recovery log and refuses to overwrite validation retry evidence. A
fifth argument supplies a distinct `validation_...` label for an additional
validation attempt; the prior validation JSON is archived under that label.
`restore` runs just the authenticated materialization-restoration step on an
already registered copy.

Recovery verifies every file in the initial copy against the original bundle,
reconstructs registry identity, opens the copy through scWorkbench, and imports
the unchanged 23 gene sets from
`reports/signature_report_20260921T221910Z/final_accepted_feature_sets.tsv`.
It writes and verifies an equivalent CSV for scWorkbench's documented table
file parser. A bounded ephemeral session is sufficient for this cheap import;
the staged import must contain all 23 sets before its evidence is promoted.
The script checks all cell sets, cell metadata, unaffected registry-table
identities, and exact gene memberships. Full recipient validation must pass,
including content verification of
`cell_vector_family_20260921_155957_5n7edbnw`. Before validation, recovery restores
the five missing materialization files declared by the interrupted 2026-09-21
commit. The retained checkpoint is under
`sessions/parse_signature_ucell_20260921T152017/`. The restoration script
authenticates it through its pointer and metadata hashes, writes the saved
primitives using the engine's original RDS format, and requires both original
materialization identities and original raw-file SHA-256 values to match before
placing any file in the copy. It uses the normal scWorkbench writer lease and
leaves registry identity and the saved checkpoint unchanged. No scoring or
statistical analysis is rerun. After validation, the runner verifies that every
original bundle file still has its initial SHA-256.

Within `finalization/`, `completion.json` records engine validation,
`materialization_restoration.json` records the restored files and identities,
and `source_preservation.json` records the final original-file check. All are
required before accepting the copy in `config/bundles.yml` and
`config/parse_cytokines.yml`.
The registry must also record the required `runtime.r_library` and
`runtime.scworkbench_version`. Parse processing/report runners resolve these
through `src/parse_runtime.sh`; R registry validation loads the same runtime
before opening the bundle. The unmodified image still has the original limit.

The existing stage-09 prevention change remains in place: future scoring runs
keep per-cell score families in durable sessions and promote compact summaries.
That change is independent of recovery of the already committed score family.
