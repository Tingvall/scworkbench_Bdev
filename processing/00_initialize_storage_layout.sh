#!/usr/bin/env bash
set -euo pipefail

# Create dataset-centered mirrored-storage directories. Canonical objects live
# directly under intermediate_objects/<dataset_id>/; no index layer or
# large-file symlinks are used.

STORAGE_ROOT="${1:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev}"
INTERMEDIATE_ROOT="${STORAGE_ROOT}/intermediate_objects"

if [[ ! -d "${STORAGE_ROOT}/datasets" || ! -d "${INTERMEDIATE_ROOT}" ]]; then
  echo "Expected datasets/ and intermediate_objects/ under ${STORAGE_ROOT}" >&2
  exit 1
fi

dataset_ids=(
  bdev_reference
  parse_cytokines
  scmulti_flt3ko_il7ko
  bdev_integrated
)

for dataset_id in "${dataset_ids[@]}"; do
  mkdir -p -- \
    "${INTERMEDIATE_ROOT}/${dataset_id}/prepared" \
    "${INTERMEDIATE_ROOT}/${dataset_id}/reports" \
    "${INTERMEDIATE_ROOT}/${dataset_id}/work"
done

cat > "${INTERMEDIATE_ROOT}/README.md" <<'EOF'
# Intermediate objects

Each dataset has one directory directly under `intermediate_objects/`:

| Dataset | Role | Current key object |
|---|---|---|
| `bdev_reference` | Reference | `bdev_reference/bundle/` |
| `parse_cytokines` | Parse cytokine analysis | `parse_cytokines/prepared/parse_cytokines_unfiltered_seurat.rds` |
| `scmulti_flt3ko_il7ko` | KO multiome source analysis | Prepared object pending |
| `bdev_integrated` | Bdev-centered integration | Bundle pending |

Within each dataset directory:

- `prepared/` contains outside-bundle source objects and validation evidence.
- `bundle/` is the canonical scWorkbench bundle and exists only after creation.
- `reports/` contains generated reports and run-specific validation evidence.
- `work/` contains temporary, failed, or resumable working state.

Run `processing/00_initialize_storage_layout.sh` from the repository to refresh
the directory skeleton after new datasets are added.
EOF

for dataset_id in "${dataset_ids[@]}"; do
  cat > "${INTERMEDIATE_ROOT}/${dataset_id}/README.md" <<EOF
# ${dataset_id}

- Prepared source objects: \`prepared/\`
- scWorkbench bundle: \`bundle/\`
- Reports: \`reports/\`
- Working files: \`work/\`

Raw delivered inputs remain under \`../../datasets/\`. See \`../README.md\` and
the repository's \`config/${dataset_id}.yml\` for status and provenance.
EOF
done

echo "Dataset directories initialized at ${INTERMEDIATE_ROOT}"
