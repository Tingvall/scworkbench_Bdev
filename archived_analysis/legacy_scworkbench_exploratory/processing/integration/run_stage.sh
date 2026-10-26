#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
stage="${1:?Usage: run_stage.sh STAGE RUN_ID}"
run_id="${2:?A stable run ID is required}"
[[ "${run_id}" =~ ^[A-Za-z0-9_]+$ ]] || exit 2
case "${stage}" in
  preflight) script="processing/integration/01_preflight_inputs.R" ;;
  create) script="processing/integration/02_create_bdev_integrated.R" ;;
  score) script="processing/integration/03_score_parse_signatures.R" ;;
  dynamics) script="processing/integration/04_cytokine_pseudotime.R" ;;
  validate) script="processing/integration/05_validate_cytokine_bundle.R" ;;
  figures) script="figure_scripts/bdev_integrated/03_sliding_window_cytokine_report.R" ;;
  figures_5pct) script="figure_scripts/bdev_integrated/04_add_five_percent_iqr.R" ;;
  *) echo "Unknown stage: ${stage}" >&2; exit 2 ;;
esac
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
OUTPUT_ROOT="${STORAGE}/bdev_integrated"
review_status="$(python3 - "${REPO_ROOT}/config/bdev_integrated.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle).get('review_status', ''))
PY
)"
if [[ "${review_status}" == awaiting_user_review && "${stage}" != figures && "${stage}" != figures_5pct ]]; then
  echo "Integration is awaiting figure-report review; only the read-only figures stage is enabled." >&2
  exit 2
fi
library="$(python3 - "${REPO_ROOT}/config/bundles.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    bundles = yaml.safe_load(handle)['bundles']
runtime = bundles['bdev_integrated'].get('runtime') or bundles['parse_cytokines']['runtime']
print(runtime['r_library'])
PY
)"
[[ -d "${library}/scworkbench" ]] || { echo "Registered runtime missing: ${library}" >&2; exit 1; }
SCWORKBENCH_RUNTIME_ARGS=(--bind "${library}:${library}:ro" --env "R_LIBS_USER=${library}")
if [[ "${stage}" == figures || "${stage}" == figures_5pct ]]; then
  mkdir -p "${OUTPUT_ROOT}/reports/${run_id}"
  OUTPUT_BIND=(--bind "${OUTPUT_ROOT}/reports/${run_id}:${OUTPUT_ROOT}/reports/${run_id}:rw")
else
  mkdir -p "${OUTPUT_ROOT}"
  OUTPUT_BIND=(--bind "${OUTPUT_ROOT}:${OUTPUT_ROOT}:rw")
fi
ulimit -v $((64 * 1024 * 1024))
exec singularity exec --cleanenv --no-home \
  "${SCWORKBENCH_RUNTIME_ARGS[@]}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  "${OUTPUT_BIND[@]}" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4 \
  --pwd "${REPO_ROOT}" "${STORAGE}/scworkbench-bdev-signac.sif" \
  Rscript --vanilla "${REPO_ROOT}/${script}" "${run_id}"
