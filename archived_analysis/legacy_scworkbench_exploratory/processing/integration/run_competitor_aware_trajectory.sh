#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
stage="${1:?Usage: run_competitor_aware_trajectory.sh score|figures RUN_ID}"
run_id="${2:?A stable run ID is required}"
[[ "${run_id}" =~ ^[A-Za-z0-9_]+$ ]] || exit 2

case "${stage}" in
  score) script="processing/integration/06_score_competitor_aware_signatures.R" ;;
  figures) script="figure_scripts/bdev_integrated/05_competitor_aware_sliding_window_report.R" ;;
  *) echo "Unknown stage: ${stage}" >&2; exit 2 ;;
esac

STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
OUTPUT_ROOT="${STORAGE}/bdev_integrated"
library="$(python3 - "${REPO_ROOT}/config/bundles.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    bundles = yaml.safe_load(handle)['bundles']
print(bundles['bdev_integrated']['runtime']['r_library'])
PY
)"
[[ -d "${library}/scworkbench" ]] || { echo "Registered runtime missing: ${library}" >&2; exit 1; }
mkdir -p "${OUTPUT_ROOT}/reports/${run_id}" "${OUTPUT_ROOT}/sessions/${run_id}"

ulimit -v $((64 * 1024 * 1024))
exec singularity exec --cleanenv --no-home \
  --bind "${library}:${library}:ro" --env "R_LIBS_USER=${library}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${OUTPUT_ROOT}:${OUTPUT_ROOT}:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4 \
  --pwd "${REPO_ROOT}" "${STORAGE}/scworkbench-bdev-signac.sif" \
  Rscript --vanilla "${REPO_ROOT}/${script}" "${run_id}"
