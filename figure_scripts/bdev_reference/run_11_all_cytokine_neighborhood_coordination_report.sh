#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
BDEV_ROOT="${STORAGE}/bdev_reference"
BUNDLE="${BDEV_ROOT}/bundle"
UCELL_REPORT="${BDEV_ROOT}/reports/all_cytokine_refined_ucell_v1_20260925"
BALANCED_REPORT="${BDEV_ROOT}/reports/all_cytokine_three_method_comparison_v1_20260928"
DRIVER_REPORT="${BDEV_ROOT}/reports/all_cytokine_refined_drivers_v5_mean_internal_scale_bb_20260928"
RUN_TAG="${NEIGHBOR_COORDINATION_RUN_TAG:-all_cytokine_neighborhood_coordination_v1_$(date -u +%Y%m%dT%H%M%SZ)}"
OUTPUT_DIR="${1:-${BDEV_ROOT}/reports/${RUN_TAG}}"

EXPECTED_FINGERPRINT="$(python3 - "${REPO_ROOT}/config/bdev_reference.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle)['validation']['bundle_state_fingerprint'])
PY
)"

[[ -f "${SIF}" ]] || { echo "Container image does not exist: ${SIF}" >&2; exit 1; }
[[ -d "${BUNDLE}" ]] || { echo "Bdev bundle does not exist: ${BUNDLE}" >&2; exit 1; }
[[ -f "${UCELL_REPORT}/completion.json" ]] || { echo "UCell report missing" >&2; exit 1; }
[[ -f "${BALANCED_REPORT}/completion.json" ]] || { echo "Balanced report missing" >&2; exit 1; }
[[ -f "${DRIVER_REPORT}/completion.json" ]] || { echo "Driver report missing" >&2; exit 1; }
case "${OUTPUT_DIR}" in
  "${BDEV_ROOT}"/*) ;;
  *) echo "Output must be inside ${BDEV_ROOT}: ${OUTPUT_DIR}" >&2; exit 1 ;;
esac
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output already exists: ${OUTPUT_DIR}" >&2; exit 1; }

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="singularity"
else
  echo "Apptainer or Singularity is required." >&2
  exit 1
fi

ulimit -v $((64 * 1024 * 1024))
exec "${CONTAINER_RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:/repo:ro" \
  --bind "${BDEV_ROOT}:${BDEV_ROOT}:rw" \
  --pwd /repo "${SIF}" \
  Rscript --vanilla \
  "/repo/figure_scripts/bdev_reference/11_all_cytokine_neighborhood_coordination_report.R" \
  "${BUNDLE}" "${UCELL_REPORT}" "${BALANCED_REPORT}" "${DRIVER_REPORT}" \
  "${OUTPUT_DIR}" "${EXPECTED_FINGERPRINT}"
