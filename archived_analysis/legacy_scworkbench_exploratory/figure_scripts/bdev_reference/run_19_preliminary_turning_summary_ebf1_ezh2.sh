#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
BDEV_ROOT="${STORAGE}/bdev_reference"
BUNDLE="${BDEV_ROOT}/bundle"
WINDOW_REPORT="${BDEV_ROOT}/reports/rna_atac_latent_progress_w30_s5_v1_20261008"
SMOOTHED_REPORT="${BDEV_ROOT}/reports/preliminary_smoothed_slope_expression_w30_s5_span008_20261008"
TURNING_REPORT="${BDEV_ROOT}/reports/preliminary_slope_turning_regions_w30_s5_span008_20261008"
RUN_TAG="${TURNING_SUMMARY_RUN_TAG:-preliminary_turning_summary_ebf1_ezh2_w30_s5_$(date -u +%Y%m%dT%H%M%SZ)}"
OUTPUT_DIR="${1:-${BDEV_ROOT}/reports/${RUN_TAG}}"

EXPECTED_FINGERPRINT="$(python3 - "${REPO_ROOT}/config/bdev_reference.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle)['validation']['bundle_state_fingerprint'])
PY
)"

[[ -f "${SIF}" ]] || { echo "Container image missing: ${SIF}" >&2; exit 1; }
[[ -d "${BUNDLE}" ]] || { echo "Bdev bundle missing: ${BUNDLE}" >&2; exit 1; }
[[ -f "${TURNING_REPORT}/selected_turning_regions.tsv" ]] || {
  echo "Turning report missing: ${TURNING_REPORT}" >&2; exit 1; }
case "${OUTPUT_DIR}" in
  "${BDEV_ROOT}"/*) ;;
  *) echo "Output must be inside ${BDEV_ROOT}: ${OUTPUT_DIR}" >&2; exit 1 ;;
esac
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output exists: ${OUTPUT_DIR}" >&2; exit 1; }

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
  "/repo/figure_scripts/bdev_reference/19_preliminary_turning_summary_ebf1_ezh2.R" \
  "${BUNDLE}" "${WINDOW_REPORT}" "${SMOOTHED_REPORT}" "${TURNING_REPORT}" \
  "${OUTPUT_DIR}" "${EXPECTED_FINGERPRINT}"
