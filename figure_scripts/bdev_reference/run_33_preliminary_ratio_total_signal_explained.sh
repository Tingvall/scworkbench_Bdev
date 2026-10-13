#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
BDEV_ROOT="${STORAGE}/bdev_reference"
SMOOTHED_REPORT="${BDEV_ROOT}/reports/preliminary_smoothed_slope_expression_w30_s5_span008_20261008"
TURN_PATTERN_REPORT="${BDEV_ROOT}/reports/preliminary_modality_turn_pattern_screen_v2_20261009"
DEPTH_REPORT="${BDEV_ROOT}/reports/preliminary_inverse_hits_depth_diagnostic_v1_20261009"
RUN_TAG="${RATIO_DEPTH_RUN_TAG:-preliminary_ratio_total_signal_explained_$(date -u +%Y%m%dT%H%M%SZ)}"
OUTPUT_DIR="${1:-${BDEV_ROOT}/reports/${RUN_TAG}}"

[[ -f "${SIF}" ]] || { echo "Container image missing: ${SIF}" >&2; exit 1; }
[[ -f "${DEPTH_REPORT}/smoothed_window_total_rna_atac.tsv" ]] || {
  echo "Depth report missing: ${DEPTH_REPORT}" >&2; exit 1; }
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

exec "${CONTAINER_RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:/repo:ro" \
  --bind "${BDEV_ROOT}:${BDEV_ROOT}:rw" \
  --pwd /repo "${SIF}" \
  Rscript --vanilla \
  "/repo/figure_scripts/bdev_reference/33_preliminary_ratio_total_signal_explained.R" \
  "${SMOOTHED_REPORT}" "${TURN_PATTERN_REPORT}" "${DEPTH_REPORT}" \
  "${OUTPUT_DIR}"
