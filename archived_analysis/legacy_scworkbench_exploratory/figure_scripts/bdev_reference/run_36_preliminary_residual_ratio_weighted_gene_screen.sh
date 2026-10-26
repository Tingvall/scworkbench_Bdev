#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
BDEV_ROOT="${STORAGE}/bdev_reference"
BUNDLE="${BDEV_ROOT}/bundle"
WINDOW_REPORT="${BDEV_ROOT}/reports/rna_atac_latent_progress_w30_s5_v1_20261008"
RATIO_DEPTH_REPORT="${BDEV_ROOT}/reports/preliminary_ratio_total_signal_explained_v1_20261009"
RUN_TAG="${RESIDUAL_GENE_RUN_TAG:-preliminary_residual_ratio_weighted_gene_screen_$(date -u +%Y%m%dT%H%M%SZ)}"
OUTPUT_DIR="${1:-${BDEV_ROOT}/reports/${RUN_TAG}}"
BATCH_SIZE="${RESIDUAL_GENE_BATCH_SIZE:-256}"

EXPECTED_FINGERPRINT="$(python3 - "${REPO_ROOT}/config/bdev_reference.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle)['validation']['bundle_state_fingerprint'])
PY
)"

[[ -f "${SIF}" ]] || { echo "Container image missing: ${SIF}" >&2; exit 1; }
[[ -d "${BUNDLE}" ]] || { echo "Bdev bundle missing: ${BUNDLE}" >&2; exit 1; }
[[ -f "${RATIO_DEPTH_REPORT}/ratio_total_signal_fitted_and_residual_windows.tsv" ]] || {
  echo "Ratio-depth report missing: ${RATIO_DEPTH_REPORT}" >&2; exit 1; }
case "${OUTPUT_DIR}" in
  "${BDEV_ROOT}"/*) ;;
  *) echo "Output must be inside ${BDEV_ROOT}: ${OUTPUT_DIR}" >&2; exit 1 ;;
esac
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output exists: ${OUTPUT_DIR}" >&2; exit 1; }

LOCAL_APPTAINER="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
if [[ -x "${LOCAL_APPTAINER}" ]]; then
  CONTAINER_RUNTIME="${LOCAL_APPTAINER}"
elif command -v apptainer >/dev/null 2>&1; then
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
  "/repo/figure_scripts/bdev_reference/36_preliminary_residual_ratio_weighted_gene_screen.R" \
  "${BUNDLE}" "${WINDOW_REPORT}" "${RATIO_DEPTH_REPORT}" "${OUTPUT_DIR}" \
  "${EXPECTED_FINGERPRINT}" "${BATCH_SIZE}"
