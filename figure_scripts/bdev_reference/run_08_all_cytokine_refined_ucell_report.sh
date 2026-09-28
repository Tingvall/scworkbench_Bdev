#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
BDEV_ROOT="${STORAGE}/bdev_reference"
BUNDLE="${BDEV_ROOT}/bundle"
SOURCE_DIR="${1:?Usage: run_08_all_cytokine_refined_ucell_report.sh <parse_refinement_report> [output_dir]}"
RUN_TAG="${ALL_CYTOKINE_BDEV_RUN_TAG:-all_cytokine_refined_ucell_v1_$(date -u +%Y%m%dT%H%M%SZ)}"
OUTPUT_DIR="${2:-${BDEV_ROOT}/reports/${RUN_TAG}}"
BDEV_THREADS="${BDEV_THREADS:-4}"

EXPECTED_FINGERPRINT="$(python3 - "${REPO_ROOT}/config/bdev_reference.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle)['validation']['bundle_state_fingerprint'])
PY
)"

[[ -f "${SIF}" ]] || { echo "Container image does not exist: ${SIF}" >&2; exit 1; }
[[ -d "${BUNDLE}" ]] || { echo "Bdev bundle does not exist: ${BUNDLE}" >&2; exit 1; }
[[ "${BDEV_THREADS}" =~ ^[1-9][0-9]*$ ]] || {
  echo "BDEV_THREADS must be a positive integer." >&2
  exit 1
}
[[ -f "${SOURCE_DIR}/completion.json" && -f "${SOURCE_DIR}/signature_summary.tsv" && -f "${SOURCE_DIR}/final_gene_membership.tsv" ]] || {
  echo "Incomplete Parse refinement report: ${SOURCE_DIR}" >&2
  exit 1
}
SOURCE_DIR="$(cd -- "${SOURCE_DIR}" && pwd)"
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
  --bind "${SOURCE_DIR}:/source:ro" \
  --env "BDEV_THREADS=${BDEV_THREADS}" \
  --pwd /repo "${SIF}" \
  Rscript --vanilla \
  "/repo/figure_scripts/bdev_reference/08_all_cytokine_refined_ucell_report.R" \
  "${BUNDLE}" /source "${OUTPUT_DIR}" "${EXPECTED_FINGERPRINT}"
