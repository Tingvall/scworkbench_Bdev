#!/usr/bin/env bash
set -euo pipefail

# Render one treatment-by-sample violin page for every scored discovery gene
# set and document every exclusion. Optional arg:
#   1. Existing parse_cytokines bundle directory
#   2. Existing completed signature report directory

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

IMAGE="${SCWORKBENCH_IMAGE:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif}"
PARSE_ROOT="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines"
source "${REPO_ROOT}/src/parse_runtime.sh"
configure_parse_runtime "${1:-}"
REPORT_DIR="${2:-${PARSE_ROOT}/reports/signature_report_20260921T221910Z}"
PARSE_THREADS="${PARSE_THREADS:-2}"
PARSE_MAX_MEMORY_GIB="${PARSE_MAX_MEMORY_GIB:-8}"

for value_name in PARSE_THREADS PARSE_MAX_MEMORY_GIB; do
  value="${!value_name}"
  if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
    echo "${value_name} must be a positive integer." >&2
    exit 1
  fi
done

TOTAL_MEMORY_KIB="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
MEMORY_LIMIT_KIB="$((PARSE_MAX_MEMORY_GIB * 1024 * 1024))"
MAX_SHARED_SERVER_KIB="$((TOTAL_MEMORY_KIB * 40 / 100))"
if (( MEMORY_LIMIT_KIB > MAX_SHARED_SERVER_KIB )); then
  echo "Requested memory ceiling exceeds 40% of server RAM; refusing to start." >&2
  exit 1
fi

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="singularity"
else
  echo "Apptainer or Singularity is required." >&2
  exit 1
fi
if [[ ! -f "${IMAGE}" ]]; then
  echo "Container image does not exist: ${IMAGE}" >&2
  exit 1
fi
if [[ ! -d "${BUNDLE_DIR}" ]]; then
  echo "Bundle directory does not exist: ${BUNDLE_DIR}" >&2
  exit 1
fi
if [[ ! -f "${REPORT_DIR}/all_signature_status.tsv" || ! -f "${REPORT_DIR}/ucell_sample_medians.tsv" ]]; then
  echo "A completed signature report with status and sample-median tables is required: ${REPORT_DIR}" >&2
  exit 1
fi

ulimit -v "${MEMORY_LIMIT_KIB}"

exec "${CONTAINER_RUNTIME}" exec --cleanenv --no-home \
  "${SCWORKBENCH_RUNTIME_ARGS[@]}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${PARSE_ROOT}:${PARSE_ROOT}:rw" \
  --env "OMP_NUM_THREADS=${PARSE_THREADS}" \
  --env "OPENBLAS_NUM_THREADS=${PARSE_THREADS}" \
  --env "MKL_NUM_THREADS=${PARSE_THREADS}" \
  --pwd "${REPO_ROOT}" \
  "${IMAGE}" \
  Rscript --vanilla \
  "${SCRIPT_DIR}/03_signature_sample_treatment_violins.R" \
  "${BUNDLE_DIR}" \
  "${REPORT_DIR}"
