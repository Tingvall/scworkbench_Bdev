#!/usr/bin/env bash
set -euo pipefail

# Create the parse_cytokines scWorkbench bundle from the accepted pre-filtered
# Seurat RDS. Optional positional arguments:
#   1. Pre-filtered input Seurat RDS
#   2. New bundle directory
#   3. New report directory

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

IMAGE="${SCWORKBENCH_IMAGE:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif}"
PARSE_ROOT="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines"
INPUT_RDS="${1:-${PARSE_ROOT}/prepared/parse_cytokines_prefiltered_min500counts_seurat.rds}"
BUNDLE_DIR="${2:-${PARSE_ROOT}/bundle}"
RUN_ID="bundle_import_$(date -u +%Y%m%dT%H%M%SZ)"
REPORT_DIR="${3:-${PARSE_ROOT}/reports/${RUN_ID}}"
PARSE_THREADS="${PARSE_THREADS:-4}"
PARSE_MAX_MEMORY_GIB="${PARSE_MAX_MEMORY_GIB:-160}"
SCWORKBENCH_CREATED_BY="${SCWORKBENCH_CREATED_BY:-jtg}"

if [[ ! "${PARSE_THREADS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "PARSE_THREADS must be a positive integer." >&2
  exit 1
fi
if [[ ! "${PARSE_MAX_MEMORY_GIB}" =~ ^[1-9][0-9]*$ ]]; then
  echo "PARSE_MAX_MEMORY_GIB must be a positive integer." >&2
  exit 1
fi

TOTAL_MEMORY_KIB="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
MEMORY_LIMIT_KIB="$((PARSE_MAX_MEMORY_GIB * 1024 * 1024))"
MAX_SHARED_SERVER_KIB="$((TOTAL_MEMORY_KIB * 40 / 100))"
if (( MEMORY_LIMIT_KIB > MAX_SHARED_SERVER_KIB )); then
  echo "Requested memory ceiling (${PARSE_MAX_MEMORY_GIB} GiB) exceeds 40% of server RAM; refusing to start." >&2
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
if [[ ! -f "${INPUT_RDS}" ]]; then
  echo "Input Seurat RDS does not exist: ${INPUT_RDS}" >&2
  exit 1
fi
if [[ -e "${BUNDLE_DIR}" ]]; then
  echo "Bundle destination already exists; refusing to overwrite: ${BUNDLE_DIR}" >&2
  exit 1
fi
if [[ -e "${REPORT_DIR}" ]]; then
  echo "Report destination already exists; refusing to overwrite: ${REPORT_DIR}" >&2
  exit 1
fi

BUNDLE_PARENT="$(dirname -- "${BUNDLE_DIR}")"
REPORT_PARENT="$(dirname -- "${REPORT_DIR}")"
if [[ ! -d "${BUNDLE_PARENT}" ]]; then
  echo "Bundle parent directory does not exist: ${BUNDLE_PARENT}" >&2
  exit 1
fi
if [[ ! -d "${REPORT_PARENT}" ]]; then
  echo "Report parent directory does not exist: ${REPORT_PARENT}" >&2
  exit 1
fi

IMAGE_SHA256="$(sha256sum -- "${IMAGE}" | awk '{print $1}')"
ulimit -v "${MEMORY_LIMIT_KIB}"

exec "${CONTAINER_RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${BUNDLE_PARENT}:${BUNDLE_PARENT}:rw" \
  --bind "${REPORT_PARENT}:${REPORT_PARENT}:rw" \
  --bind "${INPUT_RDS}:${INPUT_RDS}:ro" \
  --env "SCWORKBENCH_IMAGE=${IMAGE}" \
  --env "SCWORKBENCH_IMAGE_SHA256=${IMAGE_SHA256}" \
  --env "SCWORKBENCH_CREATED_BY=${SCWORKBENCH_CREATED_BY}" \
  --env "PARSE_THREADS=${PARSE_THREADS}" \
  --env "PARSE_MAX_MEMORY_GIB=${PARSE_MAX_MEMORY_GIB}" \
  --env "OMP_NUM_THREADS=${PARSE_THREADS}" \
  --env "OPENBLAS_NUM_THREADS=${PARSE_THREADS}" \
  --env "MKL_NUM_THREADS=${PARSE_THREADS}" \
  --env "VECLIB_MAXIMUM_THREADS=${PARSE_THREADS}" \
  --env "NUMEXPR_NUM_THREADS=${PARSE_THREADS}" \
  --pwd "${REPO_ROOT}" \
  "${IMAGE}" \
  Rscript --vanilla \
  "${SCRIPT_DIR}/03_create_parse_cytokines.R" \
  "${INPUT_RDS}" \
  "${BUNDLE_DIR}" \
  "${REPORT_DIR}"

