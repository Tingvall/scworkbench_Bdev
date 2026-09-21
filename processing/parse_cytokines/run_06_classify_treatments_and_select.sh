#!/usr/bin/env bash
set -euo pipefail

# Classify every Parse sample/cell by treatment category and derive the
# QC-accepted control/single downstream analysis population.
# Optional positional arguments:
#   1. Existing parse_cytokines bundle directory
#   2. New run-specific report directory
#   3. Reviewed treatment-category mapping TSV

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

IMAGE="${SCWORKBENCH_IMAGE:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif}"
PARSE_ROOT="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines"
BUNDLE_DIR="${1:-${PARSE_ROOT}/bundle}"
RUN_ID="treatment_selection_$(date -u +%Y%m%dT%H%M%SZ)"
REPORT_DIR="${2:-${PARSE_ROOT}/reports/${RUN_ID}}"
TREATMENT_MAP="${3:-${REPO_ROOT}/config/parse_cytokines_treatment_categories.tsv}"
PARSE_THREADS="${PARSE_THREADS:-4}"
PARSE_MAX_MEMORY_GIB="${PARSE_MAX_MEMORY_GIB:-96}"
SCWORKBENCH_CREATED_BY="${SCWORKBENCH_CREATED_BY:-jtg}"

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
if [[ ! -d "${BUNDLE_DIR}" ]]; then
  echo "Bundle directory does not exist: ${BUNDLE_DIR}" >&2
  exit 1
fi
if [[ ! -f "${TREATMENT_MAP}" ]]; then
  echo "Treatment mapping does not exist: ${TREATMENT_MAP}" >&2
  exit 1
fi
if [[ -e "${REPORT_DIR}" ]]; then
  echo "Report destination already exists; refusing to overwrite: ${REPORT_DIR}" >&2
  exit 1
fi
if [[ ! -d "$(dirname -- "${REPORT_DIR}")" ]]; then
  echo "Report parent directory does not exist: $(dirname -- "${REPORT_DIR}")" >&2
  exit 1
fi

ulimit -v "${MEMORY_LIMIT_KIB}"

exec "${CONTAINER_RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${PARSE_ROOT}:${PARSE_ROOT}:rw" \
  --env "SCWORKBENCH_CREATED_BY=${SCWORKBENCH_CREATED_BY}" \
  --env "OMP_NUM_THREADS=${PARSE_THREADS}" \
  --env "OPENBLAS_NUM_THREADS=${PARSE_THREADS}" \
  --env "MKL_NUM_THREADS=${PARSE_THREADS}" \
  --env "VECLIB_MAXIMUM_THREADS=${PARSE_THREADS}" \
  --env "NUMEXPR_NUM_THREADS=${PARSE_THREADS}" \
  --pwd "${REPO_ROOT}" \
  "${IMAGE}" \
  Rscript --vanilla \
  "${SCRIPT_DIR}/06_classify_treatments_and_select.R" \
  "${BUNDLE_DIR}" \
  "${REPORT_DIR}" \
  "${TREATMENT_MAP}"
