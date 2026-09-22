#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

IMAGE="${SCWORKBENCH_IMAGE:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif}"
PARSE_ROOT="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines"
BUNDLE_DIR="${1:-${PARSE_ROOT}/bundle}"
RUN_ID="signature_states_$(date -u +%Y%m%dT%H%M%SZ)"
REPORT_DIR="${2:-${PARSE_ROOT}/reports/${RUN_ID}}"
SESSION_ROOT="${3:-${PARSE_ROOT}/sessions}"
RESUME_WORKSPACE="${4:-}"
PARSE_THREADS="${PARSE_THREADS:-8}"
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

for required in "${IMAGE}" "${BUNDLE_DIR}"; do
  if [[ ! -e "${required}" ]]; then
    echo "Required input does not exist: ${required}" >&2
    exit 1
  fi
done
if [[ -e "${REPORT_DIR}" ]]; then
  echo "Report destination already exists; refusing to overwrite: ${REPORT_DIR}" >&2
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
  "${SCRIPT_DIR}/07_build_signature_states.R" \
  "${BUNDLE_DIR}" \
  "${REPORT_DIR}" \
  "${SESSION_ROOT}" \
  ${RESUME_WORKSPACE:+"${RESUME_WORKSPACE}"}
