#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

IMAGE="${SCWORKBENCH_IMAGE:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif}"
PARSE_ROOT="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines"
BUNDLE_DIR="${1:-${PARSE_ROOT}/bundle}"
STATE_REPORT_DIR="${2:-}"
RUN_ID="signature_discovery_$(date -u +%Y%m%dT%H%M%SZ)"
OUTPUT_DIR="${3:-${PARSE_ROOT}/reports/${RUN_ID}}"
PARSE_THREADS="${PARSE_THREADS:-8}"
PARSE_MAX_MEMORY_GIB="${PARSE_MAX_MEMORY_GIB:-96}"

if [[ -z "${STATE_REPORT_DIR}" ]]; then
  STATE_REPORT_DIR="$(find "${PARSE_ROOT}/reports" -maxdepth 1 -type d -name 'signature_states_*' -print | sort | tail -n 1)"
fi
if [[ -z "${STATE_REPORT_DIR}" || ! -f "${STATE_REPORT_DIR}/completion.json" ]]; then
  echo "A completed signature_states report directory is required." >&2
  exit 1
fi
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
if [[ -e "${OUTPUT_DIR}" ]]; then
  echo "Output destination already exists; refusing to overwrite: ${OUTPUT_DIR}" >&2
  exit 1
fi

ulimit -v "${MEMORY_LIMIT_KIB}"

exec "${CONTAINER_RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${PARSE_ROOT}:${PARSE_ROOT}:rw" \
  --env "PARSE_THREADS=${PARSE_THREADS}" \
  --env "OMP_NUM_THREADS=${PARSE_THREADS}" \
  --env "OPENBLAS_NUM_THREADS=${PARSE_THREADS}" \
  --env "MKL_NUM_THREADS=${PARSE_THREADS}" \
  --pwd "${REPO_ROOT}" \
  "${IMAGE}" \
  Rscript --vanilla \
  "${SCRIPT_DIR}/08_discover_specific_signatures.R" \
  "${BUNDLE_DIR}" \
  "${STATE_REPORT_DIR}" \
  "${OUTPUT_DIR}"
