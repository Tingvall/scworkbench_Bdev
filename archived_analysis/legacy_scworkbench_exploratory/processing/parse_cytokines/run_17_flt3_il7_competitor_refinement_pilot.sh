#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
REST_ROOT="${STORAGE}/parse_cytokines_rest_only"
IMAGE="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
OUTPUT_DIR="${1:-${REST_ROOT}/reports/rest_only_20260924/flt3_il7_competitor_refinement_top100_v1}"
PARSE_THREADS="${PARSE_THREADS:-4}"
PARSE_MAX_MEMORY_GIB="${PARSE_MAX_MEMORY_GIB:-48}"

source "${REPO_ROOT}/src/parse_runtime.sh"
configure_parse_runtime

for value_name in PARSE_THREADS PARSE_MAX_MEMORY_GIB; do
  value="${!value_name}"
  [[ "${value}" =~ ^[1-9][0-9]*$ ]] || { echo "${value_name} must be a positive integer." >&2; exit 1; }
done

TOTAL_MEMORY_KIB="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
MEMORY_LIMIT_KIB="$((PARSE_MAX_MEMORY_GIB * 1024 * 1024))"
MAX_SHARED_SERVER_KIB="$((TOTAL_MEMORY_KIB * 40 / 100))"
(( MEMORY_LIMIT_KIB <= MAX_SHARED_SERVER_KIB )) || {
  echo "Requested memory ceiling exceeds 40% of server RAM; refusing to start." >&2
  exit 1
}

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="singularity"
else
  echo "Apptainer or Singularity is required." >&2
  exit 1
fi
[[ -f "${IMAGE}" ]] || { echo "Container image does not exist: ${IMAGE}" >&2; exit 1; }
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output already exists: ${OUTPUT_DIR}" >&2; exit 1; }

ulimit -v "${MEMORY_LIMIT_KIB}"
exec "${CONTAINER_RUNTIME}" exec --cleanenv --no-home \
  "${SCWORKBENCH_RUNTIME_ARGS[@]}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${REST_ROOT}:${REST_ROOT}:rw" \
  --env "OMP_NUM_THREADS=${PARSE_THREADS}" \
  --env "OPENBLAS_NUM_THREADS=${PARSE_THREADS}" \
  --env "MKL_NUM_THREADS=${PARSE_THREADS}" \
  --env "PARSE_THREADS=${PARSE_THREADS}" \
  --pwd "${REPO_ROOT}" \
  "${IMAGE}" \
  Rscript --vanilla \
  "${SCRIPT_DIR}/17_flt3_il7_competitor_refinement_pilot.R" \
  "${OUTPUT_DIR}" 100
