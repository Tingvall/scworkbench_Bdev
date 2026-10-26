#!/usr/bin/env bash
set -euo pipefail

# Run the outside-scWorkbench Seurat prefilter in the supplied image.
# Optional positional arguments:
#   1. Unfiltered input Seurat RDS
#   2. New prefiltered output Seurat RDS

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

IMAGE="${SCWORKBENCH_IMAGE:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif}"
PREPARED_DIR="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines/prepared"
INPUT_RDS="${1:-${PREPARED_DIR}/parse_cytokines_unfiltered_seurat.rds}"
OUTPUT_RDS="${2:-${PREPARED_DIR}/parse_cytokines_prefiltered_min500counts_seurat.rds}"
PARSE_THREADS="${PARSE_THREADS:-4}"
PARSE_MAX_MEMORY_GIB="${PARSE_MAX_MEMORY_GIB:-160}"

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
if [[ "$(realpath -- "${INPUT_RDS}")" == "$(realpath -m -- "${OUTPUT_RDS}")" ]]; then
  echo "Input and output RDS paths must differ." >&2
  exit 1
fi
if [[ -e "${OUTPUT_RDS}" ]]; then
  echo "Output RDS already exists; refusing to overwrite: ${OUTPUT_RDS}" >&2
  exit 1
fi

OUTPUT_DIR="$(dirname -- "${OUTPUT_RDS}")"
if [[ ! -d "${OUTPUT_DIR}" ]]; then
  echo "Output directory does not exist: ${OUTPUT_DIR}" >&2
  exit 1
fi

IMAGE_SHA256="$(sha256sum -- "${IMAGE}" | awk '{print $1}')"
ulimit -v "${MEMORY_LIMIT_KIB}"

exec "${CONTAINER_RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${OUTPUT_DIR}:${OUTPUT_DIR}:rw" \
  --bind "${INPUT_RDS}:${INPUT_RDS}:ro" \
  --env "SCWORKBENCH_IMAGE=${IMAGE}" \
  --env "SCWORKBENCH_IMAGE_SHA256=${IMAGE_SHA256}" \
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
  "${SCRIPT_DIR}/02_prefilter_min500counts.R" \
  "${INPUT_RDS}" \
  "${OUTPUT_RDS}"
