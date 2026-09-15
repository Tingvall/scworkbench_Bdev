#!/usr/bin/env bash
set -euo pipefail

# Run the outside-bundle Parse preparation in the supplied, versioned image.
# Optional positional arguments:
#   1. Parse all-sample directory (or DGE_unfiltered directory)
#   2. Empty output directory

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

IMAGE="${SCWORKBENCH_IMAGE:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif}"
INPUT_ROOT="${1:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/datasets/parse_cytokine_stimulation/parse_cytokine_stimulation_tailmaker_pipeline_output/all-sample}"
OUTPUT_DIR="${2:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines/prepared}"
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
if [[ ! -d "${INPUT_ROOT}" ]]; then
  echo "Parse input directory does not exist: ${INPUT_ROOT}" >&2
  exit 1
fi

OUTPUT_PARENT="$(dirname -- "${OUTPUT_DIR}")"
mkdir -p -- "${OUTPUT_PARENT}"
IMAGE_SHA256="$(sha256sum -- "${IMAGE}" | awk '{print $1}')"
ulimit -v "${MEMORY_LIMIT_KIB}"

exec "${CONTAINER_RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${INPUT_ROOT}:${INPUT_ROOT}:ro" \
  --bind "${OUTPUT_PARENT}:${OUTPUT_PARENT}:rw" \
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
  "${SCRIPT_DIR}/01_readparsebio_unfiltered.R" \
  "${INPUT_ROOT}" \
  "${OUTPUT_DIR}"
