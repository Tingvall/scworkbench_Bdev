#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
IMAGE="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
OUTPUT="${1:-${REPO_ROOT}/figure_outputs/parse_cytokines/signal_inventory_$(date -u +%Y%m%dT%H%M%SZ)}"
source "${REPO_ROOT}/src/parse_runtime.sh"
configure_parse_runtime
mkdir -p "${OUTPUT}"
OUTPUT="$(cd -- "${OUTPUT}" && pwd)"
if [[ -n "$(ls -A "${OUTPUT}")" ]]; then
  echo "Output directory must be empty: ${OUTPUT}" >&2
  exit 2
fi
ulimit -v $((16 * 1024 * 1024))
exec singularity exec --cleanenv --no-home \
  "${SCWORKBENCH_RUNTIME_ARGS[@]}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${OUTPUT}:${OUTPUT}:rw" \
  --env OMP_NUM_THREADS=2 --env OPENBLAS_NUM_THREADS=2 --env MKL_NUM_THREADS=2 \
  --pwd "${REPO_ROOT}" "${IMAGE}" \
  Rscript --vanilla "${SCRIPT_DIR}/04_cytokine_signal_inventory.R" "${OUTPUT}"
