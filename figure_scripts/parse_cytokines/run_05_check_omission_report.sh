#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
OUTPUT="${1:-${REPO_ROOT}/figure_outputs/parse_cytokines/check_omissions_$(date -u +%Y%m%dT%H%M%SZ)}"
source "${REPO_ROOT}/src/parse_runtime.sh"
configure_parse_runtime
mkdir -p "${OUTPUT}"
OUTPUT="$(cd -- "${OUTPUT}" && pwd)"
[[ -z "$(ls -A "${OUTPUT}")" ]] || { echo "Output must be empty: ${OUTPUT}" >&2; exit 2; }
ulimit -v $((8 * 1024 * 1024))
exec singularity exec --cleanenv --no-home "${SCWORKBENCH_RUNTIME_ARGS[@]}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${OUTPUT}:${OUTPUT}:rw" --env OMP_NUM_THREADS=2 --env OPENBLAS_NUM_THREADS=2 \
  --pwd "${REPO_ROOT}" "${STORAGE}/scworkbench-bdev-signac.sif" \
  Rscript --vanilla "${SCRIPT_DIR}/05_check_omission_report.R" "${OUTPUT}"
