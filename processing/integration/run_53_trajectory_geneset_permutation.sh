#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BUNDLE="${STORAGE}/bdev_integrated/bundle"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
MSIGDB_LIBRARY="${MSIGDB_LIBRARY:-/tmp/scworkbench_bdev_gsea_rlib}"
MSIGDB_CACHE="${MSIGDB_CACHE:-/tmp/scworkbench_bdev_gsea_cache}"
RUN_ID="${1:-trajectory_geneset_permutation_20261011}"
GROUP="${2:-all}"
N_PERMUTATIONS="${N_PERMUTATIONS:-1000}"
MAX_SETS="${MAX_SETS:-}"
OUTPUT_ROOT="${STORAGE}/bdev_integrated/reports/${RUN_ID}"

[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ ]] || { echo "Invalid run ID" >&2; exit 2; }
[[ "${GROUP}" =~ ^(all|cytokines|bp|mf|cc)$ ]] || { echo "Invalid group: ${GROUP}" >&2; exit 2; }
[[ "${N_PERMUTATIONS}" =~ ^[1-9][0-9]*$ ]] || { echo "N_PERMUTATIONS must be positive" >&2; exit 2; }
[[ -x "${RUNTIME}" && -f "${SIF}" && -d "${BUNDLE}" ]] || exit 1
[[ -d "${REGISTERED_LIBRARY}/scworkbench" ]] || { echo "Registered scWorkbench runtime is missing" >&2; exit 1; }
[[ -d "${MSIGDB_LIBRARY}/msigdbr" ]] || { echo "msigdbr runtime is missing: ${MSIGDB_LIBRARY}" >&2; exit 1; }
mkdir -p "${OUTPUT_ROOT}"

run_group() {
  local analysis_group="$1"
  local target="${OUTPUT_ROOT}/${analysis_group}"
  local optional_limit=()
  if [[ -n "${MAX_SETS}" ]]; then optional_limit=("${MAX_SETS}"); fi
  "${RUNTIME}" exec --cleanenv --no-home \
    --env "R_LIBS_USER=${MSIGDB_LIBRARY}:${REGISTERED_LIBRARY}" \
    --env "XDG_CACHE_HOME=${MSIGDB_CACHE}" \
    --env OMP_NUM_THREADS=2 --env OPENBLAS_NUM_THREADS=2 --env MKL_NUM_THREADS=2 \
    --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
    --bind "${STORAGE}:${STORAGE}:ro" \
    --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
    --bind "${OUTPUT_ROOT}:${OUTPUT_ROOT}:rw" \
    --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
    "${REPO_ROOT}/processing/integration/53_trajectory_geneset_permutation_group.R" \
    "${BUNDLE}" "${target}" "${analysis_group}" "${N_PERMUTATIONS}" \
    "${optional_limit[@]}"
}

if [[ "${GROUP}" == "all" ]]; then
  for analysis_group in cytokines bp mf cc; do
    run_group "${analysis_group}" 2>&1 | tee "${OUTPUT_ROOT}/${analysis_group}.log"
  done
else
  run_group "${GROUP}" 2>&1 | tee "${OUTPUT_ROOT}/${GROUP}.log"
fi

python3 "${REPO_ROOT}/processing/integration/summarize_trajectory_geneset_permutation.py" "${OUTPUT_ROOT}"
