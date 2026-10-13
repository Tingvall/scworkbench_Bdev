#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BUNDLE="${STORAGE}/bdev_integrated/bundle"
SOURCE="${STORAGE}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925"
RUN_ID="${1:-cytokine_ucell_rna_feature_residual_dynamic_range_perm1000_v1_20261012}"
OUTPUT="${STORAGE}/bdev_integrated/reports/${RUN_ID}"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
BASE_SEED="${BASE_SEED:-314159}"
N_WORKERS="${N_WORKERS:-8}"

[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ ]] || exit 2
[[ "${BASE_SEED}" =~ ^[0-9]+$ ]] || exit 2
[[ "${N_WORKERS}" =~ ^[1-9][0-9]*$ ]] || exit 2
[[ -d "${BUNDLE}" && -f "${SOURCE}/completion.json" ]] || exit 1
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.incomplete" ]] || exit 1

"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${REGISTERED_LIBRARY}" \
  --env "N_WORKERS=${N_WORKERS}" \
  --env "CYTOKINE_UCELL_COMPLEXITY_RESIDUALIZE=true" \
  --env "CYTOKINE_UCELL_BUNDLE=${BUNDLE}" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/processing/integration/66_cytokine_ucell_dynamic_range.R" \
  "${SOURCE}" "${OUTPUT}" "${BASE_SEED}" 1000 \
  2>&1 | tee "${STORAGE}/bdev_integrated/reports/${RUN_ID}.log"
