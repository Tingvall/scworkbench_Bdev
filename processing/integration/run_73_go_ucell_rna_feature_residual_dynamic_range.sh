#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BUNDLE="${STORAGE}/bdev_integrated/bundle"
PREFLIGHT="${STORAGE}/bdev_integrated/reports/trajectory_go_depth_adjusted_preflight_v1_20261011"
CYTOKINE="${STORAGE}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925"
PRIOR="${STORAGE}/bdev_integrated/reports/go_ucell_median_perm1000_v1_20261012"
RANKINGS="${PRIOR}/trajectory_cell_rankings_maxrank1500.rds"
RUN_ID="${1:-go_ucell_rna_feature_residual_dynamic_range_perm1000_v1_20261012}"
OUTPUT="${STORAGE}/bdev_integrated/reports/${RUN_ID}"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
BASE_SEED="${BASE_SEED:-314159}"
N_WORKERS="${N_WORKERS:-8}"

[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ ]] || exit 2
[[ "${BASE_SEED}" =~ ^[0-9]+$ && "${N_WORKERS}" =~ ^[1-9][0-9]*$ ]] || exit 2
[[ -d "${BUNDLE}" && -f "${PREFLIGHT}/completion.json" && -f "${RANKINGS}" ]] || exit 1
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.incomplete" ]] || exit 1

"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${REGISTERED_LIBRARY}" \
  --env "N_WORKERS=${N_WORKERS}" \
  --env "GO_UCELL_INFERENCE_MODE=dynamic_range" \
  --env "GO_UCELL_COMPLEXITY_RESIDUALIZE=true" \
  --env "GO_UCELL_RANKINGS_PATH=${RANKINGS}" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/processing/integration/65_go_ucell_median_permutation.R" \
  "${BUNDLE}" "${PREFLIGHT}/all_group_set_inventory.tsv.gz" \
  "${PREFLIGHT}/all_group_eligible_mapped_memberships.tsv.gz" \
  "${CYTOKINE}" "${OUTPUT}" "${BASE_SEED}" 1000 \
  2>&1 | tee "${STORAGE}/bdev_integrated/reports/${RUN_ID}.log"
