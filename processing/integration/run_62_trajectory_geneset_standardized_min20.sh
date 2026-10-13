#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BUNDLE="${STORAGE}/bdev_integrated/bundle"
PREFLIGHT="${STORAGE}/bdev_integrated/reports/trajectory_go_depth_adjusted_preflight_v1_20261011"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
RUN_ID="${1:-trajectory_go_gene_standardized_min20_perm1000_v1_20261011}"
GROUP="${2:-all}"
BASE_SEED="${BASE_SEED:-314159}"
N_WORKERS="${N_WORKERS:-8}"
OUTPUT_ROOT="${STORAGE}/bdev_integrated/reports/${RUN_ID}"

[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ ]] || exit 2
[[ "${GROUP}" =~ ^(all|bp|mf|cc|cytokines)$ ]] || exit 2
[[ "${BASE_SEED}" =~ ^[0-9]+$ ]] || exit 2
[[ "${N_WORKERS}" =~ ^[1-9][0-9]*$ ]] || exit 2
[[ -f "${PREFLIGHT}/completion.json" && -d "${BUNDLE}" ]] || exit 1
mkdir -p "${OUTPUT_ROOT}"

run_group() {
  local analysis_group="$1"
  local target="${OUTPUT_ROOT}/${analysis_group}"
  [[ ! -e "${target}" && ! -e "${target}.incomplete" ]] || {
    echo "Output or incomplete output exists: ${target}" >&2; return 1;
  }
  "${RUNTIME}" exec --cleanenv --no-home \
    --env "R_LIBS_USER=${REGISTERED_LIBRARY}" \
    --env "N_WORKERS=${N_WORKERS}" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
    --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
    --bind "${STORAGE}:${STORAGE}:ro" \
    --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
    --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
    "${REPO_ROOT}/processing/integration/62_trajectory_geneset_standardized_min20.R" \
    "${BUNDLE}" "${PREFLIGHT}/all_group_set_inventory.tsv.gz" \
    "${PREFLIGHT}/all_group_eligible_mapped_memberships.tsv.gz" \
    "${analysis_group}" "${target}" "${BASE_SEED}" 1000 \
    2>&1 | tee "${OUTPUT_ROOT}/${analysis_group}.log"
}

if [[ "${GROUP}" == "all" ]]; then
  for analysis_group in cytokines cc mf bp; do run_group "${analysis_group}"; done
else
  run_group "${GROUP}"
fi

