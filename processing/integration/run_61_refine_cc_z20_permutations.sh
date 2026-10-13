#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BUNDLE="${STORAGE}/bdev_integrated/bundle"
PREFLIGHT="${STORAGE}/bdev_integrated/reports/trajectory_go_depth_adjusted_preflight_v1_20261011"
SOURCE_ROOT="${STORAGE}/bdev_integrated/reports/trajectory_go_depth_adjusted_adaptive_v1_20261011"
OUTPUT="${SOURCE_ROOT}/targeted_permutation_pilots/cc_z20_100k"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
BASE_SEED="${BASE_SEED:-314159}"
N_WORKERS="${N_WORKERS:-8}"

[[ "${BASE_SEED}" =~ ^[0-9]+$ ]] || exit 2
[[ "${N_WORKERS}" =~ ^[1-9][0-9]*$ ]] || exit 2
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.incomplete" ]] || {
  echo "Output or incomplete output exists: ${OUTPUT}" >&2; exit 1;
}
mkdir -p "$(dirname -- "${OUTPUT}")"

"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${REGISTERED_LIBRARY}" \
  --env "N_WORKERS=${N_WORKERS}" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/processing/integration/61_refine_cc_z20_permutations.R" \
  "${BUNDLE}" "${PREFLIGHT}/all_group_set_inventory.tsv.gz" \
  "${PREFLIGHT}/all_group_eligible_mapped_memberships.tsv.gz" \
  "${SOURCE_ROOT}/cc" "${OUTPUT}" "${BASE_SEED}" 100000 \
  2>&1 | tee "${SOURCE_ROOT}/targeted_permutation_pilots/cc_z20_100k.log"
