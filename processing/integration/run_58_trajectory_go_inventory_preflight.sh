#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BUNDLE="${STORAGE}/bdev_integrated/bundle"
CYTOKINE_SOURCE="${STORAGE}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
MSIGDB_LIBRARY="${MSIGDB_LIBRARY:-/tmp/scworkbench_bdev_gsea_rlib}"
MSIGDB_CACHE="${MSIGDB_CACHE:-/tmp/scworkbench_bdev_gsea_cache}"
RUN_ID="${1:-trajectory_go_depth_adjusted_preflight_v1_20261011}"
OUTPUT_DIR="${STORAGE}/bdev_integrated/reports/${RUN_ID}"

[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ ]] || exit 2
[[ -x "${RUNTIME}" && -f "${SIF}" && -d "${BUNDLE}" ]] || exit 1
[[ -d "${REGISTERED_LIBRARY}/scworkbench" ]] || exit 1
[[ -d "${MSIGDB_LIBRARY}/msigdbr" ]] || exit 1
[[ -f "${CYTOKINE_SOURCE}/passing_signature_inventory.tsv" ]] || exit 1
[[ -f "${CYTOKINE_SOURCE}/passing_gene_membership.tsv" ]] || exit 1
[[ ! -e "${OUTPUT_DIR}" && ! -e "${OUTPUT_DIR}.incomplete" ]] || {
  echo "Output already exists: ${OUTPUT_DIR}" >&2; exit 1;
}

"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${MSIGDB_LIBRARY}:${REGISTERED_LIBRARY}" \
  --env "XDG_CACHE_HOME=${MSIGDB_CACHE}" \
  --env OMP_NUM_THREADS=2 --env OPENBLAS_NUM_THREADS=2 --env MKL_NUM_THREADS=2 \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/processing/integration/58_trajectory_go_inventory_preflight.R" \
  "${BUNDLE}" \
  "${CYTOKINE_SOURCE}/passing_signature_inventory.tsv" \
  "${CYTOKINE_SOURCE}/passing_gene_membership.tsv" \
  "${OUTPUT_DIR}"
