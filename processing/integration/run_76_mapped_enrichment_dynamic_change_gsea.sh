#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
INPUT="${STORAGE}/bdev_integrated/reports/flt3l_il7_ko_vs_wt_pseudotime_matched_gsea_v1_20261012"
RUN_ID="${1:-flt3l_il7_ko_vs_wt_dynamic_change_gsea_v1_20261012}"
OUTPUT="${STORAGE}/bdev_integrated/reports/${RUN_ID}"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
R_LIB="${GSEA_R_LIB:-/tmp/scworkbench_bdev_gsea_rlib}"
CACHE="${GSEA_CACHE_DIR:-/tmp/scworkbench_bdev_gsea_cache}"
SEED="${GSEA_SEED:-1729}"

[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ && "${SEED}" =~ ^[0-9]+$ ]] || exit 2
[[ -d "${INPUT}" && ! -e "${OUTPUT}" && ! -e "${OUTPUT}.incomplete" ]] || exit 1
mkdir -p "${R_LIB}" "${CACHE}"

"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${R_LIB}" --env "XDG_CACHE_HOME=${CACHE}" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/processing/integration/76_mapped_enrichment_dynamic_change_gsea.R" \
  "${INPUT}" "${OUTPUT}" "${SEED}" \
  2>&1 | tee "${STORAGE}/bdev_integrated/reports/${RUN_ID}.log"
