#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
DYNAMIC="${STORAGE}/bdev_integrated/reports/flt3l_il7_ko_vs_wt_dynamic_change_gsea_v1_20261012"
TOP_GO="${STORAGE}/bdev_integrated/reports/flt3l_il7_ko_vs_wt_dynamic_change_top_go_v1_20261012"
RUN_ID="${1:-flt3l_il7_ko_vs_wt_fdr_supported_go_curves_v1_20261012}"
OUTPUT="${STORAGE}/bdev_integrated/reports/${RUN_ID}"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
R_LIB="${GSEA_R_LIB:-/tmp/scworkbench_bdev_gsea_rlib}"

[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ ]] || exit 2
[[ -d "${DYNAMIC}" && -d "${TOP_GO}" && ! -e "${OUTPUT}" && \
  ! -e "${OUTPUT}.incomplete" ]] || exit 1

"${RUNTIME}" exec --cleanenv --no-home --env "R_LIBS_USER=${R_LIB}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/figure_scripts/bdev_integrated/40_fdr_supported_go_dynamic_curves.R" \
  "${DYNAMIC}" "${TOP_GO}" "${OUTPUT}"
