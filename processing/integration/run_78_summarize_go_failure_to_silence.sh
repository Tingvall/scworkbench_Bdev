#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
INPUT="${STORAGE}/bdev_integrated/reports/flt3l_il7_ko_vs_wt_dynamic_change_gsea_v1_20261012"
RUN_ID="${1:-flt3l_il7_ko_vs_wt_go_failure_to_silence_v1_20261012}"
OUTPUT="${STORAGE}/bdev_integrated/reports/${RUN_ID}"
TOP_N="${TOP_GO_N:-5}"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
R_LIB="${GSEA_R_LIB:-/tmp/scworkbench_bdev_gsea_rlib}"

[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ && "${TOP_N}" =~ ^[1-9][0-9]*$ ]] || exit 2
[[ -d "${INPUT}" && ! -e "${OUTPUT}" && ! -e "${OUTPUT}.incomplete" ]] || exit 1

"${RUNTIME}" exec --cleanenv --no-home --env "R_LIBS_USER=${R_LIB}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/processing/integration/78_summarize_go_failure_to_silence.R" \
  "${INPUT}" "${OUTPUT}" "${TOP_N}"
