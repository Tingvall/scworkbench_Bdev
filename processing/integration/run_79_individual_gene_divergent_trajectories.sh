#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
INPUT="${STORAGE}/bdev_integrated/reports/flt3l_il7_ko_vs_wt_dynamic_change_gsea_v1_20261012"
SESSION_STAGE="${STORAGE}/bdev_integrated/sessions/ko_wt_enriched_pt_peak_de_v1/stage"
PROJECTION="${STORAGE}/bdev_integrated/reports/ko_mapping_reference_umap_v1_20261001/query_reference_umap_projection.tsv.gz"
RUN_ID="${1:-flt3l_il7_ko_vs_wt_individual_gene_divergence_v1_20261012}"
OUTPUT="${STORAGE}/bdev_integrated/reports/${RUN_ID}"
TOP_REPORT_N="${TOP_GENE_REPORT_N:-20}"
TOP_CURVE_N="${TOP_GENE_CURVE_N:-12}"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
R_LIB="${GSEA_R_LIB:-/tmp/scworkbench_bdev_gsea_rlib}"

[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ && "${TOP_REPORT_N}" =~ ^[1-9][0-9]*$ \
  && "${TOP_CURVE_N}" =~ ^[1-9][0-9]*$ ]] || exit 2
(( TOP_CURVE_N <= TOP_REPORT_N )) || exit 2
[[ -d "${INPUT}" && -d "${SESSION_STAGE}" && -f "${PROJECTION}" && \
  ! -e "${OUTPUT}" && ! -e "${OUTPUT}.incomplete" ]] || exit 1

"${RUNTIME}" exec --cleanenv --no-home --env "R_LIBS_USER=${R_LIB}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/processing/integration/79_individual_gene_divergent_trajectories.R" \
  "${INPUT}" "${SESSION_STAGE}" "${PROJECTION}" "${OUTPUT}" \
  "${TOP_REPORT_N}" "${TOP_CURVE_N}"
