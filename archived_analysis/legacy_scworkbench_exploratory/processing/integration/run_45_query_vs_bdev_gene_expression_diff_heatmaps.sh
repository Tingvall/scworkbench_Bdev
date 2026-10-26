#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
SOURCE_REPORT="${STORAGE}/bdev_integrated/reports/bdev_variable_gene_expression_pt_heatmaps_v1_20261009_run4"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
OUTPUT_DIR="${1:-${STORAGE}/bdev_integrated/reports/query_vs_bdev_gene_expression_diff_heatmaps_$(date -u +%Y%m%dT%H%M%SZ)}"
LOW_CELL_THRESHOLD="${GENE_EXPRESSION_LOW_CELL_THRESHOLD:-20}"
DISPLAY_LIMIT="${GENE_EXPRESSION_DIFF_DISPLAY_LIMIT:-3}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
[[ -x "${RUNTIME}" ]] || { echo "Apptainer missing: ${RUNTIME}" >&2; exit 1; }
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output exists: ${OUTPUT_DIR}" >&2; exit 1; }

exec "${RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:/repo:ro" \
  --bind "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw" \
  --pwd /repo "${SIF}" Rscript --vanilla \
  /repo/processing/integration/45_query_vs_bdev_gene_expression_diff_heatmaps.R \
  "${SOURCE_REPORT}" "${OUTPUT_DIR}" "${LOW_CELL_THRESHOLD}" "${DISPLAY_LIMIT}"
