#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
REPORT_DIR="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/bdev_integrated/reports/go_ucell_rna_feature_residual_dynamic_range_perm1000_v1_20261012/summary_figures/significant_region_residual_ucell_kmeans_k25_20_20/figure_report"
SOURCE_HTML="${REPORT_DIR}/FIGURE_REPORT.html"
OUTPUT_PDF="${REPORT_DIR}/FIGURE_REPORT.pdf"
OUTPUT_STAGE="${OUTPUT_PDF}.incomplete"
CHROME="/usr/bin/google-chrome"

[[ -f "${SOURCE_HTML}" && -x "${CHROME}" ]] || exit 1
[[ ! -e "${OUTPUT_STAGE}" ]] || exit 1

PRINT_HTML="$(mktemp --suffix=.html)"
cleanup() { rm -f -- "${PRINT_HTML}"; }
trap cleanup EXIT

Rscript --vanilla \
  "${REPO_ROOT}/figure_scripts/bdev_integrated/38_prepare_cluster_figure_report_print_html.R" \
  "${SOURCE_HTML}" "${PRINT_HTML}"

"${CHROME}" --headless --no-sandbox --disable-gpu \
  --allow-file-access-from-files --run-all-compositor-stages-before-draw \
  --virtual-time-budget=5000 --no-pdf-header-footer \
  --print-to-pdf="${OUTPUT_STAGE}" "file://${PRINT_HTML}"

[[ -s "${OUTPUT_STAGE}" ]] || exit 1
mv -f -- "${OUTPUT_STAGE}" "${OUTPUT_PDF}"
