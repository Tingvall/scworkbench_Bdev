#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BDEV_BUNDLE="${STORAGE}/bdev_reference/bundle"
QUERY_BUNDLE="${STORAGE}/scmulti_flt3ko_il7ko/bundle"
REGION_REPORT="${STORAGE}/bdev_integrated/reports/ko_wt_enriched_peak_rna_de_v1_20261001"
QUERY_SESSION="${STORAGE}/bdev_integrated/sessions/ko_wt_enriched_pt_peak_de_v1"
BDEV_SESSION="${STORAGE}/bdev_integrated/sessions/bdev_variable_gene_expression_norm_v1"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
OUTPUT_DIR="${1:-${STORAGE}/bdev_integrated/reports/ko_specific_gene_expression_diff_heatmaps_$(date -u +%Y%m%dT%H%M%SZ)}"
TOP_N="${KO_SPECIFIC_TOP_N:-30}"
LOW_CELL_THRESHOLD="${GENE_EXPRESSION_LOW_CELL_THRESHOLD:-20}"
DISPLAY_LIMIT="${KO_SPECIFIC_DISPLAY_LIMIT:-2}"
BDEV_VARIABILITY_QUANTILE="${BDEV_VARIABILITY_QUANTILE:-0.75}"
BDEV_FP="$(python3 - "${REPO_ROOT}/config/bdev_reference.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle)['validation']['bundle_state_fingerprint'])
PY
)"
QUERY_FP="$(python3 - "${REPO_ROOT}/config/scmulti_flt3ko_il7ko.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle)['validation']['bundle_state_fingerprint'])
PY
)"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
[[ -x "${RUNTIME}" ]] || { echo "Apptainer missing: ${RUNTIME}" >&2; exit 1; }
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output exists: ${OUTPUT_DIR}" >&2; exit 1; }

exec "${RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:/repo:ro" \
  --bind "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw" \
  --pwd /repo "${SIF}" Rscript --vanilla \
  /repo/processing/integration/46_ko_specific_gene_expression_diff_heatmaps.R \
  "${BDEV_BUNDLE}" "${QUERY_BUNDLE}" "${REGION_REPORT}" \
  "${QUERY_SESSION}" "${BDEV_SESSION}" "${OUTPUT_DIR}" \
  "${BDEV_FP}" "${QUERY_FP}" "${TOP_N}" "${LOW_CELL_THRESHOLD}" \
  "${DISPLAY_LIMIT}" "${BDEV_VARIABILITY_QUANTILE}"
