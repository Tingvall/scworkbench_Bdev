#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BDEV_BUNDLE="${STORAGE}/bdev_reference/bundle"
QUERY_BUNDLE="${STORAGE}/scmulti_flt3ko_il7ko/bundle"
REGION_REPORT="${STORAGE}/bdev_integrated/reports/ko_wt_enriched_peak_rna_de_v1_20261001"
QUERY_SESSION="${STORAGE}/bdev_integrated/sessions/ko_wt_enriched_pt_peak_de_v1"
BDEV_SESSION_ROOT="${STORAGE}/bdev_integrated/sessions"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
OUTPUT_DIR="${1:-${STORAGE}/bdev_integrated/reports/bdev_variable_gene_expression_pt_heatmaps_$(date -u +%Y%m%dT%H%M%SZ)}"
TOP_N="${GENE_EXPRESSION_TOP_N:-30}"
MINIMUM_CELLS="${GENE_EXPRESSION_MINIMUM_CELLS:-20}"
SEED="${GENE_EXPRESSION_SEED:-1729}"
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
  /repo/processing/integration/44_bdev_variable_gene_expression_pt_heatmaps.R \
  "${BDEV_BUNDLE}" "${QUERY_BUNDLE}" "${REGION_REPORT}" \
  "${QUERY_SESSION}" "${BDEV_SESSION_ROOT}" "${OUTPUT_DIR}" \
  "${BDEV_FP}" "${QUERY_FP}" "${TOP_N}" "${MINIMUM_CELLS}" "${SEED}"
