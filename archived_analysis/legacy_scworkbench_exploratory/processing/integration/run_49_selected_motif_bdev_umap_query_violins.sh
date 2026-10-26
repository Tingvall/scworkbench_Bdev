#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BDEV_ROOT="${STORAGE}/bdev_reference"
INTEGRATED_ROOT="${STORAGE}/bdev_integrated"
BDEV_BUNDLE="${BDEV_ROOT}/bundle"
MAPPING_BUNDLE="${INTEGRATED_ROOT}/ko_mapping_v1/bundle"
JOINT_MATRIX="${MAPPING_BUNDLE}/mat/mat_frag_ATAC_bdev_mm10_peaks_20260930_122757"
REGION_REPORT="${INTEGRATED_ROOT}/reports/ko_wt_enriched_peak_rna_de_v1_20261001"
BDEV_RDATA="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/datasets/scmulti_bdev/Bdev.Rdata"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
OUTPUT_DIR="${1:-${INTEGRATED_ROOT}/reports/selected_motif_bdev_umap_query_violins_$(date -u +%Y%m%dT%H%M%SZ)}"
SEED="${MOTIF_SEED:-1729}"
BDEV_FP="$(python3 - "${REPO_ROOT}/config/bdev_reference.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle)['validation']['bundle_state_fingerprint'])
PY
)"
MAPPING_FP="$(python3 - "${MAPPING_BUNDLE}/manifest.json" <<'PY'
import json, sys
with open(sys.argv[1]) as handle:
    print(json.load(handle)['current_state']['fingerprint'])
PY
)"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
[[ -x "${RUNTIME}" ]] || { echo "Apptainer missing: ${RUNTIME}" >&2; exit 1; }
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output exists: ${OUTPUT_DIR}" >&2; exit 1; }

exec "${RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:/repo:ro" \
  --bind "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw" \
  --pwd /repo "${SIF}" Rscript --vanilla \
  /repo/processing/integration/49_selected_motif_bdev_umap_query_violins.R \
  "${BDEV_BUNDLE}" "${MAPPING_BUNDLE}" "${JOINT_MATRIX}" \
  "${REGION_REPORT}" "${BDEV_RDATA}" "${OUTPUT_DIR}" \
  "${BDEV_FP}" "${MAPPING_FP}" "${SEED}" \
  "snai_tcf3_runx1_srebf2_stat5"
