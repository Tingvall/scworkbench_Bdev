#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
INTEGRATED_ROOT="${STORAGE}/bdev_integrated"
SOURCE_REPORT="${INTEGRATED_ROOT}/reports/bdev_selected_motif_pt_heatmaps_v1_20261010_run5"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
OUTPUT_DIR="${1:-${INTEGRATED_ROOT}/reports/bdev_up_down_motif_trajectory_adjusted_wt_ko_$(date -u +%Y%m%dT%H%M%SZ)}"
MINIMUM_BDEV_CHANGE="${MOTIF_MINIMUM_BDEV_CHANGE_SD:-1.5}"
MINIMUM_ABSOLUTE_SPEARMAN="${MOTIF_MINIMUM_ABSOLUTE_SPEARMAN:-0.6}"
MINIMUM_EFFECT="${MOTIF_ADJUSTED_MINIMUM_EFFECT:-0.5}"
MINIMUM_SIGN_CONSISTENCY="${MOTIF_ADJUSTED_MINIMUM_SIGN_CONSISTENCY:-0.7}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
[[ -x "${RUNTIME}" ]] || { echo "Apptainer missing: ${RUNTIME}" >&2; exit 1; }
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output exists: ${OUTPUT_DIR}" >&2; exit 1; }

exec "${RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:/repo:ro" \
  --bind "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw" \
  --pwd /repo "${SIF}" Rscript --vanilla \
  /repo/processing/integration/47_bdev_up_down_motif_trajectory_adjusted_wt_ko.R \
  "${SOURCE_REPORT}" "${OUTPUT_DIR}" "${MINIMUM_BDEV_CHANGE}" \
  "${MINIMUM_ABSOLUTE_SPEARMAN}" "${MINIMUM_EFFECT}" \
  "${MINIMUM_SIGN_CONSISTENCY}"
