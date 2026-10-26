#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
INTEGRATED_ROOT="${STORAGE}/bdev_integrated"
MOTIF_REPORT="${INTEGRATED_ROOT}/reports/bdev_selected_motif_pt_heatmaps_v1_20261010_run5"
SCORE_REPORT="${INTEGRATED_ROOT}/reports/bdev_strong_directional_motif_trajectory_adjusted_wt_ko_v1_20261010"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
OUTPUT_DIR="${1:-${INTEGRATED_ROOT}/reports/top_adjusted_wt_ko_motif_heatmaps_$(date -u +%Y%m%dT%H%M%SZ)}"
TOP_N="${TOP_ADJUSTED_MOTIFS_PER_DIRECTION:-20}"
ACTIVITY_LIMIT="${MOTIF_ACTIVITY_DISPLAY_LIMIT:-2.5}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
[[ -x "${RUNTIME}" ]] || { echo "Apptainer missing: ${RUNTIME}" >&2; exit 1; }
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output exists: ${OUTPUT_DIR}" >&2; exit 1; }

exec "${RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:/repo:ro" \
  --bind "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw" \
  --pwd /repo "${SIF}" Rscript --vanilla \
  /repo/processing/integration/48_top_adjusted_wt_ko_motif_heatmaps.R \
  "${MOTIF_REPORT}" "${SCORE_REPORT}" "${OUTPUT_DIR}" \
  "${TOP_N}" "${ACTIVITY_LIMIT}"
