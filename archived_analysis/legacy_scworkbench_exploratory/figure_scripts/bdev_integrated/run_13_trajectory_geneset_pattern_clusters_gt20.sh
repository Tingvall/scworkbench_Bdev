#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
REPORT_ROOT="${1:-${STORAGE}/bdev_integrated/reports/trajectory_go_depth_adjusted_adaptive_v1_20261011}"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
BB_CHANGE_POINTS="${STORAGE}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925/bb_up_change_points.tsv"

[[ -f "${REPORT_ROOT}/completion.json" && -f "${BB_CHANGE_POINTS}" ]] || exit 1
"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${REGISTERED_LIBRARY}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}/bdev_reference:${STORAGE}/bdev_reference:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/figure_scripts/bdev_integrated/13_trajectory_geneset_pattern_clusters_gt20.R" \
  "${REPORT_ROOT}" "${BB_CHANGE_POINTS}"
