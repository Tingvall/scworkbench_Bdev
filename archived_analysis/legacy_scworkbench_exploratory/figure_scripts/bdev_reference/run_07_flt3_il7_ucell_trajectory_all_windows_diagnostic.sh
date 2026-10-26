#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
RUN_ID="${1:-flt3_il7_refined_ucell_trajectory_all_windows_20260925}"
[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ ]] || exit 2

STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
SIF="${STORAGE}/scworkbench-bdev-signac.sif"
SOURCE_DIR="${STORAGE}/bdev_reference/reports/flt3_il7_refined_ucell_trajectory_5pct_v4_20260925"
OUTPUT_DIR="${STORAGE}/bdev_reference/reports/${RUN_ID}"

[[ -f "${SIF}" && -f "${SOURCE_DIR}/sliding_window_summary.tsv" ]]
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output already exists: ${OUTPUT_DIR}" >&2; exit 1; }
mkdir -p "${OUTPUT_DIR}"

exec singularity exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:/repo:ro" \
  --bind "${SOURCE_DIR}:/source:ro" \
  --bind "${OUTPUT_DIR}:/output:rw" \
  --pwd /repo "${SIF}" \
  Rscript --vanilla \
  "/repo/figure_scripts/bdev_reference/07_flt3_il7_ucell_trajectory_all_windows_diagnostic.R" \
  /source /output
