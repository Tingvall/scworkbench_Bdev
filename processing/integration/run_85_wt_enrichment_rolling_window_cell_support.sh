#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
source_report="${storage}/bdev_integrated/reports/bdev_cutoff025_mapping_report_v4_20261012"
output_dir="${1:-${storage}/bdev_integrated/reports/wt_enrichment_rolling_3bin_cell_support_v2_20261013}"
runtime="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"

[[ -x "${runtime}" ]] || runtime="$(command -v apptainer || command -v singularity)"
[[ -x "${runtime}" ]] || { echo "No Apptainer/Singularity runtime found" >&2; exit 1; }
[[ -f "${sif}" ]] || { echo "Missing image: ${sif}" >&2; exit 1; }
[[ -d "${source_report}" ]] || { echo "Missing source report: ${source_report}" >&2; exit 1; }
[[ ! -e "${output_dir}" ]] || { echo "Output already exists: ${output_dir}" >&2; exit 1; }
case "${output_dir}" in
  "${storage}/bdev_integrated/reports/"*) ;;
  *) echo "Output must be below ${storage}/bdev_integrated/reports" >&2; exit 1 ;;
esac

"${runtime}" exec --cleanenv --no-home \
  -B "${repo_root}:/repo:ro" -B "${storage}:${storage}:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 \
  --env MKL_NUM_THREADS=4 --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/85_wt_enrichment_rolling_window_cell_support.R \
  "${source_report}" "${output_dir}" 10 3 3

printf 'Rolling-window support report: %s\n' "${output_dir}/REPORT.md"
