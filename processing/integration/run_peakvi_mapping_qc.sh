#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
phase_root="${storage}/bdev_integrated/ko_mapping_v1"
bundle_dir="${phase_root}/bundle"
stage_dir="${phase_root}/sessions/peakvi_full_mapping_v1/stage"
mapping_json="${storage}/bdev_integrated/reports/ko_mapping_peakvi_full_v1_20260930/mapping.json"
report_dir="${storage}/bdev_integrated/reports/ko_mapping_qc_v1_20260930"

if command -v apptainer >/dev/null 2>&1; then container_runtime="apptainer"; else container_runtime="singularity"; fi
[[ -d "${library}/scworkbench" ]] || { echo "Registered scWorkbench runtime missing" >&2; exit 1; }
[[ -f "${mapping_json}" ]] || { echo "Completed full mapping report missing" >&2; exit 1; }
[[ -d "${stage_dir}" ]] || { echo "Completed full mapping stage bundle missing" >&2; exit 1; }
[[ ! -e "${report_dir}" ]] || { echo "QC report already exists: ${report_dir}" >&2; exit 1; }
mkdir -p "${storage}/bdev_integrated/reports"
ulimit -v $((64 * 1024 * 1024))

"${container_runtime}" exec --cleanenv --no-home \
  -B "${library}:${library}:ro" --env "R_LIBS_USER=${library}" \
  -B "${repo_root}:/repo:ro" -B "${storage}:${storage}:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4 \
  --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/11_peakvi_mapping_qc.R \
  "${bundle_dir}" "${stage_dir}" "${mapping_json}" "${report_dir}"

printf 'Mapping QC report: %s\n' "${report_dir}/REPORT.md"
