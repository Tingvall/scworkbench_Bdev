#!/usr/bin/env bash
set -euo pipefail

stage="${1:?Usage: run_ko_mapping.sh create|smoke}"
if [[ "${stage}" != create && "${stage}" != smoke ]]; then
  echo "Usage: $0 create|smoke" >&2
  exit 2
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
bdev_bundle="${storage}/bdev_reference/bundle"
ko_root="${storage}/scmulti_flt3ko_il7ko"
ko_bundle="${ko_root}/bundle"
mapping="${ko_root}/reports/bdev_peak_liftover_v1_20260930/accepted_bdev_peak_mapping.tsv.gz"
ko_matrix="${ko_root}/prepared/bdev_peak_quantification_v1_20260930/ATAC_bdev_mm10_peaks_ko_counts"
phase_root="${storage}/bdev_integrated/ko_mapping_v1"
bundle_dir="${phase_root}/bundle"
create_report="${storage}/bdev_integrated/reports/ko_mapping_bundle_v1_20260930"
smoke_report="${storage}/bdev_integrated/reports/ko_mapping_peakvi_smoke_v1_20260930"
session_root="${phase_root}/sessions"

[[ -d "${library}/scworkbench" ]] || { echo "Registered scWorkbench runtime missing: ${library}" >&2; exit 1; }

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
else
  container_runtime="singularity"
fi

mkdir -p "${phase_root}" "${storage}/bdev_integrated/reports"
ulimit -v $((64 * 1024 * 1024))

common=(
  exec --cleanenv --no-home
  -B "${library}:${library}:ro"
  --env "R_LIBS_USER=${library}"
  --env "SCWORKBENCH_PYBENCH_PYTHON=/opt/pyenv/bin/python"
  -B "${repo_root}:/repo:ro"
  -B "${storage}:${storage}:rw"
  --env OMP_NUM_THREADS=4
  --env OPENBLAS_NUM_THREADS=4
  --env MKL_NUM_THREADS=4
  --pwd /repo
  "${sif}"
  Rscript --vanilla
)

if [[ "${stage}" == create ]]; then
  "${container_runtime}" "${common[@]}" \
    /repo/processing/integration/07_create_ko_mapping_bundle.R \
    "${bdev_bundle}" "${ko_bundle}" "${mapping}" "${ko_matrix}" \
    "${bundle_dir}" "${create_report}" jtg
else
  "${container_runtime}" "${common[@]}" \
    /repo/processing/integration/08_peakvi_mapping_smoke.R \
    "${bundle_dir}" "${smoke_report}" "${session_root}"
fi

printf 'Mapping bundle: %s\n' "${bundle_dir}"
printf 'Creation report: %s\n' "${create_report}/creation.json"
printf 'Smoke report: %s\n' "${smoke_report}/smoke.json"

