#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
report_dir="${storage}/bdev_integrated/reports/ko_shared_enriched_peak_rna_de_v1_20261002"

for path in "${sif}" "${library}/scworkbench" \
  "${report_dir}/completion.json"; do
  [[ -e "${path}" ]] || { echo "Required input missing: ${path}" >&2; exit 1; }
done

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
elif [[ -x /home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer ]]; then
  container_runtime="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
elif command -v singularity >/dev/null 2>&1; then
  container_runtime="singularity"
else
  echo "No Apptainer or Singularity runtime found" >&2
  exit 1
fi

"${container_runtime}" exec --cleanenv --no-home \
  -B "${library}:${library}:ro" --env "R_LIBS_USER=${library}" \
  -B "${repo_root}:/repo:ro" -B "${storage}:${storage}:rw" \
  --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/21_add_shared_region_flt3_il7r_plots.R \
  "${report_dir}" "v5_flt3_il7r_target_plots_20261002" execute
