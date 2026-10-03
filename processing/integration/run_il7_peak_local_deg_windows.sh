#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
peak_report="${storage}/bdev_integrated/reports/il7_responsive_peak_annotation_v1_20261002"
deg_report="${storage}/bdev_integrated/reports/ko_shared_enriched_peak_rna_de_v1_20261002"
bdev_rdata="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/datasets/scmulti_bdev/Bdev.Rdata"
output="${storage}/bdev_integrated/reports/il7_responsive_peak_local_deg_v1_20261002"

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
else
  container_runtime="singularity"
fi
common=(
  exec --cleanenv --no-home
  -B "${library}:${library}:ro"
  --env "R_LIBS_USER=${library}"
  -B "${repo_root}:/repo:ro"
  -B "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw"
  --pwd /repo
  "${sif}"
  Rscript --vanilla
)
peak_sha="$(${container_runtime} "${common[@]}" -e \
  'cat(digest::digest(commandArgs(TRUE)[1], algo="sha256", file=TRUE))' \
  "${peak_report}/completion.json")"
deg_sha="$(${container_runtime} "${common[@]}" -e \
  'cat(digest::digest(commandArgs(TRUE)[1], algo="sha256", file=TRUE))' \
  "${deg_report}/completion.json")"
"${container_runtime}" "${common[@]}" \
  /repo/processing/integration/26_il7_peak_local_deg_windows.R \
  "${peak_report}" "${deg_report}" "${bdev_rdata}" "${output}" \
  "${peak_sha}" "${deg_sha}" 26_il7_peak_local_deg_windows.R:v1 \
  3000000 0.05
printf 'Local DEG report: %s\n' "${output}/REPORT.md"
