#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
global_report="${storage}/bdev_integrated/reports/il7_responsive_peak_global_v1_20261002"
bdev_rdata="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/datasets/scmulti_bdev/Bdev.Rdata"
output="${storage}/bdev_integrated/reports/il7_responsive_peak_annotation_v1_20261002"

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
completion_sha="$(${container_runtime} "${common[@]}" -e \
  'cat(digest::digest(commandArgs(TRUE)[1], algo="sha256", file=TRUE))' \
  "${global_report}/completion.json")"
"${container_runtime}" "${common[@]}" \
  /repo/processing/integration/25_annotate_il7_validated_peaks.R \
  "${global_report}" "${bdev_rdata}" "${output}" "${completion_sha}" \
  25_annotate_il7_validated_peaks.R:v1 2000
printf 'Annotation report: %s\n' "${output}/REPORT.md"
