#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
local_report="${storage}/bdev_integrated/reports/il7_responsive_peak_local_deg_v1_20261002"
output="${storage}/bdev_integrated/reports/il7_responsive_peak_linked_gene_volcano_v1_20261002"

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
  -B "${storage}:${storage}:rw"
  --pwd /repo
  "${sif}"
  Rscript --vanilla
)
local_sha="$(${container_runtime} "${common[@]}" -e \
  'cat(digest::digest(commandArgs(TRUE)[1], algo="sha256", file=TRUE))' \
  "${local_report}/completion.json")"
"${container_runtime}" "${common[@]}" \
  /repo/processing/integration/27_plot_il7_peak_linked_gene_volcano.R \
  "${local_report}" "${output}" "${local_sha}" \
  27_plot_il7_peak_linked_gene_volcano.R:v3_top50 50
printf 'Volcano report: %s\n' "${output}/REPORT.md"
