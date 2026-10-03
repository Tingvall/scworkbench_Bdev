#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
runtime="${SCWORKBENCH_RUNTIME:-/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
local_report="${storage}/bdev_integrated/reports/flt3l_responsive_peak_local_deg_v1_20261002"
output="${storage}/bdev_integrated/reports/flt3l_responsive_peak_linked_gene_volcano_v1_20261002"

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
local_sha="$(${runtime} "${common[@]}" -e \
  'cat(digest::digest(commandArgs(TRUE)[1], algo="sha256", file=TRUE))' \
  "${local_report}/completion.json")"
"${runtime}" "${common[@]}" \
  /repo/processing/integration/33_plot_flt3l_peak_linked_gene_volcano.R \
  "${local_report}" "${output}" "${local_sha}" \
  33_plot_flt3l_peak_linked_gene_volcano.R:v1_top50 50
printf 'FLT3L target-gene volcano report: %s\n' "${output}/REPORT.md"
