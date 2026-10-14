#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
image="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
source_report="${storage}/bdev_integrated/reports/tcf3_ebf1_coverage_placseq_tracks_v4_20261005"
matched_report="${storage}/bdev_integrated/reports/ebf1_pm2mb_depth_complexity_matched_da_v1_20261013"
output_dir="${1:-${storage}/bdev_integrated/reports/ebf1_locus_matched_da_e2a_placseq_v1_20261013}"

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
else
  container_runtime="singularity"
fi

"${container_runtime}" exec --cleanenv \
  -B /mnt/sata2:/mnt/sata2 \
  -B "${repo_root}:/repo" \
  "${image}" Rscript --vanilla \
  /repo/processing/integration/94_ebf1_locus_matched_da_e2a_placseq.R \
  "${source_report}" "${matched_report}" "${output_dir}" \
  "94_ebf1_locus_matched_da_e2a_placseq.R:v1"

printf 'Matched Ebf1 locus figure: %s/REPORT.md\n' "${output_dir}"
