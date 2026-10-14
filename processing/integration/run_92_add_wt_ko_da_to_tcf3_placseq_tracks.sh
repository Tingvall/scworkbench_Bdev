#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
image="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
source_report="${storage}/bdev_integrated/reports/tcf3_ebf1_coverage_placseq_tracks_v4_20261005"
da_report="${storage}/bdev_integrated/reports/ebf1_pm2mb_shared_support_da_v1_20261013"
output_dir="${1:-${storage}/bdev_integrated/reports/e2a_tcf3_placseq_with_shared_support_wt_ko_da_v2_20261013}"

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
else
  container_runtime="singularity"
fi

"${container_runtime}" exec --cleanenv \
  -B /mnt/sata2:/mnt/sata2 \
  -B "${repo_root}:/repo" \
  "${image}" Rscript --vanilla \
  /repo/processing/integration/92_add_wt_ko_da_to_tcf3_placseq_tracks.R \
  "${source_report}" \
  "${da_report}/all_targets_ebf1_pm2mb_shared_min10_wt_vs_ko_da.tsv" \
  "${output_dir}" \
  "92_add_wt_ko_da_to_tcf3_placseq_tracks.R:v2"

printf 'E2A/TCF3 + WT/KO DA + PLAC-seq report: %s/REPORT.md\n' "${output_dir}"
