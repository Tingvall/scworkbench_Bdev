#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
image="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
tcf_root="${storage}/gse259355_tcf3_prob"
three_way_report="${storage}/bdev_integrated/reports/ebf1_putative_enhancers_three_way_v2_20261005"
output_dir="${1:-${storage}/bdev_integrated/reports/tcf3_ebf1_coverage_placseq_tracks_v4_20261005}"

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
else
  container_runtime="singularity"
fi

"${container_runtime}" exec --cleanenv \
  -B /mnt/sata2:/mnt/sata2 \
  -B "${repo_root}:/repo" \
  "${image}" Rscript --vanilla \
  /repo/processing/integration/36_tcf3_ebf1_coverage_placseq_tracks.R \
  "${tcf_root}/prepared/sample_manifest.tsv" \
  "${tcf_root}/prepared/bam" \
  "${tcf_root}/prepared/qc" \
  "${three_way_report}" "${output_dir}" \
  "36_tcf3_ebf1_coverage_placseq_tracks.R:v4.1-inplace" 44617317 314159

printf 'TCF3 Ebf1 coverage report: %s/REPORT.html\n' "${output_dir}"
