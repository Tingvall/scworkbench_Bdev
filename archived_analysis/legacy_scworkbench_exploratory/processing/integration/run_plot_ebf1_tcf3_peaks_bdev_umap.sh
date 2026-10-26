#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
image="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
bdev_bundle="${storage}/bdev_reference/bundle"
session_workspace="${storage}/bdev_integrated/sessions/ebf1_putative_enhancers_three_way_v1"
three_way_report="${storage}/bdev_integrated/reports/ebf1_putative_enhancers_three_way_v2_20261005"
output_dir="${1:-${storage}/bdev_integrated/reports/ebf1_tcf3_peaks_bdev_umap_v3_20261005}"
expected_bdev_fp="1027f6390c6b853acc2bf40d3fae9abb59393f4256cb39285685eb2b2973f7ec"

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
else
  container_runtime="singularity"
fi

"${container_runtime}" exec --cleanenv \
  -B /mnt/sata2:/mnt/sata2 \
  -B "${repo_root}:/repo" \
  "${image}" Rscript --vanilla \
  /repo/processing/integration/35_plot_ebf1_tcf3_peaks_bdev_umap.R \
  "${bdev_bundle}" "${session_workspace}" \
  "${three_way_report}/all_259_bdev_ebf1_window_peaks_three_way_tests.tsv" \
  "${output_dir}" "${expected_bdev_fp}" \
  "35_plot_ebf1_tcf3_peaks_bdev_umap.R:v3" 314159

printf 'Bdev UMAP report: %s/REPORT.html\n' "${output_dir}"
