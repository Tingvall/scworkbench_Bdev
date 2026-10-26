#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
image="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
bdev_bundle="${storage}/bdev_reference/bundle"
global_report="${storage}/bdev_integrated/reports/il7_responsive_peak_global_v1_20261002"
ebf1_da_report="${storage}/bdev_integrated/reports/ebf1_enhancer_il7_enriched_pt_da_v1_20261001"
tcf3_report="${storage}/gse259355_tcf3_prob/reports/tcf3_depleted_vs_control_da_v1"
placseq="${storage}/bdev_integrated/prepared/gse162858_placseq/GSE162858_230238_H3K4me3_FitHiChIP.interactions_FitHiC_Q0.05_MergeNearContacts.bed.txt.gz"
session_root="${storage}/bdev_integrated/sessions"
output_dir="${1:-${storage}/bdev_integrated/reports/ebf1_putative_enhancers_three_way_v2_20261005}"
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
  /repo/processing/integration/34_ebf1_putative_enhancers_three_way.R \
  "${bdev_bundle}" \
  "${global_report}/bdev_il7_peak2_cells.tsv" \
  "${ebf1_da_report}/ebf1_pm3mb_bdev_peak_universe.tsv" \
  "${tcf3_report}/tcf3_depleted_vs_control_differential_accessibility.tsv" \
  "${placseq}" "${session_root}" "${output_dir}" "${expected_bdev_fp}" \
  "34_ebf1_putative_enhancers_three_way.R:v2" 44617317 314159

printf 'Report: %s/REPORT.html\n' "${output_dir}"
