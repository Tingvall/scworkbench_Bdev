#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev"
dataset_id="gse259355_tcf3_prob"
raw_root="${storage}/datasets/${dataset_id}"
analysis_root="${storage}/intermediate_objects/${dataset_id}"
work_root="${analysis_root}/work"
prepared_root="${analysis_root}/prepared"
report_root="${analysis_root}/reports/tcf3_depleted_vs_control_da_v1"
runtime="${analysis_root}/runtime"
manifest="${repo_root}/processing/${dataset_id}/gse259355_tcf3_prob_samples.tsv"
threads="${THREADS:-12}"

mkdir -p "${raw_root}" "${work_root}" "${prepared_root}" "${analysis_root}/reports"

"${repo_root}/processing/${dataset_id}/01_prepare_atac_counts.sh" \
  "${manifest}" "${raw_root}" "${work_root}" "${prepared_root}" "${runtime}" "${threads}"

if [[ -n "${ONLY_SAMPLE:-}" ]]; then
  echo "Completed requested sample ${ONLY_SAMPLE}; skipping cohort-level DESeq2"
  exit 0
fi

if [[ ! -e "${report_root}" ]]; then
  "${runtime}/bin/Rscript" --vanilla \
    "${repo_root}/processing/${dataset_id}/02_deseq2_differential_accessibility.R" \
    "${prepared_root}/counts/gse259355_tcf3_prob_featurecounts.tsv" \
    "${manifest}" "${report_root}" \
    "02_deseq2_differential_accessibility.R:v1"
fi

echo "Report: ${report_root}/REPORT.md"
