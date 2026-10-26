#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
image="${SCWORKBENCH_IMAGE:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif}"
root="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/gse259355_tcf3_prob"
counts="${root}/prepared/counts/gse259355_tcf3_prob_featurecounts.tsv"
peaks="${root}/prepared/counts/gse259355_tcf3_prob_consensus_peaks.mm10.bed"
manifest="${repo_root}/processing/gse259355_tcf3_prob/gse259355_tcf3_prob_samples.tsv"
bundle="${root}/bundle"
report="${root}/reports/source_bundle_v1"

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
else
  container_runtime="singularity"
fi

"${container_runtime}" exec --cleanenv \
  -B /mnt/sata2:/mnt/sata2 \
  -B "${repo_root}:/repo" \
  "${image}" \
  Rscript --vanilla /repo/processing/gse259355_tcf3_prob/03_create_source_bundle.R \
  "${counts}" "${peaks}" "${manifest}" "${bundle}" "${report}"

printf 'Bundle: %s\n' "${bundle}"
printf 'Validation: %s\n' "${report}/validation.json"
