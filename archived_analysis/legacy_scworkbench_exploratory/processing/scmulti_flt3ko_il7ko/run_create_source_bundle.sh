#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
sif="${SCWORKBENCH_IMAGE:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif}"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scmulti_flt3ko_il7ko"
outs="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/datasets/scmulti_flt3ko_il7ko/FLT3LKO_IL7KO_BM/outs"
qc_selection="${storage}/reports/qc_selection_v1_20260929/cell_qc_selection.tsv.gz"
inventory="${storage}/reports/unfiltered_qc_v1_20260929/input_inventory.tsv"
bundle_dir="${storage}/bundle"
report_dir="${storage}/reports/source_bundle_v1_20260929"

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
else
  container_runtime="singularity"
fi

"${container_runtime}" exec --cleanenv \
  -B /mnt/sata2:/mnt/sata2 \
  -B /mnt/synology:/mnt/synology \
  -B /home/jtg/resources:/home/jtg/resources \
  -B "${repo_root}:/repo" \
  "${sif}" \
  Rscript --vanilla /repo/processing/scmulti_flt3ko_il7ko/05_create_source_bundle.R \
  "${outs}" "${qc_selection}" "${bundle_dir}" "${report_dir}" "${inventory}"

printf 'Bundle: %s\n' "${bundle_dir}"
printf 'Creation report: %s\n' "${report_dir}/completion.json"
