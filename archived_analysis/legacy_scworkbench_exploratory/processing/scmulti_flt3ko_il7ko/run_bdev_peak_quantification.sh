#!/usr/bin/env bash
set -euo pipefail

stage="${1:-all}"
if [[ "${stage}" != audit && "${stage}" != quantify && "${stage}" != validate && "${stage}" != all ]]; then
  echo "Usage: $0 audit|quantify|validate|all" >&2
  exit 2
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
sif="${SCWORKBENCH_IMAGE:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif}"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scmulti_flt3ko_il7ko"
bdev_matrix="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/bdev_reference/bundle/mat/mat_frag_peaks_counts_20260914_162209"
inputs="${storage}/reference_mapping_inputs"
audit_dir="${storage}/reports/bdev_peak_liftover_v1_20260930"
quant_dir="${storage}/prepared/bdev_peak_quantification_v1_20260930"
validation_dir="${storage}/reports/bdev_peak_quantification_v1_20260930"
fragments="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/datasets/scmulti_flt3ko_il7ko/FLT3LKO_IL7KO_BM/outs/atac_fragments.tsv.gz"
cells="${storage}/reports/qc_selection_v1_20260929/downstream_cell_ids.txt"
qc_selection="${storage}/reports/qc_selection_v1_20260929/cell_qc_selection.tsv.gz"
bdev_bundle="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/bdev_reference/bundle"
ko_bundle="${storage}/bundle"

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
else
  container_runtime="singularity"
fi

run_r() {
  "${container_runtime}" exec --cleanenv \
    -B /mnt/sata2:/mnt/sata2 \
    -B "${repo_root}:/repo" \
    "${sif}" Rscript --vanilla "$@"
}

if [[ "${stage}" == audit || "${stage}" == all ]]; then
  run_r /repo/processing/scmulti_flt3ko_il7ko/06_prepare_bdev_peak_universe.R \
    "${bdev_matrix}" \
    "${inputs}/mm10ToMm39.over.chain" \
    "${inputs}/mm39ToMm10.over.chain" \
    "${audit_dir}" 0.90
fi

if [[ "${stage}" == quantify || "${stage}" == all ]]; then
  run_r /repo/processing/scmulti_flt3ko_il7ko/07_quantify_bdev_peaks_in_ko.R \
    "${audit_dir}/accepted_bdev_peak_mapping.tsv.gz" \
    "${fragments}" "${cells}" "${quant_dir}" 2000
fi

if [[ "${stage}" == validate || "${stage}" == all ]]; then
  run_r /repo/processing/scmulti_flt3ko_il7ko/08_validate_bdev_peak_quantification.R \
    "${quant_dir}" "${qc_selection}" "${bdev_matrix}" \
    "${bdev_bundle}" "${ko_bundle}" "${validation_dir}" \
    1027f6390c6b853acc2bf40d3fae9abb59393f4256cb39285685eb2b2973f7ec
fi

printf 'LiftOver audit: %s\n' "${audit_dir}/completion.json"
printf 'Quantification: %s\n' "${quant_dir}/completion.json"
printf 'Validation: %s\n' "${validation_dir}/validation.json"
