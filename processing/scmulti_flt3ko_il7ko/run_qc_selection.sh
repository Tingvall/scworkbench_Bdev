#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
sif="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif"
qc_root="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scmulti_flt3ko_il7ko/reports/unfiltered_qc_v1_20260929"
selection_root="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scmulti_flt3ko_il7ko/reports/qc_selection_v1_20260929"

if [[ ! -f "${selection_root}/completion.json" ]]; then
  singularity exec --cleanenv -B /mnt/sata2:/mnt/sata2 -B "${repo_root}:/repo" "${sif}" \
    Rscript --vanilla /repo/processing/scmulti_flt3ko_il7ko/03_mark_qc_analysis_set.R \
    "${qc_root}" "${selection_root}"
fi

singularity exec --cleanenv -B /mnt/sata2:/mnt/sata2 -B "${repo_root}:/repo" "${sif}" \
  Rscript --vanilla /repo/processing/scmulti_flt3ko_il7ko/04_render_qc_selection_report.R \
  "${selection_root}" /repo/figure_scripts/scmulti_flt3ko_il7ko/02_qc_selection_report.Rmd \
  "${selection_root}/REPORT.md"

printf 'QC selection report: %s\n' "${selection_root}/REPORT.md"
