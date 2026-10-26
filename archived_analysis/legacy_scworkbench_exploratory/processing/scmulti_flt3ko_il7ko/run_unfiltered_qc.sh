#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
sif="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif"
outs="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/datasets/scmulti_flt3ko_il7ko/FLT3LKO_IL7KO_BM/outs"
per_barcode="/mnt/synology/Naive/Personal_Storage/jtg/data/scmulti/counts/flt3lko_il7ko/per_barcode_metrics"
reference="/home/jtg/resources/refs/refdata-cellranger-arc-GRCm39-2024-A"
report_root="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scmulti_flt3ko_il7ko/reports/unfiltered_qc_v1_20260929"

if [[ ! -f "${report_root}/completion.json" ]]; then
  singularity exec --cleanenv \
    -B /mnt/sata2:/mnt/sata2 \
    -B /mnt/synology:/mnt/synology \
    -B /home/jtg/resources:/home/jtg/resources \
    -B "${repo_root}:/repo" \
    "${sif}" \
    Rscript --vanilla /repo/processing/scmulti_flt3ko_il7ko/01_calculate_unfiltered_qc.R \
    "${outs}" "${per_barcode}" "${reference}" "${report_root}"
fi

singularity exec --cleanenv \
  -B /mnt/sata2:/mnt/sata2 \
  -B "${repo_root}:/repo" \
  "${sif}" \
  Rscript --vanilla /repo/processing/scmulti_flt3ko_il7ko/02_render_unfiltered_qc_report.R \
  "${report_root}" \
  /repo/figure_scripts/scmulti_flt3ko_il7ko/01_unfiltered_qc_report.Rmd \
  "${report_root}/REPORT.md"

printf 'QC report: %s\n' "${report_root}/REPORT.md"
