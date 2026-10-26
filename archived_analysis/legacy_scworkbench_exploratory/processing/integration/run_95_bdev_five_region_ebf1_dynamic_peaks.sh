#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
bundle="${storage}/bdev_integrated/ko_mapping_v1/bundle"
source_report="${storage}/bdev_integrated/reports/ebf1_pm2mb_shared_support_da_v1_20261013"
tcf3_coverage="${storage}/bdev_integrated/reports/tcf3_ebf1_coverage_placseq_tracks_v4_20261005/condition_mean_normalized_50bp_coverage.tsv.gz"
density="${storage}/bdev_integrated/reports/bdev_cutoff025_mapping_report_v4_20261012/query_density_cutoff025.tsv.gz"
output_dir="${1:-${storage}/bdev_integrated/reports/bdev_five_region_ebf1_dynamic_peaks_v3_20261013}"
session_root="${storage}/bdev_integrated/sessions"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
runtime="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"

bundle_fp="$(python3 - "${bundle}/manifest.json" <<'PY'
import json, sys
with open(sys.argv[1]) as handle:
    print(json.load(handle)['current_state']['fingerprint'])
PY
)"

"${runtime}" exec --cleanenv --no-home \
  --bind "${repo_root}:/repo:ro" \
  --bind "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 \
  --env MKL_NUM_THREADS=4 --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/95_bdev_five_region_ebf1_dynamic_peaks.R \
  "${bundle}" "${source_report}" "${density}" "${output_dir}" \
  "${session_root}" "${bundle_fp}" "v3_20261013" "${tcf3_coverage}"

printf 'Five-state Bdev Ebf1 dynamic peak report: %s/REPORT.md\n' "${output_dir}"
