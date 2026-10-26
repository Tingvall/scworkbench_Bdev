#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
bundle="${storage}/bdev_integrated/ko_mapping_v1/bundle"
support="${storage}/bdev_integrated/reports/ebf1_pm2mb_shared_support_da_v1_20261013"
complexity="${storage}/bdev_integrated/reports/ebf1_pm2mb_depth_complexity_matched_da_v1_20261013"
tcf3="${storage}/bdev_integrated/reports/ebf1_putative_enhancers_three_way_v2_20261005/all_259_bdev_ebf1_window_peaks_three_way_tests.tsv"
output_dir="${1:-${storage}/bdev_integrated/reports/tcf3_peaks_detection_normalized_wt_ko_ratio_v2_20261013}"
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
  --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/97_tcf3_peaks_detection_normalized_wt_ko_ratio.R \
  "${bundle}" "${support}" "${complexity}" "${tcf3}" "${output_dir}" \
  "${bundle_fp}" "v2_20261013"

printf 'TCF3 normalized WT/KO ratio report: %s/REPORT.md\n' "${output_dir}"
