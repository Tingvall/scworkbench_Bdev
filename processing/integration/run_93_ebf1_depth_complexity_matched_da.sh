#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
bundle="${storage}/bdev_integrated/ko_mapping_v1/bundle"
source_report="${storage}/bdev_integrated/reports/ebf1_pm2mb_shared_support_da_v1_20261013"
output_dir="${1:-${storage}/bdev_integrated/reports/ebf1_pm2mb_depth_complexity_matched_da_v1_20261013}"
session_root="${storage}/bdev_integrated/sessions"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
runtime="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
caliper_sd="${MATCH_CALIPER_SD:-0.2}"

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
  /repo/processing/integration/93_ebf1_depth_complexity_matched_da.R \
  "${bundle}" "${source_report}" "${output_dir}" "${session_root}" \
  "${bundle_fp}" "${caliper_sd}" "v1_20261013"

printf 'Ebf1 depth/complexity-matched DA report: %s/REPORT.md\n' "${output_dir}"
