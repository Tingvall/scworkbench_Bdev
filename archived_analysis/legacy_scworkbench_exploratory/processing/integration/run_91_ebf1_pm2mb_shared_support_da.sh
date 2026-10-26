#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
bundle="${storage}/bdev_integrated/ko_mapping_v1/bundle"
peak_mapping="${storage}/scmulti_flt3ko_il7ko/prepared/bdev_peak_quantification_v1_20260930/accepted_bdev_peak_mapping.tsv.gz"
support_report="${storage}/bdev_integrated/reports/wt_enrichment_rolling_3bin_cell_support_v2_20261013"
query_support="${storage}/bdev_integrated/reports/bdev_trajectory_support_mapping_comparison_v2_20261012/query_trajectory_support_and_acceptance.tsv.gz"
output_dir="${1:-${storage}/bdev_integrated/reports/ebf1_pm2mb_shared_support_da_v1_20261013}"
session_root="${storage}/bdev_integrated/sessions"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
runtime="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
seed="${DA_SEED:-1729}"

bundle_fp="$(python3 - "${bundle}/manifest.json" <<'PY'
import json, sys
with open(sys.argv[1]) as handle:
    print(json.load(handle)['current_state']['fingerprint'])
PY
)"
[[ -x "${runtime}" ]] || { echo "Apptainer missing: ${runtime}" >&2; exit 1; }
for path in "${bundle}" "${peak_mapping}" "${support_report}" "${query_support}"; do
  [[ -e "${path}" ]] || { echo "Required input missing: ${path}" >&2; exit 1; }
done
[[ ! -e "${output_dir}" ]] || { echo "Output exists: ${output_dir}" >&2; exit 1; }
[[ ! -e "${output_dir}.incomplete" ]] || { echo "Incomplete output exists: ${output_dir}.incomplete" >&2; exit 1; }
case "${output_dir}" in
  "${storage}/bdev_integrated/reports/"*) ;;
  *) echo "Output must be below ${storage}/bdev_integrated/reports" >&2; exit 1 ;;
esac

"${runtime}" exec --cleanenv --no-home \
  --bind "${repo_root}:/repo:ro" \
  --bind "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 \
  --env MKL_NUM_THREADS=4 --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/91_ebf1_pm2mb_shared_support_da.R \
  "${bundle}" "${peak_mapping}" "${support_report}" "${query_support}" \
  "${output_dir}" "${session_root}" "${bundle_fp}" "${seed}" "v1_20261013" &
analysis_pid=$!
while kill -0 "${analysis_pid}" 2>/dev/null; do
  sleep 20
  if kill -0 "${analysis_pid}" 2>/dev/null; then
    printf 'Ebf1 ±2 Mb shared-support DA still running (%s)\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  fi
done
set +e
wait "${analysis_pid}"
status=$?
set -e
if [[ "${status}" -ne 0 ]]; then
  printf 'Ebf1 ±2 Mb shared-support DA failed with exit status %d\n' "${status}" >&2
  exit "${status}"
fi
printf 'Ebf1 ±2 Mb shared-support DA report: %s\n' "${output_dir}/REPORT.md"
