#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
bdev_bundle="${storage}/bdev_reference/bundle"
mapping_bundle="${storage}/bdev_integrated/ko_mapping_v1/bundle"
joint_matrix="${mapping_bundle}/mat/mat_frag_ATAC_bdev_mm10_peaks_20260930_122757"
support_report="${storage}/bdev_integrated/reports/wt_enrichment_rolling_3bin_cell_support_v2_20261013"
query_support="${storage}/bdev_integrated/reports/bdev_trajectory_support_mapping_comparison_v2_20261012/query_trajectory_support_and_acceptance.tsv.gz"
bdev_rdata="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/datasets/scmulti_bdev/Bdev.Rdata"
selection_report="${storage}/bdev_integrated/reports/fine_rolling_motif_trajectory_heatmaps_v2_20261013"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
output_dir="${1:-${storage}/bdev_integrated/reports/single_cell_axis_15nn_motif_heatmaps_v2_20261013}"
neighbor_count="${MOTIF_PT_NEIGHBOR_COUNT:-15}"
maximum_radius="${MOTIF_PT_MAX_NEIGHBOR_RADIUS:-0.25}"
seed="${MOTIF_SEED:-1729}"
runtime="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"

bdev_fp="$(python3 - "${repo_root}/config/bdev_reference.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle)['validation']['bundle_state_fingerprint'])
PY
)"
mapping_fp="$(python3 - "${mapping_bundle}/manifest.json" <<'PY'
import json, sys
with open(sys.argv[1]) as handle:
    print(json.load(handle)['current_state']['fingerprint'])
PY
)"
[[ -x "${runtime}" ]] || { echo "Apptainer missing: ${runtime}" >&2; exit 1; }
[[ ! -e "${output_dir}" ]] || { echo "Output exists: ${output_dir}" >&2; exit 1; }
case "${output_dir}" in
  "${storage}/bdev_integrated/reports/"*) ;;
  *) echo "Output must be below ${storage}/bdev_integrated/reports" >&2; exit 1 ;;
esac

"${runtime}" exec --cleanenv --no-home \
  --bind "${repo_root}:/repo:ro" \
  --bind "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 \
  --env MKL_NUM_THREADS=4 --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/89_single_cell_axis_15nn_motif_heatmaps.R \
  "${bdev_bundle}" "${mapping_bundle}" "${joint_matrix}" \
  "${support_report}" "${query_support}" "${bdev_rdata}" \
  "${selection_report}" "${output_dir}" "${bdev_fp}" "${mapping_fp}" \
  "${neighbor_count}" "${maximum_radius}" "${seed}" &
analysis_pid=$!
while kill -0 "${analysis_pid}" 2>/dev/null; do
  sleep 20
  if kill -0 "${analysis_pid}" 2>/dev/null; then
    printf 'Single-cell-axis %s-NN motif analysis still running (%s)\n' "${neighbor_count}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  fi
done
set +e
wait "${analysis_pid}"
status=$?
set -e
if [[ "${status}" -ne 0 ]]; then
  printf 'Single-cell-axis %s-NN motif analysis failed with exit status %d\n' "${neighbor_count}" "${status}" >&2
  exit "${status}"
fi
printf 'Single-cell-axis %s-NN motif report: %s\n' "${neighbor_count}" "${output_dir}/REPORT.md"
