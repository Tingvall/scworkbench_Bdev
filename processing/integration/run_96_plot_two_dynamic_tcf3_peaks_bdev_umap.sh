#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
bundle="${storage}/bdev_reference/bundle"
source_values="${storage}/bdev_integrated/reports/ebf1_tcf3_peaks_bdev_umap_v3_20261005/all_cells_eight_peak_log1p_cp10k_accessibility.tsv.gz"
three_way="${storage}/bdev_integrated/reports/ebf1_putative_enhancers_three_way_v2_20261005/all_259_bdev_ebf1_window_peaks_three_way_tests.tsv"
output_dir="${1:-${storage}/bdev_integrated/reports/two_dynamic_tcf3_peaks_bdev_umap_v1_20261013}"
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
  /repo/processing/integration/96_plot_two_dynamic_tcf3_peaks_bdev_umap.R \
  "${bundle}" "${source_values}" "${three_way}" "${output_dir}" \
  "${bundle_fp}" "v1_20261013"

printf 'Two-peak Bdev UMAP report: %s/REPORT.md\n' "${output_dir}"
