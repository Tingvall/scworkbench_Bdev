#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
output_dir="${1:-${storage}/bdev_integrated/reports/single_cell_axis_30nn_motif_heatmaps_v2_20261013}"
export MOTIF_PT_NEIGHBOR_COUNT=30
export MOTIF_PT_MAX_NEIGHBOR_RADIUS=Inf
exec bash "${script_dir}/run_89_single_cell_axis_15nn_motif_heatmaps.sh" "${output_dir}"
