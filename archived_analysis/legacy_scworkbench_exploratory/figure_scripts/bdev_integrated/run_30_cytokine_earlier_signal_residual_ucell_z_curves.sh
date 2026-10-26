#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
DYNAMIC="${STORAGE}/bdev_integrated/reports/cytokine_ucell_dynamic_range_perm1000_v1_20261012"
COMPARISON="${DYNAMIC}/diagnostics/rna_feature_residualization_curve_comparison"
CURVES="${COMPARISON}/cytokine_raw_residual_ucell_original_z_scaled_0_1_curves.tsv.gz"
EARLIER="${STORAGE}/bdev_integrated/reports/reviewed_cytokine_38_positive_peak_windows_depth_adjusted_v2_20261011"
MEMBERSHIP="${STORAGE}/bdev_integrated/reports/trajectory_go_gene_standardized_min20_perm1000_v1_20261011/summary_figures/cytokine_ucell_standalone_comparison/cytokine_ucell_kmeans_membership.tsv"
BB="${STORAGE}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925/bb_up_change_points.tsv"
OUTPUT="${COMPARISON}/earlier_depth_adjusted_signal_residual_ucell_z_comparison"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"

[[ -f "${CURVES}" && -f "${EARLIER}/window_statistics.tsv.gz" ]] || exit 1
[[ -f "${MEMBERSHIP}" && -f "${BB}" ]] || exit 1
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.incomplete" ]] || exit 1

"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${REGISTERED_LIBRARY}" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}/bdev_reference:${STORAGE}/bdev_reference:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/figure_scripts/bdev_integrated/30_cytokine_earlier_signal_residual_ucell_z_curves.R" \
  "${CURVES}" "${EARLIER}/window_statistics.tsv.gz" "${MEMBERSHIP}" \
  "${BB}" "${EARLIER}/completion.json" "${OUTPUT}"
