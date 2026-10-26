#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
Z_ROOT="${STORAGE}/bdev_integrated/reports/trajectory_go_gene_standardized_min20_perm1000_v1_20261011"
COMPARISON="${Z_ROOT}/summary_figures/cytokine_ucell_standalone_comparison"
PERMUTATION="${STORAGE}/bdev_integrated/reports/cytokine_ucell_median_perm1000_v1_20261011"
MEMBERSHIP="${COMPARISON}/cytokine_ucell_kmeans_membership.tsv"
BB="${STORAGE}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925/bb_up_change_points.tsv"
OUTPUT="${PERMUTATION}/summary_figures/raw_ucell_vs_original_z_passing_comparison"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"

[[ -f "${COMPARISON}/completion.json" && -f "${PERMUTATION}/completion.json" ]] || exit 1
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.incomplete" ]] || exit 1

"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${REGISTERED_LIBRARY}" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}/bdev_reference:${STORAGE}/bdev_reference:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/figure_scripts/bdev_integrated/23_cytokine_raw_ucell_original_z_passing_heatmap.R" \
  "${COMPARISON}" "${PERMUTATION}" "${MEMBERSHIP}" "${BB}" "${OUTPUT}"
