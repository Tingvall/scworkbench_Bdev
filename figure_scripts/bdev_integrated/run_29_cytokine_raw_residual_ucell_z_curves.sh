#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BUNDLE="${STORAGE}/bdev_integrated/bundle"
NATIVE="${STORAGE}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925"
DYNAMIC="${STORAGE}/bdev_integrated/reports/cytokine_ucell_dynamic_range_perm1000_v1_20261012"
Z_DIR="${STORAGE}/bdev_integrated/reports/trajectory_go_gene_standardized_min20_perm1000_v1_20261011/cytokines"
BB="${NATIVE}/bb_up_change_points.tsv"
OUTPUT="${DYNAMIC}/diagnostics/rna_feature_residualization_curve_comparison"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"

[[ -d "${BUNDLE}" && -f "${NATIVE}/ucell_cell_scores.tsv.gz" ]] || exit 1
[[ -f "${DYNAMIC}/window_curves.tsv.gz" && -f "${Z_DIR}/window_curves.tsv.gz" ]] || exit 1
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.incomplete" ]] || exit 1

"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${REGISTERED_LIBRARY}" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}/bdev_reference:${STORAGE}/bdev_reference:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/figure_scripts/bdev_integrated/29_cytokine_raw_residual_ucell_z_curves.R" \
  "${BUNDLE}" "${NATIVE}" "${DYNAMIC}" "${Z_DIR}" "${BB}" "${OUTPUT}"
