#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
UCELL="${STORAGE}/bdev_integrated/reports/go_ucell_dynamic_range_perm1000_v1_20261012"
DIAGNOSTIC="${UCELL}/diagnostics/epigenetic_boundary_peak_alignment_rna_feature_residualization"
Z_ROOT="${STORAGE}/bdev_integrated/reports/trajectory_go_gene_standardized_min20_perm1000_v1_20261011"
SELECTION="${UCELL}/summary_tables/top3_epigenetic_dynamic_gene_sets/top3_epigenetic_dynamic_gene_sets_by_group.tsv"
BB="${STORAGE}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925/bb_up_change_points.tsv"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"

[[ -d "${DIAGNOSTIC}" && -d "${Z_ROOT}" && -f "${SELECTION}" ]] || exit 1

"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${REGISTERED_LIBRARY}" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/figure_scripts/bdev_integrated/28_epigenetic_raw_residual_ucell_z_curves.R" \
  "${DIAGNOSTIC}" "${Z_ROOT}" "${SELECTION}" "${BB}" "${DIAGNOSTIC}"
