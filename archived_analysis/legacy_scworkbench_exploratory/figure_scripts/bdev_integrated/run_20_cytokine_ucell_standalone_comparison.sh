#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
UCELL_ROOT="${STORAGE}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925"
Z_ROOT="${STORAGE}/bdev_integrated/reports/trajectory_go_gene_standardized_min20_perm1000_v1_20261011"
BB="${UCELL_ROOT}/bb_up_change_points.tsv"
OUTPUT_DIR="${1:-${Z_ROOT}/summary_figures/cytokine_ucell_standalone_comparison}"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"

[[ -f "${UCELL_ROOT}/completion.json" && -f "${Z_ROOT}/cytokines/completion.json" ]] || exit 1
[[ ! -e "${OUTPUT_DIR}" && ! -e "${OUTPUT_DIR}.incomplete" ]] || exit 1
"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${REGISTERED_LIBRARY}" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}/bdev_reference:${STORAGE}/bdev_reference:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/figure_scripts/bdev_integrated/20_cytokine_ucell_standalone_comparison.R" \
  "${UCELL_ROOT}" "${Z_ROOT}" "${BB}" "${OUTPUT_DIR}"
