#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BUNDLE="${STORAGE}/bdev_integrated/bundle"
SOURCE_REPORT="${STORAGE}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925"
INVENTORY="${SOURCE_REPORT}/passing_signature_inventory.tsv"
MEMBERSHIP="${SOURCE_REPORT}/passing_gene_membership.tsv"
VALIDATED_REPORT="${STORAGE}/bdev_integrated/reports/reviewed_cytokine_38_permutation_validation_v1_20261011"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
RUN_ID="${1:-reviewed_cytokine_38_positive_peak_windows_depth_adjusted_v2_20261011}"
N_PERMUTATIONS="${N_PERMUTATIONS:-1000}"
PERMUTATION_SEED="${PERMUTATION_SEED:-2718}"
CLUSTER_ALPHA="${CLUSTER_ALPHA:-0.01}"
REGION_Q_CUTOFF="${REGION_Q_CUTOFF:-0.01}"
OUTPUT_DIR="${STORAGE}/bdev_integrated/reports/${RUN_ID}"

[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ ]] || exit 2
[[ "${N_PERMUTATIONS}" =~ ^[1-9][0-9]*$ ]] || exit 2
[[ "${PERMUTATION_SEED}" =~ ^[0-9]+$ ]] || exit 2
[[ -f "${INVENTORY}" && -f "${MEMBERSHIP}" && -d "${BUNDLE}" ]] || exit 1
[[ -f "${VALIDATED_REPORT}/completion.json" ]] || exit 1
[[ ! -e "${OUTPUT_DIR}" && ! -e "${OUTPUT_DIR}.incomplete" ]] || {
  echo "Output or incomplete output already exists: ${OUTPUT_DIR}" >&2; exit 1;
}

"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${REGISTERED_LIBRARY}" \
  --env OMP_NUM_THREADS=2 --env OPENBLAS_NUM_THREADS=2 --env MKL_NUM_THREADS=2 \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/processing/integration/55_reviewed_cytokine_multiregion_report.R" \
  "${BUNDLE}" "${INVENTORY}" "${MEMBERSHIP}" "${VALIDATED_REPORT}" \
  "${OUTPUT_DIR}" "${N_PERMUTATIONS}" "${PERMUTATION_SEED}" \
  positive_only "${CLUSTER_ALPHA}" "${REGION_Q_CUTOFF}" peak_window
