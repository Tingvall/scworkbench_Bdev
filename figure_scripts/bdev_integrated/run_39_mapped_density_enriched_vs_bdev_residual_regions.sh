#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
DENSITY="${STORAGE}/bdev_integrated/reports/ko_cytokine_gene_expression_bdev_normalized_v4_20261001/query_mapped_pt_cell_distribution.tsv"
SHARED="${STORAGE}/bdev_integrated/reports/ko_shared_enriched_peak_rna_de_v1_20261002/shared_enriched_pt_regions.tsv"
RESIDUAL="${STORAGE}/bdev_integrated/reports/cytokine_ucell_rna_feature_residual_dynamic_range_perm1000_v1_20261012/significant_positive_regions.tsv"
OUTPUT="${STORAGE}/bdev_integrated/reports/cytokine_ucell_rna_feature_residual_dynamic_range_perm1000_v1_20261012/summary_figures/mapped_density_enriched_vs_bdev_residual_regions"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
REGISTERED_LIBRARY="${STORAGE}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"

[[ -f "${DENSITY}" && -f "${SHARED}" && -f "${RESIDUAL}" ]] || exit 1
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.incomplete" ]] || exit 1

"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${REGISTERED_LIBRARY}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/figure_scripts/bdev_integrated/39_mapped_density_enriched_vs_bdev_residual_regions.R" \
  "${DENSITY}" "${SHARED}" "${RESIDUAL}" "${OUTPUT}"
