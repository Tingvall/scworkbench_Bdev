#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BUNDLE="${STORAGE}/scmulti_flt3ko_il7ko/bundle"
SESSION="${STORAGE}/bdev_integrated/sessions/ko_wt_enriched_pt_peak_de_v1"
PROJECTION="${STORAGE}/bdev_integrated/reports/ko_mapping_reference_umap_v1_20261001/query_reference_umap_projection.tsv.gz"
DENSITY="${STORAGE}/bdev_integrated/reports/ko_cytokine_gene_expression_bdev_normalized_v4_20261001/query_mapped_pt_cell_distribution.tsv"
REGIONS="${STORAGE}/bdev_integrated/reports/ko_shared_enriched_peak_rna_de_v1_20261002/shared_enriched_pt_regions.tsv"
GENES="${STORAGE}/bdev_integrated/reports/ko_shared_enriched_peak_rna_de_v1_20261002/flt3l_shared_peak1_all_genes_normalized_counts_and_de.tsv"
CYTOKINES="${STORAGE}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925/passing_gene_membership.tsv"
RUN_ID="${1:-flt3l_il7_ko_vs_wt_pseudotime_matched_gsea_v1_20261012}"
OUTPUT="${STORAGE}/bdev_integrated/reports/${RUN_ID}"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
R_LIB="${GSEA_R_LIB:-/tmp/scworkbench_bdev_gsea_rlib}"
CACHE="${GSEA_CACHE_DIR:-/tmp/scworkbench_bdev_gsea_cache}"
SEED="${GSEA_SEED:-1729}"

[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ && "${SEED}" =~ ^[0-9]+$ ]] || exit 2
for path in "${BUNDLE}" "${SESSION}" "${PROJECTION}" "${DENSITY}" \
  "${REGIONS}" "${GENES}" "${CYTOKINES}"; do [[ -e "${path}" ]] || exit 1; done
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.incomplete" ]] || exit 1
mkdir -p "${R_LIB}" "${CACHE}"

"${RUNTIME}" exec --cleanenv --no-home \
  --env "R_LIBS_USER=${R_LIB}" --env "XDG_CACHE_HOME=${CACHE}" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${STORAGE}/bdev_integrated:${STORAGE}/bdev_integrated:rw" \
  --pwd "${REPO_ROOT}" "${SIF}" Rscript --vanilla \
  "${REPO_ROOT}/processing/integration/75_mapped_enrichment_pseudotime_matched_gsea.R" \
  "${BUNDLE}" "${SESSION}" "${PROJECTION}" "${DENSITY}" "${REGIONS}" \
  "${GENES}" "${CYTOKINES}" "${OUTPUT}" "${SEED}" \
  2>&1 | tee "${STORAGE}/bdev_integrated/reports/${RUN_ID}.log"
