#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
DE_REPORT="${STORAGE}/bdev_integrated/reports/ko_shared_enriched_peak_rna_de_v1_20261002"
FLT3L_DE="${DE_REPORT}/flt3l_shared_peak1_all_genes_normalized_counts_and_de.tsv"
IL7_DE="${DE_REPORT}/il7_shared_peak2_all_genes_normalized_counts_and_de.tsv"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
OUTPUT_DIR="${1:-${STORAGE}/bdev_integrated/reports/flt3l_il7_ko_vs_wt_hallmark_gsea_$(date -u +%Y%m%dT%H%M%SZ)}"
SEED="${GSEA_SEED:-1729}"
GSEA_R_LIB="${GSEA_R_LIB:-/tmp/scworkbench_bdev_gsea_rlib}"
GSEA_CACHE_DIR="${GSEA_CACHE_DIR:-/tmp/scworkbench_bdev_gsea_cache}"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"

[[ -x "${RUNTIME}" ]] || { echo "Apptainer missing: ${RUNTIME}" >&2; exit 1; }
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output exists: ${OUTPUT_DIR}" >&2; exit 1; }
mkdir -p "${GSEA_R_LIB}" "${GSEA_CACHE_DIR}"

if ! "${RUNTIME}" exec --cleanenv --no-home \
    --env R_LIBS_USER="${GSEA_R_LIB}" \
    --env XDG_CACHE_HOME="${GSEA_CACHE_DIR}" \
    "${SIF}" Rscript --vanilla -e \
    'quit(status=ifelse(requireNamespace("fgsea",quietly=TRUE) && requireNamespace("msigdbr",quietly=TRUE),0,1))'; then
  "${RUNTIME}" exec --cleanenv --no-home \
    --env R_LIBS_USER="${GSEA_R_LIB}" \
    --env XDG_CACHE_HOME="${GSEA_CACHE_DIR}" \
    "${SIF}" Rscript --vanilla -e \
    'options(repos=c(CRAN="https://cloud.r-project.org")); if (!requireNamespace("msigdbr",quietly=TRUE)) install.packages("msigdbr",lib=Sys.getenv("R_LIBS_USER")); if (!requireNamespace("fgsea",quietly=TRUE)) BiocManager::install("fgsea",lib=Sys.getenv("R_LIBS_USER"),ask=FALSE,update=FALSE)'
fi

exec "${RUNTIME}" exec --cleanenv --no-home \
  --env R_LIBS_USER="${GSEA_R_LIB}" \
  --env XDG_CACHE_HOME="${GSEA_CACHE_DIR}" \
  --bind "${REPO_ROOT}:/repo:ro" \
  --bind "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw" \
  --pwd /repo "${SIF}" Rscript --vanilla \
  /repo/processing/integration/52_flt3l_il7_ko_vs_wt_hallmark_gsea.R \
  "${FLT3L_DE}" "${IL7_DE}" "${OUTPUT_DIR}" "${SEED}"
