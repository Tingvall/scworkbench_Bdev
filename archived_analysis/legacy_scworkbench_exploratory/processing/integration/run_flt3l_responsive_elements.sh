#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
project_storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
runtime="${SCWORKBENCH_RUNTIME:-/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
bdev_bundle="${storage}/bdev_reference/bundle"
query_bundle="${storage}/bdev_integrated/ko_mapping_v1/bundle"
peak_mapping="${storage}/scmulti_flt3ko_il7ko/prepared/bdev_peak_quantification_v1_20260930/accepted_bdev_peak_mapping.tsv.gz"
shared_report="${storage}/bdev_integrated/reports/ko_shared_enriched_peak_rna_de_v1_20261002"
signatures="${storage}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925/passing_gene_membership.tsv"
bdev_rdata="${project_storage}/datasets/scmulti_bdev/Bdev.Rdata"
global_report="${storage}/bdev_integrated/reports/flt3l_responsive_peak_global_v1_20261002"
annotation_report="${storage}/bdev_integrated/reports/flt3l_responsive_peak_annotation_v1_20261002"
local_report="${storage}/bdev_integrated/reports/flt3l_responsive_peak_local_deg_v1_20261002"
umap_report="${storage}/bdev_integrated/reports/flt3l_responsive_peak_linked_top20_bdev_umap_v1_20261002"
bdev_sessions="${storage}/bdev_integrated/sessions"
query_sessions="${storage}/bdev_integrated/ko_mapping_v1/sessions"
bdev_fingerprint="1027f6390c6b853acc2bf40d3fae9abb59393f4256cb39285685eb2b2973f7ec"
query_fingerprint="68507e76b152839a5f93283a3e227d98ab3ec4eaf7804337090cfa9401aacef1"

mkdir -p "${storage}/bdev_integrated/reports" "${bdev_sessions}" "${query_sessions}"
ulimit -v $((64 * 1024 * 1024))
common=(
  exec --cleanenv --no-home
  -B "${library}:${library}:ro"
  --env "R_LIBS_USER=${library}"
  -B "${repo_root}:/repo:ro"
  -B "${project_storage}:${project_storage}:rw"
  --env OMP_NUM_THREADS=4
  --env OPENBLAS_NUM_THREADS=4
  --env MKL_NUM_THREADS=4
  --pwd /repo
  "${sif}"
  Rscript --vanilla
)

file_sha256() {
  "${runtime}" "${common[@]}" -e \
    'cat(digest::digest(commandArgs(TRUE)[1], algo="sha256", file=TRUE))' "$1"
}

if [[ ! -f "${global_report}/completion.json" ]]; then
  "${runtime}" "${common[@]}" \
    /repo/processing/integration/29_flt3l_responsive_peak_global.R \
    "${bdev_bundle}" "${query_bundle}" "${peak_mapping}" "${shared_report}" \
    "${signatures}" "${bdev_rdata}" "${global_report}" "${bdev_sessions}" \
    "${query_sessions}" "${bdev_fingerprint}" "${query_fingerprint}" \
    29_flt3l_responsive_peak_global.R:v1 2718
fi

global_sha="$(file_sha256 "${global_report}/completion.json")"
if [[ ! -f "${annotation_report}/completion.json" ]]; then
  "${runtime}" "${common[@]}" \
    /repo/processing/integration/30_annotate_flt3l_validated_peaks.R \
    "${global_report}" "${bdev_rdata}" "${annotation_report}" \
    "${global_sha}" 30_annotate_flt3l_validated_peaks.R:v1 2000
fi

annotation_sha="$(file_sha256 "${annotation_report}/completion.json")"
deg_sha="$(file_sha256 "${shared_report}/completion.json")"
if [[ ! -f "${local_report}/completion.json" ]]; then
  "${runtime}" "${common[@]}" \
    /repo/processing/integration/31_flt3l_peak_local_deg_windows.R \
    "${annotation_report}" "${shared_report}" "${bdev_rdata}" \
    "${local_report}" "${annotation_sha}" "${deg_sha}" \
    31_flt3l_peak_local_deg_windows.R:v1 3000000 0.05
fi

local_sha="$(file_sha256 "${local_report}/completion.json")"
if [[ ! -f "${umap_report}/completion.json" ]]; then
  "${runtime}" "${common[@]}" \
    /repo/processing/integration/32_plot_flt3l_peak_linked_top20_bdev_umap.R \
    "${local_report}" "${bdev_bundle}" "${umap_report}" "${local_sha}" \
    "${bdev_fingerprint}" \
    32_plot_flt3l_peak_linked_top20_bdev_umap.R:v1 20
fi

printf 'Global report: %s\n' "${global_report}/REPORT.md"
printf 'Annotation report: %s\n' "${annotation_report}/REPORT.md"
printf 'Putative target-gene report: %s\n' "${local_report}/REPORT.md"
printf 'Bdev UMAP report: %s\n' "${umap_report}/REPORT.md"
