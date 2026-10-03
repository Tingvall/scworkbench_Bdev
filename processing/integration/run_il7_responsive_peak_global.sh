#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
bdev_bundle="${storage}/bdev_reference/bundle"
query_bundle="${storage}/bdev_integrated/ko_mapping_v1/bundle"
peak_mapping="${storage}/scmulti_flt3ko_il7ko/prepared/bdev_peak_quantification_v1_20260930/accepted_bdev_peak_mapping.tsv.gz"
shared_report="${storage}/bdev_integrated/reports/ko_shared_enriched_peak_rna_de_v1_20261002"
signatures="${storage}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925/passing_gene_membership.tsv"
bdev_rdata="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/datasets/scmulti_bdev/Bdev.Rdata"
global_report="${storage}/bdev_integrated/reports/il7_responsive_peak_global_v1_20261002"
ebf1_report="${storage}/bdev_integrated/reports/ebf1_il7_responsive_enhancers_v1_20261002"
bdev_sessions="${storage}/bdev_integrated/sessions"
query_sessions="${storage}/bdev_integrated/ko_mapping_v1/sessions"

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
else
  container_runtime="singularity"
fi

mkdir -p "${storage}/bdev_integrated/reports" "${bdev_sessions}" "${query_sessions}"
ulimit -v $((64 * 1024 * 1024))
common=(
  exec --cleanenv --no-home
  -B "${library}:${library}:ro"
  --env "R_LIBS_USER=${library}"
  -B "${repo_root}:/repo:ro"
  -B "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw"
  --env OMP_NUM_THREADS=4
  --env OPENBLAS_NUM_THREADS=4
  --env MKL_NUM_THREADS=4
  --pwd /repo
  "${sif}"
  Rscript --vanilla
)

if [[ ! -f "${global_report}/completion.json" ]]; then
  "${container_runtime}" "${common[@]}" \
    /repo/processing/integration/23_il7_responsive_peak_global.R \
    "${bdev_bundle}" "${query_bundle}" "${peak_mapping}" "${shared_report}" \
    "${signatures}" "${bdev_rdata}" "${global_report}" "${bdev_sessions}" \
    "${query_sessions}" \
    1027f6390c6b853acc2bf40d3fae9abb59393f4256cb39285685eb2b2973f7ec \
    68507e76b152839a5f93283a3e227d98ab3ec4eaf7804337090cfa9401aacef1 \
    23_il7_responsive_peak_global.R:v1 1729
fi

completion_sha="$(${container_runtime} "${common[@]}" -e \
  'cat(digest::digest(commandArgs(TRUE)[1], algo="sha256", file=TRUE))' \
  "${global_report}/completion.json")"
"${container_runtime}" "${common[@]}" \
  /repo/processing/integration/24_ebf1_il7_responsive_enhancers.R \
  "${global_report}" "${ebf1_report}" \
  24_ebf1_il7_responsive_enhancers.R:v1 "${completion_sha}"

printf 'Global report: %s\n' "${global_report}/REPORT.md"
printf 'Ebf1 report: %s\n' "${ebf1_report}/REPORT.md"
