#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
runtime="${SCWORKBENCH_RUNTIME:-/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
volcano_report="${storage}/bdev_integrated/reports/il7_responsive_peak_linked_gene_volcano_v1_20261002"
bdev_bundle="${storage}/bdev_reference/bundle"
output="${storage}/bdev_integrated/reports/il7_responsive_peak_linked_top20_bdev_umap_v1_20261002"
bdev_fingerprint="1027f6390c6b853acc2bf40d3fae9abb59393f4256cb39285685eb2b2973f7ec"

common=(
  exec --cleanenv --no-home
  -B "${library}:${library}:ro"
  --env "R_LIBS_USER=${library}"
  -B "${repo_root}:/repo:ro"
  -B "${storage}:${storage}:rw"
  --pwd /repo
  "${sif}"
  Rscript --vanilla
)
volcano_sha="$(${runtime} "${common[@]}" -e \
  'cat(digest::digest(commandArgs(TRUE)[1], algo="sha256", file=TRUE))' \
  "${volcano_report}/completion.json")"
"${runtime}" "${common[@]}" \
  /repo/processing/integration/28_plot_il7_peak_linked_top20_bdev_umap.R \
  "${volcano_report}" "${bdev_bundle}" "${output}" "${volcano_sha}" \
  "${bdev_fingerprint}" 28_plot_il7_peak_linked_top20_bdev_umap.R:v1 20
printf 'Bdev UMAP report: %s\n' "${output}/REPORT.md"
