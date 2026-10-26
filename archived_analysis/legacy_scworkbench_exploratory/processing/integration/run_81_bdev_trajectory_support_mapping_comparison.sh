#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
output_dir="${1:-${storage}/bdev_integrated/reports/bdev_trajectory_support_mapping_comparison_v2_20261012}"
runtime="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"

[[ -x "${runtime}" ]] || runtime="$(command -v apptainer || command -v singularity)"
[[ -x "${runtime}" ]] || { echo "No Apptainer/Singularity runtime found" >&2; exit 1; }
[[ -f "${sif}" ]] || { echo "Missing image: ${sif}" >&2; exit 1; }
[[ -d "${library}/scworkbench" ]] || { echo "Missing scWorkbench runtime library" >&2; exit 1; }
[[ ! -e "${output_dir}" ]] || { echo "Output already exists: ${output_dir}" >&2; exit 1; }
case "${output_dir}" in
  "${storage}/bdev_integrated/reports/"*) ;;
  *) echo "Output must be below ${storage}/bdev_integrated/reports" >&2; exit 1 ;;
esac

"${runtime}" exec --cleanenv --no-home \
  -B "${library}:${library}:ro" --env "R_LIBS_USER=${library}" \
  -B "${repo_root}:/repo:ro" -B "${storage}:${storage}:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4 \
  --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/81_bdev_trajectory_support_mapping_comparison.R \
  "${storage}/bdev_integrated/ko_mapping_v1/bundle" \
  "${storage}/bdev_integrated/ko_mapping_v1/sessions/peakvi_full_mapping_v1/stage" \
  "${storage}/bdev_reference/bundle" \
  "${storage}/bdev_integrated/reports/ko_mapping_reference_umap_v1_20261001/query_reference_umap_projection.tsv.gz" \
  "${storage}/bdev_integrated/reports/ko_mapping_qc_v1_20260930/query_bdev_neighbors_k50.tsv.gz" \
  "${storage}/bdev_integrated/reports/ko_wt_enriched_peak_rna_de_v1_20261001/flt3l_il7_parse_signal_aggregate_cell_scores.tsv.gz" \
  "${storage}/bdev_integrated/reports/ko_wt_enriched_peak_rna_de_v1_20261001/wt_enriched_pt_regions.tsv" \
  "${storage}/bdev_integrated/reports/ko_mapping_peakvi_full_v1_20260930/mapping.json" \
  "${output_dir}" "v2_20261012"

printf 'Bdev trajectory-support mapping comparison: %s\n' "${output_dir}/REPORT.md"
