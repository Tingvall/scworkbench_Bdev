#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
bdev_bundle="${storage}/bdev_reference/bundle"
qc_dir="${storage}/bdev_integrated/reports/ko_mapping_qc_v1_20260930"
report_dir="${1:-${storage}/bdev_integrated/reports/ko_mapping_reference_umap_v1_20261001}"

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
elif [[ -x /home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer ]]; then
  container_runtime="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
elif command -v singularity >/dev/null 2>&1; then
  container_runtime="singularity"
else
  echo "No Apptainer or Singularity runtime found" >&2
  exit 1
fi
[[ -d "${library}/scworkbench" ]] || { echo "Registered scWorkbench runtime missing" >&2; exit 1; }
[[ -d "${bdev_bundle}" ]] || { echo "Bdev reference bundle missing" >&2; exit 1; }
[[ -f "${qc_dir}/query_mapping.tsv.gz" ]] || { echo "Query mapping table missing" >&2; exit 1; }
[[ -f "${qc_dir}/query_bdev_neighbors_k50.tsv.gz" ]] || { echo "Bdev neighbor table missing" >&2; exit 1; }
[[ -f "${qc_dir}/metrics.json" ]] || { echo "Mapping metrics missing" >&2; exit 1; }
[[ ! -e "${report_dir}" ]] || { echo "Projection report already exists: ${report_dir}" >&2; exit 1; }

case "${report_dir}" in
  "${storage}/bdev_integrated/reports/"*) ;;
  *) echo "Report must be inside ${storage}/bdev_integrated/reports: ${report_dir}" >&2; exit 1 ;;
esac

mkdir -p "$(dirname "${report_dir}")"
ulimit -v $((64 * 1024 * 1024))

"${container_runtime}" exec --cleanenv --no-home \
  -B "${library}:${library}:ro" --env "R_LIBS_USER=${library}" \
  -B "${repo_root}:/repo:ro" -B "${storage}:${storage}:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4 \
  --pwd /repo "${sif}" Rscript --vanilla \
  /repo/figure_scripts/bdev_integrated/06_peakvi_reference_umap_projection.R \
  "${bdev_bundle}" "${qc_dir}/query_mapping.tsv.gz" \
  "${qc_dir}/query_bdev_neighbors_k50.tsv.gz" "${qc_dir}/metrics.json" \
  "${report_dir}"

printf 'Reference-anchored projection report: %s\n' "${report_dir}/REPORT.md"
