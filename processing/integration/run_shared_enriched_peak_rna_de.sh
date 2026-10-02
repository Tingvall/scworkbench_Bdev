#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
ko_bundle="${storage}/scmulti_flt3ko_il7ko/bundle"
distribution_report="${storage}/bdev_integrated/reports/ko_cytokine_gene_expression_bdev_normalized_v4_20261001"
projection_report="${storage}/bdev_integrated/reports/ko_mapping_reference_umap_v1_20261001"
original_region_report="${storage}/bdev_integrated/reports/ko_wt_enriched_peak_rna_de_v1_20261001"
session_workspace="${storage}/bdev_integrated/sessions/ko_wt_enriched_pt_peak_de_v1"
output_dir="${1:-${storage}/bdev_integrated/reports/ko_shared_enriched_peak_rna_de_v1_20261002}"

if (( $# > 1 )); then
  echo "Usage: $0 [OUTPUT_DIR]" >&2
  exit 1
fi
for path in "${sif}" "${library}/scworkbench" "${ko_bundle}" \
  "${distribution_report}/query_mapped_pt_cell_distribution.tsv" \
  "${projection_report}/query_reference_umap_projection.tsv.gz" \
  "${original_region_report}/wt_enriched_pt_regions.tsv" \
  "${session_workspace}/CURRENT.json"; do
  [[ -e "${path}" ]] || { echo "Required input missing: ${path}" >&2; exit 1; }
done
[[ ! -e "${output_dir}" ]] || { echo "Output already exists: ${output_dir}" >&2; exit 1; }
[[ ! -e "${output_dir}.incomplete" ]] || { echo "Incomplete output exists: ${output_dir}.incomplete" >&2; exit 1; }
case "${output_dir}" in
  "${storage}/bdev_integrated/reports/"*) ;;
  *) echo "Output must be inside ${storage}/bdev_integrated/reports" >&2; exit 1 ;;
esac

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
ko_fp="$(python3 - "${repo_root}/config/bundles.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle)["bundles"]["scmulti_flt3ko_il7ko"]["accepted_state"]["bundle_state_fingerprint"])
PY
)"
mkdir -p "$(dirname "${output_dir}")"
ulimit -v $((64 * 1024 * 1024))
"${container_runtime}" exec --cleanenv --no-home \
  -B "${library}:${library}:ro" --env "R_LIBS_USER=${library}" \
  -B "${repo_root}:/repo:ro" -B "${storage}:${storage}:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4 \
  --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/17_shared_enriched_peak_rna_de.R \
  "${ko_bundle}" "${distribution_report}" "${projection_report}" \
  "${original_region_report}" "${session_workspace}" "${output_dir}" \
  "${ko_fp}" "v1_20261002" execute
printf 'Shared-enriched peak RNA DE report: %s\n' "${output_dir}/REPORT.html"
