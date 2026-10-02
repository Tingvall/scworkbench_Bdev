#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
ko_bundle="${storage}/scmulti_flt3ko_il7ko/bundle"
peak2_report="${storage}/bdev_integrated/reports/ebf1_enhancer_il7_peak2_matrix_aligned_atac_tracks_v3_20261002"
session_workspace="${storage}/bdev_integrated/sessions/il7_peak2_il7r_expression_v1"
output_dir="${1:-${storage}/bdev_integrated/reports/il7_peak2_top_de_genes_v1_20261002}"

if (( $# > 1 )); then
  echo "Usage: $0 [OUTPUT_DIR]" >&2
  exit 1
fi
for path in "${sif}" "${library}/scworkbench" "${ko_bundle}" \
  "${peak2_report}/all_four_samples_il7_peak2_pt_cell_assignments.tsv.gz" \
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
  /repo/figure_scripts/bdev_integrated/09_il7_peak2_top_de_genes.R \
  "${ko_bundle}" "${peak2_report}" "${session_workspace}" "${output_dir}" \
  "${ko_fp}" "v1_20261002"

printf 'IL-7 peak-2 top differential genes report: %s\n' "${output_dir}/REPORT.html"
