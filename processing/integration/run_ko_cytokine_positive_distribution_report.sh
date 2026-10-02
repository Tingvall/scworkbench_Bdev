#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
bdev_bundle="${storage}/bdev_reference/bundle"
ucell_report="${storage}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925"
score_report="${storage}/bdev_reference/reports/all_cytokine_breadth_intensity_v2_20260929"
score_promotion="${storage}/bdev_reference/reports/all_cytokine_breadth_intensity_promotion_v1_20260929"
projection_report="${storage}/bdev_integrated/reports/ko_mapping_reference_umap_v1_20261001"
output_dir="${1:-${storage}/bdev_integrated/reports/ko_cytokine_positive_distribution_v2_20261001}"
replace_mode=""
extra_args=()
if (( $# > 1 )); then
  for flag in "${@:2}"; do
    case "${flag}" in
      replace) replace_mode="replace" ;;
      bdev_normalized|gene_expression_panels) ;;
      *) echo "Optional flags must be 'replace', 'bdev_normalized', or 'gene_expression_panels': ${flag}" >&2; exit 1 ;;
    esac
    extra_args+=("${flag}")
  done
fi

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

for path in "${bdev_bundle}" "${ucell_report}" "${score_report}" \
  "${score_promotion}" "${projection_report}"; do
  [[ -e "${path}" ]] || { echo "Required input missing: ${path}" >&2; exit 1; }
done
[[ -d "${library}/scworkbench" ]] || { echo "Registered scWorkbench runtime missing" >&2; exit 1; }
if [[ -e "${output_dir}" && "${replace_mode}" != "replace" ]]; then
  echo "Report already exists: ${output_dir}" >&2
  exit 1
fi
[[ ! -e "${output_dir}.incomplete" ]] || { echo "Incomplete report exists: ${output_dir}.incomplete" >&2; exit 1; }
[[ ! -e "${output_dir}.replacement_backup" ]] || { echo "Replacement backup exists: ${output_dir}.replacement_backup" >&2; exit 1; }
case "${output_dir}" in
  "${storage}/bdev_integrated/reports/"*) ;;
  *) echo "Output must be inside ${storage}/bdev_integrated/reports: ${output_dir}" >&2; exit 1 ;;
esac

bdev_fp="$(python3 - "${repo_root}/config/bundles.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle)["bundles"]["bdev_reference"]["accepted_state"]["bundle_state_fingerprint"])
PY
)"

mkdir -p "$(dirname "${output_dir}")"
ulimit -v $((64 * 1024 * 1024))
"${container_runtime}" exec --cleanenv --no-home \
  -B "${library}:${library}:ro" --env "R_LIBS_USER=${library}" \
  -B "${repo_root}:/repo:ro" -B "${storage}:${storage}:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4 \
  --pwd /repo "${sif}" Rscript --vanilla \
  /repo/figure_scripts/bdev_integrated/08_ko_cytokine_positive_distribution_report.R \
  "${bdev_bundle}" "${ucell_report}" "${score_report}" "${score_promotion}" \
  "${projection_report}" "${output_dir}" "${bdev_fp}" "${extra_args[@]}"

printf 'Revised cytokine distribution report: %s\n' "${output_dir}/REPORT.html"
