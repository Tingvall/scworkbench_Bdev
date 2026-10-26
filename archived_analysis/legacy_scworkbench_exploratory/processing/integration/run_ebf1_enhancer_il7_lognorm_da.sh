#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
bundle="${storage}/bdev_integrated/ko_mapping_v1/bundle"
peak_mapping="${storage}/scmulti_flt3ko_il7ko/prepared/bdev_peak_quantification_v1_20260930/accepted_bdev_peak_mapping.tsv.gz"
pt_report="${storage}/bdev_integrated/reports/ko_wt_enriched_peak_rna_de_v1_20261001"
raw_report="${storage}/bdev_integrated/reports/ebf1_enhancer_il7_enriched_pt_da_v1_20261001"
session_root="${storage}/bdev_integrated/sessions"
output_dir="${1:-${storage}/bdev_integrated/reports/ebf1_enhancer_il7_enriched_pt_lognorm_da_v1_20261001}"

if (( $# > 1 )); then echo "Usage: $0 [OUTPUT_DIR]" >&2; exit 1; fi
for path in "${sif}" "${library}/scworkbench" "${bundle}" "${peak_mapping}" \
  "${pt_report}" "${raw_report}"; do
  [[ -e "${path}" ]] || { echo "Required input missing: ${path}" >&2; exit 1; }
done
[[ ! -e "${output_dir}" ]] || { echo "Output already exists: ${output_dir}" >&2; exit 1; }
[[ ! -e "${output_dir}.incomplete" ]] || { echo "Incomplete output exists: ${output_dir}.incomplete" >&2; exit 1; }
case "${output_dir}" in
  "${storage}/bdev_integrated/reports/"*) ;;
  *) echo "Output must be inside ${storage}/bdev_integrated/reports" >&2; exit 1 ;;
esac

runtime="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
[[ -x "${runtime}" ]] || runtime="$(command -v apptainer || command -v singularity)"
bundle_fp="$(python3 - "${repo_root}/config/bundles.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    x = yaml.safe_load(handle)
print(x["bundles"]["bdev_integrated"]["staged_derivatives"]["ko_mapping_v1"]["bundle_state_fingerprint"])
PY
)"

mkdir -p "$(dirname "${output_dir}")" "${session_root}"
ulimit -v $((64 * 1024 * 1024))
"${runtime}" exec --cleanenv --no-home \
  -B "${library}:${library}:ro" --env "R_LIBS_USER=${library}" \
  -B "${repo_root}:/repo:ro" -B "${storage}:${storage}:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4 \
  --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/15_ebf1_enhancer_il7_lognorm_da.R \
  "${bundle}" "${peak_mapping}" "${pt_report}" "${raw_report}" \
  "${output_dir}" "${session_root}" "${bundle_fp}" "v1_20261001"

printf 'Normalized Ebf1 differential-accessibility report: %s\n' "${output_dir}/REPORT.html"
