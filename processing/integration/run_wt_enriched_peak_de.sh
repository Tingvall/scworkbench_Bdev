#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
ko_bundle="${storage}/scmulti_flt3ko_il7ko/bundle"
distribution_report="${storage}/bdev_integrated/reports/ko_cytokine_gene_expression_bdev_normalized_v4_20261001"
projection_report="${storage}/bdev_integrated/reports/ko_mapping_reference_umap_v1_20261001"
parse_signal_report="${storage}/bdev_reference/reports/all_cytokine_refined_ucell_v1_20260925"
session_root="${storage}/bdev_integrated/sessions"
lineage_workbook="${storage}/bdev_integrated/external_inputs/S2589004225011162_mmc2.xlsx"
lineage_tsv="${storage}/bdev_integrated/external_inputs/S2589004225011162_Table_S1B_lineage_genes.tsv"
lineage_workbook_sha="1190244c4718f8af3cb91f8d882c4936f200c452222a8b7e9db6ad6584d62fb3"
output_dir="${1:-${storage}/bdev_integrated/reports/ko_wt_enriched_peak_rna_de_v1_20261001}"

if (( $# > 1 )); then
  echo "Usage: $0 [OUTPUT_DIR]" >&2
  exit 1
fi
for path in "${sif}" "${library}/scworkbench" "${ko_bundle}" \
  "${distribution_report}" "${projection_report}"; do
  [[ -e "${path}" ]] || { echo "Required input missing: ${path}" >&2; exit 1; }
done
[[ -f "${parse_signal_report}/passing_gene_membership.tsv" ]] || {
  echo "Accepted Parse-derived signal membership missing" >&2; exit 1;
}
[[ -f "${parse_signal_report}/completion.json" ]] || {
  echo "Accepted Parse-derived signal completion missing" >&2; exit 1;
}
[[ -f "${lineage_workbook}" ]] || { echo "Required lineage workbook missing: ${lineage_workbook}" >&2; exit 1; }
observed_lineage_sha="$(sha256sum "${lineage_workbook}" | awk '{print $1}')"
[[ "${observed_lineage_sha}" == "${lineage_workbook_sha}" ]] || {
  echo "Lineage workbook checksum mismatch" >&2; exit 1;
}
python3 "${repo_root}/processing/integration/extract_lineage_gene_sets.py" \
  "${lineage_workbook}" "${lineage_tsv}"
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

mkdir -p "$(dirname "${output_dir}")" "${session_root}"
ulimit -v $((64 * 1024 * 1024))
"${container_runtime}" exec --cleanenv --no-home \
  -B "${library}:${library}:ro" --env "R_LIBS_USER=${library}" \
  -B "${repo_root}:/repo:ro" -B "${storage}:${storage}:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4 \
  --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/12_wt_enriched_peak_de.R \
  "${ko_bundle}" "${distribution_report}" "${projection_report}" \
  "${session_root}" "${output_dir}" "${ko_fp}" "${parse_signal_report}" \
  "${lineage_workbook}" "${lineage_tsv}" "${lineage_workbook_sha}" \
  "v1_20261001"

printf 'WT-enriched peak DE report: %s\n' "${output_dir}/REPORT.html"
