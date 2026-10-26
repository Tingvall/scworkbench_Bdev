#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
da_report="${storage}/bdev_integrated/reports/ebf1_enhancer_il7_enriched_pt_da_v1_20261001"
pt_report="${storage}/bdev_integrated/reports/ko_wt_enriched_peak_rna_de_v1_20261001"
projection_report="${storage}/bdev_integrated/reports/ko_mapping_reference_umap_v1_20261001"
peak_mapping="${storage}/scmulti_flt3ko_il7ko/prepared/bdev_peak_quantification_v1_20260930/accepted_bdev_peak_mapping.tsv.gz"
fragments="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/datasets/scmulti_flt3ko_il7ko/FLT3LKO_IL7KO_BM/outs/atac_fragments.tsv.gz"
cell_qc="${storage}/scmulti_flt3ko_il7ko/reports/unfiltered_qc_v1_20260929/cell_qc_metrics.tsv.gz"
output_dir="${1:-${storage}/bdev_integrated/reports/ebf1_enhancer_il7_peak2_atac_tracks_v1_20261001}"

if (( $# > 1 )); then echo "Usage: $0 [OUTPUT_DIR]" >&2; exit 1; fi
for path in "${sif}" "${da_report}" "${pt_report}" "${projection_report}" "${peak_mapping}" \
  "${fragments}" "${fragments}.tbi" "${cell_qc}"; do
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
mkdir -p "$(dirname "${output_dir}")"
ulimit -v $((64 * 1024 * 1024))
"${runtime}" exec --cleanenv --no-home \
  -B "${repo_root}:/repo:ro" -B "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4 \
  --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/14_ebf1_peak2_atac_tracks.R \
  "${da_report}" "${pt_report}" "${projection_report}" "${peak_mapping}" "${fragments}" \
  "${cell_qc}" "${output_dir}" "v1_20261001"

printf 'Ebf1 peak-2 ATAC track report: %s\n' "${output_dir}/REPORT.html"
