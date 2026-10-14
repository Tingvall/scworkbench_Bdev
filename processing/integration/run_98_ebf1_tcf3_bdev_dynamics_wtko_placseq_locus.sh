#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
bundle="${storage}/bdev_integrated/ko_mapping_v1/bundle"
five="${storage}/bdev_integrated/reports/bdev_five_region_ebf1_dynamic_peaks_v3_20261013"
ratio="${storage}/bdev_integrated/reports/tcf3_peaks_detection_normalized_wt_ko_ratio_v2_20261013"
tcf3="${storage}/bdev_integrated/reports/ebf1_putative_enhancers_three_way_v2_20261005/all_259_bdev_ebf1_window_peaks_three_way_tests.tsv"
plac="${storage}/bdev_integrated/reports/tcf3_ebf1_coverage_placseq_tracks_v4_20261005/ebf1_promoter_placseq_links_in_evaluated_region.tsv"
output_dir="${1:-${storage}/bdev_integrated/reports/ebf1_tcf3_bdev_dynamics_wtko_placseq_locus_v1_20261013}"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
runtime="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"

bundle_fp="$(python3 - "${bundle}/manifest.json" <<'PY'
import json, sys
with open(sys.argv[1]) as handle:
    print(json.load(handle)['current_state']['fingerprint'])
PY
)"

"${runtime}" exec --cleanenv --no-home \
  --bind "${repo_root}:/repo:ro" \
  --bind "/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev:rw" \
  --pwd /repo "${sif}" Rscript --vanilla \
  /repo/processing/integration/98_ebf1_tcf3_bdev_dynamics_wtko_placseq_locus.R \
  "${bundle}" "${five}" "${ratio}" "${tcf3}" "${plac}" \
  "${output_dir}" "${bundle_fp}" "v1_20261013"

printf 'Integrated Ebf1 locus report: %s/REPORT.md\n' "${output_dir}"
