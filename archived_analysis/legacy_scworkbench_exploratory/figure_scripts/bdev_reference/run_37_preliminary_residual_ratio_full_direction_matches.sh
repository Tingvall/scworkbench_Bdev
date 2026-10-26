#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
BDEV_ROOT="${STORAGE}/bdev_reference"
SIF="${SCWORKBENCH_IMAGE:-${STORAGE}/scworkbench-bdev-signac.sif}"
OUTPUT_DIR="${1:-${BDEV_ROOT}/reports/preliminary_residual_ratio_full_direction_matches_$(date -u +%Y%m%dT%H%M%SZ)}"
FINGERPRINT="$(python3 - "${REPO_ROOT}/config/bdev_reference.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle)['validation']['bundle_state_fingerprint'])
PY
)"
RUNTIME="/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer"
[[ -x "${RUNTIME}" ]] || { echo "Apptainer missing: ${RUNTIME}" >&2; exit 1; }
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output exists: ${OUTPUT_DIR}" >&2; exit 1; }

exec "${RUNTIME}" exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:/repo:ro" \
  --bind "${BDEV_ROOT}:${BDEV_ROOT}:rw" --pwd /repo "${SIF}" \
  Rscript --vanilla \
  /repo/figure_scripts/bdev_reference/37_preliminary_residual_ratio_full_direction_matches.R \
  "${BDEV_ROOT}/bundle" \
  "${BDEV_ROOT}/reports/rna_atac_latent_progress_w30_s5_v1_20261008" \
  "${BDEV_ROOT}/reports/preliminary_ratio_total_signal_explained_v1_20261009" \
  "${BDEV_ROOT}/reports/preliminary_residual_ratio_direction_gene_screen_v4_20261009" \
  "${OUTPUT_DIR}" "${FINGERPRINT}"
