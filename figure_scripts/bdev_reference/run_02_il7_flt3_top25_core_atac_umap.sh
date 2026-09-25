#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
RUN_ID="${1:-il7_flt3_top25_core_atac_umap_20260925}"
[[ "${RUN_ID}" =~ ^[A-Za-z0-9_]+$ ]] || exit 2

STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
SIF="${STORAGE}/scworkbench-bdev-signac.sif"
BUNDLE="${STORAGE}/bdev_reference/bundle"
SOURCE_DIR="${STORAGE}/parse_cytokines_rest_only/reports/rest_only_20260924/competitor_aware_bdev25_dual_v1"
OUTPUT_DIR="${STORAGE}/bdev_reference/reports/${RUN_ID}"
SCRIPT="${REPO_ROOT}/figure_scripts/bdev_reference/02_il7_flt3_top25_core_atac_umap.R"

EXPECTED_FINGERPRINT="$(python3 - "${REPO_ROOT}/config/bdev_reference.yml" <<'PY'
import sys, yaml
with open(sys.argv[1]) as handle:
    print(yaml.safe_load(handle)['validation']['bundle_state_fingerprint'])
PY
)"

[[ -f "${SIF}" && -d "${BUNDLE}" && -f "${SOURCE_DIR}/top25_core_gene_membership.tsv" ]]
[[ ! -e "${OUTPUT_DIR}" ]] || { echo "Output already exists: ${OUTPUT_DIR}" >&2; exit 1; }
mkdir -p "${OUTPUT_DIR}"

ulimit -v $((64 * 1024 * 1024))
exec singularity exec --cleanenv --no-home \
  --bind "${REPO_ROOT}:/repo:ro" \
  --bind "${BUNDLE}:/bundle:ro" \
  --bind "${SOURCE_DIR}:/source:ro" \
  --bind "${OUTPUT_DIR}:/output:rw" \
  --pwd /repo "${SIF}" \
  Rscript --vanilla \
  "/repo/figure_scripts/bdev_reference/02_il7_flt3_top25_core_atac_umap.R" \
  /bundle /source/top25_core_gene_membership.tsv /output \
  "${EXPECTED_FINGERPRINT}"
