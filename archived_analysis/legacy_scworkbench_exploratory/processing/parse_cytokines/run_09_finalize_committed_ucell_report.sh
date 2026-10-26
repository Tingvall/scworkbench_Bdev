#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

IMAGE="${SCWORKBENCH_IMAGE:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif}"
PARSE_ROOT="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines"
source "${REPO_ROOT}/src/parse_runtime.sh"
configure_parse_runtime "${1:-}"
DISCOVERY_DIR="${2:-${PARSE_ROOT}/reports/signature_discovery_20260921T221102Z}"
REPORT_DIR="${3:-${PARSE_ROOT}/reports/signature_report_20260921T221910Z}"
SESSION_ROOT="${4:-${PARSE_ROOT}/sessions}"
SCORE_ARTIFACT_ID="${5:-cell_vector_family_20260921_155957_5n7edbnw}"
SCORE_SUMMARY_ARTIFACT_ID="${6:-obs_table_20260921_164902_2mk1awzd}"
QC_KEY_ARTIFACT_ID="${7:-obs_table_20260921_164903_d0ndemln}"
QC_PAIRWISE_ARTIFACT_ID="${8:-obs_table_20260921_164903_vtghjfz9}"
COMMIT_ID="${9:-commit_20260921_155948_s32blwqt}"
PARSE_MAX_MEMORY_GIB="${PARSE_MAX_MEMORY_GIB:-96}"
SCWORKBENCH_CREATED_BY="${SCWORKBENCH_CREATED_BY:-jtg}"

if [[ ! "${PARSE_MAX_MEMORY_GIB}" =~ ^[1-9][0-9]*$ ]]; then
  echo "PARSE_MAX_MEMORY_GIB must be a positive integer." >&2
  exit 1
fi
TOTAL_MEMORY_KIB="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
MEMORY_LIMIT_KIB="$((PARSE_MAX_MEMORY_GIB * 1024 * 1024))"
MAX_SHARED_SERVER_KIB="$((TOTAL_MEMORY_KIB * 40 / 100))"
if (( MEMORY_LIMIT_KIB > MAX_SHARED_SERVER_KIB )); then
  echo "Requested memory ceiling exceeds 40% of server RAM; refusing to start." >&2
  exit 1
fi
if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="singularity"
else
  echo "Apptainer or Singularity is required." >&2
  exit 1
fi
if [[ -e "${REPORT_DIR}" ]]; then
  echo "Report destination already exists; refusing to overwrite: ${REPORT_DIR}" >&2
  exit 1
fi

ulimit -v "${MEMORY_LIMIT_KIB}"

exec "${CONTAINER_RUNTIME}" exec --cleanenv --no-home \
  "${SCWORKBENCH_RUNTIME_ARGS[@]}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" \
  --bind "${PARSE_ROOT}:${PARSE_ROOT}:rw" \
  --env "SCWORKBENCH_CREATED_BY=${SCWORKBENCH_CREATED_BY}" \
  --pwd "${REPO_ROOT}" \
  "${IMAGE}" \
  Rscript --vanilla \
  "${SCRIPT_DIR}/09_finalize_committed_ucell_report.R" \
  "${BUNDLE_DIR}" \
  "${DISCOVERY_DIR}" \
  "${REPORT_DIR}" \
  "${SESSION_ROOT}" \
  "${SCORE_ARTIFACT_ID}" \
  "${SCORE_SUMMARY_ARTIFACT_ID}" \
  "${QC_KEY_ARTIFACT_ID}" \
  "${QC_PAIRWISE_ARTIFACT_ID}" \
  "${COMMIT_ID}"
