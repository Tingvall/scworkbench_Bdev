#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
STORAGE="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
stage="${1:?Usage: run_rest_only.sh discover|prepare|score|report}"
case "${stage}" in
  discover) script="processing/parse_cytokines/11_discover_rest_only_signatures.R" ;;
  prepare) script="processing/parse_cytokines/12_prepare_rest_only_bundle.R" ;;
  score) script="processing/parse_cytokines/13_score_rest_only_signatures.R" ;;
  report) script="figure_scripts/parse_cytokines/06_rest_only_signal_report.R" ;;
  *) echo "Unknown stage: ${stage}" >&2; exit 2 ;;
esac
source "${REPO_ROOT}/src/parse_runtime.sh"
configure_parse_runtime
OUTPUT="$(python3 - "${REPO_ROOT}/config/parse_cytokines_rest_only.yml" <<'PY'
import sys,yaml
with open(sys.argv[1]) as h:print(yaml.safe_load(h)['storage_root'])
PY
)"
mkdir -p "${OUTPUT}"
ulimit -v $((48 * 1024 * 1024))
exec singularity exec --cleanenv --no-home "${SCWORKBENCH_RUNTIME_ARGS[@]}" \
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" --bind "${STORAGE}:${STORAGE}:ro" \
  --bind "${OUTPUT}:${OUTPUT}:rw" \
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4 \
  --pwd "${REPO_ROOT}" "${STORAGE}/scworkbench-bdev-signac.sif" \
  Rscript --vanilla "${REPO_ROOT}/${script}"
