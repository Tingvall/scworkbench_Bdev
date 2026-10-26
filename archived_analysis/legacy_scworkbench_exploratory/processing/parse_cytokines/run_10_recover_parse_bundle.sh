#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
PARSE_ROOT="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/parse_cytokines"
IMAGE="${SCWORKBENCH_IMAGE:-/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects/scworkbench-bdev-signac.sif}"
RUN_ID="${1:?Supply a unique recovery run ID}"
PHASE="${2:-prepare}"
[[ "${RUN_ID}" =~ ^signature_recovery_[0-9]{8}T[0-9]{6}Z$ ]] || { echo 'Invalid run ID' >&2; exit 1; }
WORK="${PARSE_ROOT}/work/${RUN_ID}"
REPORT_ROOT="${PARSE_ROOT}/reports/${RUN_ID}"
RUNTIME="${3:-${WORK}/runtime}"
BUNDLE="${PARSE_ROOT}/bundle_recovered_${RUN_ID#signature_recovery_}"
ACCEPTED="${PARSE_ROOT}/reports/signature_report_20260921T221910Z/final_accepted_feature_sets.tsv"
if command -v apptainer >/dev/null 2>&1; then CONTAINER_RUNTIME=apptainer; else CONTAINER_RUNTIME=singularity; fi
ulimit -v 100663296
CONTAINER=("${CONTAINER_RUNTIME}" exec --cleanenv --no-home
  --bind "${REPO_ROOT}:${REPO_ROOT}:ro" --bind "${PARSE_ROOT}:${PARSE_ROOT}:rw"
  --bind "${PARSE_ROOT}/bundle:${PARSE_ROOT}/bundle:ro"
  --env OMP_NUM_THREADS=4 --env OPENBLAS_NUM_THREADS=4 --env MKL_NUM_THREADS=4)

case "${PHASE}" in
  prepare)
    [[ ! -e "${WORK}" && ! -e "${REPORT_ROOT}" ]] || { echo 'Recovery destination exists' >&2; exit 1; }
    mkdir -p "${WORK}" "${REPORT_ROOT}"
    "${CONTAINER[@]}" "${IMAGE}" bash "${SCRIPT_DIR}/recovery/build_runtime.sh" \
      "${REPO_ROOT}" "${RUNTIME}" "${REPORT_ROOT}/runtime_tests" \
      2>&1 | tee "${REPORT_ROOT}/runtime_build_and_tests.log"
    ;;
  recover)
    TEST_REPORT="${4:-${REPORT_ROOT}/runtime_tests}"
    [[ -f "${TEST_REPORT}/runtime_tests_passed" ]] || { echo 'Runtime tests have not passed' >&2; exit 1; }
    mkdir -p "${WORK}" "${REPORT_ROOT}"
    REPORT="${REPORT_ROOT}/finalization"
    [[ ! -e "${REPORT}" ]] || { echo 'Finalization report already exists; inspect its checkpoint before retrying' >&2; exit 1; }
    mkdir "${REPORT}"
    sha256sum "${IMAGE}" "${SCRIPT_DIR}/10_recover_parse_bundle.R" \
      "${SCRIPT_DIR}/recovery/copy_bundle.py" "${REPO_ROOT}/src/scworkbench_registry_json.R" \
      "${REPO_ROOT}/src/scworkbench_object_row_keys.R" > "${REPORT}/recovery_input_sha256.txt"
    cp "${SCRIPT_DIR}/10_recover_parse_bundle.R" "${REPORT}/executed_recovery.R"
    cp "${SCRIPT_DIR}/run_10_recover_parse_bundle.sh" "${REPORT}/executed_runner.sh"
    cp "${SCRIPT_DIR}/recovery/validate_recovered_bundle.R" "${REPORT}/executed_validator.R"
    cp "${SCRIPT_DIR}/recovery/restore_materializations.R" "${REPORT}/executed_materialization_restoration.R"
    python3 "${SCRIPT_DIR}/recovery/copy_bundle.py" copy "${PARSE_ROOT}/bundle" "${BUNDLE}" "${REPORT}"
    "${CONTAINER[@]}" --env "R_LIBS_USER=${RUNTIME}/library" "${IMAGE}" \
      Rscript --vanilla "${SCRIPT_DIR}/10_recover_parse_bundle.R" \
      "${PARSE_ROOT}/bundle" "${BUNDLE}" "${REPORT}" "${ACCEPTED}" \
      2>&1 | tee "${REPORT}/recovery.log"
    python3 "${SCRIPT_DIR}/recovery/copy_bundle.py" verify-source "${PARSE_ROOT}/bundle" "${BUNDLE}" "${REPORT}"
    ;;
  restore)
    REPORT="${REPORT_ROOT}/finalization"
    [[ -f "${REPORT}/scientific_state_checks.json" && ! -e "${REPORT}/materialization_restoration.json" ]] || { echo 'Completed registration and a new restoration report required' >&2; exit 1; }
    cp "${SCRIPT_DIR}/recovery/restore_materializations.R" "${REPORT}/executed_materialization_restoration.R"
    "${CONTAINER[@]}" --env "R_LIBS_USER=${RUNTIME}/library" "${IMAGE}" \
      Rscript --vanilla "${SCRIPT_DIR}/recovery/restore_materializations.R" \
      "${BUNDLE}" "${REPORT}" "${PARSE_ROOT}/sessions/parse_signature_ucell_20260921T152017" \
      2>&1 | tee "${REPORT}/materialization_restoration.log"
    ;;
  validate)
    TEST_REPORT="${4:-${REPORT_ROOT}/runtime_tests}"
    REPORT="${REPORT_ROOT}/finalization"
    [[ -f "${TEST_REPORT}/runtime_tests_passed" && -f "${REPORT}/scientific_state_checks.json" ]] || { echo 'Tested runtime and completed registration checks required' >&2; exit 1; }
    VALIDATION_LABEL="${5:-validation_retry}"
    [[ "${VALIDATION_LABEL}" =~ ^validation_[a-z0-9_]+$ ]] || { echo "Invalid validation label" >&2; exit 1; }
    [[ ! -e "${REPORT}/completion.json" && ! -e "${REPORT}/${VALIDATION_LABEL}.log" ]] || { echo 'Validation evidence already exists; inspect it before retrying' >&2; exit 1; }
    cp "${SCRIPT_DIR}/recovery/validate_recovered_bundle.R" "${REPORT}/${VALIDATION_LABEL}_validator.R"
    if [[ -f "${REPORT}/validation.json" ]]; then cp "${REPORT}/validation.json" "${REPORT}/${VALIDATION_LABEL}_prior_validation.json"; fi
    sha256sum "${SCRIPT_DIR}/recovery/validate_recovered_bundle.R" "${RUNTIME}/build_provenance.json" > "${REPORT}/${VALIDATION_LABEL}_input_sha256.txt"
    "${CONTAINER[@]}" --env "R_LIBS_USER=${RUNTIME}/library" "${IMAGE}" \
      Rscript --vanilla "${SCRIPT_DIR}/recovery/validate_recovered_bundle.R" \
      "${PARSE_ROOT}/bundle" "${BUNDLE}" "${REPORT}" \
      2>&1 | tee "${REPORT}/${VALIDATION_LABEL}.log"
    python3 "${SCRIPT_DIR}/recovery/copy_bundle.py" verify-source "${PARSE_ROOT}/bundle" "${BUNDLE}" "${REPORT}"
    ;;
  *) echo "Unknown recovery phase: ${PHASE}" >&2; exit 1 ;;
esac
