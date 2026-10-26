#!/usr/bin/env bash
# Run inside the supplied image. All destinations are new, separate directories.
set -euo pipefail
repo="${1:?Repository root required}"
runtime="${2:?New runtime directory required}"
report="${3:?New runtime test report directory required}"
[[ ! -e "${report}" ]] || { echo 'Runtime test report already exists' >&2; exit 1; }
mkdir -p "${report}"
python3 "${repo}/processing/parse_cytokines/recovery/prepare_runtime.py" "${repo}" "${runtime}"
R CMD INSTALL --library="${runtime}/library" "${runtime}/source"
export R_LIBS_USER="${runtime}/library"
Rscript --vanilla "${repo}/processing/parse_cytokines/recovery/test_registry_json.R"
Rscript --vanilla "${repo}/processing/parse_cytokines/recovery/test_object_row_keys.R"
Rscript --vanilla "${repo}/processing/parse_cytokines/recovery/test_feature_set_csv_registration.R"
Rscript --vanilla "${repo}/processing/parse_cytokines/recovery/test_artifact_listing.R"
Rscript --vanilla "${repo}/processing/parse_cytokines/recovery/test_engine_recovery.R" "${runtime}/source" "${report}"
touch "${report}/runtime_tests_passed"
