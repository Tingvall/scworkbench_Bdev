# Resolve the active Parse bundle and optional isolated R library from the
# central registry. Source this after REPO_ROOT is defined in a Parse runner.
# Host dependency: Python 3 with PyYAML (also used for project configuration).
configure_parse_runtime() {
  local settings
  settings="$(python3 - "${REPO_ROOT}/config/bundles.yml" <<'PY'
import sys
from pathlib import Path
import yaml
entry = yaml.safe_load(Path(sys.argv[1]).read_text())['bundles']['parse_cytokines']
print(entry['bundle_dir'])
print((entry.get('runtime') or {}).get('r_library') or '')
PY
  )" || return
  local -a values
  mapfile -t values <<< "${settings}"
  BUNDLE_DIR="${1:-${values[0]}}"
  SCWORKBENCH_RUNTIME_ARGS=()
  local library="${SCWORKBENCH_R_LIBRARY:-${values[1]:-}}"
  if [[ -n "${library}" ]]; then
    [[ -d "${library}/scworkbench" ]] || { echo "Registered scWorkbench runtime is missing: ${library}" >&2; return 1; }
    SCWORKBENCH_RUNTIME_ARGS=(--bind "${library}:${library}:ro" --env "R_LIBS_USER=${library}")
  fi
}
