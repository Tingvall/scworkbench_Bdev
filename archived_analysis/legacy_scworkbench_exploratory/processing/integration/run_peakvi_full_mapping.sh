#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
library="${storage}/parse_cytokines/work/signature_recovery_20260923T200125Z/runtime_v4/library"
phase_root="${storage}/bdev_integrated/ko_mapping_v1"
bundle_dir="${phase_root}/bundle"
report_dir="${storage}/bdev_integrated/reports/ko_mapping_peakvi_full_v1_20260930"
session_root="${phase_root}/sessions"
monitor_dir="${phase_root}/monitoring/peakvi_full_mapping_v1"

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
else
  container_runtime="singularity"
fi

[[ -d "${library}/scworkbench" ]] || { echo "Registered scWorkbench runtime missing" >&2; exit 1; }
[[ ! -e "${report_dir}" ]] || { echo "Full mapping report already exists: ${report_dir}" >&2; exit 1; }
mkdir -p "${monitor_dir}" "${storage}/bdev_integrated/reports"

metrics_log="${monitor_dir}/gpu_metrics.tsv"
events_log="${monitor_dir}/courtesy_events.log"
training_log="${monitor_dir}/training.log"
printf 'timestamp\tindex\tname\tmemory_total_mib\tmemory_used_mib\tmemory_free_mib\tgpu_util_percent\tmemory_util_percent\ttemperature_c\tpower_w\n' > "${metrics_log}"
printf '%s\tmonitor_started\n' "$(date --iso-8601=seconds)" > "${events_log}"

common=(
  exec --nv --cleanenv --no-home
  -B "${library}:${library}:ro"
  --env "R_LIBS_USER=${library}"
  --env "SCWORKBENCH_PYBENCH_PYTHON=/opt/pyenv/bin/python"
  -B "${repo_root}:/repo:ro"
  -B "${storage}:${storage}:rw"
  --env OMP_NUM_THREADS=4
  --env OPENBLAS_NUM_THREADS=4
  --env MKL_NUM_THREADS=4
  --pwd /repo
  "${sif}"
  Rscript --vanilla
  /repo/processing/integration/10_peakvi_full_mapping.R
  "${bundle_dir}" "${report_dir}" "${session_root}"
)

setsid "${container_runtime}" "${common[@]}" > "${training_log}" 2>&1 &
training_pid=$!
owner_uid="$(id -u)"
paused=0

cleanup() {
  if [[ "${paused}" -eq 1 ]]; then
    kill -CONT -- "-${training_pid}" 2>/dev/null || true
  fi
}
trap cleanup EXIT INT TERM

while kill -0 "${training_pid}" 2>/dev/null; do
  timestamp="$(date --iso-8601=seconds)"
  gpu_line="$(${container_runtime} exec --nv --cleanenv --no-home "${sif}" nvidia-smi \
    --query-gpu=index,name,memory.total,memory.used,memory.free,utilization.gpu,utilization.memory,temperature.gpu,power.draw \
    --format=csv,noheader,nounits 2>/dev/null || true)"
  if [[ -n "${gpu_line}" ]]; then
    printf '%s\t%s\n' "${timestamp}" "$(tr -d ' ' <<< "${gpu_line}" | tr ',' '\t')" >> "${metrics_log}"
  fi

  foreign=0
  while IFS=, read -r gpu_pid process_name used_memory; do
    gpu_pid="$(tr -d ' ' <<< "${gpu_pid}")"
    [[ "${gpu_pid}" =~ ^[0-9]+$ ]] || continue
    process_uid="$(ps -o uid= -p "${gpu_pid}" 2>/dev/null | tr -d ' ' || true)"
    if [[ -n "${process_uid}" && "${process_uid}" != "${owner_uid}" ]]; then
      foreign=1
      printf '%s\tforeign_gpu_process\tuid=%s\tpid=%s\tname=%s\tmemory=%s\n' \
        "${timestamp}" "${process_uid}" "${gpu_pid}" "${process_name}" "${used_memory}" >> "${events_log}"
    fi
  done < <(${container_runtime} exec --nv --cleanenv --no-home "${sif}" nvidia-smi \
    --query-compute-apps=pid,process_name,used_memory --format=csv,noheader,nounits 2>/dev/null || true)

  if [[ "${foreign}" -eq 1 && "${paused}" -eq 0 ]]; then
    kill -STOP -- "-${training_pid}"
    paused=1
    printf '%s\ttraining_paused_for_foreign_user\n' "${timestamp}" >> "${events_log}"
  elif [[ "${foreign}" -eq 0 && "${paused}" -eq 1 ]]; then
    kill -CONT -- "-${training_pid}"
    paused=0
    printf '%s\ttraining_resumed\n' "${timestamp}" >> "${events_log}"
  fi
  sleep 10
done

set +e
wait "${training_pid}"
status=$?
set -e
printf '%s\tmonitor_finished\texit_status=%s\n' "$(date --iso-8601=seconds)" "${status}" >> "${events_log}"
exit "${status}"
