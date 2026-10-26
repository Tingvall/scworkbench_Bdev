#!/usr/bin/env bash
set -euo pipefail

# Cooperative controller for an already-running full PeakVI mapping process.
# While another user's CUDA process is present, allow this analysis to run for
# 15 seconds and yield for 45 seconds. With no foreign process, keep it running.

storage="/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/intermediate_objects"
sif="${SCWORKBENCH_IMAGE:-${storage}/scworkbench-bdev-signac.sif}"
monitor_dir="${storage}/bdev_integrated/ko_mapping_v1/monitoring/peakvi_full_mapping_v1"
events_log="${monitor_dir}/low_duty_cycle_events.log"
lock_dir="${monitor_dir}/.low_duty_controller.lock"
owner_uid="$(id -u)"

if command -v apptainer >/dev/null 2>&1; then
  container_runtime="apptainer"
else
  container_runtime="singularity"
fi

mkdir -p "${monitor_dir}"
if ! mkdir "${lock_dir}" 2>/dev/null; then
  echo "A low-duty controller is already active: ${lock_dir}" >&2
  exit 1
fi
trap 'rmdir "${lock_dir}" 2>/dev/null || true' EXIT

target_pid="$(pgrep -u "${owner_uid}" -f '/repo/processing/integration/10_peakvi_full_mapping.R' | head -n 1 || true)"
if [[ ! "${target_pid}" =~ ^[0-9]+$ ]]; then
  echo "No active full PeakVI mapping R process found" >&2
  exit 1
fi
target_pgid="$(ps -o pgid= -p "${target_pid}" | tr -d ' ')"
if [[ ! "${target_pgid}" =~ ^[0-9]+$ ]]; then
  echo "Could not resolve the PeakVI process group" >&2
  exit 1
fi

foreign_gpu_present() {
  local gpu_pid process_name used_memory process_uid
  while IFS=, read -r gpu_pid process_name used_memory; do
    gpu_pid="$(tr -d ' ' <<< "${gpu_pid}")"
    [[ "${gpu_pid}" =~ ^[0-9]+$ ]] || continue
    process_uid="$(ps -o uid= -p "${gpu_pid}" 2>/dev/null | tr -d ' ' || true)"
    if [[ -n "${process_uid}" && "${process_uid}" != "${owner_uid}" ]]; then
      return 0
    fi
  done < <(${container_runtime} exec --nv --cleanenv --no-home "${sif}" nvidia-smi \
    --query-compute-apps=pid,process_name,used_memory --format=csv,noheader,nounits 2>/dev/null || true)
  return 1
}

restore_expected_state() {
  if foreign_gpu_present; then
    kill -STOP -- "-${target_pgid}" 2>/dev/null || true
  else
    kill -CONT -- "-${target_pgid}" 2>/dev/null || true
  fi
}
trap restore_expected_state INT TERM

printf '%s\tcontroller_started\ttarget_pid=%s\ttarget_pgid=%s\trun_seconds=15\tyield_seconds=45\n' \
  "$(date --iso-8601=seconds)" "${target_pid}" "${target_pgid}" >> "${events_log}"

while kill -0 "${target_pid}" 2>/dev/null; do
  if foreign_gpu_present; then
    kill -CONT -- "-${target_pgid}"
    printf '%s\tshared_run_window_started\n' "$(date --iso-8601=seconds)" >> "${events_log}"
    sleep 15
    if kill -0 "${target_pid}" 2>/dev/null; then
      kill -STOP -- "-${target_pgid}"
      printf '%s\tshared_yield_window_started\n' "$(date --iso-8601=seconds)" >> "${events_log}"
      sleep 45
    fi
  else
    kill -CONT -- "-${target_pgid}" 2>/dev/null || true
    sleep 10
  fi
done

printf '%s\tcontroller_finished\n' "$(date --iso-8601=seconds)" >> "${events_log}"
