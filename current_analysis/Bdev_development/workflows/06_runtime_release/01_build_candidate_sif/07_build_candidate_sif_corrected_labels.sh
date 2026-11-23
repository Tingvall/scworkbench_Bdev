#!/usr/bin/env bash
set -euo pipefail

repo_root=/home/jtg/msi_projects/scworkbench_Bdev
definition="$repo_root/current_analysis/Bdev_development/workflows/06_runtime_release/01_build_candidate_sif/07_bdev_candidate_sif_corrected_labels.def"
output=/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/containers/06_runtime_release/01_build_candidate_sif/bdev_candidate_sif_corrected_labels_20261118T022000Z/05_objects/01_scworkbench_bdev_candidate_amd64_20261118T022000Z.sif
log=/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/containers/06_runtime_release/01_build_candidate_sif/bdev_candidate_sif_corrected_labels_20261118T022000Z/04_logs/01_apptainer_build.log
apptainer=/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer

test ! -e "$output"
mkdir -p "$(dirname "$output")" "$(dirname "$log")"
build_tmp=$(mktemp -d /tmp/bdev-candidate-sif-build.XXXXXX)
trap 'rm -rf "$build_tmp"' EXIT
APPTAINER_TMPDIR="$build_tmp" "$apptainer" build "$output" "$definition" 2>&1 | tee "$log"

