#!/usr/bin/env bash
set -euo pipefail

repo=/home/jtg/msi_projects/scworkbench_Bdev
config="$repo/current_analysis/Bdev_development/config/05_trajectory_group_contrasts/04_flt3lko_il7ko_full_mapping_fixed_umap.yaml"
sif=/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif
library=/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/intermediate_objects/05_trajectory_group_contrasts/04_map_flt3lko_il7ko_full/flt3lko_il7ko_full_mapping_fixed_umap_20261118T022951Z/05_objects/01_r_library
apptainer=/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer

test -f "$config"
test ! -e "$library"
mkdir -p "$library"
"$apptainer" exec --bind /mnt/sata2:/mnt/sata2 "$sif" \
  R CMD INSTALL -l "$library" /home/jtg/msi_projects/scworkbench_code/scworkbench
"$apptainer" exec --bind /mnt/sata2:/mnt/sata2 --env R_LIBS_USER="$library" "$sif" \
  R CMD INSTALL -l "$library" /home/jtg/msi_projects/scworkbench_code/scwbRepresentation
"$apptainer" exec --bind /mnt/sata2:/mnt/sata2 --env R_LIBS_USER="$library" "$sif" \
  R CMD INSTALL -l "$library" /home/jtg/msi_projects/scworkbench_code/scwbLatent
"$apptainer" exec --bind /mnt/sata2:/mnt/sata2 --env R_LIBS_USER="$library" "$sif" \
  R CMD INSTALL -l "$library" /home/jtg/msi_projects/scworkbench_code/scwbReferenceMap

"$apptainer" exec --bind /mnt/sata2:/mnt/sata2 --env R_LIBS_USER="$library" "$sif" \
  Rscript -e 'stopifnot(as.character(packageVersion("scworkbench")) == "0.1.2", as.character(packageVersion("scwbRepresentation")) == "0.1.1", as.character(packageVersion("scwbLatent")) == "0.1.3", as.character(packageVersion("scwbReferenceMap")) == "0.1.0"); writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("R_LIBS_USER"), "RUNTIME_SESSION_INFO.txt"))'
