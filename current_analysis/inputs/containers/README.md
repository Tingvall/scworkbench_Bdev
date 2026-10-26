# scWorkbench containers

Container binaries live in the paired mirrored-storage directory, not in Git:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/
```

After the updated image has finished transferring, record its final filename,
cryptographic digest, scworkbench version or source revision, architecture, and
documentation entry point in `../../config/`. Do not configure workflows
against temporary or partially transferred filenames.
