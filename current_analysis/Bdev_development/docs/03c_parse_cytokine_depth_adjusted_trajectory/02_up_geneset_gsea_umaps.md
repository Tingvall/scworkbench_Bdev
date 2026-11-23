# Up-signature local-GSEA UMAPs

Run `bdev_cytokine_up_gsea_umaps_20261106T210417Z` adds a matched ATAC-UMAP
view for every one of the 21 Parse38 `up` signatures from the completed 03c
qualification. This is a separate post-qualification visualization run; it did
not alter the parent run or promote any output to implementation.

Each map uses all Bdev cells as gray context. Color is restricted to cells in
the authenticated 1,810-cell fitted-L4 fixture that fall inside a supported
merged local-GSEA interval. Within each interval, a cell receives the NES of
the nearest constituent supported-window center. Every panel uses the same
diverging local-NES scale from -3 to +3.

Ten signatures have at least one supported interval. The other eleven retain a
complete gray context map with an explicit `No supported local-GSEA region`
annotation rather than being omitted. The figure index records the supported
window and interval counts, colored-cell count, observed NES range, and PNG/PDF
filename for every signature.

The result root is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/03c_parse_cytokine_depth_adjusted_trajectory/02_plot_up_geneset_gsea_umaps/bdev_cytokine_up_gsea_umaps_20261106T210417Z/
```

The final inventory contains 47 files, including 21 PNG/PDF figure pairs. An
independent audit reproduced every recorded size and SHA-256. The inventory is
`06_manifests/02_run_file_inventory.json`, with SHA-256
`0b5ec2e7eab8e622f222a8c879aca5eeb76381a578568c15eb298ee502fd8e38`.

## Combined supported-signature figure

Run `bdev_cytokine_supported_up_gsea_grid_20261106T211733Z` combines the ten
`up` signatures having at least one supported interval into one 2-by-5 panel.
It reuses the same gray background, colored-cell assignments, and shared NES
scale as the individual maps. Its result root is:

```text
/mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/03c_parse_cytokine_depth_adjusted_trajectory/03_combine_supported_up_geneset_gsea_umaps/bdev_cytokine_supported_up_gsea_grid_20261106T211733Z/
```

The six inventoried outputs rehashed successfully. The combined-run inventory
is `06_manifests/02_run_file_inventory.json`, with SHA-256
`132dbb7e8aeeb4f4c74dbc19ae21f98854cc3e9bb3402b258ab9f5923992585a`.
