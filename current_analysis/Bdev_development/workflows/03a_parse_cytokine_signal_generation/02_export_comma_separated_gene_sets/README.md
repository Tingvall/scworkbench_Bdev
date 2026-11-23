# 02 Export comma-separated Parse38 gene sets

This derived development workflow converts the frozen, checksummed Parse38 GMT
into a tab-separated table with one row per set and a comma-separated gene list
in `genes_comma_separated`. It preserves source set order, gene order, names,
descriptions, and capitalization. It does not alter the completed reconstruction
run or select a subset of cytokine programs.

Run from the repository root:

```bash
python3 current_analysis/Bdev_development/workflows/03a_parse_cytokine_signal_generation/02_export_comma_separated_gene_sets/01_export_comma_separated_gene_sets.py \
  --gmt /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/03a_parse_cytokine_signal_generation/01_recreate_all_cytokine_genesets/parse38_standalone_reconstruction_20261103T052911Z/01_tables/05_passing_gene_sets.gmt \
  --output /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/Bdev_development/results/03a_parse_cytokine_signal_generation/02_export_comma_separated_gene_sets/parse38_comma_separated_20261118T180640Z \
  --run-identity parse38_comma_separated_20261118T180640Z
```

The workflow requires the pinned GMT checksum and validates 38 sets and 934
memberships before writing the export and its manifest.
