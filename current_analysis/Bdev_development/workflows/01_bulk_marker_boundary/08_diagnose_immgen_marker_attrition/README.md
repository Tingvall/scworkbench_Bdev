# ImmGen marker attrition diagnostic

This workflow reads the frozen outputs of step 07 and reports where genes are
lost for each marker-set target. It does not refit a model, change an approved
threshold, use Bdev single-cell data, or create a gene-set import artifact.

Run from the repository root:

```bash
python current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/08_diagnose_immgen_marker_attrition/01_diagnose_marker_attrition.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/08_immgen_marker_attrition_diagnostic.yaml
```

Generated tables, the machine-readable report, and the closure manifest are
written to the configured paired-storage workflow path.

