# Simple mean-contrast bootstrap stability

This diagnostic workflow resamples ImmGen libraries with replacement within
each of the 64 populations for 30 seeded rounds. Each round regenerates the
population-balanced mean-contrast top-20 sets and applies the selected native
scWorkbench raw-UCell argmax rule with margin 0.01 to the fixed 1,746-cell Bdev
fixture. Bdev is used only to measure transfer stability and never to select
genes. For efficiency, all round-specific signatures are scored in one native
UCell operation. The operation's documented argmax, tie, score, and margin
contract is then applied to each round in memory to avoid repeatedly loading the
195 MB score artifact. Its round-1 output must exactly match a completed native
reducer operation. UCell scores each signature independently against the same
fixed cell-rank matrix, so the batched score is equivalent to separate passes.

Singleton populations are explicitly marked structurally unassessable. No
artifacts or labels are promoted, the selected default is not changed, and the
source bundle fingerprint must remain unchanged.

```bash
python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/28_simple_mean_contrast_bootstrap_stability/01_generate_bootstrap_sets.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/29_simple_mean_contrast_bootstrap_stability.yaml

/home/jtg/.local/opt/apptainer/1.5.4/bin/apptainer exec \
  --bind /mnt/sata2:/mnt/sata2 \
  /mnt/sata2/jtg_storage/msi_projects/scworkbench_Bdev/current_analysis/inputs/containers/scworkbench-generic-collaborator-amd64-20261029.sif \
  Rscript current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/28_simple_mean_contrast_bootstrap_stability/02_score_bootstrap_sets.R \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/29_simple_mean_contrast_bootstrap_stability.yaml

python3 current_analysis/Bdev_development/workflows/01_bulk_marker_boundary/28_simple_mean_contrast_bootstrap_stability/03_summarize_and_plot.py \
  --config current_analysis/Bdev_development/config/01_bulk_marker_boundary/29_simple_mean_contrast_bootstrap_stability.yaml
```
