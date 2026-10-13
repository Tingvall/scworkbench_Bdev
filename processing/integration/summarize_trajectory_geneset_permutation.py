#!/usr/bin/env python3

import csv
import json
import pathlib
import sys


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: summarize_trajectory_geneset_permutation.py OUTPUT_ROOT")
    root = pathlib.Path(sys.argv[1]).resolve()
    groups = ("cytokines", "bp", "mf", "cc")
    summaries = []
    combined = []
    fieldnames = None
    for group in groups:
        completion_path = root / group / "completion.json"
        result_path = root / group / "significant_gene_sets.tsv"
        if not completion_path.exists() or not result_path.exists():
            continue
        completion = json.loads(completion_path.read_text())
        summaries.append(completion)
        with result_path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames:
                fieldnames = fieldnames or reader.fieldnames
            combined.extend(reader)

    if not summaries:
        raise SystemExit("No completed group outputs found")
    if fieldnames:
        with (root / "significant_gene_sets_all_collections.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                                    extrasaction="ignore", quoting=csv.QUOTE_MINIMAL)
            writer.writeheader()
            writer.writerows(combined)

    labels = {"cytokines": "Cytokines", "bp": "GO Biological Process",
              "mf": "GO Molecular Function", "cc": "GO Cellular Component"}
    lines = ["# Trajectory gene-set permutation screen", "",
             "All screens used exact Bdev `pt`, `RNA_log1p_cp10k`, mean member-gene expression, 30 equal-width bins, no smoothing, and BH correction separately within each collection.", "",
             "| Collection | Tested | Significant (q < 0.05) | Permutations |",
             "|---|---:|---:|---:|"]
    for summary in summaries:
        lines.append(f"| {labels[summary['group']]} | {summary['tested_sets']} | "
                     f"{summary['significant_q_lt_0_05']} | {summary['n_permutations']} |")
    lines += ["", "A significant result is significant for overall trajectory dynamics. `peak_pseudotime` or `trough_pseudotime` localizes that pattern but is not a bin-specific hypothesis test.", "",
              "Full significant results are in `significant_gene_sets_all_collections.tsv`; each collection directory also contains all tested terms, curves, overlap QC, fit diagnostics, and its durable unpromoted scWorkbench session."]
    (root / "REPORT.md").write_text("\n".join(lines) + "\n")
    (root / "completion.json").write_text(json.dumps({
        "status": "complete", "groups": [x["group"] for x in summaries],
        "total_tested": sum(x["tested_sets"] for x in summaries),
        "total_significant_within_collection_q_lt_0_05": sum(
            x["significant_q_lt_0_05"] for x in summaries),
        "group_summaries": summaries,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
