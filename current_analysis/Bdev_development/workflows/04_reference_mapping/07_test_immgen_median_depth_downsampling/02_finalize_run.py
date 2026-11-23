#!/usr/bin/env python3
"""Inventory the median-depth ImmGen pseudo-cell diagnostic."""

import argparse
import hashlib
import json
from pathlib import Path

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    root = Path(config["storage"]["results_root"])
    output = root / "06_manifests" / "02_run_file_inventory.json"
    files = [{"path": str(path.relative_to(root)), "bytes": path.stat().st_size,
              "sha256": sha256(path)}
             for path in sorted(root.rglob("*")) if path.is_file() and path != output]
    required = {
        "01_tables/01_pseudocell_design.tsv",
        "01_tables/02_pseudocell_neighbor_ledger.tsv",
        "01_tables/03_pseudocell_l4_predictions.tsv",
        "01_tables/04_library_consensus_neighbors.tsv",
        "01_tables/05_library_downsample_comparison.tsv",
        "01_tables/06_neighbor_reuse_comparison.tsv",
        "02_figures/01_population_downsample_l4.png",
        "02_figures/01_population_downsample_l4.pdf",
        "02_figures/02_full_vs_downsample_l4.png",
        "02_figures/02_full_vs_downsample_l4.pdf",
        "02_figures/03_neighbor_reuse_comparison.png",
        "02_figures/03_neighbor_reuse_comparison.pdf",
        "03_reports/01_median_depth_test_summary.json",
        "04_logs/01_median_depth_test.log",
        "06_manifests/01_execution_receipt.json",
    }
    missing = sorted(required - {item["path"] for item in files})
    if missing:
        raise RuntimeError(f"Required outputs missing: {missing}")
    inventory = {
        "schema_version": "bdev.immgen_median_depth_test_inventory.v1",
        "status": "complete_diagnostic", "run_identity": config["run_identity"],
        "config_path": str(config_path), "config_sha256": sha256(config_path),
        "file_count": len(files), "files": files, "diagnostic_only": True,
        "source_mutated": False, "paper_analysis": False,
        "implementation_unlocked": False, "sif_built": False,
    }
    output.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": inventory["status"], "files": len(files),
                      "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()

