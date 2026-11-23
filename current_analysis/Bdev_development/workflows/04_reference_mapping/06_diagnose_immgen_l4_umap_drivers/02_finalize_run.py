#!/usr/bin/env python3
"""Inventory the ImmGen L4 and fixed-UMAP driver diagnostic."""

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
    files = [
        {"path": str(path.relative_to(root)), "bytes": path.stat().st_size,
         "sha256": sha256(path)}
        for path in sorted(root.rglob("*")) if path.is_file() and path != output
    ]
    required = {
        "01_tables/01_neighbor_l4_umap_driver_ledger.tsv",
        "01_tables/02_library_l4_umap_reconstruction.tsv",
        "01_tables/03_population_l4_distribution_summary.tsv",
        "01_tables/04_population_reference_cluster_drivers.tsv",
        "02_figures/01_population_neighbor_l4_distributions.png",
        "02_figures/01_population_neighbor_l4_distributions.pdf",
        "02_figures/02_population_fixed_umap_drivers.png",
        "02_figures/02_population_fixed_umap_drivers.pdf",
        "03_reports/01_diagnostic_summary.json",
        "04_logs/01_build_diagnostic.log",
        "06_manifests/01_execution_receipt.json",
    }
    missing = sorted(required - {entry["path"] for entry in files})
    if missing:
        raise RuntimeError(f"Required outputs missing: {missing}")
    inventory = {
        "schema_version": "bdev.immgen_l4_umap_driver_inventory.v1",
        "status": "complete_diagnostic",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "file_count": len(files),
        "files": files,
        "diagnostic_only": True,
        "source_mutated": False,
        "paper_analysis": False,
        "implementation_unlocked": False,
        "sif_built": False,
    }
    output.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": inventory["status"], "files": len(files),
                      "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()

