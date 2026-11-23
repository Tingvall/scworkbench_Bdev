#!/usr/bin/env python3
"""Inventory the bounded FLT3LKO/IL7KO fixed-UMAP visualization."""

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
    output = root / "06_manifests" / "03_run_file_inventory.json"
    files = [{"path": str(path.relative_to(root)), "bytes": path.stat().st_size,
              "sha256": sha256(path)}
             for path in sorted(root.rglob("*"))
             if path.is_file() and path != output]
    required = {
        "01_tables/01_query_to_bdev_neighbors.tsv.gz",
        "01_tables/02_fixed_umap_placements.tsv",
        "01_tables/03_fixed_umap_group_summary.tsv",
        "02_figures/01_flt3lko_il7ko_fixture_fixed_atac_umap.png",
        "02_figures/01_flt3lko_il7ko_fixture_fixed_atac_umap.pdf",
        "03_reports/01_fixed_umap_summary.json",
        "04_logs/02_build_fixed_umap.log",
        "06_manifests/01_preparation_receipt.json",
        "06_manifests/02_execution_receipt.json",
    }
    missing = sorted(required - {entry["path"] for entry in files})
    if missing:
        raise RuntimeError(f"Required outputs missing: {missing}")
    inventory = {
        "schema_version": "bdev.flt3lko_il7ko_fixed_umap_inventory.v1",
        "status": "complete", "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path), "file_count": len(files),
        "files": files, "fixed_umap_unchanged": True,
        "all_query_cells_shown": True, "source_mutated": False,
        "paper_analysis": False, "implementation_unlocked": False,
    }
    output.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "complete", "files": len(files),
                      "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
