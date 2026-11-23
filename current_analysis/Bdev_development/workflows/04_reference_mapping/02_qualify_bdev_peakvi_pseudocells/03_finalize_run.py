#!/usr/bin/env python3
"""Inventory and rehash the bounded PeakVI pseudocell qualification."""

from __future__ import annotations

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


def merge(base: dict, override: dict) -> dict:
    result = dict(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(path: Path) -> dict:
    current = yaml.safe_load(path.read_text())
    base_name = current.pop("base_config", None)
    if base_name is None:
        return current
    return merge(load_config((path.parent / base_name).resolve()), current)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = load_config(config_path)
    root = Path(config["storage"]["results_root"])
    output = root / "06_manifests" / "04_run_file_inventory.json"
    files = [
        {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(root.rglob("*")) if path.is_file() and path != output
    ]
    required = {
        "01_tables/01_reference_and_holdout_cells.tsv",
        "01_tables/02_synthetic_query_design.tsv",
        "01_tables/03_query_transfer_predictions.tsv",
        "01_tables/04_metrics_by_query_type.tsv",
        "03_reports/01_peakvi_pseudocell_qualification_summary.json",
        "06_manifests/01_runtime_install.json",
        "06_manifests/02_execution_receipt.json",
        "06_manifests/03_bundle_inventory.json",
    }
    missing = sorted(required - {entry["path"] for entry in files})
    if missing:
        raise RuntimeError(f"Required outputs are missing: {missing}")
    manifest = {
        "schema_version": "bdev.peakvi_pseudocell_run_inventory.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "files": files,
        "file_count": len(files),
        "paper_analysis": False,
        "real_immgen_projection": False,
        "implementation_unlocked": False,
        "sif_built": False,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "complete", "files": len(files), "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
