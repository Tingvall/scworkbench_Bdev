#!/usr/bin/env python3
"""Write a checksummed inventory for the completed four-model comparison."""

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


def load_config(path: Path) -> dict:
    current = yaml.safe_load(path.read_text())
    base_name = current.pop("base_config", None)
    if base_name is None:
        return current
    base = load_config((path.parent / base_name).resolve())
    for key, value in current.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            base[key].update(value)
        else:
            base[key] = value
    return base


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = load_config(config_path)
    root = Path(config["storage"]["results_root"])
    summary_path = root / "03_reports" / "01_external_four_model_benchmark_summary.json"
    receipt_path = root / "06_manifests" / "03_external_execution_receipt.json"
    if json.loads(summary_path.read_text())["status"] != "complete_external_comparison_no_automatic_promotion":
        raise RuntimeError("external comparison summary is not complete")
    if json.loads(receipt_path.read_text())["status"] != "complete":
        raise RuntimeError("external execution receipt is not complete")

    inventory_path = root / "06_manifests" / "04_run_file_inventory.json"
    records = []
    for path in sorted(p for p in root.rglob("*") if p.is_file() and p != inventory_path):
        records.append({
            "path": str(path.relative_to(root)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    output = {
        "schema_version": "bdev.frozen_peakvi_immgen_four_model_inventory.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "external_config_path": str(config_path),
        "external_config_sha256": sha256(config_path),
        "files": records,
    }
    inventory_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
