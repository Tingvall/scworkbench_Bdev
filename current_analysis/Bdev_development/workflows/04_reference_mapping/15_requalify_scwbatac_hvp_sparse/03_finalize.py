#!/usr/bin/env python3
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    config_path = parser.parse_args().config.resolve()
    config = yaml.safe_load(config_path.read_text())
    root = Path(config["storage"]["results_root"])
    report = root / "03_reports" / "01_hvp_sparse_requalification_summary.json"
    execution = root / "06_manifests" / "02_execution_receipt.json"
    if not report.exists() or not execution.exists():
        raise RuntimeError("qualification outputs are incomplete")
    summary = json.loads(report.read_text())
    if not summary.get("gate_passed"):
        raise RuntimeError("qualification gate did not pass")
    inventory_path = root / "06_manifests" / "03_run_inventory.json"
    files = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path != inventory_path:
            files.append({
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path)
            })
    inventory = {
        "schema_version": "bdev.scwbatac_hvp_sparse_inventory.v1",
        "status": "complete_passed",
        "run_identity": config["run_identity"],
        "config_sha256": sha256(config_path),
        "implementation_unlocked": False,
        "sif_built": False,
        "files": files
    }
    inventory_path.write_text(
        json.dumps(inventory, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
