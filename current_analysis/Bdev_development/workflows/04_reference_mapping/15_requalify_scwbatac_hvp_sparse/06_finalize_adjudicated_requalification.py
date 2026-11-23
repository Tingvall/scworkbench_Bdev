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
    root = Path(config["execution"]["result_root"])
    report = root / "03_reports" / config["outputs"]["report"]
    receipt = root / "06_manifests" / config["outputs"]["receipt"]
    for path in (report, receipt):
        if not path.exists():
            raise RuntimeError(f"missing adjudication artifact: {path}")
        if not json.loads(path.read_text()).get("gate_passed"):
            raise RuntimeError(f"adjudication did not pass: {path}")

    inventory_path = root / "06_manifests" / "05_adjudicated_run_inventory.json"
    files = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path != inventory_path:
            files.append({
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
    inventory = {
        "schema_version": "bdev.scwbatac_hvp_sparse_adjudicated_inventory.v1",
        "status": "complete_passed",
        "execution_run_identity": config["execution"]["run_identity"],
        "adjudication_identity": config["adjudication_identity"],
        "adjudication_config_sha256": sha256(config_path),
        "adjudication_report_sha256": sha256(report),
        "adjudication_receipt_sha256": sha256(receipt),
        "implementation_unlocked": False,
        "sif_built": False,
        "files": files,
    }
    inventory_path.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
