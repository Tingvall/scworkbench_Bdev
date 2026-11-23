#!/usr/bin/env python3
"""Inventory the bounded Step 4 qualification outputs."""

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
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    root = Path(config["storage"]["results_root"])
    output = root / "06_manifests" / "04_run_file_inventory.json"
    files = [
        {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(root.rglob("*")) if path.is_file() and path != output
    ]
    required = {
        "01_tables/03_held_out_transfer_predictions.tsv",
        "01_tables/04_fold_metrics.tsv",
        "01_tables/09_selected_immgen_atac_populations.tsv",
        "01_tables/10_selected_immgen_atac_libraries.tsv",
        "03_reports/01_reference_mapping_qualification_summary.json",
        "06_manifests/01_preparation_and_immgen_preflight.json",
        "06_manifests/02_runtime_install.json",
        "06_manifests/03_execution_receipt.json",
    }
    observed = {entry["path"] for entry in files}
    missing = sorted(required - observed)
    if missing:
        raise RuntimeError(f"Required Step 4 outputs are missing: {missing}")
    manifest = {
        "schema_version": "bdev.reference_mapping_run_inventory.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "files": files,
        "file_count": len(files),
        "paper_analysis": False,
        "implementation_unlocked": False,
        "sif_built": False,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "complete", "files": len(files), "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
