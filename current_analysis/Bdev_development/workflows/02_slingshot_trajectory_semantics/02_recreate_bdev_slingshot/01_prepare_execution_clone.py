#!/usr/bin/env python3
"""Prepare an isolated, reflink-capable Bdev bundle clone for trajectory validation."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
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
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    config = yaml.safe_load(config_path.read_text())

    source = Path(config["inputs"]["source_bundle"]["path"])
    destination = Path(config["storage"]["execution_clone"])
    results = Path(config["storage"]["results_root"])
    expected_sha = config["inputs"]["source_bundle"]["manifest_sha256"]
    observed_sha = sha256(source / "manifest.json")
    if observed_sha != expected_sha:
        raise RuntimeError("Archived source manifest checksum changed")
    if destination.exists():
        raise RuntimeError(f"Refusing to overwrite execution clone: {destination}")

    for name in ("01_tables", "02_figures", "03_reports", "04_logs", "06_manifests"):
        (results / name).mkdir(parents=True, exist_ok=True)
    destination.parent.mkdir(parents=True, exist_ok=True)
    Path(config["storage"]["sessions"]).mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["cp", "--archive", "--reflink=auto", "--", str(source), str(destination)],
        check=True,
    )
    clone_sha = sha256(destination / "manifest.json")
    if clone_sha != observed_sha:
        raise RuntimeError("Execution-clone manifest differs from the source")

    receipt = {
        "schema_version": "bdev.slingshot_execution_clone_receipt.v1",
        "status": "copied_not_yet_migrated",
        "run_identity": config["run_identity"],
        "source": str(source),
        "destination": str(destination),
        "source_manifest_sha256": observed_sha,
        "clone_manifest_sha256_before_migration": clone_sha,
        "archive_mutated": False,
        "copy_command": ["cp", "--archive", "--reflink=auto"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
    }
    output = results / "06_manifests" / "01_execution_clone_receipt.json"
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
