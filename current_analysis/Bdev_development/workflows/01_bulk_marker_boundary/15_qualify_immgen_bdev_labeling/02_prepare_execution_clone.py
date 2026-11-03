#!/usr/bin/env python3
"""Create an active-storage writable clone while preserving the archive."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import yaml


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).resolve().read_text())
    inputs = yaml.safe_load(Path(config["inputs_manifest"]).read_text())
    source = Path(inputs["bdev_reference"]["bundle_path"])
    destination = Path(config["storage"]["execution_clone"])
    results_root = Path(config["storage"]["results_root"])
    receipt_path = results_root / "06_manifests" / "02_execution_clone_receipt.json"
    if destination.exists():
        raise RuntimeError(f"Execution clone already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["cp", "--archive", "--reflink=auto", "--", str(source), str(destination)],
        check=True,
    )
    source_manifest_sha = sha256_file(source / "manifest.json")
    destination_manifest_sha = sha256_file(destination / "manifest.json")
    if source_manifest_sha != destination_manifest_sha:
        raise RuntimeError("Copied manifest differs from the archived source")
    receipt = {
        "schema_version": "bdev.execution_clone_receipt.v1",
        "status": "copied_not_yet_migrated",
        "source": str(source),
        "destination": str(destination),
        "source_manifest_sha256": source_manifest_sha,
        "destination_manifest_sha256_before_migration": destination_manifest_sha,
        "archive_mutated": False,
        "copy_command": ["cp", "--archive", "--reflink=auto"],
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

