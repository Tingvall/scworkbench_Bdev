#!/usr/bin/env python3
"""Verify inputs and clone the mapping bundle for ATAC-UMAP placement."""

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
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    results = Path(config["storage"]["results_root"])
    for category in ("01_tables", "02_figures", "03_reports", "04_logs", "05_objects", "06_manifests"):
        (results / category).mkdir(parents=True, exist_ok=True)
    Path(config["storage"]["sessions"]).mkdir(parents=True, exist_ok=True)

    mapping = Path(config["inputs"]["mapping_bundle"]["path"])
    bdev = Path(config["inputs"]["bdev_bundle"]["path"])
    checks = (
        (mapping / "manifest.json", config["inputs"]["mapping_bundle"]["manifest_sha256"]),
        (bdev / "manifest.json", config["inputs"]["bdev_bundle"]["manifest_sha256"]),
        (Path(config["inputs"]["mapping_table"]["path"]), config["inputs"]["mapping_table"]["sha256"]),
        (Path(config["inputs"]["authenticated_l4"]["path"]), config["inputs"]["authenticated_l4"]["sha256"]),
    )
    for path, expected in checks:
        if sha256(path) != expected:
            raise RuntimeError(f"Input checksum changed: {path}")

    destination = Path(config["storage"]["execution_bundle"])
    if destination.exists():
        raise RuntimeError(f"Refusing to overwrite execution bundle: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["cp", "--recursive", "--reflink=auto", "--", str(mapping), str(destination)], check=True)
    if sha256(destination / "manifest.json") != config["inputs"]["mapping_bundle"]["manifest_sha256"]:
        raise RuntimeError("Execution clone differs before neighbor calculation")

    receipt = {
        "schema_version": "bdev.immgen_atac_umap_preparation.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "mapping_manifest_sha256": sha256(mapping / "manifest.json"),
        "bdev_manifest_sha256": sha256(bdev / "manifest.json"),
        "execution_bundle": str(destination),
        "source_mutated": False,
    }
    output = results / "06_manifests" / "01_preparation_receipt.json"
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

