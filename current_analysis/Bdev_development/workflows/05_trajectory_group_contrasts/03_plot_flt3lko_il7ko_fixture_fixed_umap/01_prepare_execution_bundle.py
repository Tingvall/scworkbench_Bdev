#!/usr/bin/env python3
"""Validate fixed-UMAP inputs and reflink-clone the projection bundle."""

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

    inputs = config["inputs"]
    checks = (
        (Path(config["runtime"]["sif"]), config["runtime"]["sif_sha256"]),
        (Path(inputs["projection_bundle"]["path"]) / "manifest.json",
         inputs["projection_bundle"]["manifest_sha256"]),
        (Path(inputs["fixed_bdev_bundle"]["path"]) / "manifest.json",
         inputs["fixed_bdev_bundle"]["manifest_sha256"]),
        (Path(inputs["authenticated_l4"]["model_path"]),
         inputs["authenticated_l4"]["model_sha256"]),
        (Path(inputs["authenticated_l4"]["pseudotime_path"]),
         inputs["authenticated_l4"]["pseudotime_sha256"]),
        (Path(inputs["slot_map"]["path"]), inputs["slot_map"]["sha256"]),
        (Path(inputs["mapping_outcomes"]["path"]),
         inputs["mapping_outcomes"]["sha256"]),
    )
    for path, expected in checks:
        observed = sha256(path)
        if observed != expected:
            raise RuntimeError(f"Input checksum changed: {path}: {observed}")

    source = Path(inputs["projection_bundle"]["path"])
    destination = Path(config["storage"]["execution_bundle"])
    if destination.exists():
        raise RuntimeError(f"Refusing to overwrite execution bundle: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["cp", "--recursive", "--reflink=auto", "--",
                    str(source), str(destination)], check=True)
    if sha256(destination / "manifest.json") != inputs["projection_bundle"]["manifest_sha256"]:
        raise RuntimeError("Execution clone differs before neighbor calculation")

    receipt = {
        "schema_version": "bdev.flt3lko_il7ko_fixed_umap_preparation.v1",
        "status": "complete",
        "run_identity": config["run_identity"],
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
        "execution_bundle": str(destination),
        "source_projection_manifest_sha256": sha256(source / "manifest.json"),
        "fixed_bdev_manifest_sha256": sha256(
            Path(inputs["fixed_bdev_bundle"]["path"]) / "manifest.json"),
        "source_mutated": False,
    }
    output = results / "06_manifests" / "01_preparation_receipt.json"
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
